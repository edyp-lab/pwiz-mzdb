﻿//
// $Id: UnifiData.cpp 11392 2017-09-18 18:41:40Z chambm $
//
//
// Original author: Matt Chambers <matt.chambers .@. vanderbilt.edu>
//
// Copyright 2009 Vanderbilt University - Nashville, TN 37232
//
// Licensed under Creative Commons 3.0 United States License, which requires:
//  - Attribution
//  - Noncommercial
//  - No Derivative Works
//
// http://creativecommons.org/licenses/by-nc-nd/3.0/us/
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//


#define PWIZ_SOURCE

#ifndef PWIZ_READER_UNIFI
#error compiler is not MSVC or DLL not available
#else // PWIZ_READER_UNIFI

#pragma unmanaged
#include "UnifiData.hpp"
#include "pwiz/utility/misc/DateTime.hpp"
#include "pwiz/utility/misc/String.hpp"
#include "pwiz/utility/misc/Container.hpp"
#include "pwiz/utility/misc/Exception.hpp"

#pragma warning(push)
#pragma warning(disable: 4400 4538 4564)

#pragma managed
#include "pwiz/utility/misc/cpp_cli_utilities.hpp"
#using <System.dll>
#using <System.Core.dll>
#using <System.Xml.dll>
#using <System.Net.dll>
#using <System.Web.dll>
#using <System.Net.Http.dll>
#using <System.Net.Http.WebRequest.dll>
using namespace pwiz::util;
using namespace System;
using namespace System::Web;
using namespace System::Net;
using namespace System::Net::Http::Headers;
using namespace Newtonsoft::Json;
using namespace Newtonsoft::Json::Linq;
using namespace System::Runtime::Caching::Generic;
using namespace System::Collections::Generic;
using namespace System::Collections::Concurrent;
using System::Threading::Tasks::Task;
using System::Threading::Tasks::TaskScheduler;
using System::Threading::Tasks::Schedulers::QueuedTaskScheduler;
using System::Net::Http::HttpClient;
using System::Uri;
using IdentityModel::Client::TokenClient;
using IdentityModel::Client::TokenResponse;

namespace pwiz {
namespace vendor_api {
namespace UNIFI {

[ProtoBuf::ProtoContract]
enum class ProtoEnergyLevel
{
    Unknown = 0,
    Low = 1,
    High = 2
};

[ProtoBuf::ProtoContract]
enum class ProtoPolarity
{
    Unknown = 0,
    Negative = 1,
    Positive = 2
};

ref class MSeMassSpectrum;
ref class MassSpectrum;

[ProtoBuf::ProtoContract]
[ProtoBuf::ProtoInclude(200, MassSpectrum::typeid)]
public ref class Spectrum abstract
{
public:

    [ProtoBuf::ProtoMember(1)]
    property cli::array<double>^ Intensities;
};

[ProtoBuf::ProtoContract]
[ProtoBuf::ProtoInclude(100, MSeMassSpectrum::typeid)]
public ref class MassSpectrum abstract : Spectrum
{
public:

    [ProtoBuf::ProtoMember(1)]
    property cli::array<double>^ Masses;

    [ProtoBuf::ProtoMember(2)]
    property cli::array<int>^ ScanSize;
};

[ProtoBuf::ProtoContract]
public ref class MSeMassSpectrum : MassSpectrum
{
public:

    [ProtoBuf::ProtoMember(1)]
    property double RetentionTime;

    [ProtoBuf::ProtoMember(2)]
    property ProtoEnergyLevel EnergyLevel;

    [ProtoBuf::ProtoMember(3)]
    property ProtoPolarity IonizationPolarity;
};


ref class ParallelDownloadQueue
{
    generic<typename T, typename TResult>
    ref class Bind1
    {
        initonly T arg;
        Func<T, TResult>^ const f;
        TResult _() { return f(arg); }

        public:
        initonly Func<TResult>^ binder;
        Bind1(Func<T, TResult>^ f, T arg) : f(f), arg(arg) {
            binder = gcnew Func<TResult>(this, &Bind1::_);
        }
    };

    ref class Binder abstract sealed // static
    {
        public:
        generic<typename T, typename TResult>
            static Func<TResult>^ Create(Func<T, TResult>^ f, T arg) {
                return (gcnew Bind1<T, TResult>(f, arg))->binder;
            }
    };

    /// returns a JSON array or protobuf stream of the spectral intensities and masses (if the HTTP Accept header specifies 'application/octet-stream')
    System::String^ spectrumEndpoint(size_t skip, size_t top) { return _sampleResultUrl + "/spectra/mass.mse?$skip=" + skip + "&$top=" + top; }

    Uri^ _sampleResultUrl;
    System::String^ _accessToken;
    HttpClient^ _httpClient;
    int _numSpectra;
    IMemoryCache<int, MSeMassSpectrum^>^ _cache;
    IDictionary<int, Task^>^ _tasksByIndex;
    System::Threading::CancellationTokenSource^ _cancelTokenSource;
    int _chunkSize;
    int _concurrentTasks;
    QueuedTaskScheduler^ _queueScheduler;
    TaskScheduler^ _primaryScheduler;
    TaskScheduler^ _readaheadScheduler;
    System::Threading::EventWaitHandle^ _waitForStart;
    DateTime _startTime;

    //ConcurrentQueue<int>^ _chunkQueue;
    ConcurrentQueue<HttpClient^>^ _httpClients;

    public:
    ParallelDownloadQueue(Uri^ url, System::String^ token, HttpClient^ client, int numSpectra, IMemoryCache<int, MSeMassSpectrum^>^ cache, IDictionary<int, Task^>^ tasksByIndex, int chunkSize, int concurrentTasks)
    {
        _sampleResultUrl = url;
        _accessToken = token;
        _httpClient = client;
        _numSpectra = numSpectra;
        _cache = cache;
        _tasksByIndex = tasksByIndex;
        _httpClients = gcnew System::Collections::Concurrent::ConcurrentQueue<HttpClient^>();
        _cancelTokenSource = gcnew System::Threading::CancellationTokenSource();
        _chunkSize = chunkSize;
        _concurrentTasks = concurrentTasks;
        for (int i = 0; i < concurrentTasks+2; ++i)
        {
            auto webRequestHandler = gcnew System::Net::Http::WebRequestHandler();
            webRequestHandler->UnsafeAuthenticatedConnectionSharing = true;
            webRequestHandler->PreAuthenticate = true;

            auto httpClient = gcnew HttpClient(webRequestHandler);
            httpClient->BaseAddress = gcnew Uri(_sampleResultUrl->GetLeftPart(System::UriPartial::Authority));
            httpClient->DefaultRequestHeaders->Authorization = gcnew AuthenticationHeaderValue("Bearer", _accessToken);
            httpClient->DefaultRequestHeaders->Accept->Add(gcnew MediaTypeWithQualityHeaderValue("application/octet-stream"));
            httpClient->Timeout = System::TimeSpan::FromSeconds(60);

            _httpClients->Enqueue(httpClient);
        }

        _queueScheduler = gcnew QueuedTaskScheduler();
        _primaryScheduler = _queueScheduler->ActivateNewQueue(0);
        _readaheadScheduler = _queueScheduler->ActivateNewQueue(1);

        _waitForStart = gcnew System::Threading::EventWaitHandle(false, System::Threading::EventResetMode::AutoReset);
        _startTime = DateTime::UtcNow;

        //_chunkQueue = gcnew ConcurrentQueue<int>();
        //for (int i = 0; i < _numSpectra; i += chunkSize)
        //    _chunkQueue->Enqueue(i);
    }

    ~ParallelDownloadQueue()
    {
        Console::WriteLine("Disposing queue and cancelling requests.");
        _cancelTokenSource->Cancel();
        //for each (Task^ task in _tasksByIndex->Values)
        //    task->Wait();
    }

    int runChunkTask(size_t taskIndex)
    {
        int currentThreadId = System::Threading::Thread::CurrentThread->ManagedThreadId;
#ifdef _WIN32 //DEBUG
        Console::WriteLine((DateTime::UtcNow - _startTime).ToString("\\[h\\:mm\\:ss\\]\\ ") + "Requesting chunk {0} on thread {1}", taskIndex, currentThreadId);
#endif
        HttpClient^ httpClient = nullptr;
        while (!_httpClients->TryDequeue(httpClient)) {}
        /*if (!_httpClientByThread->ContainsKey(currentThreadId))
        {
            auto webRequestHandler = gcnew System::Net::Http::WebRequestHandler();
            webRequestHandler->UnsafeAuthenticatedConnectionSharing = true;
            webRequestHandler->PreAuthenticate = true;

            httpClient = gcnew HttpClient(webRequestHandler);
            httpClient->BaseAddress = gcnew Uri(_sampleResultUrl->GetLeftPart(System::UriPartial::Authority));
            httpClient->DefaultRequestHeaders->Authorization = gcnew AuthenticationHeaderValue("Bearer", _accessToken);
            httpClient->DefaultRequestHeaders->Accept->Add(gcnew MediaTypeWithQualityHeaderValue("application/octet-stream"));
            httpClient->Timeout = System::TimeSpan::FromSeconds(60);

            _httpClientByThread[currentThreadId] = httpClient;
            Console::WriteLine("CREATING NEW HTTPCLIENT FOR THREAD {0}: {1} clients total", currentThreadId, _httpClientByThread->Count);
        }
        else
            httpClient = _httpClientByThread[currentThreadId];*/

        DateTime start = DateTime::UtcNow;
        System::Net::Http::HttpRequestMessage^ request;
        System::Net::Http::HttpResponseMessage^ response;
        long bytesDownloaded = 0;
        for (int streamRetryCount = 1, streamMaxRetryCount = 15; streamRetryCount <= streamMaxRetryCount; ++streamRetryCount)
        {
            try
            {
                auto requestStart = start;
                for (int requestRetryCount = 1, requestMaxRetryCount = 15; requestRetryCount <= requestMaxRetryCount; ++requestRetryCount)
                {
                    try
                    {
                        requestStart = DateTime::UtcNow;
                        request = gcnew System::Net::Http::HttpRequestMessage(System::Net::Http::HttpMethod::Get, spectrumEndpoint(taskIndex, min(_numSpectra, taskIndex + _chunkSize)));
                        response = httpClient->SendAsync(request, System::Net::Http::HttpCompletionOption::ResponseHeadersRead)->Result;
                        if (response->IsSuccessStatusCode)
                            break;
                    }
                    catch (Exception^ e)
                    {
                        if (requestRetryCount < requestMaxRetryCount)
                        {
                            // try again
#ifdef _WIN32 //DEBUG
                            Console::WriteLine((DateTime::UtcNow - _startTime).ToString("\\[h\\:mm\\:ss\\]\\ ") + System::String::Format("Retrying spectra chunk request {0} on thread {1} (attempt #{3}) due to error ({2})", taskIndex, currentThreadId, e->ToString()->Replace("\r", "")->Split(L'\n')[0], requestRetryCount));
#endif
                            System::Threading::Thread::Sleep(2000 * Math::Pow(2, requestRetryCount));
                        }
                        else
                            throw gcnew Exception(System::String::Format("error requesting spectra chunk {0} on thread {1} ({2})", taskIndex, currentThreadId, e->ToString()->Replace("\r", "")->Split(L'\n')[0]));
                    }
                }
                if (!response->IsSuccessStatusCode)
                    throw gcnew Exception(System::String::Format("error requesting spectra chunk {0} (HTTP {1})", taskIndex, response->StatusCode));

                DateTime stop = DateTime::UtcNow;

                _waitForStart->Set();

#ifdef _WIN32 //DEBUG
                //if (streamRetryCount == 1)
                    Console::WriteLine((DateTime::UtcNow - _startTime).ToString("\\[h\\:mm\\:ss\\]\\ ") + "Starting chunk {0} ({1}ms to send request and read response headers)", taskIndex, (stop - requestStart).TotalMilliseconds);
#endif

                start = DateTime::UtcNow;
                auto lastSpectrum = start;
                for (int i = 0; i < _chunkSize && (taskIndex + i) < _numSpectra; ++i)
                {
                    //if(!_cache->Contains(taskIndex + i))// ((i % 50) == 0)
                    {
                        auto stream = response->Content->ReadAsStreamAsync()->Result;
                        auto spectrum = ProtoBuf::Serializer::DeserializeWithLengthPrefix<MSeMassSpectrum^>(stream, ProtoBuf::PrefixStyle::Base128);
                        if (spectrum == nullptr)
                            throw gcnew Exception(System::String::Format("deserialized null spectrum for index {0} (spectrum {1})", i, taskIndex + i));

                        bytesDownloaded += sizeof(double) * spectrum->Masses->Length * 2;

                        if (bytesDownloaded == 0)
                            throw gcnew Exception(System::String::Format("empty array for index {0} (spectrum {1})", i, taskIndex + i));

                        // if a spectrum takes more than 5 seconds, something's wrong
                        if ((DateTime::UtcNow - lastSpectrum).TotalSeconds > 10.0)
                            throw gcnew Exception(System::String::Format("download is going too slow in chunk {0} on thread {1}: one spectrum ({2}) took {3}s", taskIndex, currentThreadId, taskIndex+i, (DateTime::UtcNow - lastSpectrum).TotalSeconds));
                        lastSpectrum = DateTime::UtcNow;

                        //spectrum->Masses = gcnew cli::array<double>(10);
                        //spectrum->Intensities = gcnew cli::array<double>(10);
                        if (!_cache->Contains(taskIndex + i))
                            _cache->Add(taskIndex + i, spectrum);
                    }
                    /*else
                    {
                        auto spectrum = gcnew MSeMassSpectrum();
                        spectrum->Masses = gcnew cli::array<double>(10);
                        spectrum->Intensities = gcnew cli::array<double>(10);
                        _cache->Add(taskIndex + i, spectrum);
                    }*/
                }
                break; // successfully streamed chunk
            }
            catch (Exception^ e)
            {
                if (streamRetryCount < streamMaxRetryCount)
                {
                    // try again
#ifdef _WIN32 //DEBUG
                    Console::WriteLine((DateTime::UtcNow - _startTime).ToString("\\[h\\:mm\\:ss\\]\\ ") + System::String::Format("Retrying spectra chunk download {0} on thread {1} (attempt #{3}) due to error ({2})", taskIndex, currentThreadId, e->ToString()->Replace("\r", "")->Split(L'\n')[0], streamRetryCount));
#endif
                    System::Threading::Thread::Sleep(2000 * Math::Pow(2, streamRetryCount));
                    bytesDownloaded = 0;
                }
                else
                    throw gcnew Exception(System::String::Format("error deserializing spectra chunk {0} protobuf: {1}", taskIndex, e->ToString()->Replace("\r", "")->Split(L'\n')[0]), e);
            }
            finally
            {
                delete request;

                if (response != nullptr)
                    delete response;
            }
        }
#ifdef _WIN32 //DEBUG
        DateTime stop = DateTime::UtcNow;
        Console::WriteLine((DateTime::UtcNow - _startTime).ToString("\\[h\\:mm\\:ss\\]\\ ") + "FINISHED chunk {0} on thread {1} ({2} bytes in {3}s)", taskIndex, currentThreadId, bytesDownloaded, (stop - start).TotalSeconds);
#endif
        _tasksByIndex->Remove(taskIndex); // remove the task
        _httpClients->Enqueue(httpClient); // add client back to queue

        return 0;
    }

    Task^ getChunkTask(size_t taskIndex, bool primary, bool waitForStart)
    {
        // if the taskIndex is greater than the number of spectra, do nothing
        if (taskIndex >= _numSpectra)
            return nullptr;

        Task^ chunkTask;
        if (!_tasksByIndex->TryGetValue(taskIndex, chunkTask))
        {
            auto f = gcnew Func<size_t, int>(this, &ParallelDownloadQueue::runChunkTask);
            chunkTask = Task::Factory->StartNew(Binder::Create(f, taskIndex), System::Threading::CancellationToken::None, System::Threading::Tasks::TaskCreationOptions::PreferFairness, primary ? _primaryScheduler : _readaheadScheduler);
            if (waitForStart)
                _waitForStart->WaitOne();
            _tasksByIndex->Add(taskIndex, chunkTask);
        }
        return chunkTask;
    }
};

class UnifiData::Impl
{
    public:
    Impl(const std::string& sampleResultUrl) : _acquisitionStartTime(blt::not_a_date_time)
    {
        try
        {
            if (bal::starts_with(sampleResultUrl, "http://") || bal::starts_with(sampleResultUrl, "https://"))
                _sampleResultUrl = gcnew Uri(ToSystemString(sampleResultUrl));
            else
                _sampleResultUrl = gcnew Uri(ToSystemString("https://" + sampleResultUrl));

            auto webRequestHandler = gcnew System::Net::Http::WebRequestHandler();
            webRequestHandler->UnsafeAuthenticatedConnectionSharing = true;
            webRequestHandler->PreAuthenticate = true;

            _httpClient = gcnew HttpClient(webRequestHandler);

            getAccessToken();
            getHttpClient();

            getSampleMetadata();
            getNumberOfSpectra();

            _chunkSize = (int) std::ceil(_numSpectra/100.0);//200;
            _chunkReadahead = 3;
            _cacheSize = _chunkSize * _chunkReadahead * 2;

            _cache = gcnew MemoryCache<int, MSeMassSpectrum^>(_cacheSize);
            _cache->SetPolicy(LruEvictionPolicy<int, MSeMassSpectrum^>::typeid->GetGenericTypeDefinition());

            _tasksByIndex = gcnew ConcurrentDictionary<int, Task^>();
            _queue = gcnew ParallelDownloadQueue(_sampleResultUrl, _accessToken, _httpClient, _numSpectra, _cache, _tasksByIndex, _chunkSize, _chunkReadahead+1);
        }
        CATCH_AND_FORWARD
    }

    friend class UnifiData;

    private:
    const string identityServerBaseAddress = "https://unifiapi.waters.com:50333/identity";
    const string tokenEndpoint = identityServerBaseAddress + "/connect/token";

    /// returns JSON describing the sampleResult, for example:
    //{
    //  "id" : "d686e598-7621-4ad6-9efb-321c015ff59a",
    //  "name" : "MM  5.0 ug/kg",
    //  "description" : "Matrix matched standard 5.0 ug/kg",
    //  "sample" : {
    //    "description": "",
    //    "originalSampleId" : null,
    //    "id" : "f0f33a91-f668-4f9e-a153-4734109cdf41",
    //    "gender" : "",
    //    "name" : "MM2",
    //    "assayConditions" : "",
    //    "bracketGroup" : "",
    //    "dose" : "NaN",
    //    "dosingRoute" : "",
    //    "day" : "",
    //    "eCordId" : "",
    //    "experimentalConcentration" : "",
    //    "groupId" : "",
    //    "injectionId" : "",
    //    "matrix" : "",
    //    "solventDelay" : "NaN",
    //    "species" : "",
    //    "studyId" : "",
    //    "subjectId" : "",
    //    "molForm" : "",
    //    "preparation" : "",
    //    "sampleType" : "Standard",
    //    "batchId" : "",
    //    "studyName" : "",
    //    "sampleLevel" : "Unspecified",
    //    "sampleWeight" : 1,
    //    "dilution" : 1,
    //    "replicateNumber" : 1,
    //    "wellPosition" : "1:C,4",
    //    "injectionVolume" : 10,
    //    "acquisitionRunTime" : 17,
    //    "acquisitionStartTime" : "0001-01-01T00:00:00Z",
    //    "time" : "NaN",
    //    "processingOptions" : "QuantitationStd",
    //    "processingFunction" : "",
    //    "processingSequenceNumber" : 0
    //    }
    //}
    System::String^ sampleResultMetadataEndpoint() { return _sampleResultUrl->AbsoluteUri; }

    /// returns a JSON array describing the 'functions' (aka scan events), for example:
    //{
    //  "value": [
    //    {
    //      "id": "2f52bf3f-d4b5-472d-ab6d-4f55d6b7665e",
    //      "name": "1: TOF MSᴱ (50-1200) 4eV ESI+ - Low CE",
    //      "isCentroidData": false,
    //      "isRetentionData": true,
    //      "isIonMobilityData": false,
    //      "hasCCSCalibration": false,
    //      "detectorType": "MS",
    //      "analyticalTechnique": {
    //        "hardwareName": "",
    //        "scanningMethod": "MS",
    //        "massAnalyser": "QTOF",
    //        "ionisationMode": "+",
    //        "ionisationType": "ESI",
    //        "lowMass": 50,
    //        "highMass": 1200,
    //        "adcGroup": {
    //          "acquisitionMode": "ADC_PD",
    //          "acquisitionFrequency": "NaN",
    //          "ionResponses": [
    //            {
    //              "ionType": "LEU_ENK",
    //              "charge": 1,
    //              "averageIonArea": 33.887678
    //            }
    //          ]
    //        },
    //        "tofGroup": {
    //          "nominalResolution": 7000,
    //          "mseLevel": "Low",
    //          "pusherFrequency": "NaN"
    //        },
    //        "quadGroup": null
    //      },
    //      "axisX": {
    //        "label": "Observed mass",
    //        "unit": "m/z",
    //        "lowerBound": 50,
    //        "upperBound": 1200
    //      },
    //      "axisY": {
    //        "label": "Intensity",
    //        "unit": "Counts",
    //        "lowerBound": "NaN",
    //        "upperBound": "NaN"
    //      }
    //    },
    //    {
    //      "id": "bbcfa906-4545-49ba-aa84-d020eb41e518",
    //      "name": "Reference TOF MSᴱ (50-1200) 6eV ESI+",
    //      "isCentroidData": false,
    //      "isRetentionData": true,
    //      "isIonMobilityData": false,
    //      "hasCCSCalibration": false,
    //      "detectorType": "MS",
    //      "analyticalTechnique": {
    //        "hardwareName": "",
    //        "scanningMethod": "MS",
    //        "massAnalyser": "QTOF",
    //        "ionisationMode": "+",
    //        "ionisationType": "ESI",
    //        "lowMass": 50,
    //        "highMass": 1200,
    //        "adcGroup": {
    //          "acquisitionMode": "ADC_PD",
    //          "acquisitionFrequency": "NaN",
    //          "ionResponses": [
    //            {
    //              "ionType": "LEU_ENK",
    //              "charge": 1,
    //              "averageIonArea": 33.887678
    //            }
    //          ]
    //        },
    //        "tofGroup": {
    //          "nominalResolution": 7000,
    //          "mseLevel": "Unknown",
    //          "pusherFrequency": "NaN"
    //        },
    //        "quadGroup": null
    //      },
    //      "axisX": {
    //        "label": "Observed mass",
    //        "unit": "m/z",
    //        "lowerBound": 50,
    //        "upperBound": 1200
    //      },
    //      "axisY": {
    //        "label": "Intensity",
    //        "unit": "Counts",
    //        "lowerBound": "NaN",
    //        "upperBound": "NaN"
    //      }
    //    },
    //    {
    //      "id": "290805b1-faf2-44ca-af1b-9f8c5d3873ac",
    //      "name": "2: TOF MSᴱ (50-1200) 10-45eV ESI+ - High CE",
    //      "isCentroidData": false,
    //      "isRetentionData": true,
    //      "isIonMobilityData": false,
    //      "hasCCSCalibration": false,
    //      "detectorType": "MS",
    //      "analyticalTechnique": {
    //        "hardwareName": "",
    //        "scanningMethod": "MS",
    //        "massAnalyser": "QTOF",
    //        "ionisationMode": "+",
    //        "ionisationType": "ESI",
    //        "lowMass": 50,
    //        "highMass": 1200,
    //        "adcGroup": {
    //          "acquisitionMode": "ADC_PD",
    //          "acquisitionFrequency": "NaN",
    //          "ionResponses": [
    //            {
    //              "ionType": "LEU_ENK",
    //              "charge": 1,
    //              "averageIonArea": 33.887678
    //            }
    //          ]
    //        },
    //        "tofGroup": {
    //          "nominalResolution": 7000,
    //          "mseLevel": "High",
    //          "pusherFrequency": "NaN"
    //        },
    //        "quadGroup": null
    //      },
    //      "axisX": {
    //        "label": "Observed mass",
    //        "unit": "m/z",
    //        "lowerBound": 50,
    //        "upperBound": 1200
    //      },
    //      "axisY": {
    //        "label": "Intensity",
    //        "unit": "Counts",
    //        "lowerBound": "NaN",
    //        "upperBound": "NaN"
    //      }
    //    }
    //  ]
    //}
    System::String^ functionInfoEndpoint() { return _sampleResultUrl + "/spectrumInfo"; }

    /// returns a simple integer count of the number of spectra in this sampleResult
    System::String^ spectrumCountEndpoint() { return _sampleResultUrl + "/spectra/mass.mse/$count"; }

    /// returns a JSON array or protobuf stream of the spectral intensities and masses (if the HTTP Accept header specifies 'application/octet-stream')
    System::String^ spectrumEndpoint(size_t skip, size_t top) { return _sampleResultUrl + "/spectra/mass.mse?$skip=" + skip + "&$top=" + top; }


    gcroot<Uri^> _sampleResultUrl;
    gcroot<System::String^> _accessToken;
    gcroot<HttpClient^> _httpClient;
    gcroot<ParallelDownloadQueue^> _queue;

    int _numSpectra;
    string _sampleName;
    string _sampleDescription;
    blt::local_date_time _acquisitionStartTime; // UTC

    gcroot<MemoryCache<int, MSeMassSpectrum^>^> _cache;
    gcroot<IDictionary<int, Task^>^> _tasksByIndex;
    int _chunkSize;
    int _chunkReadahead;
    int _cacheSize;

    size_t taskIndexFromSpectrumIndex(size_t index) { return (index / _chunkSize) * _chunkSize; }

    void getAccessToken()
    {
        if (String::IsNullOrEmpty(_sampleResultUrl->UserInfo))
            throw user_error("UserInfo null; username and password must be specified in the sample result URL (e.g. username:password@unifiserver.com:{port}/{sampleResultPath})");
        auto userPassPair = _sampleResultUrl->UserInfo->Split(':');
        if (userPassPair->Length != 2 || String::IsNullOrEmpty(userPassPair[0]) || String::IsNullOrEmpty(userPassPair[1]))
            throw user_error("UserInfo not a pair of values; username and password must be specified in the sample result URL (e.g. username:password@unifiserver.com:{port}/{sampleResultPath})");
        auto username = userPassPair[0];
        auto password = userPassPair[1];

        auto fields = gcnew System::Collections::Generic::Dictionary<System::String^, System::String^>();
        fields->Add(IdentityModel::OidcConstants::TokenRequest::GrantType, IdentityModel::OidcConstants::GrantTypes::Password);
        fields->Add(IdentityModel::OidcConstants::TokenRequest::UserName, username);
        fields->Add(IdentityModel::OidcConstants::TokenRequest::Password, password);
        fields->Add(IdentityModel::OidcConstants::TokenRequest::Scope, "unifi");

        auto tokenClient = gcnew TokenClient(ToSystemString(tokenEndpoint), "resourceownerclient", "secret", IdentityModel::Client::AuthenticationStyle::BasicAuthentication);
        TokenResponse^ response = tokenClient->RequestAsync(fields, System::Threading::CancellationToken::None)->Result;
        if (response->IsError)
            throw user_error("authentication error: incorrect username or password? (" + ToStdString(response->Error) + ")");

        _accessToken = response->AccessToken;
    }

    void getHttpClient()
    {
        _httpClient->DefaultRequestHeaders->Authorization = gcnew AuthenticationHeaderValue("Bearer", _accessToken);
        _httpClient->BaseAddress = gcnew Uri(_sampleResultUrl->GetLeftPart(System::UriPartial::Authority));
    }

    void getSampleMetadata()
    {
        String^ json;
        try
        {
            auto response = _httpClient->GetAsync(sampleResultMetadataEndpoint())->Result;
            if (!response->IsSuccessStatusCode)
                throw std::runtime_error("error getting sample result metadata (" + ToStdString(response->StatusCode.ToString()) + ")");

            json = response->Content->ReadAsStringAsync()->Result;
        }
        catch (Exception^ e)
        {
            throw std::runtime_error("error getting sample result metadata: " + ToStdString(e->ToString()->Split(L'\n')[0]));
        }

        try
        {
            auto o = JObject::Parse(json);
            _sampleName = ToStdString(o->SelectToken("$.name")->ToString());
            _sampleDescription = ToStdString(o->SelectToken("$.description")->ToString());

            auto acquisitionTime = (System::DateTime) o->SelectToken("$.sample.acquisitionStartTime");

            // these are Boost.DateTime restrictions enforced because one of the test files had a corrupt date
            if (acquisitionTime.Year > 10000)
                acquisitionTime = acquisitionTime.AddYears(10000 - acquisitionTime.Year);
            else if (acquisitionTime.Year < 1400)
                acquisitionTime = acquisitionTime.AddYears(1400 - acquisitionTime.Year);

            bpt::ptime pt(boost::gregorian::date(acquisitionTime.Year, boost::gregorian::greg_month(acquisitionTime.Month), acquisitionTime.Day),
                          bpt::time_duration(acquisitionTime.Hour, acquisitionTime.Minute, acquisitionTime.Second, bpt::millisec(acquisitionTime.Millisecond).fractional_seconds()));
            _acquisitionStartTime = blt::local_date_time(pt, blt::time_zone_ptr()); // UTC
        }
        catch (Exception^ e)
        {
            throw std::runtime_error("error parsing sample result metadata: " + ToStdString(e->Message));
        }
        catch (std::exception& e)
        {
            throw e;
        }
    }

    void getNumberOfSpectra()
    {
        auto response = _httpClient->GetAsync(spectrumCountEndpoint())->Result;
        if (!response->IsSuccessStatusCode)
            throw std::runtime_error("error getting number of spectra (" + ToStdString(response->StatusCode.ToString()) + ")");
        if (!Int32::TryParse(response->Content->ReadAsStringAsync()->Result, _numSpectra))
            throw std::runtime_error("error parsing number of spectra \"" + ToStdString(response->Content->ReadAsStringAsync()->Result) + "\"");
    }

    void convertWatersToPwizSpectrum(MSeMassSpectrum^ spectrum, UnifiSpectrum& result, bool getBinaryData)
    {
        result.retentionTime = spectrum->RetentionTime;
        result.scanPolarity = (Polarity)spectrum->IonizationPolarity;
        result.arrayLength = spectrum->Masses->Length;

        if (getBinaryData)
        {
            ToStdVector(spectrum->Masses, result.mzArray);
            ToStdVector(spectrum->Intensities, result.intensityArray);
        }
    }

    void getSpectrum(size_t index, UnifiSpectrum& result, bool getBinaryData)
    {
        try
        {
            //Console::WriteLine("\ngetSpectrum {0}", index);
            MSeMassSpectrum^ spectrum = nullptr;
            int taskIndex = taskIndexFromSpectrumIndex(index);

            if (_cache->TryGet(index, spectrum))
            {
                convertWatersToPwizSpectrum(spectrum, result, getBinaryData);

                // also queue the next 5 chunks
                for (int i = 1; i <= _chunkReadahead; ++i)
                {
                    if ((taskIndex + _chunkSize * i) > _numSpectra)
                        break;
                    if (!_cache->Contains(taskIndex + ((_chunkSize+1) * i - 1))) // if cache contains last index for chunk, don't requeue it
                        _queue->getChunkTask(taskIndex + _chunkSize * i, false, false);
                }

                // remove earlier spectra from the cache
                /*if (index > _chunkSize && (index % (_chunkSize/2)) == 0)
                {
                    for (size_t i = index - _chunkSize/2; i > 0; --i)
                        _cache->Remove(i - 1);
#ifdef _WIN32 //DEBUG
                    Console::WriteLine(System::String::Format("Removed tasks {0}-{1} from cache", "0", index - _chunkSize));
#endif
                }*/

                return;
            }

            // get the primary chunkTask;
            // if cache is empty, wait for primary chunk download to start before starting other chunks
            auto chunkTask = _queue->getChunkTask(taskIndex, true, _cache->Count == 0);

            // also queue the next 5 chunks
            for (int i = 1; i <= _chunkReadahead; ++i)
            {
                if ((taskIndexFromSpectrumIndex(index) + _chunkSize * i) > _numSpectra)
                    break;
                if (!_cache->Contains(taskIndex + ((_chunkSize + 1) * i - 1))) // if cache contains last index for chunk, don't requeue it
                    _queue->getChunkTask(taskIndex + _chunkSize * i, false, false);
            }
#ifdef _WIN32 //DEBUG
            Console::WriteLine("WAITING for chunk {0}", taskIndex);
#endif
            chunkTask->Wait(); // wait for the task to finish

            spectrum = _cache->Get(index);
            if (spectrum == nullptr)
                throw std::runtime_error("spectrum still null after task finished: " + lexical_cast<string>(index));

            convertWatersToPwizSpectrum(spectrum, result, getBinaryData);

        } CATCH_AND_FORWARD_EX(index)
    }

    /*virtual InstrumentModel getInstrumentModel() const;
    virtual IonSourceType getIonSourceType() const;
    virtual blt::local_date_time getSampleAcquisitionTime(int sample, bool adjustToHostTime) const;

    virtual ExperimentPtr getExperiment(int sample, int period, int experiment) const;
    virtual SpectrumPtr getSpectrum(int sample, int period, int experiment, int cycle) const;
    virtual SpectrumPtr getSpectrum(ExperimentPtr experiment, int cycle) const;*/
};


/*struct ExperimentImpl : public Experiment
{
    ExperimentImpl(const WiffFileImpl* wifffile, int sample, int period, int experiment);

    virtual int getSampleNumber() const {return sample;}
    virtual int getPeriodNumber() const {return period;}
    virtual int getExperimentNumber() const {return experiment;}

    virtual size_t getSRMSize() const;
    virtual void getSRM(size_t index, Target& target) const;
    virtual void getSIC(size_t index, std::vector<double>& times, std::vector<double>& intensities) const;
    virtual void getSIC(size_t index, std::vector<double>& times, std::vector<double>& intensities,
                        double& basePeakX, double& basePeakY) const;

    virtual bool getHasIsolationInfo() const;
    virtual void getIsolationInfo(double& centerMz, double& lowerLimit, double& upperLimit) const;

    virtual void getAcquisitionMassRange(double& startMz, double& stopMz) const;
    virtual ScanType getScanType() const;
    virtual ExperimentType getExperimentType() const;
    virtual Polarity getPolarity() const;

    virtual double convertCycleToRetentionTime(int cycle) const;
    virtual double convertRetentionTimeToCycle(double rt) const;

    virtual void getTIC(std::vector<double>& times, std::vector<double>& intensities) const;
    virtual void getBPC(std::vector<double>& times, std::vector<double>& intensities) const;

    const WiffFileImpl* wifffile_;
    gcroot<MSExperiment^> msExperiment;
    int sample, period, experiment;

    size_t transitionCount;
    typedef map<pair<double, double>, pair<int, int> > TransitionParametersMap;
    TransitionParametersMap transitionParametersMap;

    const vector<double>& cycleTimes() const {initializeTIC(); return cycleTimes_;}
    const vector<double>& cycleIntensities() const {initializeTIC(); return cycleIntensities_;}
    const vector<double>& basePeakMZs() const {initializeBPC(); return basePeakMZs_;}
    const vector<double>& basePeakIntensities() const {initializeBPC(); return basePeakIntensities_;}

    private:
    void initializeTIC() const;
    void initializeBPC() const;
    mutable bool initializedTIC_;
    mutable bool initializedBPC_;
    mutable vector<double> cycleTimes_;
    mutable vector<double> cycleIntensities_;
    mutable vector<double> basePeakMZs_;
    mutable vector<double> basePeakIntensities_;
};

typedef boost::shared_ptr<ExperimentImpl> ExperimentImplPtr;


struct SpectrumImpl : public Spectrum
{
    SpectrumImpl(ExperimentImplPtr experiment, int cycle);

    virtual int getSampleNumber() const {return experiment->sample;}
    virtual int getPeriodNumber() const {return experiment->period;}
    virtual int getExperimentNumber() const {return experiment->experiment;}
    virtual int getCycleNumber() const {return cycle;}

    virtual int getMSLevel() const;

    virtual bool getHasIsolationInfo() const;
    virtual void getIsolationInfo(double& centerMz, double& lowerLimit, double& upperLimit) const;

    virtual bool getHasPrecursorInfo() const;
    virtual void getPrecursorInfo(double& selectedMz, double& intensity, int& charge) const;

    virtual double getStartTime() const;

    virtual bool getDataIsContinuous() const {return pointsAreContinuous;}
    size_t getDataSize(bool doCentroid, bool ignoreZeroIntensityPoints) const;
    virtual void getData(bool doCentroid, std::vector<double>& mz, std::vector<double>& intensities, bool ignoreZeroIntensityPoints) const;

    virtual double getSumY() const {return sumY;}
    virtual double getBasePeakX() const {initializeBasePeak(); return bpX;}
    virtual double getBasePeakY() const {initializeBasePeak(); return bpY;}
    virtual double getMinX() const {return minX;}
    virtual double getMaxX() const {return maxX;}

    ExperimentImplPtr experiment;
    gcroot<MassSpectrumInfo^> spectrumInfo;
    mutable gcroot<MassSpectrum^> spectrum;

#if __CLR_VER > 40000000 // .NET 4
    mutable gcroot<cli::array<PeakClass^>^> peakList;
#endif

    int cycle;

    // data points
    double sumY, minX, maxX;
    vector<double> x, y;
    bool pointsAreContinuous;

    // precursor info
    double selectedMz, intensity;
    int charge;

    private:

    mutable double bpX, bpY;
    void initializeBasePeak() const
    {
        if (bpY == -1)
        {
            bpY = experiment->basePeakIntensities()[cycle - 1];
            bpX = bpY > 0 ? experiment->basePeakMZs()[cycle - 1] : 0;
        }
    }
};

typedef boost::shared_ptr<SpectrumImpl> SpectrumImplPtr;*/


/*InstrumentModel WiffFileImpl::getInstrumentModel() const
{
    try
    {
        String^ modelName = sample->Details->InstrumentName->ToUpperInvariant()->Replace(" ", "")->Replace("API", "");
        if (modelName == "UNKNOWN")                 return InstrumentModel_Unknown;
        if (modelName->Contains("2000QTRAP"))       return API2000QTrap; // predicted
        if (modelName->Contains("2000"))            return API2000;
        if (modelName->Contains("2500QTRAP"))       return API2000QTrap; // predicted
        if (modelName->Contains("3000"))            return API3000; // predicted
        if (modelName->Contains("3200QTRAP"))       return API3200QTrap;
        if (modelName->Contains("3200"))            return API3200; // predicted
        if (modelName->Contains("3500QTRAP"))       return API3500QTrap; // predicted
        if (modelName->Contains("4000QTRAP"))       return API4000QTrap;
        if (modelName->Contains("4000"))            return API4000; // predicted
        if (modelName->Contains("QTRAP4500"))       return API4500QTrap;
        if (modelName->Contains("4500"))            return API4500;
        if (modelName->Contains("5000"))            return API5000; // predicted
        if (modelName->Contains("QTRAP5500"))       return API5500QTrap; // predicted
        if (modelName->Contains("5500"))            return API5500; // predicted
        if (modelName->Contains("QTRAP6500"))       return API6500QTrap; // predicted
        if (modelName->Contains("6500"))            return API6500; // predicted
        if (modelName->Contains("QSTARPULSAR"))     return QStarPulsarI; // also covers variants like "API QStar Pulsar i, 0, Qstar"
        if (modelName->Contains("QSTARXL"))         return QStarXL;
        if (modelName->Contains("QSTARELITE"))      return QStarElite;
        if (modelName->Contains("QSTAR"))           return QStar; // predicted
        if (modelName->Contains("TRIPLETOF4600"))   return API4600TripleTOF; // predicted
        if (modelName->Contains("TRIPLETOF5600"))   return API5600TripleTOF;
        if (modelName->Contains("TRIPLETOF6600"))   return API6600TripleTOF; // predicted
        if (modelName->Contains("NLXTOF"))          return NlxTof; // predicted
        if (modelName->Contains("100LC"))           return API100LC; // predicted
        if (modelName->Contains("100"))             return API100; // predicted
        if (modelName->Contains("150MCA"))          return API150MCA; // predicted
        if (modelName->Contains("150EX"))           return API150EX; // predicted
        if (modelName->Contains("165"))             return API165; // predicted
        if (modelName->Contains("300"))             return API300; // predicted
        if (modelName->Contains("350"))             return API350; // predicted
        if (modelName->Contains("365"))             return API365; // predicted
        if (modelName->Contains("X500QTOF"))        return X500QTOF;
        throw gcnew Exception("unknown instrument type: " + sample->Details->InstrumentName);
    }
    CATCH_AND_FORWARD
}*/

/*IonSourceType WiffFileImpl::getIonSourceType() const
{
    try {return (IonSourceType) 0;} CATCH_AND_FORWARD
}

blt::local_date_time WiffFileImpl::getSampleAcquisitionTime(int sample, bool adjustToHostTime) const
{
    try
    {
        setSample(sample);

        System::DateTime acquisitionTime = this->sample->Details->AcquisitionDateTime;
        bpt::ptime pt(boost::gregorian::date(acquisitionTime.Year, boost::gregorian::greg_month(acquisitionTime.Month), acquisitionTime.Day),
            bpt::time_duration(acquisitionTime.Hour, acquisitionTime.Minute, acquisitionTime.Second, bpt::millisec(acquisitionTime.Millisecond).fractional_seconds()));

        if (adjustToHostTime)
        {
            bpt::time_duration tzOffset = bpt::second_clock::universal_time() - bpt::second_clock::local_time();
            return blt::local_date_time(pt + tzOffset, blt::time_zone_ptr()); // treat time as if it came from host's time zone; actual time zone may not be provided by Sciex
        }
        else
            return blt::local_date_time(pt, blt::time_zone_ptr());
    }
    CATCH_AND_FORWARD
}*/


/*ExperimentPtr WiffFileImpl::getExperiment(int sample, int period, int experiment) const
{
    setExperiment(sample, period, experiment);
    ExperimentImplPtr msExperiment(new ExperimentImpl(this, sample, period, experiment));
    return msExperiment;
}


ExperimentImpl::ExperimentImpl(const WiffFileImpl* wifffile, int sample, int period, int experiment)
: wifffile_(wifffile), sample(sample), period(period), experiment(experiment), transitionCount(0), initializedTIC_(false), initializedBPC_(false)
{
    try
    {
        wifffile_->setExperiment(sample, period, experiment);
        msExperiment = wifffile_->msSample->GetMSExperiment(experiment-1);

        if ((ExperimentType) msExperiment->Details->ExperimentType == MRM)
            transitionCount = msExperiment->Details->MassRangeInfo->Length;
    }
    CATCH_AND_FORWARD
}

void ExperimentImpl::initializeTIC() const
{
    if (initializedTIC_)
        return;

    try
    {
        TotalIonChromatogram^ tic = msExperiment->GetTotalIonChromatogram();
        ToStdVector(tic->GetActualXValues(), cycleTimes_);
        ToStdVector(tic->GetActualYValues(), cycleIntensities_);

        initializedTIC_ = true;
    }
    CATCH_AND_FORWARD
}

void ExperimentImpl::initializeBPC() const
{
    if (initializedBPC_)
        return;

    try
    {
        BasePeakChromatogramSettings^ bpcs = gcnew BasePeakChromatogramSettings(0, 0, gcnew array<double>(0), gcnew array<double>(0));
        BasePeakChromatogram^ bpc = msExperiment->GetBasePeakChromatogram(bpcs);
        BasePeakChromatogramInfo^ bpci = bpc->Info;
        ToStdVector(bpc->GetActualYValues(), basePeakIntensities_);

        basePeakMZs_.resize(cycleTimes_.size());
        for (size_t i = 0; i < cycleTimes_.size(); ++i)
            basePeakMZs_[i] = bpci->GetBasePeakMass(i);

        initializedBPC_ = true;
    }
    CATCH_AND_FORWARD
}

size_t ExperimentImpl::getSRMSize() const
{
    try {return transitionCount;} CATCH_AND_FORWARD
}

void ExperimentImpl::getSRM(size_t index, Target& target) const
{
    try
    {
        if (index >= transitionCount)
            throw std::out_of_range("[Experiment::getSRM()] index out of range");

        MRMMassRange^ transition = (MRMMassRange^) msExperiment->Details->MassRangeInfo[index];
        //const pair<int, int>& e = transitionParametersMap.find(make_pair(transition->Q1Mass->MassAsDouble, transition->Q3Mass->MassAsDouble))->second;

        target.type = TargetType_SRM;
        target.Q1 = transition->Q1Mass;
        target.Q3 = transition->Q3Mass;
        target.dwellTime = transition->DwellTime;
        // TODO: store RTWindow?

        {
            // TODO: use NaN to indicate these values should be considered missing?
            target.collisionEnergy = 0;
            target.declusteringPotential = 0;
        }
    }
    CATCH_AND_FORWARD
}

void ExperimentImpl::getSIC(size_t index, std::vector<double>& times, std::vector<double>& intensities) const
{
    double x, y;
    getSIC(index, times, intensities, x, y);
}

void ExperimentImpl::getSIC(size_t index, std::vector<double>& times, std::vector<double>& intensities,
                            double& basePeakX, double& basePeakY) const
{
    try
    {
        if (index >= transitionCount)
            throw std::out_of_range("[Experiment::getSIC()] index out of range");

        ExtractedIonChromatogramSettings^ option = gcnew ExtractedIonChromatogramSettings(index);
        ExtractedIonChromatogram^ xic = msExperiment->GetExtractedIonChromatogram(option);

        ToStdVector(xic->GetActualXValues(), times);
        ToStdVector(xic->GetActualYValues(), intensities);

        basePeakY = xic->MaximumYValue;
        basePeakX = 0;
        for (size_t i=0; i < intensities.size(); ++i)
            if (intensities[i] == basePeakY)
            {
                basePeakX = times[i];
                break;
            }
    }
    CATCH_AND_FORWARD
}

void ExperimentImpl::getAcquisitionMassRange(double& startMz, double& stopMz) const
{
    try
    {
        if ((ExperimentType) msExperiment->Details->ExperimentType != MRM)
        {
            startMz = msExperiment->Details->StartMass;
            stopMz = msExperiment->Details->EndMass;
        }
        else
            startMz = stopMz = 0;
    }
    CATCH_AND_FORWARD
}

ScanType ExperimentImpl::getScanType() const
{
    try {return (ScanType) msExperiment->Details->SpectrumType;} CATCH_AND_FORWARD
}

ExperimentType ExperimentImpl::getExperimentType() const
{
    try {return (ExperimentType) msExperiment->Details->ExperimentType;} CATCH_AND_FORWARD
}

Polarity ExperimentImpl::getPolarity() const
{
    try {return (Polarity) msExperiment->Details->Polarity;} CATCH_AND_FORWARD
}

double ExperimentImpl::convertCycleToRetentionTime(int cycle) const
{
    try {return msExperiment->GetRTFromExperimentScanIndex(cycle);} CATCH_AND_FORWARD
}

double ExperimentImpl::convertRetentionTimeToCycle(double rt) const
{
    try {return msExperiment->RetentionTimeToExperimentScan(rt);} CATCH_AND_FORWARD
}

void ExperimentImpl::getTIC(std::vector<double>& times, std::vector<double>& intensities) const
{
    try
    {
        times = cycleTimes();
        intensities = cycleIntensities();
    }
    CATCH_AND_FORWARD
}

void ExperimentImpl::getBPC(std::vector<double>& times, std::vector<double>& intensities) const
{
    try
    {
        times = cycleTimes();
        intensities = basePeakIntensities();
    }
    CATCH_AND_FORWARD
}

bool ExperimentImpl::getHasIsolationInfo() const
{
    return (ExperimentType)msExperiment->Details->ExperimentType == Product &&
           msExperiment->Details->MassRangeInfo->Length > 0;
}

void ExperimentImpl::getIsolationInfo(double& centerMz, double& lowerLimit, double& upperLimit) const
{
    if (!getHasIsolationInfo())
        return;

    try
    {
        double isolationWidth = ((FragmentBasedScanMassRange^)(msExperiment->Details->MassRangeInfo[0]))->IsolationWindow;

        centerMz = (double)((FragmentBasedScanMassRange^)(msExperiment->Details->MassRangeInfo[0]))->FixedMasses[0];
        lowerLimit = centerMz - isolationWidth / 2;
        upperLimit = centerMz + isolationWidth / 2;
    }
    CATCH_AND_FORWARD
}


SpectrumImpl::SpectrumImpl(ExperimentImplPtr experiment, int cycle)
: experiment(experiment), cycle(cycle), selectedMz(0), bpY(-1), bpX(-1)
{
    try
    {
        spectrumInfo = experiment->msExperiment->GetMassSpectrumInfo(cycle-1);

        pointsAreContinuous = !spectrumInfo->CentroidMode;

        sumY = experiment->cycleIntensities()[cycle-1];
        //minX = experiment->; // TODO Mass range?
        //maxX = spectrum->MaximumXValue;

        if (spectrumInfo->IsProductSpectrum)
        {
            selectedMz = spectrumInfo->ParentMZ;
            intensity = 0;
            charge = spectrumInfo->ParentChargeState;
        }
    }
    CATCH_AND_FORWARD
}

int SpectrumImpl::getMSLevel() const
{
    try {return spectrumInfo->MSLevel == 0 ? 1 : spectrumInfo->MSLevel;} CATCH_AND_FORWARD
}

bool SpectrumImpl::getHasIsolationInfo() const { return experiment->getHasIsolationInfo(); }

void SpectrumImpl::getIsolationInfo(double& centerMz, double& lowerLimit, double& upperLimit) const
{
    if (!getHasIsolationInfo())
        return;

    try
    {
        double isolationWidth = ((FragmentBasedScanMassRange^)(experiment->msExperiment->Details->MassRangeInfo[0]))->IsolationWindow;
        centerMz = getHasPrecursorInfo() ? selectedMz : (double)((FragmentBasedScanMassRange^)(experiment->msExperiment->Details->MassRangeInfo[0]))->FixedMasses[0];
        lowerLimit = centerMz - isolationWidth / 2;
        upperLimit = centerMz + isolationWidth / 2;
    }
    CATCH_AND_FORWARD
}

bool SpectrumImpl::getHasPrecursorInfo() const
{
    return selectedMz != 0;
}

void SpectrumImpl::getPrecursorInfo(double& selectedMz, double& intensity, int& charge) const
{
    selectedMz = this->selectedMz;
    intensity = this->intensity;
    charge = this->charge;
}

double SpectrumImpl::getStartTime() const
{
    return spectrumInfo->StartRT;
}

size_t SpectrumImpl::getDataSize(bool doCentroid, bool ignoreZeroIntensityPoints) const
{
    try
    {
#if __CLR_VER > 40000000 // .NET 4
        if (doCentroid)
        {
            if ((cli::array<PeakClass^>^) peakList == nullptr) peakList = experiment->msExperiment->GetPeakArray(cycle-1);
            return (size_t) peakList->Length;
        }
        else
#endif
        {
            if ((MassSpectrum^) spectrum == nullptr)
            {
                spectrum = experiment->msExperiment->GetMassSpectrum(cycle-1);
#if __CLR_VER > 40000000 // the .NET 4 version has an efficient way to add zeros
                if (!ignoreZeroIntensityPoints && pointsAreContinuous)
                {
                    ExperimentType experimentType = experiment->getExperimentType();
                    if (experimentType != MRM && experimentType != SIM)
                        experiment->msExperiment->AddZeros((MassSpectrum^) spectrum, 1);
                }
#endif
            }
            return (size_t) spectrum->NumDataPoints;
        }
    }
    CATCH_AND_FORWARD
}

void SpectrumImpl::getData(bool doCentroid, std::vector<double>& mz, std::vector<double>& intensities, bool ignoreZeroIntensityPoints) const
{
    try
    {
#if __CLR_VER > 40000000 // .NET 4
        if (doCentroid && pointsAreContinuous)
        {
            if ((cli::array<PeakClass^>^) peakList == nullptr) peakList = experiment->msExperiment->GetPeakArray(cycle-1);
            size_t numPoints = peakList->Length;
            mz.resize(numPoints);
            intensities.resize(numPoints);
            for (size_t i=0; i < numPoints; ++i)
            {
                PeakClass^ peak = peakList[(int)i];
                mz[i] = peak->xValue;
                intensities[i] = peak->area;
            }
        }
        else
#endif
        {
            if ((MassSpectrum^) spectrum == nullptr)
            {
                spectrum = experiment->msExperiment->GetMassSpectrum(cycle-1);
#if __CLR_VER > 40000000 // the .NET 4 version has an efficient way to add zeros
                if (!ignoreZeroIntensityPoints && pointsAreContinuous)
                {
                    ExperimentType experimentType = experiment->getExperimentType();
                    if (experimentType != MRM && experimentType != SIM)
                        experiment->msExperiment->AddZeros((MassSpectrum^) spectrum, 1);
                }
#endif
            }
            ToStdVector(spectrum->GetActualXValues(), mz);
            ToStdVector(spectrum->GetActualYValues(), intensities);
        }
    }
    CATCH_AND_FORWARD
}


SpectrumPtr WiffFileImpl::getSpectrum(int sample, int period, int experiment, int cycle) const
{
    try
    {
        ExperimentPtr msExperiment = getExperiment(sample, period, experiment);
        return getSpectrum(msExperiment, cycle);
    }
    CATCH_AND_FORWARD
}

SpectrumPtr WiffFileImpl::getSpectrum(ExperimentPtr experiment, int cycle) const
{
    SpectrumImplPtr spectrum(new SpectrumImpl(boost::static_pointer_cast<ExperimentImpl>(experiment), cycle));
    return spectrum;
}*/




PWIZ_API_DECL
UnifiData::UnifiData(const std::string& sampleResultUrl)
    : _impl(new Impl(sampleResultUrl))
{
}

PWIZ_API_DECL
UnifiData::~UnifiData()
{
}

PWIZ_API_DECL size_t UnifiData::numberOfSpectra() const { return (size_t) _impl->_numSpectra; }

PWIZ_API_DECL void UnifiData::getSpectrum(size_t index, UnifiSpectrum& spectrum, bool getBinaryData) const { _impl->getSpectrum(index, spectrum, getBinaryData); }

PWIZ_API_DECL const boost::local_time::local_date_time& UnifiData::getAcquisitionStartTime() const { return _impl->_acquisitionStartTime; }
PWIZ_API_DECL const std::string& UnifiData::getSampleName() const { return _impl->_sampleName; }
PWIZ_API_DECL const std::string& UnifiData::getSampleDescription() const { return _impl->_sampleDescription; }

} // ABI
} // vendor_api
} // pwiz

#pragma warning(pop)

#endif