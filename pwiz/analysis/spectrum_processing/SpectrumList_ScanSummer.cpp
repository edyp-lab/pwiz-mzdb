//
// $Id: SpectrumList_ScanSummer.cpp 6567 2014-08-01 20:14:32Z pcbrefugee $
//
//
// Original author: William French <william.french .@. vanderbilt.edu>
//
// Copyright 2008 Spielberg Family Center for Applied Proteomics
//   Cedars-Sinai Medical Center, Los Angeles, California  90048
// Copyright 2008 Vanderbilt University - Nashville, TN 37232
//
// Licensed under the Apache License, Version 2.0 (the "License"); 
// you may not use this file except in compliance with the License. 
// You may obtain a copy of the License at 
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software 
// distributed under the License is distributed on an "AS IS" BASIS, 
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
// See the License for the specific language governing permissions and 
// limitations under the License.
//

#define PWIZ_SOURCE

#include "SpectrumList_ScanSummer.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/data/vendor_readers/Waters/SpectrumList_Waters.hpp"

bool sortFunc (parentIon i, parentIon j) { return (i.mz<j.mz); } // comparator for sorting of parentIon by m/z 
bool sortRtime (precursorGroup i, precursorGroup j) { return (i.rTimeStart<j.rTimeStart); } // comparator for sorting of precursor List by rtime 
bool pIonCompare (parentIon i, double mz) { return (i.mz<mz); } // comparator for searching parentIon by m/z
bool pGroupCompare (precursorGroup i, double mz) { return ( i.precursorMZ < mz ); } // comparator for searching precursorGroups by m/z


namespace pwiz {
namespace analysis {


using namespace msdata;

void SpectrumList_ScanSummer::pushSpectrum(const SpectrumIdentity& spectrumIdentity)
{
    indexMap.push_back(spectrumIdentity.index);
    spectrumIdentities.push_back(spectrumIdentity);
    spectrumIdentities.back().index = spectrumIdentities.size()-1;
}

double SpectrumList_ScanSummer::getPrecursorMz(const msdata::Spectrum& spectrum) const
{
    for (size_t i=0; i<spectrum.precursors.size(); i++)
    {
        for (size_t j=0; j<spectrum.precursors[i].selectedIons.size(); j++)
        {
            CVParam param = spectrum.precursors[i].selectedIons[j].cvParam(MS_selected_ion_m_z);
            if (param.cvid != CVID_Unknown)
                return lexical_cast<double>(param.value);
        }




    }
    return 0;
}

//void SpectrumList_ScanSummer::sumSubScansResample( vector<double> & x, vector<double> & y, size_t refIndex, DetailLevel detailLevel ) const
//{
//
//    if (x.size() != y.size())
//        throw runtime_error("[SpectrumList_ScanSummer::sumSubScans()] x and y arrays must be the same size");
//
//    int nPointsPerDalton = 800;
//
//    // get the index for the precursorGroupList
//    SpectrumPtr s = inner_->spectrum(refIndex, detailLevel ); // get binary data
//    double precursorMZ = getPrecursorMz(*s); // get the precursor value
//
//    vector<precursorGroup>::const_iterator pGroupIt = lower_bound(precursorList.begin(),precursorList.end(),precursorMZ,pGroupCompare); // returns iterator to first element in precursorList where mz >= precursorMZ
//    while ( pGroupIt->indexList[0] != (int)refIndex )
//    {
//        pGroupIt++; // lower_bound returns the first element in a repeat sequence
//        if ( pGroupIt == precursorList.end() )
//            throw runtime_error("[SpectrumList_ScanSummer::sumSubScans()] Cannot find the correct precursorList element...");
//    }
//
//    // setup the intensity bins for summing over multiple sub-scans
//    double MZspacing = 1.0 / double(nPointsPerDalton);
//    vector<double> summedIntensity(int(TotalDaltons/MZspacing)+2,0.0);
//
//    for( vector<int>::const_iterator listIt = pGroupIt->indexList.begin(); listIt != pGroupIt->indexList.end(); ++listIt)
//    {
//
//        SpectrumPtr s2 = inner_->spectrum( *listIt, detailLevel );
//        vector<double>& subMZ = s2->getMZArray()->data;
//        if (subMZ.size() < 2) continue;
//        vector<double>& subIntensity = s2->getIntensityArray()->data;
//
//        int binA, binB = (subMZ[0] - lowerMZlimit) / MZspacing; // initialize 
//        for( size_t j=0, jend=subMZ.size()-1; j < jend ; ++j)
//        {
//
//            binA = binB+1; // get bucket to the right of the first data point
//            binB = (subMZ[j+1] - lowerMZlimit) / MZspacing; // get the bucket to the left of the second data point
//
//            if ( subIntensity[j] == 0 && subIntensity[j+1] == 0 ) continue; // no interpolation needed
//
//            if ( binB < binA )
//            {
//                this->warn_once("[SpectrumList_ScanSummer]: Warning, grid spacing is coarser than raw m/z spacing in at least one case" );
//            }
//
//            int k = binA;
//            while ( k <= binB ) // while loop ensures no interpolation when binB < binA
//            {
//                // get the m/z position of the current bin
//                double mzBin = double(k) * MZspacing + lowerMZlimit;
//
//                // linear interpolation
//                summedIntensity[k] += subIntensity[j] + ( subIntensity[j+1] - subIntensity[j] ) * ( ( mzBin - subMZ[j] ) / ( subMZ[j+1] - subMZ[j] ) ); 
//
//                ++k;
//            }
//
//        }
//
//    }
//
//    x.resize(summedIntensity.size());
//    y.resize(summedIntensity.size());
//    int cnt=0;
//    for (size_t i=1, iend=summedIntensity.size()-1; i < iend ; ++i)
//    {
//        // don't print zero-intensity points flanked by other zero-intensity points
//        if ( summedIntensity[i-1] != 0.0 || summedIntensity[i] != 0.0 || summedIntensity[i+1] != 0.0 )
//        {
//            x[cnt] = lowerMZlimit + double(i) * MZspacing;
//            y[cnt] = summedIntensity[i];
//            cnt++;
//        }
//    }
//    x.resize(cnt);
//    y.resize(cnt);
//
//}


//
// SpectrumList_ScanSummer
//


void SpectrumList_ScanSummer::sumSubScansNaive( vector<double> & x, vector<double> & y, size_t refIndex, DetailLevel detailLevel ) const
{

    if (x.size() != y.size())
        throw runtime_error("[SpectrumList_ScanSummer::sumSubScansNaive()] x and y arrays must be the same size");


    // get the index for the precursorGroupList
    SpectrumPtr s = inner_->spectrum(refIndex, detailLevel ); // get binary data
    double precursorMZ = getPrecursorMz(*s); // get the precursor value

    vector<precursorGroup>::const_iterator pGroupIt = lower_bound(precursorList.begin(),precursorList.end(),precursorMZ,pGroupCompare); // returns iterator to first element in precursorList where mz >= precursorMZ
    while ( pGroupIt->indexList[0] != (int)refIndex )
    {
        pGroupIt++; // lower_bound returns the first element in a repeat sequence
        if ( pGroupIt == precursorList.end() )
            throw runtime_error("[SpectrumList_ScanSummer::sumSubScans()] Cannot find the correct precursorList element...");
    }

    vector<double>& summedMZ = s->getMZArray()->data;
    vector<double>& summedIntensity = s->getIntensityArray()->data;

    vector<int>::const_iterator InitialIt = pGroupIt->indexList.begin();
    InitialIt++;
    if (InitialIt != pGroupIt->indexList.end())
    {

        for( vector<int>::const_iterator listIt = InitialIt; listIt != pGroupIt->indexList.end(); ++listIt)
        {

            SpectrumPtr s2 = inner_->spectrum( *listIt, detailLevel );
            vector<double>& subMz = s2->getMZArray()->data;
            vector<double>& subIntensity = s2->getIntensityArray()->data;

            for( size_t j=0, jend=subMz.size(); j < jend ; ++j)
            {

                // check if this m/z point was recorded from a previous sub-scan
                vector<double>::iterator pIonIt;
                pIonIt = lower_bound(summedMZ.begin(),summedMZ.end(),subMz[j]); // first element that is greater than or equal to subMz[j]
                int indexMZ = pIonIt - summedMZ.begin();
                if (pIonIt == summedMZ.end()) // first check if mzs[j] is outside search range
                {
                    if ( subMz[j] > summedMZ[indexMZ-1] ) // insert at back
                    {
                        summedMZ.push_back(subMz[j]);
                        summedIntensity.push_back(subIntensity[j]);
                    }
                    else // insert at front
                    {
                        summedMZ.insert(summedMZ.begin(),subMz[j]);
                        summedIntensity.insert(summedIntensity.begin(),subIntensity[j]);
                    }
                }
                else if (*pIonIt != subMz[j]) // if the closest value is not equal to mzs[i], start a new m/z point
                {
                    summedMZ.insert(pIonIt,subMz[j]);
                    summedIntensity.insert(summedIntensity.begin()+indexMZ,subIntensity[j]);
                }
                else // m/z value recorded from previous sub-scan for this precursor; calculate the sum
                {
                    summedIntensity[indexMZ] += subIntensity[j];
                }
            

            }

        }

    }

    x.resize(summedIntensity.size());
    y.resize(summedIntensity.size());
    for (size_t i=0, iend=summedIntensity.size(); i < iend ; ++i)
    {
            x[i] = summedMZ[i];
            y[i] = summedIntensity[i];
    }

    

}





PWIZ_API_DECL SpectrumList_ScanSummer::SpectrumList_ScanSummer(const SpectrumListPtr& original, double precursorTol, double rTimeTol)
:   SpectrumListWrapper(original), precursorTol_(precursorTol), rTimeTol_(rTimeTol)
{
    if (!inner_.get()) throw runtime_error("[SpectrumList_ScanSummer] Null pointer");

    // Some parameters
    double precursorMZ = 0.0; 
    
    ms2cnt=0;

    for (size_t i=0, end=inner_->size(); i < end; ++i )
    {
       
        const SpectrumIdentity& spectrumIdentity = inner_->spectrumIdentity(i);
        SpectrumPtr s = inner_->spectrum(i, true); 
        precursorMZ = getPrecursorMz(*s); 

        if (precursorMZ == 0.0) // ms1 scans do not need summing
        {
            pushSpectrum(spectrumIdentity);
            continue;
        }
        double rTime = s->scanList.scans[0].cvParam(MS_scan_start_time).timeInSeconds();
        ms2cnt++;

        if (ms2cnt==1) // set some parameters
        {
            lowerMZlimit = s->scanList.scans[0].scanWindows[0].cvParam(MS_scan_window_lower_limit).valueAs<double>();
            upperMZlimit = s->scanList.scans[0].scanWindows[0].cvParam(MS_scan_window_upper_limit).valueAs<double>();
            TotalDaltons = upperMZlimit - lowerMZlimit;
        }

        vector<precursorGroup>::iterator pGroupIt,prevIt;
        pGroupIt = lower_bound(precursorList.begin(),precursorList.end(),precursorMZ,pGroupCompare); // returns iterator to first element in precursorList where mz >= precursorMZ

        if ( precursorList.empty() )
        {
            pushSpectrum(spectrumIdentity); 
            precursorGroup newGroup;
            newGroup.precursorMZ = precursorMZ;
            newGroup.rTimeStart = rTime;
            newGroup.indexList.push_back(i);
            precursorList.push_back(newGroup);
        }
        else if ( pGroupIt == precursorList.end() ) 
        {

            prevIt = pGroupIt - 1;
            double lowerDiff = precursorMZ - prevIt->precursorMZ;  
            double lrTimeDiff = abs(rTime - prevIt->rTimeStart);

            if ( lowerDiff < precursorTol_ && lrTimeDiff < rTimeTol_ )
            {
                    prevIt->indexList.push_back(i);
            }
            else
            {
                pushSpectrum(spectrumIdentity); 
                precursorGroup newGroup;
                newGroup.precursorMZ = precursorMZ;
                newGroup.rTimeStart = rTime;
                newGroup.indexList.push_back(i);
                precursorList.push_back(newGroup); 
            }
        }
        else if ( pGroupIt == precursorList.begin() )
        {
            double upperDiff = pGroupIt->precursorMZ - precursorMZ; 
            double urTimeDiff = abs(rTime - pGroupIt->rTimeStart); 

            if ( upperDiff < precursorTol_ && urTimeDiff < rTimeTol_ )
            {
                    pGroupIt->indexList.push_back(i);
            }
            else
            {
                pushSpectrum(spectrumIdentity); 
                precursorGroup newGroup;
                newGroup.precursorMZ = precursorMZ;
                newGroup.rTimeStart = rTime;
                newGroup.indexList.push_back(i);
                precursorList.insert(pGroupIt,newGroup); 
            }
        }
        else
        {

            prevIt = pGroupIt - 1;
            double upperDiff = pGroupIt->precursorMZ - precursorMZ; 
            double lowerDiff = precursorMZ - prevIt->precursorMZ; 
            double urTimeDiff = abs(rTime - pGroupIt->rTimeStart); 
            double lrTimeDiff = abs(rTime - prevIt->rTimeStart);

            if ( upperDiff < precursorTol_ && upperDiff <= lowerDiff && urTimeDiff < rTimeTol_ )
            {
                    pGroupIt->indexList.push_back(i);
            }
            else if ( lowerDiff < precursorTol_ && lrTimeDiff < rTimeTol_ )
            {
                    prevIt->indexList.push_back(i);
            }
            else
            {
                pushSpectrum(spectrumIdentity); 
                precursorGroup newGroup;
                newGroup.precursorMZ = precursorMZ;
                newGroup.rTimeStart = rTime;
                newGroup.indexList.push_back(i);
                precursorList.insert(pGroupIt,newGroup); 
            }

        }

    } // end for loop over all spectra

    ms2RetentionTimes = precursorList;
    sort(ms2RetentionTimes.begin(),ms2RetentionTimes.end(),sortRtime);
    ms2cnt = 0; // reset

}

PWIZ_API_DECL size_t SpectrumList_ScanSummer::size() const
{
    return indexMap.size();
}


PWIZ_API_DECL const SpectrumIdentity& SpectrumList_ScanSummer::spectrumIdentity(size_t index) const
{
    return spectrumIdentities.at(index);
}


PWIZ_API_DECL SpectrumPtr SpectrumList_ScanSummer::spectrum(size_t index, bool getBinaryData) const
{
    return spectrum(index, getBinaryData ? DetailLevel_FullData : DetailLevel_FullMetadata);
}


PWIZ_API_DECL SpectrumPtr SpectrumList_ScanSummer::spectrum(size_t index, DetailLevel detailLevel) const
{

    size_t summedScanIndex = indexMap.at(index);
    SpectrumPtr summedSpectrum = inner_->spectrum(summedScanIndex, detailLevel);
    
    if (summedSpectrum->cvParam(MS_ms_level).valueAs<int>() > 1) // MS/MS scan
    {

        try
        {
            // output ms2 spectra by retention time, grab the appropriate spectrum
            int newIndex = ms2RetentionTimes[ms2cnt].indexList[0];
            summedSpectrum = inner_->spectrum(newIndex, detailLevel);
            ms2cnt++;

            vector<double>& mzs = summedSpectrum->getMZArray()->data;
            vector<double>& intensities = summedSpectrum->getIntensityArray()->data;
            sumSubScansNaive( mzs, intensities, newIndex, detailLevel );
            summedSpectrum->defaultArrayLength = mzs.size();

        }
        catch( exception& e )
        {
            throw runtime_error(std::string("[SpectrumList_ScanSummer::spectrum()] Error summing precursor sub-scans: ") + e.what());
        }
    
    }

    summedSpectrum->index = index; // redefine the index
    return summedSpectrum;

}


} // namespace analysis
} // namespace pwiz
