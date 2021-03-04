/*
 * Copyright 2014 CNRS.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * @file dda_producer.h
 * @brief Read spectra from file, run peakpicking algorithm and push the spectra in the cycle objects
 * @author Marc Dubois marc.dubois@ipbs.fr
 * @author Alexandre Burel alexandre.burel@unistra.fr
 */

#ifndef MZDDAPRODUCER_H
#define MZDDAPRODUCER_H

#include <unordered_map>
#include <set>

#include "windows.h"
#include "exception"

#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/utility/misc/IntegerSet.hpp"

#include "cycle_obj.hpp"

#include "../../utils/mzUtils.hpp"

using namespace std;

namespace mzdb {

/**
 * mzDDAProducer class
 * ===================
 *
 * The role of this class is to:
 *  - read spectra sequentially from pwiz object spectrumList
 *  - build some kind of cycle objects (referred as spectraCollections)
 *  - perform peak-picking
 *  - push cycle objects in a queue
 */
template< class QueueingPolicy,
          class MetadataExtractionPolicy,
          class PeakPickerPolicy,            // ability to launch peak picking process
          class SpectrumListType>          // Wiff or raw proteowizard object
class mzDDAProducer:  QueueingPolicy, MetadataExtractionPolicy {


    //several useful typedefs
    typedef typename QueueingPolicy::Obj SpectraContainer;
    typedef typename QueueingPolicy::UPtr SpectraContainerUPtr;

    typedef typename SpectraContainer::h_mz_t h_mz_t;
    typedef typename SpectraContainer::h_int_t h_int_t;
    typedef typename SpectraContainer::l_mz_t l_mz_t;
    typedef typename SpectraContainer::l_int_t l_int_t;

    typedef std::shared_ptr<mzSpectrum<h_mz_t, h_int_t> > HighResSpectrumSPtr;

private:

    /// peakpicker object
    PeakPickerPolicy m_peakPicker;
    bool m_safeMode;

    ///ms levels encountered
    set<int> m_msLevels;
    
    /// DataMode for each MS level
    map<int, DataMode> m_dataModeByMsLevel;
    
    /// true if user specifically asked to store data with double variables
    bool isNoLoss;

    /**
     * @brief _peakPickAndPush
     * performs peak-picking then push it to the queue
     *
     * @param cycle, cycle number
     * @param filetype origin file type
     * @param params for peak picking algorithm
     */
    void _peakPickAndPush(SpectraContainerUPtr& cycle, pwiz::msdata::CVID filetype,
                          mzPeakFinderUtils::PeakPickerParams& params) {

          this->m_peakPicker.start<h_mz_t, h_int_t, l_mz_t, l_int_t>(cycle, filetype, params);
          this->put(cycle);
    }

    /**
     * @brief _addSpectrum
     * Add spectrum to the current cycle object
     *
     * @param cycle
     * @param scanCount
     * @param cycleCount
     * @param isInHighRes
     * @param spec
     */
    void _addSpectrum(
        SpectraContainerUPtr& cycle,
        int scanCount,
        int cycleCount,
        bool isInHighRes,
        pwiz::msdata::SpectrumPtr spec,
        pwiz::msdata::SpectrumPtr centroidSpec,
        DataMode wantedMode
    ) {
        if (isInHighRes) {
            auto s = std::make_shared<mzSpectrum<h_mz_t, h_int_t> >(scanCount, cycleCount, spec, centroidSpec, wantedMode, m_safeMode);
            s->isInHighRes = isInHighRes;
            cycle->addHighResSpectrum(s);
        } else {
            auto s = std::make_shared<mzSpectrum<l_mz_t, l_int_t> >(scanCount, cycleCount, spec, centroidSpec, wantedMode, m_safeMode);
            s->isInHighRes = isInHighRes;
            cycle->addLowResSpectrum(s);
        }
    }

public:

    /**
     * Read all spectra from data file, performs peak-picking and push cycle objects
     * into the queue
     *
     * @param levelsToCentroid
     * @param spectrumList
     * @param nscans
     * @param bbRtWidth
     * @param filetype
     * @param params
     */
    void _produce(pwiz::util::IntegerSet& levelsToCentroid,
                  SpectrumListType* spectrumList,
                  pair<int, int>& cycleRange,
                  pair<int, int>& rtRange,
                  map<int, double>& bbRtWidth,
                  pwiz::msdata::CVID filetype,
                  mzPeakFinderUtils::PeakPickerParams& params) {

        /* initializing counter */
        int cycleCount = 0, scanCount = 1;

        unordered_map<int, SpectraContainerUPtr> spectraByMsLevel;
        HighResSpectrumSPtr currMs1(nullptr);
        size_t currMs1Id = 0;
        float lastMs1Rt = 0;
		
		std::cout << "\n--- PRODUCER START !";

        size_t size = spectrumList->size();
		boolean view = false;
        for(size_t i = 0; i < spectrumList->size(); i++) {  //JPM add && i<45000 to test performances
            scanCount = i + 1;

			if (scanCount % 1000 == 0 || scanCount >= 45220)
				std::cout << "\n--- READ " << scanCount << " scans !";
            pwiz::msdata::SpectrumPtr spectrum;
            pwiz::msdata::SpectrumPtr centroidSpectrum;
            int msLevel;
            DataMode wantedMode;
            try {
                // Retrieve the input spectrum as is
                spectrum = spectrumList->spectrum(i, true); // size_t index, bool getBinaryData
                // Retrieve the MS level
                msLevel = spectrum->cvParam(pwiz::msdata::MS_ms_level).valueAs<int>();
                
                // increment cycle number
                if (msLevel == 1 ) PwizHelper::checkCycleNumber(filetype, spectrum->id, ++cycleCount);
                float rt = static_cast<float>(spectrum->scanList.scans[0].cvParam(pwiz::msdata::MS_scan_start_time).timeInSeconds());
                if (msLevel == 1 ) lastMs1Rt = rt;
                // do not process if the cycle is not asked
                if(cycleCount < cycleRange.first) continue;
                if(cycleRange.second != 0 && cycleCount > cycleRange.second) break;
                // do not process if the rt is not asked
                if(rt < rtRange.first) continue;
                if(rtRange.second != 0 && rt > rtRange.second) break;
                // also do not process MSn spectra if their precursor has been filtered
                if(lastMs1Rt < rtRange.first) continue;
                
                m_msLevels.insert(msLevel);
                // Retrieve the effective mode
                wantedMode = m_dataModeByMsLevel[msLevel];
                // If effective mode is not profile
                if (wantedMode != PROFILE) {
                    // The content of the retrieved spectrum depends on the levelsToCentroid settings
                    // If the set is empty then we will get a PROFILE spectrum
                    // Else we will obtain a CENTROID spectrum if its msLevel belongs to levelsToCentroid
                    centroidSpectrum = mzdb::getSpectrum<SpectrumListType>(spectrumList, i, true, levelsToCentroid);
                }
            } catch (runtime_error& e) {
                std::cerr << "Runtime exception: " << e.what(); //LOG(ERROR)
            } catch (exception& e) {
				std::cerr << "Catch an exception: " << e.what();//LOG(ERROR)
				std::cerr << "Skipping scan";//LOG(ERROR)
                continue;
            } catch(...) {
				std::cerr << "\nCatch an unknown exception. Trying to recover...";//LOG(ERROR)
                continue;
            }
            
            // set new precursor
			if (scanCount % 1000 == 0 || scanCount >= 45220) std::cout << "\n--- STEP 10 ";
            if (msLevel == 1 ) {
                currMs1 = std::make_shared<mzSpectrum<h_mz_t, h_int_t> >(scanCount, cycleCount, spectrum, centroidSpectrum, wantedMode, m_safeMode);
                currMs1Id = i;
            }

            //init
			if (scanCount % 1000 == 0 || scanCount >= 45220) std::cout << "\n--- STEP 11 ";
            if (spectraByMsLevel.find(msLevel) == spectraByMsLevel.end()) {
                SpectraContainerUPtr spectraCollection(new SpectraContainer(msLevel));
                spectraCollection->parentSpectrum = currMs1;
                spectraByMsLevel[msLevel] = std::move(spectraCollection);
            }
			if (scanCount % 1000 == 0 || scanCount >= 45220) 	std::cout << "\n--- STEP 12 ";
            auto& bbRtWidthBound = bbRtWidth[msLevel];
            const float rt = PwizHelper::rtOf(spectrum);
            if(rt == 0) std::cerr << "Can not find RT for spectrum " << spectrum->id; //LOG(ERROR)
            bool isInHighRes = this->isInHighRes(spectrum, isNoLoss);
            bool added = false;

            //get a reference to the unique pointer corresponding to the current mslevel
			if (scanCount % 1000 == 0 || scanCount >= 45220) std::cout << "\n--- STEP 13 ";
            auto& container = spectraByMsLevel[msLevel];

            if (container->empty() ) {
                currMs1Id == i ? container->addHighResSpectrum(currMs1) : this->_addSpectrum(container, scanCount, cycleCount, isInHighRes, spectrum, centroidSpectrum, wantedMode);
                //this->_addSpectrum(container, scanCount, cycleCount, isInHighRes, spectrum, centroidSpectrum, wantedMode);
                added = true;
            }

            //check if the bounding box is well sized. If yes, launch peak-picking and add it to the queue
			if (scanCount % 1000 == 0 || scanCount >= 45220) std::cout << "\n--- STEP 14 ";


			if (scanCount == 45249) { //JPM
				std::cout << "\nrt:" << rt << "  container.begin:" << container->getBeginRt() << "  bbrtwidthBound:" << bbRtWidthBound;
			}

            if ( (rt - container->getBeginRt()) >= bbRtWidthBound){
				if (scanCount % 1000 == 0 || scanCount >= 45220) 	std::cout << "\n--- STEP 14a1 ";

				if (scanCount >= 45240) { //JPM
					std::cout << "\n--- Bug last time was at 45249 ";
					std::cout << "\nrt:" << rt << "  container.begin:" << container->getBeginRt() << "  bbrtwidthBound:" << bbRtWidthBound;
				}

                this->_peakPickAndPush(container, filetype, params);
				if (scanCount % 1000 == 0 || scanCount >= 45220) 	std::cout << "\n--- STEP 14a2 ";
                //---create a new container
                SpectraContainerUPtr c(new SpectraContainer(msLevel));
                //set its parent
                c->parentSpectrum = currMs1;
                //check if already been added
                if (!added) {
                    currMs1Id == i ? c->addHighResSpectrum(currMs1) : this->_addSpectrum(c, scanCount, cycleCount, isInHighRes, spectrum, centroidSpectrum, wantedMode);
                    //this->_addSpectrum(c, scanCount, cycleCount, isInHighRes, spectrum, centroidSpectrum, wantedMode);
                }
                spectraByMsLevel[msLevel] = std::move(c);
            } else {
				if (scanCount % 1000 == 0 || scanCount >= 45220) std::cout << "\n--- STEP 14b ";
                if (!added) {
                    currMs1Id == i ? container->addHighResSpectrum(currMs1) : this->_addSpectrum(container, scanCount, cycleCount, isInHighRes, spectrum, centroidSpectrum, wantedMode);
                    //this->_addSpectrum(container, scanCount, cycleCount, isInHighRes, spectrum, centroidSpectrum, wantedMode);
                }
            }
            // delete spectra objects
			if (scanCount % 1000 == 0 || scanCount >= 45220) std::cout << "\n--- STEP 15";
            spectrum.reset();
			if (scanCount % 1000 == 0 || scanCount >= 45220) std::cout << "\n--- STEP 16";
            centroidSpectrum.reset();
            
            //genericUtils::checkMemoryUsage(""); // just to avoid memory leaks
        } // end for

        //handles ending, i.e non finished cycles
		std::cout << "\n--- STEP 20";
        for(int i = 1; i <= spectraByMsLevel.size(); ++i) {
            // FIX: when there are only MS2 spectra
            if (spectraByMsLevel.find(i) != spectraByMsLevel.end()) {
                auto& container = spectraByMsLevel.at(i);
                while (container && ! container->empty()) {
                    this->_peakPickAndPush(container, filetype, params);
                }
            }
        }
		std::cout << "\n--- STEP 21"; 
        //just logging if we did not found any spectrea with mslvl = 1
        if (m_msLevels.find(1) == m_msLevels.end()) {
			std::cerr << "Did not see any msLevel 1 !";//LOG(WARNING)
        }

		std::cout << "\n--- STEP 22";
        //signify that we finished producing sending a poison pill
        SpectraContainerUPtr nullContainer(nullptr);
        this->put(nullContainer);
    }

    /**
     * @brief mzDDAProducer
     * @param queue
     * @param mzdbPath
     * @param dataModeByMsLevel
     */
    mzDDAProducer( typename QueueingPolicy::QueueType& queue,
                   MzDBFile& mzdbFile,
                   map<int, DataMode>& dataModeByMsLevel,
                   map<int, double> resolutions,
                   bool safeMode):
        QueueingPolicy(queue),
        MetadataExtractionPolicy(mzdbFile.name),
        isNoLoss(mzdbFile.isNoLoss()), 
        m_peakPicker(dataModeByMsLevel, resolutions),
        m_dataModeByMsLevel(dataModeByMsLevel),
        m_safeMode(safeMode) {
        }

    /**
     * @brief getProducerThread
     * @param levelsToCentroid
     * @param spectrumList
     * @param nscans
     * @param bbWidthManager
     * @param filetype
     * @param params
     * @return thread to join
     */
    boost::thread getProducerThread(pwiz::util::IntegerSet& levelsToCentroid,
                                    SpectrumListType* spectrumList,
                                    pair<int, int>& cycleRange,
                                    pair<int, int>& rtRange,
                                    map<int, double>& bbWidthManager,
                                    pwiz::msdata::CVID filetype,
                                    mzPeakFinderUtils::PeakPickerParams& params) {

        return boost::thread(//boost::bind(
                             &mzDDAProducer<QueueingPolicy,
                             MetadataExtractionPolicy, // TODO: create a policy which claims isInHighRes always true
                             PeakPickerPolicy, // ability to launch peak picking process
                             SpectrumListType>::_produce,
                             this,
                             std::ref(levelsToCentroid),
                             spectrumList,
                             cycleRange,
                             rtRange,
                             std::ref(bbWidthManager),
                             filetype,
                             std::ref(params)
                             );
    } // end function

};


}
#endif // MZDDACONSUMER_H
