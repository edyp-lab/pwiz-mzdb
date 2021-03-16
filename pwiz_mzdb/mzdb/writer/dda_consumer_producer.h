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
  * @file dda_consumer_producer.h
  * @brief Used for monothreading processus : Read spectra from file, run peakpicking algorithm and push the spectra in the cycle objects. Then Read cycle objects from queue, create and insert data in the mzDB file
  * @author Marc Dubois marc.dubois@ipbs.fr
  * @author Alexandre Burel alexandre.burel@unistra.fr
  */

#ifndef MZCONSUMERPRODUCER_H
#define MZCONSUMERPRODUCER_H

  //#include "../lib/sqlite3/sqlite3.h"

#include "cycle_obj.hpp"
#include "spectrum_inserter.h"
#include "bb_inserter.hpp"
#include "bb_computer.hpp"

#include <unordered_map>
#include <set>

#include "windows.h"
#include "exception"

#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/utility/misc/IntegerSet.hpp"


#include "../../utils/mzUtils.hpp"

using namespace std;

namespace mzdb {

	/**
	 * mzConsumerProducer class
	 * ====================
	 *
	 * The role of this class is to:
	 *  - read cycle objects from queue
	 *  - insert spectrum in sqlite table `spectrum`
	 *  - create and insert run slice in sqlite table `run_slice`
	 *  - create and insert bounding box in sqlite table `bounding_box`,
	 *    and sqlite virtual tables `bounding_box_rtree`, `bounding_box_msn_rtree`.
	 *
	 */
	template<class QueueingPolicy, 
		class MetadataExtractionPolicy,
		class PeakPickerPolicy,            // ability to launch peak picking process
		class SpectrumListType>
	class mzConsumerProducer : QueueingPolicy, mzSpectrumInserter, mzBBInserter, MetadataExtractionPolicy {

		//Some typedef
		typedef typename QueueingPolicy::Obj SpectraContainer;

		typedef typename QueueingPolicy::UPtr SpectraContainerUPtr;

		typedef typename SpectraContainer::h_mz_t h_mz_t;
		typedef typename SpectraContainer::h_int_t h_int_t;
		typedef typename SpectraContainer::l_mz_t l_mz_t;
		typedef typename SpectraContainer::l_int_t l_int_t;

		typedef  std::shared_ptr<mzSpectrum<h_mz_t, h_int_t> > HighResSpectrumSPtr;
		typedef  std::shared_ptr<mzSpectrum<l_mz_t, l_int_t> > LowResSpectrumSPtr;



	private:
		// list of all possible data encodings, those who are really used will be inserted on the fly in spectrum_inserter.h

		//SpectraContainerUPtr nullContainer;

		/// peakpicker object
		PeakPickerPolicy m_peakPicker;
		bool m_safeMode;

		///ms levels encountered
		set<int> m_msLevels;

		/// DataMode for each MS level
		map<int, DataMode> m_dataModeByMsLevel;

		/// true if user specifically asked to store data with double variables
		bool isNoLoss;

		int lastPercent = 0;


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
			}
			else {
				auto s = std::make_shared<mzSpectrum<l_mz_t, l_int_t> >(scanCount, cycleCount, spec, centroidSpec, wantedMode, m_safeMode);
				s->isInHighRes = isInHighRes;
				cycle->addLowResSpectrum(s);
			}
		}

	public:

		// DIA
		void _produceDIAMonoThread(
			pwiz::util::IntegerSet& levelsToCentroid,
			SpectrumListType* spectrumList,
			pair<int, int>& cycleRange,
			pair<int, int>& rtRange,
			//map<int, double>& bbRtWidth,
			pwiz::msdata::CVID filetype,
			mzPeakFinderUtils::PeakPickerParams& params,
			pwiz::msdata::MSDataPtr& msdata,
			ISerializer::xml_string_writer& serializer,
			map<int, double>& bbMzWidthByMsLevel,
			map<int, map<int, int> >& runSlices,
			int& progressionCount,
			int spectrumListSize,
			bool progressInformationEnabled
		) {

			/* initialiazing counter */
			int cycleCount = 0, scanCount = 1;
			float lastMs1Rt = 0;

			HighResSpectrumSPtr currMs1(nullptr);
			SpectraContainerUPtr cycle(nullptr);

			for (size_t i = 0; i < spectrumList->size(); i++) {
				try {
					// Retrieve the input spectrum as is
					pwiz::msdata::SpectrumPtr spectrum = spectrumList->spectrum(i, true);
					// Retrieve the MS level
					const int msLevel = spectrum->cvParam(pwiz::msdata::MS_ms_level).valueAs<int>();

					// increment cycle number
					if (msLevel == 1) PwizHelper::checkCycleNumber(filetype, spectrum->id, ++cycleCount);
					float rt = static_cast<float>(spectrum->scanList.scans[0].cvParam(pwiz::msdata::MS_scan_start_time).timeInSeconds());
					if (msLevel == 1) lastMs1Rt = rt;
					// do not process if the cycle is not asked
					if (cycleCount < cycleRange.first) continue;
					if (cycleRange.second != 0 && cycleCount > cycleRange.second) break;
					// do not process if the rt is not asked
					if (rt < rtRange.first) continue;
					if (rtRange.second != 0 && rt > rtRange.second) break;
					// also do not process MSn spectra if their precursor has been filtered
					if (lastMs1Rt < rtRange.first) continue;

					// Retrieve the effective mode
					DataMode wantedMode = m_dataModeByMsLevel[msLevel];
					// If effective mode is not profile
					pwiz::msdata::SpectrumPtr centroidSpectrum;
					if (wantedMode != PROFILE) {
						// The content of the retrieved spectrum depends on the levelsToCentroid settings
						// If the set is empty then we will get a PROFILE spectrum
						// Else we will obtain a CENTROID spectrum if its msLevel belongs to levelsToCentroid
						centroidSpectrum = mzdb::getSpectrum<SpectrumListType>(spectrumList, i, true, levelsToCentroid);
					}

					/*
					 * it should be like this:
					 * cycle = { MS1 - MS2 - MS2 - ... }
					 * but sometimes there might be no primary MS1, first spectra are MS2
					 * in this case, create a fake MS1 parent with spectrum = null
					 */
					if (!cycle && msLevel != 1) {
						// special case, there is no MS1 parent
						std::cout << "Spectrum \"" << spectrum->id << "\" has no precursor, creating a fake one";//LOG(WARNING) 
						cycleCount = 1; // new value should be 1
						// create a fake MS1 parent with spectrum = null
						currMs1 = make_shared<mzSpectrum<h_mz_t, h_int_t> >();
						currMs1->id = scanCount;
						currMs1->cycle = cycleCount;
						currMs1->originalMode = m_dataModeByMsLevel[msLevel]; // I can't know that !
						currMs1->effectiveMode = m_dataModeByMsLevel[msLevel]; // this is wantedMode and not effectiveMode (because I don't know the real originalMode)
						// initialize cycle with it
						cycle = move(SpectraContainerUPtr(new SpectraContainer(2)));
						cycle->parentSpectrum = currMs1;
						cycle->addHighResSpectrum(currMs1); // add parent spectrum to the vector to match the behavior of the dda producer
						++scanCount;
					}

					// put spectra in cycle
					if (msLevel == 1) {
						// peak picks MS2 spectra from previous cycle (if any)
						if (cycle) {
							this->_peakPickAndPush(cycle, filetype, params);
							_consumeMonoThread(msdata, serializer, bbMzWidthByMsLevel, runSlices, progressionCount, spectrumListSize, progressInformationEnabled);

						}
						// make this spectra the precursor of the new cycle
						currMs1 = std::make_shared<mzSpectrum<h_mz_t, h_int_t> >(scanCount, cycleCount, spectrum, centroidSpectrum, wantedMode, m_safeMode);
						cycle = move(SpectraContainerUPtr(new SpectraContainer(2)));
						cycle->parentSpectrum = currMs1;
						cycle->addHighResSpectrum(currMs1); // add parent spectrum to the vector to match the behavior of the dda producer
					}
					else {
						bool isInHighRes = this->isInHighRes(spectrum, isNoLoss);
						if (isInHighRes) {
							auto s = std::make_shared<mzSpectrum<h_mz_t, h_int_t> >(scanCount, cycleCount, spectrum, centroidSpectrum, wantedMode, m_safeMode);
							s->isInHighRes = isInHighRes;
							cycle->addHighResSpectrum(s);
						}
						else {
							std::cout << "Low resolution spectrum is not normal in DIA mode"; //LOG(WARNING)
							auto s = std::make_shared<mzSpectrum<l_mz_t, l_int_t> >(scanCount, cycleCount, spectrum, centroidSpectrum, wantedMode, m_safeMode);
							s->isInHighRes = isInHighRes;
							cycle->addLowResSpectrum(s);
						}
					}
					// delete spectra objects
					spectrum.reset();
					centroidSpectrum.reset();
				}
				catch (runtime_error& e) {
					std::cerr << "Runtime exception: " << e.what();//LOG(ERROR)
				}
				catch (exception& e) {
					std::cerr << "Catch an exception: " << e.what();//LOG(ERROR)
					std::cerr << "Skipping scan"; //LOG(ERROR)
					continue;
				}
				catch (...) {
					std::cerr << "\nCatch an unknown exception. Trying to recover..."; //LOG(ERROR)
					continue;
				}
				scanCount++;

			} // end for

			// ensure everyone is sent to the queue
			if (cycle && !cycle->empty()) {
				this->_peakPickAndPush(cycle, filetype, params);
				_consumeMonoThread(msdata, serializer, bbMzWidthByMsLevel, runSlices, progressionCount, spectrumListSize, progressInformationEnabled);

			}

		}

		// DDA
		void _produceDDAMonoThread(pwiz::util::IntegerSet& levelsToCentroid,
			SpectrumListType* spectrumList,
			pair<int, int>& cycleRange,
			pair<int, int>& rtRange,
			map<int, double>& bbRtWidth,
			pwiz::msdata::CVID filetype,
			mzPeakFinderUtils::PeakPickerParams& params,
			pwiz::msdata::MSDataPtr& msdata,
			ISerializer::xml_string_writer& serializer,
			map<int, double>& bbMzWidthByMsLevel,
			map<int, map<int, int> >& runSlices,
			int& progressionCount,
			int spectrumListSize,
			bool progressInformationEnabled


		) {

			/* initializing counter */
			int cycleCount = 0, scanCount = 1;

			unordered_map<int, SpectraContainerUPtr> spectraByMsLevel;
			HighResSpectrumSPtr currMs1(nullptr);
			size_t currMs1Id = 0;
			float lastMs1Rt = 0;

			//std::cout << "\n--- PRODUCER START !";

			size_t size = spectrumList->size();
			boolean view = false;
			for (size_t i = 0; i < spectrumList->size(); i++) {  //JPM add && i<45000 to test performances
				scanCount = i + 1;

				//if (scanCount % 1000 == 0)
				//	std::cout << "\n--- READ " << scanCount << " scans !";
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
					if (msLevel == 1) PwizHelper::checkCycleNumber(filetype, spectrum->id, ++cycleCount);
					float rt = static_cast<float>(spectrum->scanList.scans[0].cvParam(pwiz::msdata::MS_scan_start_time).timeInSeconds());
					if (msLevel == 1) lastMs1Rt = rt;
					// do not process if the cycle is not asked
					if (cycleCount < cycleRange.first) continue;
					if (cycleRange.second != 0 && cycleCount > cycleRange.second) break;
					// do not process if the rt is not asked
					if (rt < rtRange.first) continue;
					if (rtRange.second != 0 && rt > rtRange.second) break;
					// also do not process MSn spectra if their precursor has been filtered
					if (lastMs1Rt < rtRange.first) continue;

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
				}
				catch (runtime_error& e) {
					std::cerr << "Runtime exception: " << e.what(); //LOG(ERROR)
				}
				catch (exception& e) {
					std::cerr << "Catch an exception: " << e.what();//LOG(ERROR)
					std::cerr << "Skipping scan";//LOG(ERROR)
					continue;
				}
				catch (...) {
					std::cerr << "\nCatch an unknown exception. Trying to recover...";//LOG(ERROR)
					continue;
				}

				// set new precursor
				//if (scanCount % 1000 == 0 || scanCount >= 48460) std::cout << "\n--- STEP 10 ";
				if (msLevel == 1) {
					currMs1 = std::make_shared<mzSpectrum<h_mz_t, h_int_t> >(scanCount, cycleCount, spectrum, centroidSpectrum, wantedMode, m_safeMode);
					currMs1Id = i;
				}

				//init
				//if (scanCount % 1000 == 0 || scanCount >= 48460) std::cout << "\n--- STEP 11 ";
				if (spectraByMsLevel.find(msLevel) == spectraByMsLevel.end()) {
					SpectraContainerUPtr spectraCollection(new SpectraContainer(msLevel));
					spectraCollection->parentSpectrum = currMs1;
					spectraByMsLevel[msLevel] = std::move(spectraCollection);
				}
				//if (scanCount % 1000 == 0 || scanCount >= 48460) 	std::cout << "\n--- STEP 12 ";
				auto& bbRtWidthBound = bbRtWidth[msLevel];
				const float rt = PwizHelper::rtOf(spectrum);
				if (rt == 0) std::cerr << "Can not find RT for spectrum " << spectrum->id; //LOG(ERROR)
				bool isInHighRes = this->isInHighRes(spectrum, isNoLoss);
				bool added = false;

				//get a reference to the unique pointer corresponding to the current mslevel
				//if (scanCount % 1000 == 0 || scanCount >= 48460) std::cout << "\n--- STEP 13 ";
				auto& container = spectraByMsLevel[msLevel];

				if (container->empty()) {
					currMs1Id == i ? container->addHighResSpectrum(currMs1) : this->_addSpectrum(container, scanCount, cycleCount, isInHighRes, spectrum, centroidSpectrum, wantedMode);
					//this->_addSpectrum(container, scanCount, cycleCount, isInHighRes, spectrum, centroidSpectrum, wantedMode);
					added = true;
				}

				//check if the bounding box is well sized. If yes, launch peak-picking and add it to the queue
				//if (scanCount % 1000 == 0 || scanCount >= 48460) std::cout << "\n--- STEP 14 ";


				/*if (scanCount == 48460) { //JPM
					std::cout << "\nrt:" << rt << "  container.begin:" << container->getBeginRt() << "  bbrtwidthBound:" << bbRtWidthBound;
				}*/

				if ((rt - container->getBeginRt()) >= bbRtWidthBound) {
					//if (scanCount % 1000 == 0 || scanCount >= 48460) 	std::cout << "\n--- STEP 14a1 ";

					/*if (scanCount >= 48460) { //JPM
						std::cout << "\n--- Bug last time was at 48460 ";
						std::cout << "\nrt:" << rt << "  container.begin:" << container->getBeginRt() << "  bbrtwidthBound:" << bbRtWidthBound;
					}*/

					this->_peakPickAndPush(container, filetype, params);
					_consumeMonoThread(msdata, serializer, bbMzWidthByMsLevel, runSlices, progressionCount, spectrumListSize, progressInformationEnabled);


					//if (scanCount % 1000 == 0 || scanCount >= 48460) 	std::cout << "\n--- STEP 14a2 ";
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
				}
				else {
					//if (scanCount % 1000 == 0 || scanCount >= 48460) std::cout << "\n--- STEP 14b ";
					if (!added) {
						currMs1Id == i ? container->addHighResSpectrum(currMs1) : this->_addSpectrum(container, scanCount, cycleCount, isInHighRes, spectrum, centroidSpectrum, wantedMode);
						//this->_addSpectrum(container, scanCount, cycleCount, isInHighRes, spectrum, centroidSpectrum, wantedMode);
					}
				}
				// delete spectra objects
				//if (scanCount % 1000 == 0 || scanCount >= 48460) std::cout << "\n--- STEP 15";
				spectrum.reset();
				//if (scanCount % 1000 == 0 || scanCount >= 48460) std::cout << "\n--- STEP 16";
				centroidSpectrum.reset();

				//genericUtils::checkMemoryUsage(""); // just to avoid memory leaks
			} // end for

			//handles ending, i.e non finished cycles
			//std::cout << "\n--- STEP 20";
			for (int i = 1; i <= spectraByMsLevel.size(); ++i) {
				// FIX: when there are only MS2 spectra
				if (spectraByMsLevel.find(i) != spectraByMsLevel.end()) {
					auto& container = spectraByMsLevel.at(i);
					while (container && !container->empty()) {
						this->_peakPickAndPush(container, filetype, params);
						_consumeMonoThread(msdata, serializer, bbMzWidthByMsLevel, runSlices, progressionCount, spectrumListSize, progressInformationEnabled);

					}
				}
			}
			//std::cout << "\n--- STEP 21";
			//just logging if we did not found any spectrea with mslvl = 1
			if (m_msLevels.find(1) == m_msLevels.end()) {
				std::cerr << "Did not see any msLevel 1 !";//LOG(WARNING)
			}

			//std::cout << "\n--- END PRODUCER";


		}



		/////////////////////////////////////////////////
		void _consumeMonoThread(pwiz::msdata::MSDataPtr& msdata,
			ISerializer::xml_string_writer& serializer,
			map<int, double>& bbMzWidthByMsLevel,
			map<int, map<int, int> >& runSlices,
			int& progressionCount,
			int spectrumListSize,
			bool progressInformationEnabled) {


			SpectraContainerUPtr cycleCollection(nullptr);


			this->get(cycleCollection);


		  // TODO check if spectra are correctly inserted (spectrum values are distinguished for PROFILE and CENTROID/FITTED)
			this->insertScans<h_mz_t, h_int_t, l_mz_t, l_int_t>(cycleCollection, msdata, serializer);
			//std::cout << "\n*** Consumer STEP 12";
			vector<std::shared_ptr<mzSpectrum<h_mz_t, h_int_t> > > highResSpectra;
			vector<std::shared_ptr<mzSpectrum<l_mz_t, l_int_t> > > lowResSpectra;
			cycleCollection->getSpectra(highResSpectra, lowResSpectra); //no const since will be deleted in buildAndInsertData

			progressionCount += cycleCollection->size();

			const int& msLevel = cycleCollection->msLevel;

			// TODO check if spectra are correctly inserted (spectrum values are distinguished for PROFILE and CENTROID/FITTED)
			double bbMzWidth = (msLevel == 1 ? bbMzWidthByMsLevel[1] : bbMzWidthByMsLevel[2]);
			//std::cout << "\n*** Consumer STEP 13"; 
			this->buildAndInsertData<h_mz_t, h_int_t, l_mz_t, l_int_t>(msLevel, bbMzWidth, highResSpectra, lowResSpectra, runSlices[msLevel]);
			//std::cout << "\n*** Consumer STEP 14";
			//JPM : end part to remove to disable consumer
			if (progressInformationEnabled) {
				int newPercent = (int)(((float)progressionCount / spectrumListSize * 100.0));
				if (newPercent == lastPercent + 2.0) {
					mzdb::printProgBar(newPercent);
					lastPercent = newPercent;
				}

				if (progressionCount == spectrumListSize) {
					//LOG(INFO) << "Inserter consumer finished: reaches final progression";
					return;
				}
			}

		}

		mzConsumerProducer(typename QueueingPolicy::QueueType& queue,
			MzDBFile& mzdbFile,

			mzParamsCollecter& paramsCollecter, // consumer
			pwiz::msdata::CVID rawFileFormat,
			vector<DataEncoding> dataEncodings,

			map<int, DataMode>& dataModeByMsLevel, // producer
			map<int, double> resolutions,
			bool safeMode) :
			QueueingPolicy(queue),
			MetadataExtractionPolicy(mzdbFile.name),
			isNoLoss(mzdbFile.isNoLoss()),
			m_peakPicker(dataModeByMsLevel, resolutions),
			m_dataModeByMsLevel(dataModeByMsLevel),
			m_safeMode(safeMode),
			
			mzSpectrumInserter(mzdbFile, paramsCollecter, rawFileFormat, dataEncodings),
			mzBBInserter(mzdbFile)
		{
		}



		~mzConsumerProducer() {
			// at the end of the consuming (only one instance of it), write a description of what have been seen and done
			printGlobalInformation();
		}


	};


} // end namespace

#endif // MZCONSUMERPRODUCER_H
