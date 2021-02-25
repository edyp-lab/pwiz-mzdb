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
 * @file peak_finder_zero_bounded.hpp
 * @brief Peakpicking algorithm for Thermo spectra that must be centroided or fitted
 * @author Marc Dubois marc.dubois@ipbs.fr
 * @author Alexandre Burel alexandre.burel@unistra.fr
 */

#ifndef MZPEAKFINDERTHERMO_H
#define MZPEAKFINDERTHERMO_H

#include <fstream>

#include "ceres/ceres.h"
#include "pwiz/data/msdata/MSData.hpp"

#include "data_points_collection.hpp"
#include "ceres_problems.hpp"
#include "peak_models.hpp"

using namespace std;

namespace mzdb {


/**
 * Algorithm adapted to Thermo data
 */
namespace mzPeakFinderZeroBounded {

/**
 *
 */
template<class mz_t, class int_t>
static void findPeaks(const pwiz::msdata::SpectrumPtr& spectrum,
                      vector<std::shared_ptr<Centroid<mz_t, int_t> > >& centroids,
                      mzPeakFinderUtils::PeakPickerParams& params,
                      bool detectPeaks = false,
                      bool computeFWHM = true,
                      mzPeakFinderUtils::USE_CWT useCWT = mzPeakFinderUtils::CWT_DISABLED) {
    
    // do nothing if nothing is requested (given centroids are ok and not to be fitted)
    if(!detectPeaks && !computeFWHM) return;

    const vector<double>& mzs = spectrum->getMZArray()->data;
    const vector<double>& ints = spectrum->getIntensityArray()->data;

    if (mzs.empty() || ints.empty()) {
        // actually, it can happen in practical cases, when  a large file is processed
        //printf("Empty spectrum which is obviously unusual !\n");
        return;
    }

    if (params.adaptiveBaselineAndNoise) {
        const pair<double, double> baselineNoise = mzPeakFinderUtils::getBaselineAndNoise(ints);
        params.baseline = baselineNoise.first;
        params.noise = baselineNoise.second;
    }

    vector<mz_t> mzBuffer;
    vector<int_t> intBuffer;
    vector<std::shared_ptr<Centroid<mz_t, int_t> > > optimizedCentroids;

    size_t i = 0;
    size_t cBufferIdx = 0;
    size_t mzSize = mzs.size();
    
    /*
     * things that must be considered here:
     * - Variables detectPeaks and computeFWHM are true only if spectrum has been acquired in profile mode
     * - If detectPeaks and computeFWHM are false, it is either because the spectrum has been acquired
     *   in centroid mode, or that there is no need to redetect the centroids (and fwhm are not requested)
     *   -> In this case, there is nothing to do (hence the 'return' at the beginning of the function)
     * - This function recentroids and/or computes fwhm (and it means that the spectrum is in profile mode):
     *   > loop on all profile peaks and centroids at the same time
     *   > create a data points collections with all profile peaks between two zero-intensity peaks
     *   > if detectPeaks is true, provide the centroid peaks corresponding to this collection (setDetectedPeaks)
     *   > otherwise, use data_points_collection.detectPeaks
     *   > then use ceres optimization(computeFWHM)
     *     -> if detectPeaks is true computeFWHM is true, we will get redetected centroids with fwhm
     *     -> if detectPeaks is false computeFWHM is true, we will get original centroids with fwhm
     *     -> if detectPeaks is true computeFWHM is false, we will get redetected centroids without fwhm
     *     -> if detectPeaks is false computeFWHM is false, function would have ended already
     * 
     */
    
    mz_t lastMz = 0;
    int_t lastInt = 0;
    int_t almostNothing = (int_t)1e-3;
    while(i < mzSize) {
        /*
         * We are only looking at the current and previous lines
         * 4 possible cases: 
         * 0 -> 0 : nothing to do
         * 0 -> X : new buffers
         * X -> 0 : process buffers and clear
         * X -> Y : append to buffers
         */
        //if(lastInt == 0 && ints[i] == 0) {} else // case 1 : nothing to do
        if(lastInt == 0 && ints[i] > 0) {
            // case 2 : new buffers
            // add an empty item at the beginning of the buffer, if i == 0 create a fake mz
            mzBuffer.push_back(lastMz == 0 ? mzs[i] - almostNothing : lastMz);
            intBuffer.push_back(almostNothing);
            // append the first real point
            mzBuffer.push_back(mzs[i]);
            intBuffer.push_back(ints[i]);
        } else if(lastInt > 0 && ints[i] == 0) {
            // case 3 : process buffers
            // add an empty item at the end of the buffer
            mzBuffer.push_back(mzs[i]);
            intBuffer.push_back(almostNothing);
            // run ceres optimizer
            DataPointsCollection<mz_t, int_t> collec(mzBuffer, intBuffer, spectrum);
            if (detectPeaks) {
                // use home made peak detection (default behaviour is to use the given peaks, detectPeaks should always be false here)
                // TODO: check if home made detection is better than vendor one
                collec.detectPeaks(params, useCWT);
            } else if (!centroids.empty()) {
                // use centroids found with vendor's algorithm
                vector<std::shared_ptr<Centroid<mz_t, int_t> > > centroidBuffer;
                while(cBufferIdx < centroids.size() && centroids[cBufferIdx]->mz <= mzBuffer[mzBuffer.size() - 1]) {
                    centroidBuffer.push_back(centroids[cBufferIdx]);
                    cBufferIdx++;
                }
                // Check if centroidBuffer is empty => not expected to happen
                if (!centroidBuffer.empty()) {
                    collec.setDetectedPeaks(centroidBuffer, params, useCWT);
                    centroidBuffer.clear();
                }
            } else {
                exitOnError("Centroids must be previously computed if peak detection is disabled.");
            }
            collec.optimize(optimizedCentroids, params.optimizationOpt, computeFWHM);
            mzBuffer.clear();
            intBuffer.clear();
        } else if(lastInt > 0 && ints[i] > 0) {
            // case 4 : append to buffers
            mzBuffer.push_back(mzs[i]);
            intBuffer.push_back(ints[i]);
        }
        lastMz = mzs[i];
        lastInt = ints[i];
        i++;
    }
    
    // deal with eventual final peaks
    if(!mzBuffer.empty()) {
        // redo case 3 : process buffers (see comments above)
        mzBuffer.push_back(mzs[i-1] + almostNothing);
        intBuffer.push_back(almostNothing);
        DataPointsCollection<mz_t, int_t> collec(mzBuffer, intBuffer, spectrum);
        if (detectPeaks) {
            collec.detectPeaks(params, useCWT);
        } else if (!centroids.empty()) {
            vector<std::shared_ptr<Centroid<mz_t, int_t> > > centroidBuffer(&centroids[cBufferIdx], &centroids[centroids.size()-1]);
            if (!centroidBuffer.empty()) {
                collec.setDetectedPeaks(centroidBuffer, params, useCWT);
                centroidBuffer.clear();
            }
        } else {
            exitOnError("Centroids must be previously computed if peak detection is disabled on spectrum " + spectrum->id);
        }
        collec.optimize(optimizedCentroids, params.optimizationOpt, computeFWHM);
    }
    // clear buffers
    mzBuffer.clear();
    intBuffer.clear();
    
    // print warning if number of centroid differs, but not if centroids had to be recalculated !
    // peackpeacking is recalculated if spectrum is in PROFILE mode and number or centroids is null or equal to number of profiles data points
    if (detectPeaks && centroids.size() != mzSize && centroids.size() != optimizedCentroids.size()) {
		std::cerr << "Incoherent number of optimized centroids versus input centroids on MS" << spectrum->cvParam(pwiz::msdata::MS_ms_level).valueAs<int>() << " spectrum '" << spectrum->id << "'";//LOG(WARNING) 
    }

    // store the new centroids or fitted peaks
    centroids = optimizedCentroids;
}

}//end namespace mzThermo
}//end namespace mzdb





#endif // MZPEAKFINDERTHERMO_H


