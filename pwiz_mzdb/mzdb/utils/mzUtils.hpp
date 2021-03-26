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
 * @file mzUtils.hpp
 * @brief Global dependency regrouping functions used by many classes
 * @author Marc Dubois marc.dubois@ipbs.fr
 * @author Alexandre Burel alexandre.burel@unistra.fr
 */

#ifndef MZUTILS_HPP
#define MZUTILS_HPP

#define GLOG_NO_ABBREVIATED_SEVERITIES

#include <iostream>
#include <math.h>
#include <set>

#include <regex>

#include "pwiz/data/msdata/MSData.hpp"
#include "pwiz/analysis/peakdetect/PeakFamilyDetectorFT.hpp"
#include "pwiz/analysis/spectrum_processing/PrecursorRecalculatorDefault.hpp"

#include "libraries/SQLite/sqlite3.h"
#include "glog/logging.h"

// next two lines are required for the genericUtils namespace
//#include "windows.h"
//#include "psapi.h"

//#include "glog/logging.h"

///** START OF TYPENAME ADDON (TEMPORARY) **/
//#include <cstddef>
//#include <stdexcept>
//#include <cstring>
//#include <ostream>
//
//#ifndef _MSC_VER
//#  if __cplusplus < 201103
//#    define CONSTEXPR11_TN
//#    define CONSTEXPR14_TN
//#    define NOEXCEPT_TN
//#  elif __cplusplus < 201402
//#    define CONSTEXPR11_TN constexpr
//#    define CONSTEXPR14_TN
//#    define NOEXCEPT_TN noexcept
//#  else
//#    define CONSTEXPR11_TN constexpr
//#    define CONSTEXPR14_TN constexpr
//#    define NOEXCEPT_TN noexcept
//#  endif
//#else  // _MSC_VER
//#  if _MSC_VER < 1900
//#    define CONSTEXPR11_TN
//#    define CONSTEXPR14_TN
//#    define NOEXCEPT_TN
//#  elif _MSC_VER < 2000
//#    define CONSTEXPR11_TN constexpr
//#    define CONSTEXPR14_TN
//#    define NOEXCEPT_TN noexcept
//#  else
//#    define CONSTEXPR11_TN constexpr
//#    define CONSTEXPR14_TN constexpr
//#    define NOEXCEPT_TN noexcept
//#  endif
//#endif  // _MSC_VER
///** END OF TYPENAME ADDON (TEMPORARY) **/

using namespace std;

namespace mzdb {



#define __BEGIN_MZDB_NM namespace mzdb {

#define __END_MZDB_NM }

#define USE_STD_NM using namespace std;

//
#define MS1_BB_MZ_WIDTH_STR "ms1_bb_mz_width"

#define MSN_BB_MZ_WIDTH_STR "msn_bb_mz_width"

#define MS1_BB_TIME_WIDTH_STR "ms1_bb_time_width"

#define MSN_BB_TIME_WIDTH_STR "msn_bb_time_width"

#define IS_LOSSLESS_STR "is_lossless"

#define ORIGIN_FILE_FORMAT_STR "origin_file_format"

//xml attributes
#define XML_FLOAT "xsd:float"

#define XML_STRING "xsd:string"

#define XML_DOUBLE "xsd:double"

#define XML_INT "xsd:int"

#define XML_BOOLEAN "xsd:boolean"

#define PARAM_NAME_STR "name"

#define PARAM_VALUE_STR "value"

#define PARAM_TYPE_STR "type"

#define MS_NS_STR "MS"

#define CV_REF_STR "cvRef"

#define ACCESS_NB_STR "accession"

#define PARAMS_STR "params"

#define CV_PARAMS_STR "cvParams"

#define CV_PARAM_STR "cvParam"

#define USER_PARAMS_STR "userParams"

#define USER_PARAM_STR "userParam"

#define USER_TEXTS_STR "userTexts"

#define USER_TEXT_STR "userText"

#define TRUE_STR "true"

#define FALSE_STR "false"

#define EMPTY_STR "empty"

#define CID_STR "CID"

#define ETD_STR "ETD"

#define HCD_STR "HCD"

#define UNKNOWN_STR "unknown"

#define IN_HIGH_RES_STR "in_high_res"

#define THERMO_TRAILER "[Thermo Trailer Extra]Monoisotopic M/Z:"

#define THERMO_DATA_PROC "pwiz_mzdb_thermo_conversion"

#define ABI_DATA_PROC "pwiz_mzdb_abi_conversion"

#define ABI_SWATH_DATA_PROC "pwiz_mzdb_abi_swath_conversion"

#define XML_DATA_PROC "pwiz_mzdb_xml_conversion"


#define _64_BIT_MZ "64_bit_float_mz"

#define _64_BIT_INTENSITY "64_bit_float_intensity"

#define _32_BIT_MZ "32_bit_float_mz"

#define _32_BIT_INTENSITY "32_bit_float_intensity"

#define PROFILE_STR "profile"


/*functions*/
#define PPM2MZ(mz, ppm) mz * ppm / 1e6

#define MZ2PPM(mz, ppm) mz * 1e6 / ppm;

/* versionning */
#define BUILD_VERSION "build_version"

/* resolutions */
#define RESOLUTIONS_STR "resolutions"
#define DEFAULT_RESOLUTION 20000

/**useful typedef */
typedef unsigned char byte;
typedef unsigned int mz_uint;


///A spectrum can be stored in different ways:
///
/// -`Profile`: the spectrum has been acquired in profile and the user want to keep
/// it in profile (all data points)
///
/// -`Fitted`: the spectrum has been acquired in profile and the user want to perform a
/// peak picking on it. The `Fitted` mode is an intermediary mode between profile and
/// centroid. Each mass peak is modelized to keep most of the information: its width at half
/// maximum and m/z value plus intensity value at the apex. This mode is a good trade-off between
/// file data size and peak informations.
///
/// -`Centroid`: the spectrum could already be acquired in centroid mode. In case, it is in profile,
/// a peak picking method is performed usually preferring `vendor` algorithm
enum DataMode {
    PROFILE = 1,
    CENTROID = 12,
    FITTED = 20
};

/// a mass peak can be encoded in three different ways:
///
/// -if this a low resolution spectrum (i.e. a MS2 spectrum for exmaple), it is preferred
/// to store it with m/z and intensity values encoded in 32 bits in order to save space: enum value of
/// 8 (4 bytes + 4 bytes)
///
/// -if this a high resolution spectrum (i.e. MS1 spectrum most of the time) a high precision in
/// m/z dimension is required: it will be stored using a 64 bits (double in moder desktop architecture):
/// enum value of 12 (4 bytes +  8 bytes)
///
/// -For some reasons, you may want to encode both m/z and intensities in 64 bits: enum value of 16
/// (8 bytes + 8 bytes)
enum PeakEncoding {
    LOW_RES_PEAK = 8,
    HIGH_RES_PEAK = 12,
    NO_LOSS_PEAK = 16
};


/// A DataEncoding is the combination of the DataMode @see DataMode, the PeakEncoding @see PeakEncoding
/// and a string which indicates if data are compressed or not. This strucure is used both in reading
/// and writing tasks; so its usage can vary slightly through cases.
struct DataEncoding {

    /// sqlite DB row id of `data_encoding` table
    int id;

    ///mode @see DataMode
    DataMode mode;

    ///peakEncoding @see PeakEncoding
    PeakEncoding peakEncoding;

    ///compression
    string compression;
    
    bool hasBeenInserted;

    /**
     * DataEncoding of one or several spectrum(a)
     * @param id_
     * @param mode_ @see DataMode
     * @param pe_ @see PeakEncoding
     * @param compression default to `none`, no compression applied
     */
    DataEncoding(int id_, DataMode mode_, PeakEncoding pe_) :
        id(id_),
        mode(mode_),
        peakEncoding(pe_),
        compression("none"),
        hasBeenInserted(false) {}

    ///Default constructor
    DataEncoding(){}

    ///Create `SQL` for inserting a row in `data_encoding` table with current member values
    string buildSQL() {
        string mzPrec, intPrec;
        if (peakEncoding == NO_LOSS_PEAK) {
            mzPrec = "64"; intPrec = "64";
        } else if (peakEncoding == HIGH_RES_PEAK) {
            mzPrec = "64"; intPrec = "32";
        } else if (peakEncoding == LOW_RES_PEAK) {
            mzPrec = "32"; intPrec = "32";
        }

        string mode_str;
        if (mode == PROFILE)
            mode_str = "profile";
        else if (mode == FITTED)
            mode_str = "fitted";
        else if (mode == CENTROID)
            mode_str = "centroid";

        //return "INSERT INTO data_encoding VALUES (NULL, '" + mode_str + "', '"+ compression + "', 'little_endian', "+ mzPrec + ", " + intPrec + ", NULL);";
        return "INSERT INTO data_encoding VALUES (" + std::to_string(static_cast<long long>(id)) + ", '" + mode_str + "', '"+ compression + "', 'little_endian', "+ mzPrec + ", " + intPrec + ", NULL);";
    }
    
    void setHasBeenInserted(bool _hasBeenInserted) {
        hasBeenInserted = _hasBeenInserted;
    }
};


/// put data into a vector of bytes
/// To be more effective the vector should have some reserved space
template <typename T>
inline static void put(T data, vector<byte>& v) {
    unsigned int wpos_ = v.size();
    unsigned int s = sizeof (data);
    if (v.size() < (wpos_ + s))
        v.resize(wpos_ + s);
    memcpy(&v[wpos_], (byte*) & data, s);
}

/// getting data from the bytebuffer
template<typename T>
inline static T get(unsigned int index, vector<byte> &buffer) {
    return *((T*) & buffer[index]);
}

///Utility function to find if a value already exists in the map
template<typename T>
inline static bool isInMapKeys(int value, map<int, T>& m) {
    typename map<int, T>::iterator it;
    it = m.find(value);
    if (it == m.end())
        return false; //not found
    return true; //found
}

/// more generic version than the previous one
template <typename T> static bool isInMapKeys(typename T::key_type value, T& m) {
    typename T::iterator it;
    it = m.find(value);
    if (it == m.end())
        return false;
    return true;
}


/**
 * regroups functions taking in parameters pwiz::msdata::SpectrumPtr
 * gather some metadata on it
 */
namespace PwizHelper {

inline static float rtOf(const pwiz::msdata::SpectrumPtr& s) {
    pwiz::msdata::Scan* scan = !s->scanList.empty() ? &s->scanList.scans[0] : 0;
    if (! scan) {
        //printf("\nCan not find RT of one spectrum !\n");
        return 0;
    }
    pwiz::msdata::CVParam scanTimeParam = scan ? scan->cvParam(pwiz::msdata::MS_scan_start_time) : pwiz::msdata::CVParam();
    return static_cast<float> (scanTimeParam.timeInSeconds());
}

inline static int precursorChargeOf(const pwiz::msdata::SpectrumPtr &s) {
    const pwiz::msdata::SelectedIon& si = s->precursors.front().selectedIons.front();
    return si.cvParam(pwiz::msdata::MS_charge_state).valueAs<int>();
}

/**
 * get the precusor mz as double value
 * @param s pwiz spectrum shared pointer
 * @return precusor mz as double value
 */
inline static double precursorMzOf(const pwiz::msdata::SpectrumPtr &s) {
    const pwiz::msdata::SelectedIon& si = s->precursors.front().selectedIons.front();
    return si.cvParam(pwiz::msdata::MS_selected_ion_m_z).valueAs<double>();
}


/**
 * @brief checkCycleNumber
 * checks if cycle number matches the "real" cycle number given in the spectrum title
 *
 * @param filetype origin file type (titles are differents for each file type)
 * @param spectrumTitle in which the good cycle number may be found
 * @param cycleNumber, the computed cycle number that may be wrong
 */
inline void checkCycleNumber(pwiz::msdata::CVID filetype, string spectrumTitle, int& cycleNumber) {
    // add a special case for ABSciex files
    if(filetype == pwiz::msdata::MS_ABI_WIFF_format) {
        std::regex e ("cycle=(\\d+)");
        std::smatch match;
        if (std::regex_search(spectrumTitle, match, e) && match.size() > 0) {
            // real cycle extracted from title
            // if the spray gets lost, the spectra will not be written in the wiff file and the cycles wont be a list from 1 to the end
            int cycle = std::stoi(match.str(1));
            if(cycle != cycleNumber) {
                //printf("Cycle changed from %d to %d\n", cycleNumber, cycle);
                cycleNumber = cycle;
            }
        }
    }
}

inline int extractScanNumber(string spectrumTitle) {
    std::regex e ("scan=(\\d+)");
    std::smatch match;
    if (std::regex_search(spectrumTitle, match, e) && match.size() > 0) {
        return std::stoi(match.str(1));
    }
    return 0;
}

inline bool containsAllItems(std::vector<double> bigVector, std::vector<double> smallVector) {
    // i must know if all items in smallVector are in bigVector
    for(size_t i = 0; i < smallVector.size(); i++) {
        float item = smallVector[i];
        bool isFound = false;
        size_t j = 0;
        while(j < bigVector.size() && !isFound) {
            if(bigVector[j] == item) isFound = true;
            j++;
        }
        if(!isFound) return false;
    }
    return true;
}
inline bool contains(std::vector<double> v1, std::vector<double> v2) {
    // irrelevant if vectors are too small
    if(v1.size() + v2.size() <= 2) return false;
    if(v1.size() >= v2.size()) {
        return containsAllItems(v1, v2);
    } else {
        return containsAllItems(v2, v1);
    }
}

inline vector<double> determineIsolationWindowStarts(pwiz::msdata::SpectrumListPtr spectrumList) {
    
    int spectrumListSize = spectrumList->size();
    
    size_t nbMS1SpectraToCheck = 10; // 2 should suffice, but 10 is safer
    size_t nbMS1SpectraChecked = 0;
    size_t nbDiaLikeMS1Spectra = 0; // at the end, should be equal to nbMS1SpectraToCheck
    size_t nbSpectraWithoutExpectedCvParams = 0;
    size_t nbSpectraWithoutExpectedCvParamsToCheck = 100;
    // reference and candidate values
    vector<double> refTargets, cndTargets;
    
    for(size_t i = 0; i < spectrumList->size(); i++) {
        pwiz::msdata::SpectrumPtr spectrum = spectrumList->spectrum(i, false);
        const int& msLevel = spectrum->cvParam(pwiz::msdata::MS_ms_level).valueAs<int>();
        if(msLevel == 1) {
            // do nothing if first or second MS1 (there must be something to compare)
            if(cndTargets.size() == 0) continue;
            if(refTargets.size() == 0) refTargets = cndTargets;
            // compare target and candidates
            if(contains(refTargets, cndTargets)) nbDiaLikeMS1Spectra++;
            // keep the candidate if it is larger than the reference (a DIA reference was a subset of the candidate, if DDA it doesnt matter if the reference changes)
            if(refTargets.size() < cndTargets.size()) refTargets = cndTargets;
            // reset the candidate
            cndTargets.clear();
            // increment the counter
            nbMS1SpectraChecked++;
        } else if(msLevel == 2) {
            // look at all the MS2 spectra for the current MS1 and try to get the target value
            // there may be more than one precursor on multiplexed data (but this algorithm does not seem to work on multiplexed data anyway)
			//###VDS TimsTof: Was 
			/*pwiz::msdata::Precursor prec = spectrum->precursors[0];
			if (prec.hasCVParam(pwiz::msdata::MS_isolation_window_target_m_z)) {
				cndTargets.push_back(stod(prec.cvParam(pwiz::msdata::MS_isolation_window_target_m_z).value));
			}
			else {
				if (prec.isolationWindow.hasCVParam(pwiz::msdata::MS_isolation_window_target_m_z)) {
					cndTargets.push_back(stod(prec.isolationWindow.cvParam(pwiz::msdata::MS_isolation_window_target_m_z).value));
				}
				else {
					nbSpectraWithoutExpectedCvParams++;
				}
			}*/
			//###VDS TimsTof: End Was 
			//###VDS TimsTof: Add size Verif
			if (spectrum->precursors.size() == 0) {
				nbSpectraWithoutExpectedCvParams++;
			} else {
				pwiz::msdata::Precursor prec = spectrum->precursors[0];
				if (prec.hasCVParam(pwiz::msdata::MS_isolation_window_target_m_z)) {
					cndTargets.push_back(stod(prec.cvParam(pwiz::msdata::MS_isolation_window_target_m_z).value));
				} else {
					if (prec.isolationWindow.hasCVParam(pwiz::msdata::MS_isolation_window_target_m_z)) {
						cndTargets.push_back(stod(prec.isolationWindow.cvParam(pwiz::msdata::MS_isolation_window_target_m_z).value));
					} else {
						nbSpectraWithoutExpectedCvParams++;
					}
				}
			}
			//###VDS TimsTof: END
        }
        if(nbMS1SpectraChecked >= nbMS1SpectraToCheck)
            break;
        if(nbSpectraWithoutExpectedCvParams >= nbSpectraWithoutExpectedCvParamsToCheck)
            break;
    }
    // check what has been seen
    if(nbDiaLikeMS1Spectra != nbMS1SpectraToCheck) {
        refTargets.clear();
    }
    return refTargets;
}

} // end namespace PWIZ HELPER


/**
 * Get simple activation code as string for pwiz ``Activation`` object.
 * -HCD, CID, ETD
 * 
 * Warning: the current version of ProteoWizard does not handle well ETD metadata
 * cf. https://sourceforge.net/p/proteowizard/mailman/message/35141275/
 * "Some work was done a while ago to get better activation energy for HCD, but we never even got into ETD"
 * It means that ETD analyses will get considered like CID :(
 *
 * @param a pwiz activation object contained in _pwiz precursor_ object
 * @return string representing activation to be inserted in the database
 */
static inline string getActivationCode(const pwiz::msdata::Activation& a) {
    if (a.empty())
        return EMPTY_STR;
    /// CID (collision-induced dissociation): The dissociation of an ion after collisional excitation. The term collisional-activated dissociation is not recommended.
    if (a.hasCVParam(pwiz::msdata::MS_CID) && ! a.hasCVParam(pwiz::msdata::MS_ETD))
        return CID_STR;
    /// ETD (electron transfer dissociation): A process to fragment ions in a mass spectrometer by inducing fragmentation of cations (e.g. peptides or proteins) by transferring electrons to them.
    else if (a.hasCVParam(pwiz::msdata::MS_ETD)) //electron_transfer_dissociation))
        return ETD_STR;
    /// HCD (beam-type collision-induced dissociation): A collision-induced dissociation process that occurs in a beam-type collision cell.
    else if (a.hasCVParam(pwiz::msdata::MS_HCD)) //MS_high_energy_collision_induced_dissociation))
        return HCD_STR;
    else
        return UNKNOWN_STR;
    
}


inline string modeToString(DataMode m) {
    if(m == PROFILE)  return "PROFILE";
    if(m == CENTROID) return "CENTROID";
    if(m == FITTED)   return "FITTED";
    return "<unknown datamode>";
}

inline string toLower(std::string input) {
    std::locale locale;
    std::string output = "";
    for (std::string::size_type i = 0; i < input.length(); i++)
        output += std::tolower(input[i], locale);
    return output;
}


/**
 * Some useful functions to get first approximation of centroid
 */
namespace mzMath {


#define SIGMA_FACTOR 2.354820045

/// getting ``ymax`` with in a gaussian model
inline static double yMax( double xZero, double sigmaSquared, double x, double y) {
    double nx = x - xZero;
    return y / exp( - ( nx*nx ) /( 2 * sigmaSquared) );
}

/// compute y value given gaussian model parameters
inline static double y( double x, double xZero, double yMax, double sigmaSquared ) {
    double nx = x - xZero;
    return yMax * exp( - ( nx * nx ) / ( 2 * sigmaSquared) );
}

///compute sigma given gaussian parameters
inline static double sigma( double width, double relativeHeight ) {
    /* actually if width is null there is no problem */
    if ( width < 0 ) {
        throw std::exception("Width must >= 0");
    }
    if ( relativeHeight < 0 || relativeHeight > 1 ) {
        throw std::exception("relativeHeight must be between 0 and 1");
    }

    return width /  sqrt(- 2 * log(relativeHeight) ) ;
}

///compute width of peak in Da from gaussian parameters
inline static double width( double sigma, double relativeHeight) {
    if (! (sigma > 0) )
        throw std::exception("sigma must be > 0");
    if (! (relativeHeight <= 1) )
        throw std::exception("relative must  < 1");
    return sigma * ( 2 * sqrt(- 2 * log(relativeHeight) ) );
}

///compute the FWHM (Full Width at HALF Maximum)
inline static double fwhm( double sigma) {
    return sigma * SIGMA_FACTOR;
}

/**
 * the wizard math formula: this is used to estimate the sampling
 * when reconstructing MS signal from fitted data.
 *
 * This formula has been retrieved empirically from Thermo data.
 *
 */
inline static double mzToDeltaMz(double mz) {
    return 2E-07 * pow(mz, 1.5);
}


template< typename T>
inline bool static isFiniteNumber(T& x) {
   return (x <= DBL_MAX && x >= -DBL_MAX);
}



/**
 * MaxQuant formula (parabola) to compute m/z centroid from mass peak data points
 *
 * @param xData m/z values
 * @param yData intensities values
 * @return m/z of the centroid
 */
template<class mz_t, class int_t>
static double gaussianCentroidApex(const std::vector<mz_t>& xData, const std::vector<int_t>& yData) {
    int nb_values = xData.size();
    if (nb_values != 3)
        throw std::exception("[gaussian centroid apex] failed");

    double x_m1 = xData[0];
    double y_m1 = yData[0];
    double x_0 = xData[1];
    double y_0 = yData[1];
    double x_p1 = xData[2];
    double y_p1 = yData[2];

    double diff_log_y_0_p1 = log(y_0) - log(y_p1);
    double diff_log_y_p1_m1 = log(y_p1) - log(y_m1);
    double diff_log_y_m1_0 = log(y_m1) - log(y_0);

    double div = diff_log_y_0_p1 * x_m1 + diff_log_y_p1_m1 * x_0 + diff_log_y_m1_0 * x_p1;
    if(div != 0) {
        return 0.5 * (diff_log_y_0_p1 * pow(x_m1, 2.0) + diff_log_y_p1_m1 * pow(x_0, 2.0) + diff_log_y_m1_0 * pow(x_p1, 2.0)) / div;
    } else {
        // return average mz to avoid dividing by zero
        // this case should not happen, but it's safer !
        LOG(WARNING) << "Function gaussianCentroidApex returned an average value to avoid a division by zero !";
        return (x_m1 + x_0 + x_p1) / 3;
    }
}

/// Extension of the MaxQuant formula (parabola) to compute intensity maximum from mass peak data points
/// (Does not perform very well)
template<typename mz_t, typename int_t>
static double gaussianCentroidIntensityMax(const std::vector<mz_t>& xData, const std::vector<int_t>& yData) {
    double x1 = xData[0];
    double y1 = yData[0];
    double x2 = xData[1];
    double y2 = yData[1];
    double x3 = xData[2];
    double y3 = yData[2];

    double D = (x1 - x2) * (x1 - x3) * (x2 - x3);
    double alpha = x3 * (log(y2) - log(y1)) + x2 * (log(y1) - log(y3)) + x1 * (log(y3) - log(y2));
    double gamma = log(y1)* x2 * x3 * (x2 - x3) + log(y2)*x3*x1*(x3-x1) + log(y3) * x1 * x2 * (x1- x2);
    double beta = pow(x3, 2) * (log(y1) - log(y2)) + pow(x2, 2) * (log(y3) - log(y1)) + pow(x1, 2) * (log(y2) - log(y3));
    return exp( (gamma/D) - pow( beta/D , 2) / (4 * alpha/D) );
}

/// Extension of the MaxQuant formula (parabola) to compute sigma maximum from mass peak data points
/// (not rigorous since gaussian is a parabola only around the apex)
template<typename mz_t, typename int_t>
static double gaussianSigma(const std::vector<mz_t>& xData, const std::vector<int_t>& yData) {
    double x1 = xData[0];
    double y1 = yData[0];
    double x2 = xData[1];
    double y2 = yData[1];
    double x3 = xData[2];
    double y3 = yData[2];

    double D = (x1 - x2) * (x1 - x3) * (x2 - x3);
    double alpha = x3 * (log(y2) - log(y1)) + x2 * (log(y1) - log(y3)) + x1 * (log(y3) - log(y2));
    return sqrt( - 1. / ( 2 * (alpha/D)) );
}

}//end mzmath


///utility function to print a nice progression bar
inline static void printProgBar(int percent) {
    string bar;

    for (int i = 0; i < 50; i++) {
        if (i < (percent / 2)) {
            bar.replace(i, 1, "=");
        } else if (i == (percent / 2)) {
            bar.replace(i, 1, ">");
        } else {
            bar.replace(i, 1, " ");
        }
    }
    std::cout << "\r" "[" << bar << "] ";
    std::cout << percent << "%     " << std::flush;
    if(percent >= 100)
        std::cout << "\n";
}

//inline static string getDateAsYyyyMmDd() {
//		time_t t;
//		struct tm* tm;
//		char fdate[9];
//		time(&t);
//		tm = localtime(&t);
//		strftime(fdate, sizeof fdate, "%Y%m%d", tm);
//		return fdate;
//}

inline static void exitOnError(std::string message) {
    // print the message in such a way that user cannot think it's ok
    // the message should be really informative, for the user as well as the IT support
    // it seems complicated to get the caller function name in a portable way, so for now just set a really informative error message !
    printf("\nMZDB Failure !!!\n");
    printf("Error: %s\n\n", message.c_str());
    // delete the mzDB file if possible
    
    // exit without a runtime_error
    exit(EXIT_FAILURE);
}

}//end namespace

//namespace genericUtils {
//
//    static SIZE_T getMemoryUsage() {
//        // function based on code available here: http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
//        PROCESS_MEMORY_COUNTERS_EX pmc;
//        GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
//        SIZE_T virtualMemUsedByMe = pmc.PrivateUsage;
//        //if(virtualMemUsedByMe > 500000000) LOG(WARNING) << "MEMORY IS TOO HIGH (" << virtualMemUsedByMe/1048576 << "MB)";
//        return virtualMemUsedByMe;
//    }
//    
//    static void checkMemoryUsage(std::string message, size_t limit, size_t i) {
//        size_t memory = getMemoryUsage();
//        if(memory >= limit) {
//            //LOG(WARNING) << "Current memory is " << memory/1048576 << "MB (" << message << ")";
//            //printf("Current memory is %d MB (%s)\n", memory/1048576, message.c_str());
//            printf("%s [item %d] (%d Mo)\n", message.c_str(), i, memory/1048576);
//            if(memory >= 1200000000) {
//                //LOG(ERROR) << "Current memory is too high, exiting !";
//                printf("Current memory is too high, exiting !\n");
//                exit(EXIT_FAILURE);
//            }
//        }
//    }
//    
//    static void checkMemoryUsage(std::string message, size_t i) {
//        checkMemoryUsage(message, 500000000, i);
//    }
//    
//    static void checkMemoryUsage(std::string message) {
//        checkMemoryUsage(message, 0);
//    }
//    
//    static void checkMemoryUsage() {
//        checkMemoryUsage("");
//    }
//}

#endif // MZUTILS_HPP
