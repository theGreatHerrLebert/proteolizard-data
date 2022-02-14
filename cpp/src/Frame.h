//
// Created by david on 20.07.21.
//

#ifndef CPP_FRAME_H
#define CPP_FRAME_H

#include <tuple>
#include <string>
#include <math.h>
#include <map>
#include <algorithm>
#include "Spectrum.h"
#include "VectorizedFrame.h"
#include "../../include/eigen/Eigen/Dense"
#include "../../include/eigen/Eigen/Sparse"

/**
 * same container but trying to provide cleaner OOP interface
 */
struct TimsFramePL {
    // member data
    int frameId; 
    std::vector<int> scans;
    std::vector<double> mzs;
    std::vector<int> intensities, tofs;
    std::vector<double> inv_ion_mobs;
    
    // constructors
    TimsFramePL(){}
    TimsFramePL(int id, const std::vector<int>& scan, const std::vector<double>& mz, const std::vector<int>& intensity, const std::vector<int>& tof, const std::vector<double>& inv_ion_mob):
    frameId(id), scans(scan), mzs(mz), intensities(intensity), tofs(tof), inv_ion_mobs(inv_ion_mob){}

    // binning to a finite resolution
    TimsFramePL filterRanged(int scanMin, int scanMax, double mzMin, double mzMax, int minIntensity=1);
    TimsFramePL toResolution(const int resolution);
     friend TimsFramePL operator+(const TimsFramePL &leftFrame, const TimsFramePL &rightFrame);
    TimsFramePL fold(const int resolution, const int foldWidth);

    // vectorize to an integer mz value
    TimsFrameVectorizedPL vectorize(const int resolution);

    // get all spectra as map from scan to spectrum
    std::map<int, MzSpectrum> spectra();
};

// note: this function is not a member function!
TimsFramePL operator+(const TimsFramePL &leftFrame, const TimsFramePL &rightFrame){

    std::map<std::pair<int, double>, int> sumMap;

    // insert leftFrame values into map
    for (auto it = leftFrame.mzs.begin(); it != leftFrame.mzs.end(); ++it) {
        int i = std::distance(leftFrame.mzs.begin(), it);
        auto scan = leftFrame.scans[i];
        auto index = leftFrame.mzs[i];
        auto intensity = leftFrame.intensities[i];
        sumMap[{scan, index}] = intensity;
    }

    // insert right frame values into map or sum
    for (auto it = rightFrame.mzs.begin(); it != rightFrame.mzs.end(); ++it) {
        int i = std::distance(rightFrame.mzs.begin(), it);
        auto scan = rightFrame.scans[i];
        auto index = rightFrame.mzs[i];
        auto intensity = rightFrame.intensities[i];

        // if key is already in map, sum intensity
        if (sumMap.contains({scan, index}))
            sumMap[{scan, index}] += intensity;
        // create new key value pair otherwise
        else
            sumMap[{scan, index}] = intensity;
    }

    std::vector<int> retScans;
    retScans.reserve(sumMap.size());
    std::vector<double> retIndices;
    retIndices.reserve(sumMap.size());
    std::vector<int> retValues;
    retValues.reserve(sumMap.size());

   
   for(auto const& [key, value]: sumMap){
        retScans.push_back(key.first);
        retIndices.push_back(key.second);
        retValues.push_back(value);   
   }

    return TimsFramePL(leftFrame.frameId, retScans, retIndices, retValues, {}, {});
}

/**
 *
 * @param resolution
 * @return
 */
TimsFramePL TimsFramePL::toResolution(const int resolution) {

    std::map<std::pair<int, int>, float> intensityMap;
    double factor = pow(10.0, resolution);

    std::vector<double> resMz;
    std::vector<int> resS, resI;

    for(std::size_t i = 0; i < this->mzs.size(); i++){

        // calculate binned mz value as key_first
        int rounded = int(round(this->mzs[i] * factor));
        // get scan coordinate as key_second
        int s = this->scans[i];
        // add intensity at combined key location
        intensityMap[{s, rounded}] += this->intensities[i];
    }

    resMz.reserve(intensityMap.size());
    resI.reserve(intensityMap.size());
    resS.reserve(intensityMap.size());

    // vectorize binned spectrum
    for (const auto& [key, value] : intensityMap) {
        resMz.push_back(double(key.second) / factor);
        resS.push_back(key.first);
        resI.push_back(value);
    }

    return TimsFramePL(this->frameId, resS, resMz, resI, this->tofs, this->inv_ion_mobs);
}

/**
 *
 * @param scanMin
 * @param scanMax
 * @param mzMin
 * @param mzMax
 * @return
 */
TimsFramePL TimsFramePL::filterRanged(int scanMin, int scanMax, double mzMin, double mzMax, int minIntensity) {

    std::vector<int> retScan, retIntensity;
    std::vector<double> retMz;

    // TODO: make this more efficient (binary search?)
    for(std::size_t i = 0; i < this->mzs.size(); i++){

        int scan = this->scans[i];
        double mz = this->mzs[i];
        int intensity = this->intensities[i];

        if((scan >= scanMin) && (scan <= scanMax) && (mz >= mzMin) && (mz <= mzMax) && (intensity >= minIntensity)){
            retScan.push_back(scan);
            retMz.push_back(mz);
            retIntensity.push_back(intensity);
        }
    }

    // This guards for empty return
    if(retScan.size() > 0)
        return TimsFramePL(this->frameId, retScan, retMz, retIntensity, {1}, {0.0});

    return TimsFramePL(this->frameId, {(scanMin + scanMax) / 2}, {(mzMin + mzMax) / 2}, {0}, {1}, {1.0});
}

/**
 *
 * @param resolution
 * @param width
 * @return
 */
TimsFramePL TimsFramePL::fold(const int resolution, const int width) {
    // bin frame to same resolution along mz
    TimsFramePL frame = this->toResolution(resolution);

    // declare and reserve binned scan vector
    std::vector<int> scanFoldIds;
    scanFoldIds.reserve(frame.scans.size());

    // calculate and save scan -> scan_bin
    for(const auto& scan_index: frame.scans)
        scanFoldIds.push_back(int(floor(scan_index / width)));

    // declare summed intensity map
    std::map<std::tuple<int, double>, int> sumMap;

    // sum intensities based on index
    for(std::size_t i=0; i<frame.mzs.size(); i++)
        sumMap[{scanFoldIds[i], frame.mzs[i]}] += frame.intensities[i];

    //write back to frame
    std::vector<int> scans, intensities;
    std::vector<double> mzs;
    for(const auto& [key, intensity]: sumMap){
        const auto& [scan, mz] = key;
        scans.push_back(scan);
        mzs.push_back(mz);
        intensities.push_back(intensity);
    }

    return TimsFramePL(frameId, scans, mzs, intensities, this->tofs, this->inv_ion_mobs);
}

/**
 *
 * @param resolution
 * @return
 */
TimsFrameVectorizedPL TimsFramePL::vectorize(const int resolution) {
    
    std::map<std::pair<int, int>, int> sumMap;
    int factor = pow(10, resolution);

    std::vector<int> retScan, retIndex;
    std::vector<int> retValue;

    for (auto it = this->mzs.begin(); it != this->mzs.end(); ++it) {
        
        int i = std::distance(this->mzs.begin(), it);
        auto scan = this->scans[i];
        auto index = int(floor(round(this->mzs[i] * factor)));
        auto intensity = this->intensities[i];

        if(sumMap.contains({scan, index}))
            sumMap[{scan, index}] += intensity;
        else
            sumMap[{scan, index}] = intensity;
    }

    retScan.reserve(sumMap.size());
    retIndex.reserve(sumMap.size());
    retValue.reserve(sumMap.size());

    // vectorize binned spectrum
    for (const auto& [key, value] : sumMap) {
        retIndex.push_back(key.second);
        retScan.push_back(key.first);
        retValue.push_back(value);
    }

    return TimsFrameVectorizedPL(this->frameId, resolution, retScan, retIndex, retValue);
}

/**
 *
 * @return
 */
std::map<int, MzSpectrum> TimsFramePL::spectra() {

    std::map<int, MzSpectrum> specMap;

    for(size_t peak_id = 0; peak_id < this->scans.size(); peak_id++){
        auto scan = this->scans[peak_id];
        auto mz = this->mzs[peak_id];
        auto i = this->intensities[peak_id];
        if(specMap.contains(scan)){
            specMap[scan].mz.push_back(mz);
            specMap[scan].intensity.push_back(i);
        }
        else {
            specMap[scan] = MzSpectrum{this->frameId, static_cast<int>(scan), {mz}, {i}};
        }
    }

    return specMap;
}

#endif //CPP_FRAME_H
