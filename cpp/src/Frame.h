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
#include <utility>
#include "Spectrum.h"
#include "VectorizedFrame.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"

#include <algorithm>
#include <execution>

/**
 * create eigen sparse vector from vectorized mz spectrum
 * @param mzVector : vectorized mz spectrum to convert
 * @param numRows : dimensionality of vector
 * @return : a sparse eigen vector suited for fast vectorized operations
 */
Eigen::MatrixXd toDenseEigen(const MzVectorPL& mzVector, int numRows){

    Eigen::SparseMatrix<double> sparseVec = Eigen::SparseMatrix<double>(numRows, 1);
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(mzVector.indices.size());

    for(std::size_t i = 0; i < mzVector.indices.size(); i++)
        tripletList.emplace_back(mzVector.indices[i], 0, mzVector.values[i]);

    sparseVec.setFromTriplets(tripletList.begin(), tripletList.end());
    return sparseVec.toDense();
}

/**
 * same container but trying to provide cleaner OOP interface
 */
struct TimsFramePL {
    // member data
    int frameId{};
    std::vector<int> scans {};
    std::vector<double> mzs {};
    std::vector<int> intensities, tofs {};
    std::vector<double> inv_ion_mobs {};
    
    // constructors
    TimsFramePL()= default;
    TimsFramePL(int id, std::vector<int>  scan, std::vector<double>  mz, std::vector<int>  intensity, std::vector<int>  tof, std::vector<double>  inv_ion_mob):
    frameId(id), scans(std::move(scan)), mzs(std::move(mz)), intensities(std::move(intensity)), tofs(std::move(tof)), inv_ion_mobs(std::move(inv_ion_mob)){}

    // binning to a finite resolution
    TimsFramePL filterRanged(int scanMin, int scanMax, double mzMin, double mzMax, int minIntensity=1);
    TimsFramePL toResolution(int resolution);
     friend TimsFramePL operator+(const TimsFramePL &leftFrame, const TimsFramePL &rightFrame);
    TimsFramePL fold(int resolution, int foldWidth);

    // vectorize to an integer mz value
    TimsFrameVectorizedPL vectorize(int resolution);

    // get all spectra as map from scan to spectrum
    std::map<int, MzSpectrumPL> spectra();

    std::vector<MzSpectrumPL> exportSpectra();

    std::pair<std::vector<int>, std::pair<std::pair<std::vector<int>, std::vector<int>>,std::pair<std::vector<int>,
    std::vector<int>>>> getHashingBlocks(int resolution, int minPeaksPerWindow, int minIntensity,
                                        double windowLength, bool overlapping);

    Eigen::MatrixXd denseWindowMatrix(int resolution, int minPeaksPerWindow,
                           int minIntensity, double windowLength, bool overlapping);
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

    return {this->frameId, resS, resMz, resI, this->tofs, this->inv_ion_mobs};
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
    std::vector<double> retMz, retInv;

    // TODO: make this more efficient (binary search?)
    for(std::size_t i = 0; i < this->mzs.size(); i++){

        int scan = this->scans[i];
        double mz = this->mzs[i];
        int intensity = this->intensities[i];
        double invMob = this->inv_ion_mobs[i];

        if((scan >= scanMin) && (scan <= scanMax) && (mz >= mzMin) && (mz <= mzMax) && (intensity >= minIntensity)){
            retScan.push_back(scan);
            retMz.push_back(mz);
            retIntensity.push_back(intensity);
            retInv.push_back(invMob);
        }
    }

    // This guards for empty return
    if(!retScan.empty())
        return TimsFramePL(this->frameId, retScan, retMz, retIntensity, {1}, {retInv});

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

    return {frameId, scans, mzs, intensities, this->tofs, this->inv_ion_mobs};
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

    return {this->frameId, resolution, retScan, retIndex, retValue};
}

/**
 *
 * @return
 */
std::map<int, MzSpectrumPL> TimsFramePL::spectra() {

    std::map<int, MzSpectrumPL> specMap;

    for(size_t peak_id = 0; peak_id < this->scans.size(); peak_id++){
        auto scan = this->scans[peak_id];
        auto mz = this->mzs[peak_id];
        auto i = this->intensities[peak_id];
        if(specMap.contains(scan)){
            specMap[scan].mz.push_back(mz);
            specMap[scan].intensity.push_back(i);
        }
        else {
            specMap[scan] = MzSpectrumPL{this->frameId, static_cast<int>(scan), {mz}, {i}};
        }
    }

    return specMap;
}

std::vector<MzSpectrumPL> TimsFramePL::exportSpectra() {

    auto spectra = this->spectra();
    std::vector<MzSpectrumPL> ret;
    ret.reserve(spectra.size());

    for(auto&[key, value]: spectra)
        ret.push_back(value);

    return ret;
}

Eigen::MatrixXd TimsFramePL::denseWindowMatrix(int resolution, int minPeaksPerWindow, int minIntensity,
                                               double windowLength, bool overlapping) {

    // split frame into spectra
    auto spectra = this->spectra();
    auto numRows = int(round(2000 * pow(10, resolution)));
    std::vector<Eigen::MatrixXd> vec;

    for(auto &spectrum: spectra){
        auto windows = spectrum.second.windows(windowLength, overlapping, minPeaksPerWindow, minIntensity);
        for(auto &w: windows){
            auto v = toDenseEigen(w.second.vectorize(resolution), numRows);
            vec.push_back(v);
        }
    }

    auto numColumns = vec.size();

    Eigen::MatrixXd retMatrix(numRows, numColumns);

    for(std::size_t i = 0; i < vec.size(); i++){
        retMatrix.col(i) = vec[i];
    }

    return retMatrix;
}

std::pair<std::vector<int>, std::pair<std::pair<std::vector<int>, std::vector<int>>, std::pair<std::vector<int>,
        std::vector<int>>>> TimsFramePL::getHashingBlocks(int resolution,
                                                         int minPeaksPerWindow,
                                                         int minIntensity,
                                                         double windowLength,
                                                         bool overlapping) {

    auto spectra = this->spectra();
    std::vector<int> retScan, retBin, retWindowIndex, retIndices, retValues;

    retScan.reserve(70000);
    retBin.reserve(70000);
    retWindowIndex.reserve(70000);
    retIndices.reserve(70000);
    retValues.reserve(70000);

    int windowCounter = 0;

    for(const auto &[scan, spectrum]: spectra){

        auto windows = spectrum.windows(windowLength, overlapping,
                                        minPeaksPerWindow, minIntensity);

        for(const auto& [bin, window]: windows){

            int offset = 0;

            if (bin > 0) {
                offset = pow(10, resolution) * bin * windowLength;
            }
            else {
                offset = pow(10, resolution) * ((-bin) * windowLength - (windowLength / 2));
            }

            auto vectorizedWindow = window.vectorize(resolution);

            for(std::size_t i = 0; i < vectorizedWindow.indices.size(); i++){
                auto index = vectorizedWindow.indices[i];
                auto zeroIndex = index - offset;
                auto value = vectorizedWindow.values[i];
                retIndices.push_back(zeroIndex);
                retValues.push_back(value);
            }

            retScan.insert(retScan.end(), vectorizedWindow.indices.size(), scan);
            retBin.insert(retBin.end(),vectorizedWindow.indices.size(), bin);
            retWindowIndex.insert(retWindowIndex.end(), vectorizedWindow.indices.size(), windowCounter);
            windowCounter++;
        }

    }


    return {retWindowIndex, {{retScan, retBin}, {retIndices, retValues}}};
}

#endif //CPP_FRAME_H
