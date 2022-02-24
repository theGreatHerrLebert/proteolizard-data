//
// Created by david on 21.07.21.
//

#ifndef CPP_SPECTRUM_H
#define CPP_SPECTRUM_H

#include <math.h>

#include <utility>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "VectorizedSpectrum.h"

/**
 * container for single timsTOF mz spectra
 */
struct MzSpectrumPL {
    // data members
    int frameId{}, scanId{}; // rt coordinate
    std::vector<double> mz; // vector of mz values
    std::vector<int> intensity; // vector of intensities
    
    // constructors
    MzSpectrumPL()= default;
    MzSpectrumPL(int frame, int scan, std::vector<double> m, std::vector<int> i): frameId(frame), scanId(scan), mz(std::move(m)), intensity(std::move(i)){}

    MzSpectrumPL toResolution(int resolution);
    MzVectorPL vectorize(int resolution);

    std::map<int, MzSpectrumPL> windows(double windowLength, bool overlapping, int minPeaks, int minIntensity);
    friend MzSpectrumPL operator+(const MzSpectrumPL &leftSpec, const MzSpectrumPL &rightSpec);
};

MzSpectrumPL operator+(const MzSpectrumPL &leftSpec, const MzSpectrumPL &rightSpec){

    std::map<double, int> sumMap;

    // insert leftFrame values into map
    for (auto it = leftSpec.mz.begin(); it != leftSpec.mz.end(); ++it) {
        int i = std::distance(leftSpec.mz.begin(), it);
        auto index = leftSpec.mz[i];
        auto intensity = leftSpec.intensity[i];
        sumMap[index] = intensity;
    }

    // insert right frame values into map or sum
    for (auto it = rightSpec.mz.begin(); it != rightSpec.mz.end(); ++it) {
        int i = std::distance(rightSpec.mz.begin(), it);
        auto index = rightSpec.mz[i];
        auto intensity = rightSpec.intensity[i];

        // if key is already in map, sum intensity
        if (sumMap.contains(index))
            sumMap[index] += intensity;
        // create new key value pair otherwise
        else
            sumMap[index] = intensity;
    }

    std::vector<double> retIndices;
    retIndices.reserve(sumMap.size());
    std::vector<int> retValues;
    retValues.reserve(sumMap.size());

   
   for(auto const& [key, value]: sumMap){
        retIndices.push_back(key);
        retValues.push_back(value);   
   }
   
   return {leftSpec.frameId, leftSpec.scanId, retIndices, retValues};
}

/**
 *
 * @param resolution
 * @param normalize
 * @param sqr
 * @return
 */
MzVectorPL MzSpectrumPL::vectorize(int resolution) {

    auto tmp = this->toResolution(resolution);

    int factor = pow(10, resolution);
    std::vector<int> indices;
    indices.reserve(tmp.mz.size());

    std::vector<int> values;
    values.reserve(tmp.intensity.size());

    for(const auto m : tmp.mz)
        indices.push_back(m * factor);

    for(const auto i : tmp.intensity)
        values.push_back(i);

    return MzVectorPL{resolution, tmp.frameId, tmp.scanId, indices, values};
}

/**
 *
 * @param resolution
 * @return
 */
MzSpectrumPL MzSpectrumPL::toResolution(int resolution) {
    
    std::map<int, int> intensityMap;
    double factor = pow(10.0, resolution);

    std::vector<double> resMz;
    std::vector<int> resI;

    for(size_t i = 0; i < this->mz.size(); i++){

        // calculate binned mz value as key
        int rounded = int(roundf(this->mz[i] * factor));
        // add intensity to generated key
        intensityMap[rounded] += this->intensity[i];
    }

    resMz.reserve(intensityMap.size());
    resI.reserve(intensityMap.size());

    // get all mz values and sort them
    for (const auto& [key, value] : intensityMap) {
        resMz.push_back(double(key) / factor);
        resI.push_back(value);
    }

    return MzSpectrumPL{this->frameId, this->scanId, resMz, resI};
}

/**
 * check if a given windowId is already present in map
 * @param windowCollection collection to check
 * @param windowId id to check for
 * @return true if key is in map else false
 */
bool hasValue(std::map<int, MzSpectrumPL>& windowCollection, int windowId){
    return windowCollection.count(windowId);
}

std::map<int, MzSpectrumPL> MzSpectrumPL::windows(double windowLength, bool overlapping, int minPeaks, int minIntensity)
{
    std::map<int, MzSpectrumPL> splits;

    for(int i = 0; i < this->mz.size(); i++) {
        // get mz and intensity
        auto mmz = this->mz[i];
        auto iintensity = this->intensity[i];
        // calculate window key mz and intensity value should be sent to
        auto tmpKey = int(floor(mmz / windowLength));

        // given key is in map
        if (hasValue(splits, tmpKey)) {
            splits[tmpKey].mz.push_back(mmz);
            splits[tmpKey].intensity.push_back(iintensity);
        }
            // given key is not in map
        else {
            std::vector<double> tmpMz = {mmz};
            std::vector<int> tmpI = {iintensity};
            splits[tmpKey] = MzSpectrumPL{this->frameId, this->scanId, tmpMz, tmpI};
        }
    }

    // calculate grouping by offset
    if(overlapping) {
        std::map<int, MzSpectrumPL> splitsOffset;

        for (int i = 0; i < this->mz.size(); i++) {
            auto mmz = this->mz[i];
            auto iintensity = this->intensity[i];
            // calculate window key with offset mz and intensity value should be sent to
            auto tmpKey = -int(floor(((mmz + windowLength / 2.0) / windowLength)));

            // given key is in map
            if (hasValue(splitsOffset, tmpKey)) {
                splitsOffset[tmpKey].mz.push_back(mmz);
                splitsOffset[tmpKey].intensity.push_back(iintensity);
            }
                // given key is not in map
            else {
                std::vector<double> tmpMz = {mmz};
                std::vector<int> tmpI = {iintensity};
                splitsOffset[tmpKey] = MzSpectrumPL{this->frameId, this->scanId, tmpMz, tmpI};
            }
        }
        splits.merge(splitsOffset);
    }

    std::map<int, MzSpectrumPL> retSplits;

    for(const auto& [bin, spectrum]: splits){
        if(this->mz.size() >= minPeaks){
            auto it = *max_element(std::begin(spectrum.intensity), std::end(spectrum.intensity));
            if(it >= minIntensity)
                retSplits[bin] = spectrum;
        }
    }

    return retSplits;
}

#endif //CPP_SPECTRUM_H
