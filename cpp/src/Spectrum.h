//
// Created by david on 21.07.21.
//

#ifndef CPP_SPECTRUM_H
#define CPP_SPECTRUM_H

#include <math.h>
#include "../../include/eigen/Eigen/Dense"
#include "../../include/eigen/Eigen/Sparse"
#include "VectorizedSpectrum.h"

/**
 * container for single timsTOF mz spectra
 */
struct MzSpectrum {
    // data members
    int frameId, scanId; // rt coordinate
    std::vector<double> mz; // vector of mz values
    std::vector<int> intensity; // vector of intensities
    
    // constructors
    MzSpectrum(){}
    MzSpectrum(int frame, int scan, std::vector<double> m, std::vector<int> i): frameId(frame), scanId(scan), mz(m), intensity(i){}
    
    // 
    MzSpectrum toResolution(int resolution);
    
    // 
    MzVector vectorize(int resolution);

    std::map<int, MzSpectrum> windows(double windowLength, bool overlapping, int minPeaks, int minIntensity);
    friend MzSpectrum operator+(const MzSpectrum &leftSpec, const MzSpectrum &rightSpec);
};

MzSpectrum operator+(const MzSpectrum &leftSpec, MzSpectrum &rightSpec){

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
   
   return MzSpectrum(leftSpec.frameId, leftSpec.scanId, retIndices, retValues);
}

/**
 *
 * @param resolution
 * @param normalize
 * @param sqr
 * @return
 */
MzVector MzSpectrum::vectorize(int resolution) {

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

    return MzVector{resolution, tmp.frameId, tmp.scanId, indices, values};
}

/**
 *
 * @param resolution
 * @return
 */
MzSpectrum MzSpectrum::toResolution(int resolution) {
    
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

    return MzSpectrum{this->frameId, this->scanId, resMz, resI};
}

#endif //CPP_SPECTRUM_H
