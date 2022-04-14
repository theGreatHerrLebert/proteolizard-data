//
// Created by David Teschner on 01.02.22.
//

#ifndef CPP_VECTORIZED_SPECTRUM_H
#define CPP_VECTORIZED_SPECTRUM_H

#include <cmath>

#include <utility>
#include "Eigen/Dense"
#include "Eigen/Sparse"

/**
 * container for single vectorized timsTOF mz spectra
 */
struct MzVectorPL {
    int resolution{}, frameId{}, scanId{}; // number of decimals spectrum was binned to
    std::vector<int> indices, values; // mz indices, can be seen as vector entries
    // constructors
    MzVectorPL()= default;
    MzVectorPL(int res, int frame, int scan, std::vector<int> index, std::vector<int> value):
    resolution(res), frameId(frame), scanId(scan), indices(std::move(index)), values(std::move(value)) {}

    friend MzVectorPL operator+(const MzVectorPL &leftSpec, const MzVectorPL &rightSpec);
};

// note: this function is not a member function!
MzVectorPL operator+(const MzVectorPL &leftSpec, const MzVectorPL &rightSpec){
    
    // frames with differing resolution might not be added for now.
    if (leftSpec.resolution != rightSpec.resolution)
        return leftSpec;

    std::map<int, int> sumMap;

    // insert leftFrame values into map
    for (auto it = leftSpec.indices.begin(); it != leftSpec.indices.end(); ++it) {
        auto i = std::distance(leftSpec.indices.begin(), it);
        auto index = leftSpec.indices[i];
        auto intensity = leftSpec.values[i];
        sumMap[index] = intensity;
    }

    // insert right frame values into map or sum
    for (auto it = rightSpec.indices.begin(); it != rightSpec.indices.end(); ++it) {
        auto i = std::distance(rightSpec.indices.begin(), it);
        auto index = rightSpec.indices[i];
        auto intensity = rightSpec.values[i];

        // if key is already in map, sum intensity
        if (sumMap.contains(index))
            sumMap[index] += intensity;
        // create new key value pair otherwise
        else
            sumMap[index] = intensity;
    }

    std::vector<int> retIndices;
    retIndices.reserve(sumMap.size());
    std::vector<int> retValues;
    retValues.reserve(sumMap.size());

   
   for(auto const& [key, value]: sumMap){
        retIndices.push_back(key);
        retValues.push_back(value);   
   }

    return {leftSpec.frameId, leftSpec.scanId, leftSpec.resolution, retIndices, retValues};
}

#endif //CPP_VECTORIZED_SPECTRUM_H
