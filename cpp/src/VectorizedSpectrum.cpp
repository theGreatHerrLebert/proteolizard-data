#include "VectorizedSpectrum.h"

#include <Eigen/Sparse>

MzVectorPL::MzVectorPL(int res, int frame, int scan, std::vector<int> index, std::vector<int> value):
    resolution(res), frameId(frame), scanId(scan), indices(std::move(index)), values(std::move(value)) {}

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
