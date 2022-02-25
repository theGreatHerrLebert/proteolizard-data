//
// Created by David Teschner on 01.02.22.
//

#ifndef CPP_VECTORIZED_FRAME_H
#define CPP_VECTORIZED_FRAME_H

#include <tuple>
#include <string>
#include <math.h>
#include <map>
#include <algorithm>
#include <utility>
#include "Spectrum.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"

/**
 * a vectorized timsTOF frame
 */
struct TimsFrameVectorizedPL {
    
    int frameId{}, resolution{}; // coordinates
    std::vector<int> scans, indices, values; // coordinates

    TimsFrameVectorizedPL()= default;;
    // constructors
    TimsFrameVectorizedPL(int id, int res, std::vector<int>  scan, std::vector<int>  index, std::vector<int>  v):
    frameId(id), resolution(res), scans(std::move(scan)), indices(std::move(index)), values(std::move(v)){}

    // functions
    [[nodiscard]] TimsFrameVectorizedPL filterRanged(int scanMin, int scanMax, int indexMin, int indexMax, int intensityMin) const;
    std::map<int, MzVectorPL> spectra();
    friend TimsFrameVectorizedPL operator+(const TimsFrameVectorizedPL &leftFrame, const TimsFrameVectorizedPL &rightFrame);
};

std::map<int, MzVectorPL> TimsFrameVectorizedPL::spectra() {

    std::map<int, MzVectorPL> specMap;
    
    // insert leftFrame values into map
    for (auto it = this->indices.begin(); it != this->indices.end(); ++it) {
        int i = std::distance(this->indices.begin(), it);
        auto scan = this->scans[i];
        auto index = this->indices[i];
        auto intensity = this->values[i];

        if(specMap.contains(scan)){
            specMap[scan].indices.push_back(index);
            specMap[scan].values.push_back(intensity);
        }
        else {
            specMap[scan] = MzVectorPL{this->resolution, this->frameId, scan, {index}, {intensity}};
        }
    }
    return specMap;
}

// note: this function is not a member function!
TimsFrameVectorizedPL operator+(const TimsFrameVectorizedPL &leftFrame, const TimsFrameVectorizedPL &rightFrame){
    
    // frames with differing resolution might not be added for now.
    if (leftFrame.resolution != rightFrame.resolution)
        return leftFrame;

    std::map<std::pair<int, int>, int> sumMap;

    // insert leftFrame values into map
    for (auto it = leftFrame.indices.begin(); it != leftFrame.indices.end(); ++it) {
        int i = std::distance(leftFrame.indices.begin(), it);
        auto scan = leftFrame.scans[i];
        auto index = leftFrame.indices[i];
        auto intensity = leftFrame.values[i];
        sumMap[{scan, index}] = intensity;
    }

    // insert right frame values into map or sum
    for (auto it = rightFrame.indices.begin(); it != rightFrame.indices.end(); ++it) {
        int i = std::distance(rightFrame.indices.begin(), it);
        auto scan = rightFrame.scans[i];
        auto index = rightFrame.indices[i];
        auto intensity = rightFrame.values[i];

        // if key is already in map, sum intensity
        if (sumMap.contains({scan, index}))
            sumMap[{scan, index}] += intensity;
        // create new key value pair otherwise
        else
            sumMap[{scan, index}] = intensity;
    }

    std::vector<int> retScans;
    retScans.reserve(sumMap.size());
    std::vector<int> retIndices;
    retIndices.reserve(sumMap.size());
    std::vector<int> retValues;
    retValues.reserve(sumMap.size());

   
   for(auto const& [key, value]: sumMap){
        retScans.push_back(key.first);
        retIndices.push_back(key.second);
        retValues.push_back(value);   
   }

    return {leftFrame.frameId, leftFrame.resolution, retScans, retIndices, retValues};
}

/**
 *
 * @param scanMin
 * @param scanMax
 * @param mzMin
 * @param mzMax
 * @return
 */
TimsFrameVectorizedPL TimsFrameVectorizedPL::filterRanged(int scanMin, int scanMax, int indexMin, int indexMax, int intensityMin) const {

    std::vector<int> retIndices, retScans;
    std::vector<int> retValues;

    for(std::size_t i = 0; i < this->values.size(); i++){

        int scan = this->scans[i];
        int index = this->indices[i];
        int value = this->values[i];

        if((scan >= scanMin) && (scan <= scanMax) && (index >= indexMin) && (index <= indexMax) && (value >= intensityMin)){
            retScans.push_back(scan);
            retIndices.push_back(index);
            retValues.push_back(value);
        }
    }

    // This guards for empty return
    if(retScans.size() > 0)
        return {this->frameId, this->resolution, retScans, retIndices, retValues};

    return {this->frameId, this->resolution, {(scanMin + scanMax) / 2}, {(indexMin + indexMax) / 2}, {0}};
}

#endif //CPP_VECTORIZED_FRAME_H