//
// Created by administrator on 10.08.22.
//

#include "TimsBlockVectorized.h"

TimsBlockVectorizedPL::TimsBlockVectorizedPL(int frameStart, int frameStop,
                                             std::vector<int> frame,
                                             std::vector<int> scan,
                                             std::vector<int> mz,
                                             std::vector<int> intensity):
                                             frameIdxStart(frameStart),
                                             frameIdxStop(frameStop),
                                             frameIndices(frame),
                                             scanIndices(scan),
                                             mzIndices(mz),
                                             intensityValues(intensity) {}

std::vector<std::vector<int>> TimsBlockVectorizedPL::getData() {
    std::vector<std::vector<int>> data;
    data.push_back(frameIndices);
    data.push_back(scanIndices);
    data.push_back(mzIndices);
    data.push_back(intensityValues);
    return data;
}

TimsBlockVectorizedPL TimsBlockVectorizedPL::filterRanged(int scanMin, int scanMax, double mzMin, double mzMax,
                                                          int intensityMin) {
    std::vector<int> filteredFrame;
    std::vector<int> filteredScan;
    std::vector<int> filteredMz;
    std::vector<int> filteredIntensity;

    for (int i = 0; i < this->frameIndices.size(); i++) {
        if (this->scanIndices[i] >= scanMin && this->scanIndices[i] <= scanMax &&
            this->mzIndices[i] >= mzMin && this->mzIndices[i] <= mzMax &&
            this->intensityValues[i] >= intensityMin) {
            filteredFrame.push_back(this->frameIndices[i]);
            filteredScan.push_back(this->scanIndices[i]);
            filteredMz.push_back(this->mzIndices[i]);
            filteredIntensity.push_back(this->intensityValues[i]);
        }
    }

    return TimsBlockVectorizedPL(this->frameIdxStart, this->frameIdxStop, filteredFrame, filteredScan, filteredMz,
                                  filteredIntensity);
}