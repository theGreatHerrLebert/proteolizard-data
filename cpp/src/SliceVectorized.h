//
// Created by administrator on 12.09.22.
//

#ifndef CPP_SLICEVECTORIZED_H
#define CPP_SLICEVECTORIZED_H

#include "Frame.h"
#include "VectorizedFrame.h"

class Points3DVectorized {
public:
    Points3DVectorized(std::vector<int>& frames, std::vector<int>& scans,
             std::vector<int>& mzs, std::vector<int>& intensities);

    std::vector<int> frame;
    std::vector<int> scan;
    std::vector<int> mz;
    std::vector<int> intensity;
};

class TimsSliceVectorizedPL {
public:
    TimsSliceVectorizedPL(std::vector<TimsFrameVectorizedPL>& pf, std::vector<TimsFrameVectorizedPL>& ff);
    TimsSliceVectorizedPL filterRanged(int scanMin, int scanMax, int mzMin, int mzMax, int intensityMin);
    Points3DVectorized getPoints3D(bool precursor);

    std::vector<TimsFrameVectorizedPL> precursors;
    std::vector<TimsFrameVectorizedPL> fragments;
};

#endif //CPP_SLICEVECTORIZED_H
