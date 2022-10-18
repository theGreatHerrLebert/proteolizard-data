#ifndef CPP_TIMS_SLICE_H
#define CPP_TIMS_SLICE_H

#include "Frame.h"
#include "VectorizedFrame.h"
#include "SliceVectorized.h"

class Points3D {
public:
    Points3D(std::vector<int>& frames, std::vector<int>& scans,
        std::vector<double>& mzs, std::vector<double>& invIonMobs, std::vector<int>& intensities);

    std::vector<int> frame;
    std::vector<int> scan;
    std::vector<int> intensity;
    std::vector<double> invIonMobility;
    std::vector<double> mz;
};

class TimsSlicePL {
public:
    TimsSlicePL(std::vector<TimsFramePL>& pf, std::vector<TimsFramePL>& ff);
    TimsSlicePL filterRanged(int scanMin, int scanMax, double mzMin, double mzMax, int intensityMin, double rtMin, double rtMax);
    Points3D getPoints3D(bool precursor);
    TimsSliceVectorizedPL getVectorizedSlice(int resolution);

    std::vector<TimsFramePL> precursors;
    std::vector<TimsFramePL> fragments;
};

#endif //CPP_TIMS_SLICE_H
