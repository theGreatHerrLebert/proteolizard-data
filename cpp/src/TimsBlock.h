#ifndef CPP_TIMS_BLOCK_H
#define CPP_TIMS_BLOCK_H

#include "Frame.h"
#include "TimsBlockVectorized.h"

class TimsBlockPL {
public:
    TimsBlockPL(std::vector<int> frames, std::vector<int> scans, std::vector<int> tofs, std::vector<int> intensities,
                std::vector<double> retentionTimes, std::vector<double> invIonMobs, std::vector<double> mzs);

    std::vector<std::vector<int>> getIndices();
    std::vector<std::vector<double>> getValues();

    TimsBlockPL filterRanged(int scanMin, int scanMax, double mzMin, double mzMax, int intensityMin);
    TimsBlockVectorizedPL getBlockVectorized(int resolution);

    std::vector<int> frame;
    std::vector<int> scan;
    std::vector<int> tof;
    std::vector<int> intensity;
    std::vector<double> retentionTime;
    std::vector<double> invIonMobility;
    std::vector<double> mz;
};

#endif //CPP_TIMS_BLOCK_H
