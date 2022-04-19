#ifndef CPP_VECTORIZED_FRAME_H
#define CPP_VECTORIZED_FRAME_H

#include "Spectrum.h"
#include <map>
#include <vector>

/**
 * a vectorized timsTOF frame
 */
class TimsFrameVectorizedPL {
public:
    // constructors
    TimsFrameVectorizedPL()= default;
    TimsFrameVectorizedPL(int id, int res, std::vector<int>  scan, std::vector<int>  index, std::vector<int>  v);

    // functions
    [[nodiscard]] TimsFrameVectorizedPL filterRanged(int scanMin, int scanMax, int indexMin, int indexMax, int intensityMin) const;
    std::map<int, MzVectorPL> spectra();

    friend TimsFrameVectorizedPL operator+(const TimsFrameVectorizedPL &leftFrame, const TimsFrameVectorizedPL &rightFrame);

    int frameId{};
    int resolution{};
    std::vector<int> scans;
    std::vector<int> indices;
    std::vector<int> values;
};

TimsFrameVectorizedPL operator+(const TimsFrameVectorizedPL &leftFrame, const TimsFrameVectorizedPL &rightFrame);

#endif //CPP_VECTORIZED_FRAME_H
