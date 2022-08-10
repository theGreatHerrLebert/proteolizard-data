//
// Created by administrator on 10.08.22.
//

#ifndef CPP_TIMSBLOCKVECTORIZED_H
#define CPP_TIMSBLOCKVECTORIZED_H

#include <vector>

class TimsBlockVectorizedPL {
public:
    TimsBlockVectorizedPL(
            int frameStart,
            int frameStop,
            std::vector<int> frame,
            std::vector<int> scan,
            std::vector<int> mz,
            std::vector<int> intensity
    );

    int frameIdxStart;
    int frameIdxStop;
    std::vector<int> frameIndices;
    std::vector<int> scanIndices;
    std::vector<int> mzIndices;
    std::vector<int> intensityValues;
};

#endif //CPP_TIMSBLOCKVECTORIZED_H
