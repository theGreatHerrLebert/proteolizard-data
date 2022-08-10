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