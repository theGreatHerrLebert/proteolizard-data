//
// Created by administrator on 12.09.22.
//

#include "SliceVectorized.h"

#include <algorithm>
#include <execution>

Points3DVectorized::Points3DVectorized(
        std::vector<int>& frames,
        std::vector<int>& scans,
        std::vector<int>& mzs,
        std::vector<int>& intensities
): frame(frames), scan(scans),  mz(mzs), intensity(intensities) {}

TimsSliceVectorizedPL::TimsSliceVectorizedPL(std::vector<TimsFrameVectorizedPL>& pf, std::vector<TimsFrameVectorizedPL>& ff): precursors(pf), fragments(ff) {}

TimsSliceVectorizedPL TimsSliceVectorizedPL::filterRanged(int scanMin, int scanMax, int mzMin, int mzMax, int intensityMin){

    std::vector<TimsFrameVectorizedPL> retPrecursors;
    std::vector<TimsFrameVectorizedPL> retFragments;

    retPrecursors.resize(this->precursors.size());
    retFragments.resize(this->fragments.size());

    auto filterFrame = [&scanMin, &scanMax, &mzMin, &mzMax, &intensityMin](TimsFrameVectorizedPL& f) -> TimsFrameVectorizedPL {
        return f.filterRanged(scanMin, scanMax, mzMin, mzMax, intensityMin);
    };

    std::transform(std::execution::par_unseq, this->precursors.begin(), this->precursors.end(), retPrecursors.begin(), filterFrame);
    std::transform(std::execution::par_unseq, this->fragments.begin(), this->fragments.end(), retFragments.begin(), filterFrame);

    return {retPrecursors, retFragments};
}

Points3DVectorized TimsSliceVectorizedPL::getPoints3D(bool precursor){

    std::vector<int> retFrames, retScan, retIntensity, retMz;

    if(precursor)
    {
        for(auto& frame: this->precursors)
        {
            retFrames.insert(retFrames.end(), frame.indices.size(), frame.frameId);
            retScan.insert(retScan.end(), frame.scans.begin(), frame.scans.end());
            retMz.insert(retMz.end(), frame.indices.begin(), frame.indices.end());
            retIntensity.insert(retIntensity.end(), frame.values.begin(), frame.values.end());
        }
    }

    else
    {
        for(auto& frame: this->fragments)
        {
            retFrames.insert(retFrames.end(), frame.indices.size(), frame.frameId);
            retScan.insert(retScan.end(), frame.scans.begin(), frame.scans.end());
            retMz.insert(retMz.end(), frame.indices.begin(), frame.indices.end());
            retIntensity.insert(retIntensity.end(), frame.values.begin(), frame.values.end());
        }
    }

    return {retFrames, retScan, retMz, retIntensity};
}