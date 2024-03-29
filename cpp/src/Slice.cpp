#include "Slice.h"
#include "SliceVectorized.h"

#include <algorithm>
#include <execution>

Points3D::Points3D(
    std::vector<int>& frames,
    std::vector<int>& scans,
    std::vector<double>& mzs,
    std::vector<double>& invIonMobs,
    std::vector<int>& intensities
): frame(frames), scan(scans), intensity(intensities), invIonMobility(invIonMobs), mz(mzs) {}

TimsSlicePL::TimsSlicePL(std::vector<TimsFramePL>& pf, std::vector<TimsFramePL>& ff): precursors(pf), fragments(ff) {}

TimsSlicePL TimsSlicePL::filterRanged(int scanMin, int scanMax, double mzMin, double mzMax, int intensityMin, double rtMin, double rtMax){

    std::vector<TimsFramePL> filteredPrecursors;
    std::vector<TimsFramePL> filteredFragments;

    std::copy_if(this->precursors.begin(), this->precursors.end(), std::back_inserter(filteredPrecursors), [&rtMin, &rtMax](TimsFramePL& f){return (rtMin <= f.retentionTime) && (f.retentionTime <= rtMax) ;} );
    std::copy_if(this->fragments.begin(), this->fragments.end(), std::back_inserter(filteredFragments), [&rtMin, &rtMax](TimsFramePL& f){return (rtMin <= f.retentionTime) && (f.retentionTime <= rtMax) ;} );

    std::vector<TimsFramePL> retPrecursors;
    std::vector<TimsFramePL> retFragments;

    retPrecursors.resize(filteredPrecursors.size());
    retFragments.resize(filteredFragments.size());

    auto filterFrame = [&scanMin, &scanMax, &mzMin, &mzMax, &intensityMin](TimsFramePL& f) -> TimsFramePL {
        return f.filterRanged(scanMin, scanMax, mzMin, mzMax, intensityMin);
    };

    std::transform(std::execution::par_unseq, filteredPrecursors.begin(), filteredPrecursors.end(), retPrecursors.begin(), filterFrame);
    std::transform(std::execution::par_unseq, filteredFragments.begin(), filteredFragments.end(), retFragments.begin(), filterFrame);

    return {retPrecursors, retFragments};
}

Points3D TimsSlicePL::getPoints3D(bool precursor){

    std::vector<int> retFrames, retScan, retIntensity;
    std::vector<double> retMz, retInvIonMob;

    if(precursor)
    {
        for(auto& frame: this->precursors)
        {
            retFrames.insert(retFrames.end(), frame.mzs.size(), frame.frameId);
            retScan.insert(retScan.end(), frame.scans.begin(), frame.scans.end());
            retMz.insert(retMz.end(), frame.mzs.begin(), frame.mzs.end());
            retIntensity.insert(retIntensity.end(), frame.intensities.begin(), frame.intensities.end());
            retInvIonMob.insert(retInvIonMob.end(), frame.inv_ion_mobs.begin(), frame.inv_ion_mobs.end());
        }
    }

    else
    {
        for(auto& frame: this->fragments)
        {
            retFrames.insert(retFrames.end(), frame.mzs.size(), frame.frameId);
            retScan.insert(retScan.end(), frame.scans.begin(), frame.scans.end());
            retMz.insert(retMz.end(), frame.mzs.begin(), frame.mzs.end());
            retIntensity.insert(retIntensity.end(), frame.intensities.begin(), frame.intensities.end());
            retInvIonMob.insert(retInvIonMob.end(), frame.inv_ion_mobs.begin(), frame.inv_ion_mobs.end());
        }
    }

    return {retFrames, retScan, retMz, retInvIonMob, retIntensity};
}

TimsSliceVectorizedPL TimsSlicePL::getVectorizedSlice(int resolution) {

        std::vector<TimsFrameVectorizedPL> retPrecursors;
        std::vector<TimsFrameVectorizedPL> retFragments;

        retPrecursors.resize(this->precursors.size());
        retFragments.resize(this->fragments.size());

        auto vectorizeFrame = [&resolution](TimsFramePL& f) -> TimsFrameVectorizedPL {
            return f.vectorize(resolution);
        };

        std::transform(std::execution::par_unseq, this->precursors.begin(), this->precursors.end(), retPrecursors.begin(), vectorizeFrame);
        std::transform(std::execution::par_unseq, this->fragments.begin(), this->fragments.end(), retFragments.begin(), vectorizeFrame);

        return {retPrecursors, retFragments};
}