//
// Created by david on 20.07.21.
//

#ifndef CPP_TIMS_SLICE_H
#define CPP_TIMS_SLICE_H

#include <tuple>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>
#include <execution>
#include "Frame.h"

struct Points3D {
    std::vector<int> frame, scan, intensity;
    std::vector<double> invIonMobility, mz;
    Points3D(std::vector<int>& frames, std::vector<int>& scans,
    std::vector<double>& mzs, std::vector<double>& invIonMobs, std::vector<int>& intensities);
};

Points3D::Points3D(std::vector<int>& frames,
std::vector<int>& scans,
std::vector<double>& mzs,
std::vector<double>& invIonMobs,
std::vector<int>& intensities): frame(frames), scan(scans), mz(mzs),invIonMobility(invIonMobs), intensity(intensities){}

/**
 * @brief
 *
 */
struct TimsSlicePL {
    std::vector<TimsFramePL> precursors;
    std::vector<TimsFramePL> fragments;

    TimsSlicePL(std::vector<TimsFramePL>& pf, std::vector<TimsFramePL>& ff);
    TimsSlicePL filterRanged(int scanMin, int scanMax, double mzMin, double mzMax, int intensityMin);
    Points3D getPoints3D(bool precursor);
};

TimsSlicePL::TimsSlicePL(std::vector<TimsFramePL>& pf, std::vector<TimsFramePL>& ff): precursors(pf), fragments(ff) {}

TimsSlicePL TimsSlicePL::filterRanged(int scanMin, int scanMax, double mzMin, double mzMax, int intensityMin){

    std::vector<TimsFramePL> retPrecursors;
    std::vector<TimsFramePL> retFragments;

    retPrecursors.resize(this->precursors.size());
    retFragments.resize(this->fragments.size());

    auto filterFrame = [&scanMin, &scanMax, &mzMin, &mzMax, &intensityMin](TimsFramePL& f) -> TimsFramePL {
        return f.filterRanged(scanMin, scanMax, mzMin, mzMax, intensityMin);
    };

   std::transform(std::execution::par_unseq, this->precursors.begin(), this->precursors.end(), retPrecursors.begin(), filterFrame);
   std::transform(std::execution::par_unseq, this->fragments.begin(), this->fragments.end(), retFragments.begin(), filterFrame);

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

#endif //CPP_TIMS_SLICE_H
