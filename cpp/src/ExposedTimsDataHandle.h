//
// Created by david on 28.07.21.
//

#ifndef CPP_EXPOSEDTIMSDATAHANDLE_H
#define CPP_EXPOSEDTIMSDATAHANDLE_H

#include <vector>
#include <tuple>
#include <ostream>
#include <iostream>
#include "Frame.h"
#include <cmath>
#include <map>
#include <bits/stdc++.h>
#include "../../include/opentims/opentims++/opentims_all.cpp"
#include "Slice.h"

/**
 * helper for TDH initialization
 * @param dp data path
 * @param bp binary path
 * @return a tims data handle
 */
TimsDataHandle get_tdh(std::string dp, std::string bp){
    // need to translate TOF to MZ
    DefaultTof2MzConverterFactory::setAsDefault<BrukerTof2MzConverterFactory, const char*>(bp.c_str());
    // need to translate SCAN to 1/K0
    DefaultScan2InvIonMobilityConverterFactory::setAsDefault<BrukerScan2InvIonMobilityConverterFactory, const char*>(bp.c_str());

    return TimsDataHandle(dp);
}

/**
 *
 */
class ExposedTimsDataHandle{
private:
    TimsDataHandle handle;
    std::vector<TimsFramePL> getTimsFramesFiltered(std::vector<int> frameIds,
                                                         const int scan_min, const int scan_max,
                                                         const double mz_min, const double mz_max);
public:
    // path to dataset, path to bruker binaries
    std::string datasetPath, binaryPath;
    // handle to read from raw data
    ExposedTimsDataHandle(std::string dp, std::string bp);

    TimsFramePL getTimsFramePL(const int frameId);
    TimsSlicePL getTimsSlicePL(std::vector<int>& precursorIds, std::vector<int>& fragmentIds);
};

ExposedTimsDataHandle::ExposedTimsDataHandle(std::string dp, std::string bp) : handle(get_tdh(dp, bp)) {
    datasetPath = dp;
    binaryPath = bp;
};

TimsFramePL ExposedTimsDataHandle::getTimsFramePL(const int frameId) {
    // allocate buffer
    const size_t buffer_size_needed = handle.max_peaks_in_frame();
    std::unique_ptr<uint32_t[]> scan_ids = std::make_unique<uint32_t[]>(buffer_size_needed);
    std::unique_ptr<uint32_t[]> intens = std::make_unique<uint32_t[]>(buffer_size_needed);
    std::unique_ptr<uint32_t[]> tofs = std::make_unique<uint32_t[]>(buffer_size_needed);
    std::unique_ptr<double[]> inv_ion_mob = std::make_unique<double[]>(buffer_size_needed);
    std::unique_ptr<double[]> mz = std::make_unique<double[]>(buffer_size_needed);

    // fetch
    TimsFrame& frame = handle.get_frame(frameId);
    frame.save_to_buffs(nullptr, scan_ids.get(), tofs.get(), intens.get(), mz.get(), inv_ion_mob.get(), nullptr);

    // allocate concrete vectors
    std::vector<double> mzs, inv_ion_mobility;
    std::vector<int> scans, intensities, tof;
    
    mzs.reserve(frame.num_peaks);
    intensities.reserve(frame.num_peaks);
    inv_ion_mobility.reserve(frame.num_peaks);
    scans.reserve(frame.num_peaks);
    tof.reserve(frame.num_peaks);

    // copy
    for(size_t peak_id = 0; peak_id < frame.num_peaks; peak_id++) {
        tof.push_back(tofs[peak_id]);
        inv_ion_mobility.push_back(inv_ion_mob[peak_id]);
        intensities.push_back(intens[peak_id]);
        mzs.push_back(mz[peak_id]);
        scans.push_back(scan_ids[peak_id]);
    }

    return TimsFramePL(frameId, scans, mzs, intensities, tof, inv_ion_mobility);
}

TimsSlicePL ExposedTimsDataHandle::getTimsSlicePL(std::vector<int>& precursorIds, std::vector<int>& fragmentIds){

    std::vector<TimsFramePL> retPrecursors;
    std::vector<TimsFramePL> retFragments;

    retPrecursors.reserve(precursorIds.size());
    retFragments.reserve(fragmentIds.size());

    for(auto id: precursorIds)
    {
    retPrecursors.push_back(this->getTimsFramePL(id));
    }

    for(auto id: fragmentIds)
    {
    retFragments.push_back(this->getTimsFramePL(id));
    }

    return {retPrecursors, retFragments};
}

std::vector<TimsFramePL> ExposedTimsDataHandle::getTimsFramesFiltered(std::vector<int> frameIds,
                                                                            const int scanMin, const int scanMax,
                                                                            const double mzMin,
                                                                            const double mzMax) {
    std::sort(frameIds.begin(), frameIds.end());
    std::vector<TimsFramePL> result;
    result.reserve(frameIds.size());

    for(int id: frameIds)
        result.push_back(getTimsFramePL(id).filterRanged(scanMin, scanMax, mzMin, mzMax));

    return result;
}

#endif //CPP_EXPOSEDTIMSDATAHANDLE_H
