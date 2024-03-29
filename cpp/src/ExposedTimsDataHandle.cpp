#include "ExposedTimsDataHandle.h"

#include "opentims_all.cpp"
#include <cmath>
#include <ostream>

/**
 * helper for TDH initialization
 * @param dp data path
 * @param bp binary path
 * @return a tims data handle
 */
TimsDataHandle get_tdh(const std::string& dp, const std::string& bp){
    // need to translate TOF to MZ
    DefaultTof2MzConverterFactory::setAsDefault<BrukerTof2MzConverterFactory, const char*>(bp.c_str());
    // need to translate SCAN to 1/K0
    DefaultScan2InvIonMobilityConverterFactory::setAsDefault<BrukerScan2InvIonMobilityConverterFactory, const char*>(bp.c_str());

    return {dp};
}

class ExposedTimsDataHandle::DataHandle final {
public:
    DataHandle(const std::string& dp, const std::string& bp): handle(get_tdh(dp, bp)) {}

    TimsDataHandle handle;
};

ExposedTimsDataHandle::ExposedTimsDataHandle(const std::string& dp, const std::string& bp) :
    datasetPath(dp),
    binaryPath(bp),
    handle_(std::make_unique<DataHandle>(dp, bp)) {}

ExposedTimsDataHandle::~ExposedTimsDataHandle() = default;


TimsFramePL ExposedTimsDataHandle::getTimsFramePL(const int frameId) {
    // allocate buffer
    const size_t buffer_size_needed = handle_->handle.max_peaks_in_frame();
    std::unique_ptr<uint32_t[]> scan_ids = std::make_unique<uint32_t[]>(buffer_size_needed);
    std::unique_ptr<uint32_t[]> intens = std::make_unique<uint32_t[]>(buffer_size_needed);
    std::unique_ptr<uint32_t[]> tofs = std::make_unique<uint32_t[]>(buffer_size_needed);
    std::unique_ptr<double[]> inv_ion_mob = std::make_unique<double[]>(buffer_size_needed);
    std::unique_ptr<double[]> mz = std::make_unique<double[]>(buffer_size_needed);
    std::unique_ptr<double[]> rts = std::make_unique<double[]>(buffer_size_needed);

    // fetch
    TimsFrame& frame = handle_->handle.get_frame(frameId);
    frame.save_to_buffs(nullptr, scan_ids.get(), tofs.get(), intens.get(), mz.get(), inv_ion_mob.get(), rts.get());

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

    double retentionTime = rts[0];

    return {frameId, retentionTime, scans, mzs, intensities, tof, inv_ion_mobility};
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

/**
 * call bruker binary on first frame to generate mz-axis from tof indices
 * @return a double array of mz values, indices are tof-index
 */
std::vector<double> ExposedTimsDataHandle::getGlobalMzAxis() {

    // allocate vector
    const int buffer_size_needed = 4e5;
    std::vector<uint32_t> tofs(buffer_size_needed);
    std::vector<double> mz(buffer_size_needed);

    for(int i = 0; i < buffer_size_needed; i++)
        tofs[i] = i;

    handle_->handle.tof2mz_converter->convert(1, mz.data(), tofs.data(), buffer_size_needed);

    return mz;
}

TimsBlockPL ExposedTimsDataHandle::getTimsBlockPL(std::vector<int> &frameIds) {
    // allocate concrete vectors
    std::vector<double> rts, mzs, inv_ion_mobility;
    std::vector<int> frames, scans, intensities, tof;


    for(auto frameId: frameIds){
        // allocate buffer
        const size_t buffer_size_needed = handle_->handle.max_peaks_in_frame();
        std::unique_ptr<uint32_t[]> scan_ids = std::make_unique<uint32_t[]>(buffer_size_needed);
        std::unique_ptr<uint32_t[]> intens = std::make_unique<uint32_t[]>(buffer_size_needed);
        std::unique_ptr<uint32_t[]> tofs = std::make_unique<uint32_t[]>(buffer_size_needed);
        std::unique_ptr<double[]> rt = std::make_unique<double[]>(buffer_size_needed);
        std::unique_ptr<double[]> inv_ion_mob = std::make_unique<double[]>(buffer_size_needed);
        std::unique_ptr<double[]> mz = std::make_unique<double[]>(buffer_size_needed);

        // fetch
        TimsFrame& frame = handle_->handle.get_frame(frameId);
        frame.save_to_buffs(nullptr, scan_ids.get(), tofs.get(), intens.get(), mz.get(), inv_ion_mob.get(), rt.get());

        // copy
        for(size_t peak_id = 0; peak_id < frame.num_peaks; peak_id++) {
            tof.push_back(tofs[peak_id]);
            inv_ion_mobility.push_back(inv_ion_mob[peak_id]);
            intensities.push_back(intens[peak_id]);
            mzs.push_back(mz[peak_id]);
            scans.push_back(scan_ids[peak_id]);
            rts.push_back(rt[peak_id]);
            frames.push_back(frameId);
        }
    }

    return {frames, scans, tof, intensities, rts, inv_ion_mobility, mzs};
}