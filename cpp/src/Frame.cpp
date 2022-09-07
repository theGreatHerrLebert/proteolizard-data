#include "Frame.h"

#include "Eigen/Sparse"
#include <algorithm>
#include <execution>
#include <cmath>
#include <map>
#include <tuple>
#include <utility>

Eigen::MatrixXd toDenseEigen(const MzVectorPL& mzVector, int numRows){

    Eigen::SparseMatrix<double> sparseVec = Eigen::SparseMatrix<double>(numRows, 1);
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(mzVector.indices.size());

    for(std::size_t i = 0; i < mzVector.indices.size(); i++)
        tripletList.emplace_back(mzVector.indices[i], 0, mzVector.values[i]);

    sparseVec.setFromTriplets(tripletList.begin(), tripletList.end());
    return sparseVec.toDense();
}

// note: this function is not a member function!
TimsFramePL operator+(const TimsFramePL &leftFrame, const TimsFramePL &rightFrame){

    std::map<std::pair<int, double>, int> sumMap;

    // insert leftFrame values into map
    for (auto it = leftFrame.mzs.begin(); it != leftFrame.mzs.end(); ++it) {
        auto i = std::distance(leftFrame.mzs.begin(), it);
        auto scan = leftFrame.scans[i];
        auto index = leftFrame.mzs[i];
        auto intensity = leftFrame.intensities[i];
        sumMap[{scan, index}] = intensity;
    }

    // insert right frame values into map or sum
    for (auto it = rightFrame.mzs.begin(); it != rightFrame.mzs.end(); ++it) {
        auto i = std::distance(rightFrame.mzs.begin(), it);
        auto scan = rightFrame.scans[i];
        auto index = rightFrame.mzs[i];
        auto intensity = rightFrame.intensities[i];

        // if key is already in map, sum intensity
        if (sumMap.contains({scan, index}))
            sumMap[{scan, index}] += intensity;
            // create new key value pair otherwise
        else
            sumMap[{scan, index}] = intensity;
    }

    std::vector<int> retScans;
    retScans.reserve(sumMap.size());
    std::vector<double> retIndices;
    retIndices.reserve(sumMap.size());
    std::vector<int> retValues;
    retValues.reserve(sumMap.size());


    for(auto const& [key, value]: sumMap){
        retScans.push_back(key.first);
        retIndices.push_back(key.second);
        retValues.push_back(value);
    }

    return {leftFrame.frameId, retScans, retIndices, retValues, {}, {}};
}

HashBlock::HashBlock(int c, std::vector<int> ri, std::vector<int> s, std::vector<int> b, std::vector<int> i, std::vector<int> v):
    counter(c),
    rowIndex(std::move(ri)),
    scan(std::move(s)),
    bin(std::move(b)),
    indices(std::move(i)),
    values(std::move(v)) {}

TimsFramePL::TimsFramePL(int id, std::vector<int>  scan, std::vector<double>  mz, std::vector<int>  intensity, std::vector<int>  tof, std::vector<double>  inv_ion_mob):
    frameId(id),
    scans(std::move(scan)),
    mzs(std::move(mz)),
    intensities(std::move(intensity)),
    tofs(std::move(tof)),
    inv_ion_mobs(std::move(inv_ion_mob)) {}

TimsFramePL TimsFramePL::toResolution(const int resolution) {

    std::map<std::pair<int, int>, float> intensityMap;
    double factor = pow(10.0, resolution);

    std::vector<double> resMz;
    std::vector<int> resS, resI;

    for(std::size_t i = 0; i < this->mzs.size(); i++){

        // calculate binned mz value as key_first
        int rounded = int(round(this->mzs[i] * factor));
        // get scan coordinate as key_second
        int s = this->scans[i];
        // add intensity at combined key location
        intensityMap[{s, rounded}] += this->intensities[i];
    }

    resMz.reserve(intensityMap.size());
    resI.reserve(intensityMap.size());
    resS.reserve(intensityMap.size());

    // vectorize binned spectrum
    for (const auto& [key, value] : intensityMap) {
        resMz.push_back(double(key.second) / factor);
        resS.push_back(key.first);
        resI.push_back(value);
    }

    return {this->frameId, resS, resMz, resI, this->tofs, this->inv_ion_mobs};
}

TimsFramePL TimsFramePL::filterRanged(int scanMin, int scanMax, double mzMin, double mzMax, int minIntensity) {

    std::vector<int> retScan, retIntensity, retTof;
    std::vector<double> retMz, retInv;

    // TODO: make this more efficient (binary search?)
    for(std::size_t i = 0; i < this->mzs.size(); i++){

        int scan = this->scans[i];
        double mz = this->mzs[i];
        int tof = this->tofs[i];
        int intensity = this->intensities[i];
        double invMob = this->inv_ion_mobs[i];

        if((scan >= scanMin) && (scan <= scanMax) && (mz >= mzMin) && (mz <= mzMax) && (intensity >= minIntensity)){
            retScan.push_back(scan);
            retMz.push_back(mz);
            retIntensity.push_back(intensity);
            retInv.push_back(invMob);
            retTof.push_back(tof);
        }
    }

    // This guards for empty return
    if(!retScan.empty())
        return TimsFramePL(this->frameId, retScan, retMz, retIntensity, retTof, retInv);

    return TimsFramePL(this->frameId, {-1}, {-1}, {-1}, {-1}, {-1});
}

TimsFramePL TimsFramePL::fold(const int resolution, const int width) {
    // bin frame to same resolution along mz
    TimsFramePL frame = this->toResolution(resolution);

    // declare and reserve binned scan vector
    std::vector<int> scanFoldIds;
    scanFoldIds.reserve(frame.scans.size());

    // calculate and save scan -> scan_bin
    for(const auto& scan_index: frame.scans)
        scanFoldIds.push_back(int(floor(scan_index / width)));

    // declare summed intensity map
    std::map<std::tuple<int, double>, int> sumMap;

    // sum intensities based on index
    for(std::size_t i=0; i<frame.mzs.size(); i++)
        sumMap[{scanFoldIds[i], frame.mzs[i]}] += frame.intensities[i];

    //write back to frame
    std::vector<int> scans, intensities;
    std::vector<double> mzs;
    for(const auto& [key, intensity]: sumMap){
        const auto& [scan, mz] = key;
        scans.push_back(scan);
        mzs.push_back(mz);
        intensities.push_back(intensity);
    }

    return {frameId, scans, mzs, intensities, this->tofs, this->inv_ion_mobs};
}

TimsFrameVectorizedPL TimsFramePL::vectorize(const int resolution) {

    std::map<std::pair<int, int>, int> sumMap;
    int factor = pow(10, resolution);

    std::vector<int> retScan, retIndex;
    std::vector<int> retValue;

    for (auto it = this->mzs.begin(); it != this->mzs.end(); ++it) {

        auto i = std::distance(this->mzs.begin(), it);
        auto scan = this->scans[i];
        auto index = int(floor(round(this->mzs[i] * factor)));
        auto intensity = this->intensities[i];

        if(sumMap.contains({scan, index}))
            sumMap[{scan, index}] += intensity;
        else
            sumMap[{scan, index}] = intensity;
    }

    retScan.reserve(sumMap.size());
    retIndex.reserve(sumMap.size());
    retValue.reserve(sumMap.size());

    // vectorize binned spectrum
    for (const auto& [key, value] : sumMap) {
        retIndex.push_back(key.second);
        retScan.push_back(key.first);
        retValue.push_back(value);
    }

    return {this->frameId, resolution, retScan, retIndex, retValue};
}

std::map<int, MzSpectrumPL> TimsFramePL::spectra() {

    std::map<int, MzSpectrumPL> specMap;

    for(size_t peak_id = 0; peak_id < this->scans.size(); peak_id++){
        auto scan = this->scans[peak_id];
        auto mz = this->mzs[peak_id];
        auto i = this->intensities[peak_id];
        if(specMap.contains(scan)){
            specMap[scan].mz.push_back(mz);
            specMap[scan].intensity.push_back(i);
        }
        else {
            specMap[scan] = MzSpectrumPL{this->frameId, static_cast<int>(scan), {mz}, {i}};
        }
    }

    return specMap;
}

std::vector<MzSpectrumPL> TimsFramePL::exportSpectra() {

    auto spectra = this->spectra();
    std::vector<MzSpectrumPL> ret;
    ret.reserve(spectra.size());

    for(auto&[key, value]: spectra)
        ret.push_back(value);

    return ret;
}

Eigen::MatrixXd TimsFramePL::denseWindowMatrix(int resolution, int minPeaksPerWindow, int minIntensity,
                                               double windowLength, bool overlapping) {

    // split frame into spectra
    auto spectra = this->spectra();
    auto numRows = int(round(2000 * pow(10, resolution)));
    std::vector<Eigen::MatrixXd> vec;

    for(auto &spectrum: spectra){
        auto windows = spectrum.second.windows(windowLength, overlapping, minPeaksPerWindow, minIntensity);
        for(auto &w: windows){
            auto v = toDenseEigen(w.second.vectorize(resolution), numRows);
            vec.push_back(v);
        }
    }

    auto numColumns = vec.size();

    Eigen::MatrixXd retMatrix(numRows, numColumns);

    for(std::size_t i = 0; i < vec.size(); i++){
        retMatrix.col(i) = vec[i];
    }

    return retMatrix;
}

HashBlock TimsFramePL::getHashingBlocks(
        int resolution,
        int minPeaksPerWindow,
        int minIntensity,
        double windowLength,
        bool overlapping) {

    auto createWindows =
            [&resolution, &minPeaksPerWindow, &minIntensity, &windowLength, &overlapping]
                    (const std::pair<int, MzSpectrumPL> &p) -> LambdaReturn {

                auto windows = p.second.windows(windowLength, overlapping, minPeaksPerWindow, minIntensity);

                std::vector<int> retScan, retBin, retIndices, retValues, retScanWindow, retBinWindow;

                retScan.reserve(5000);
                retBin.reserve(5000);
                retIndices.reserve(5000);
                retValues.reserve(5000);
                retScanWindow.reserve(5000);
                retBinWindow.reserve(5000);


                for(const auto &[bin, window]: windows){
                    // first, calculate the value that has to be subtracted from indices for zero indexing
                    int offset;

                    if (bin > 0) {
                        offset = pow(10, resolution) * bin * windowLength;
                    }
                    else {
                        offset = pow(10, resolution) * ((-bin) * windowLength - (windowLength / 2));
                    }

                    auto vectorizedWindow = window.vectorize(resolution);

                    for(std::size_t i = 0; i < vectorizedWindow.indices.size(); i++){
                        auto index = vectorizedWindow.indices[i];
                        auto zeroIndex = index - offset;
                        auto value = vectorizedWindow.values[i];
                        retIndices.push_back(zeroIndex);
                        retValues.push_back(value);
                    }
                    retScan.insert(retScan.end(), vectorizedWindow.indices.size(), window.scanId);
                    retBin.insert(retBin.end(),vectorizedWindow.indices.size(), bin);
                    retScanWindow.push_back(window.scanId);
                    retBinWindow.push_back(bin);
                }

                return {retScan, retBin, retIndices, retValues, retScanWindow, retBinWindow};
            };

    auto spectra = this->spectra();

    std::vector<LambdaReturn> retWindowVec;
    retWindowVec.resize(spectra.size());

    std::transform(std::execution::par_unseq, spectra.begin(), spectra.end(), retWindowVec.begin(), createWindows);

    std::vector<int> retScan, retBin, retIndices, retValues, retScanWindows, retBinWindows;

    retScan.reserve(80000);
    retBin.reserve(80000);
    retIndices.reserve(80000);
    retValues.reserve(80000);

    retScanWindows.reserve(80000);
    retBinWindows.reserve(80000);

    for(const auto &lr: retWindowVec){
        retScan.insert(retScan.end(), lr.scans.begin(), lr.scans.end());
        retBin.insert(retBin.end(), lr.bins.begin(), lr.bins.end());
        retIndices.insert(retIndices.end(), lr.indices.begin(), lr.indices.end());
        retValues.insert(retValues.end(), lr.values.begin(), lr.values.end());
        retScanWindows.insert(retScanWindows.end(), lr.scansWindow.begin(), lr.scansWindow.end());
        retBinWindows.insert(retBinWindows.end(), lr.binsWindow.begin(), lr.binsWindow.end());
    }

    std::map<std::pair<int, int>, int> windowIndices;

    int counter = 0;
    for(std::size_t i = 0; i < retScan.size(); i++){
        std::pair<int, int> p = {retScan[i], retBin[i]};
        if(!windowIndices.contains(p)){
            windowIndices[p] = counter;
            counter++;
        }
    }

    std::vector<int> rowIndex;
    rowIndex.resize(retScan.size());

    for(std::size_t i = 0; i < retScan.size(); i++){
        std::pair<int, int> key = {retScan[i], retBin[i]};
        rowIndex[i] = windowIndices[key];
    }

    return {counter, rowIndex, retScanWindows, retBinWindows, retIndices, retValues};
}

std::pair<std::pair<std::vector<int>, std::vector<int>>, std::vector<MzVectorPL>> TimsFramePL::vectorizedWindows(
        int resolution, int minPeaksPerWindow, int minIntensity,double windowLength, bool overlapping) {


    // split frame into spectra
    auto spectra = this->spectra();
    std::vector<int> scans, bins;
    std::vector<MzVectorPL> vectorizedWindows;

    for(auto &[scan, spectrum]: spectra){
        auto windows = spectrum.windows(windowLength, overlapping, minPeaksPerWindow, minIntensity);
        for(auto &[bin, window]: windows){
            vectorizedWindows.push_back(window.vectorize(resolution));
            scans.push_back(scan);
            bins.push_back(bin);
        }
    }

    return {{scans, bins}, vectorizedWindows};

}

LambdaReturn::LambdaReturn(std::vector<int> s, std::vector<int> b, std::vector<int> i, std::vector<int> v, std::vector<int> sw, std::vector<int> bw):
    scans(std::move(s)),
    bins(std::move(b)),
    indices(std::move(i)),
    values(std::move(v)),
    scansWindow(std::move(sw)),
    binsWindow(std::move(bw)) {}
