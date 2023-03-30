#include "Spectrum.h"
#include <numeric>
#include <cmath>

MzSpectrumPL operator+(const MzSpectrumPL &leftSpec, const MzSpectrumPL &rightSpec){
    MzSpectrumPL l = leftSpec;
    return l += rightSpec;

}

MzSpectrumPL& MzSpectrumPL::operator+=(const MzSpectrumPL &rightSpec){
    std::map<double, int> sumMap;

    // insert this's values into map
    for (auto it = this->mz.begin(); it != this->mz.end(); ++it) {
        auto i = std::distance(this->mz.begin(), it);
        auto index = this->mz[i];
        auto intensity = this->intensity[i];
        sumMap[index] = intensity;
    }

    // insert right frame values into map or sum
    for (auto it = rightSpec.mz.begin(); it != rightSpec.mz.end(); ++it) {
        auto i = std::distance(rightSpec.mz.begin(), it);
        auto index = rightSpec.mz[i];
        auto intensity = rightSpec.intensity[i];

        // if key is already in map, sum intensity
        if (sumMap.contains(index))
            sumMap[index] += intensity;
            // create new key value pair otherwise
        else
            sumMap[index] = intensity;
    }

    std::vector<double> retIndices;
    retIndices.reserve(sumMap.size());
    std::vector<int> retValues;
    retValues.reserve(sumMap.size());


    for(auto const& [key, value]: sumMap){
        retIndices.push_back(key);
        retValues.push_back(value);
    }
    this->mz = retIndices;
    this->intensity = retValues;

    return *this;
}

MzSpectrumPL operator*(const MzSpectrumPL &leftSpec, const float scalar){
    std::vector<int> mult_i = leftSpec.intensity;
    for (size_t i = 0; i < mult_i.size(); i++){
        mult_i[i] *= scalar;
        }
    return MzSpectrumPL(leftSpec.frameId,leftSpec.scanId,leftSpec.mz,mult_i).filter(-1,-1, 1);
}

MzSpectrumPL operator*(const float scalar, const MzSpectrumPL &rightSpec){
    return rightSpec*scalar;
}

MzSpectrumPL::MzSpectrumPL(int frame, int scan, std::vector<double> m, std::vector<int> i):
    frameId(frame), scanId(scan), mz(std::move(m)), intensity(std::move(i)) {}

MzVectorPL MzSpectrumPL::vectorize(int resolution) const{

    // auto tmp = this->toResolution(resolution);

    std::map<int, int> intensityMap;
    double factor = pow(10.0, resolution);

    std::vector<int> indices;
    std::vector<int> values;

    for(size_t i = 0; i < this->mz.size(); i++){

        // calculate binned mz value as key
        int rounded = int(round(this->mz[i] * factor));
        // add intensity to generated key
        intensityMap[rounded] += this->intensity[i];
    }

    indices.reserve(intensityMap.size());
    values.reserve(intensityMap.size());

    // get all mz values and sort them
    for (const auto& [key, value] : intensityMap) {
        indices.push_back(key);
        values.push_back(value);
    }

    return MzVectorPL{resolution, this->frameId, this->scanId, indices, values};
}

MzSpectrumPL MzSpectrumPL::filter(double mzMin, double mzMax, int intensityMin) const {

    std::vector<int> retIntensity;
    std::vector<double> retMz;

    // TODO: make this more efficient (binary search?)
    for (std::size_t i = 0; i < this->mz.size(); i++) {

        double mz = this->mz[i];
        int intensity = this->intensity[i];

        if ( ((mz >= mzMin) || (mzMin == -1) ) &&
             ((mz <= mzMax) || (mzMax == -1) ) &&
             ((intensity >= intensityMin) || (intensityMin == -1))
           ){
            retMz.push_back(mz);
            retIntensity.push_back(intensity);
        }
    }

    return {this->frameId, this->scanId, retMz, retIntensity};
}

MzSpectrumPL MzSpectrumPL::toResolution(int resolution) const{

    std::map<int, int> intensityMap;
    double factor = pow(10.0, resolution);

    std::vector<double> resMz;
    std::vector<int> resI;

    for(size_t i = 0; i < this->mz.size(); i++){

        // calculate binned mz value as key
        int rounded = int(round(this->mz[i] * factor));
        // add intensity to generated key
        intensityMap[rounded] += this->intensity[i];
    }

    resMz.reserve(intensityMap.size());
    resI.reserve(intensityMap.size());

    // get all mz values and sort them
    for (const auto& [key, value] : intensityMap) {
        resMz.push_back(double(key) / factor);
        resI.push_back(value);
    }

    return MzSpectrumPL{this->frameId, this->scanId, resMz, resI};
}

MzSpectrumPL MzSpectrumPL::toCentroided(int baselineNoiseLevel, double sigma) const{
    // filter away noise below baseline noise level
    MzSpectrumPL filtered = this->filter(-1,-1,baselineNoiseLevel);

    // init centroided data
    std::vector<double> centMz;
    std::vector<int> centI;

    double lastMz = 0, currentMz = 0, currentIntensity = 0;
    int sumI = 0;
    double meanMz = 0;
    for(size_t i = 0; i< this->mz.size(); i++){
        currentMz = this->mz[i];
        currentIntensity = this->intensity[i];

        // if peak is too far away from last peak push centroid
        if (currentMz-lastMz > sigma && meanMz > 0){
            meanMz /= sumI;
            // push centroid
            centMz.push_back(meanMz);
            centI.push_back(sumI);

            // start new centroid
            sumI = 0;
            meanMz = 0;
        }
        meanMz += currentMz*currentIntensity;
        sumI += currentIntensity;
        lastMz = currentMz;

    }
    // push back last remaining centroid
    if (meanMz > 0){
        meanMz /= sumI;
        centMz.push_back(meanMz);
        centI.push_back(sumI);
    }



    return MzSpectrumPL{this->frameId, this->scanId, centMz, centI};
}

MzSpectrumPL& MzSpectrumPL::push(MzSpectrumPL& other){
    this->mz.insert(this->mz.end(), other.mz.begin(), other.mz.end());
    this->intensity.insert(this->intensity.end(), other.intensity.begin(), other.intensity.end());
    return *this;
}

/**
 * check if a given windowId is already present in map
 * @param windowCollection collection to check
 * @param windowId id to check for
 * @return true if key is in map else false
 */
bool hasValue(std::map<int, MzSpectrumPL>& windowCollection, int windowId){
    return windowCollection.count(windowId);
}

std::map<int, MzSpectrumPL> MzSpectrumPL::windows(double windowLength, bool overlapping, int minPeaks, int minIntensity) const
{
    std::map<int, MzSpectrumPL> splits;

    for(std::size_t i = 0; i < this->mz.size(); i++) {

        // get mz and intensity
        auto mmz = this->mz[i];
        auto iintensity = this->intensity[i];

        // calculate window key mz and intensity value should be sent to
        auto tmpKey = int(floor(mmz / windowLength));

        // given key is in map
        if (hasValue(splits, tmpKey)) {
            splits[tmpKey].mz.push_back(mmz);
            splits[tmpKey].intensity.push_back(iintensity);
        }

            // given key is not in map
        else {
            std::vector<double> tmpMz = {mmz};
            std::vector<int> tmpI = {iintensity};
            splits[tmpKey] = MzSpectrumPL{this->frameId, this->scanId, tmpMz, tmpI};
        }
    }

    // calculate grouping by offset
    if(overlapping) {
        std::map<int, MzSpectrumPL> splitsOffset;

        for (std::size_t i = 0; i < this->mz.size(); i++) {
            auto mmz = this->mz[i];
            auto iintensity = this->intensity[i];

            // calculate window key with offset mz and intensity value should be sent to
            auto tmpKey = -int(floor(((mmz + windowLength / 2.0) / windowLength)));

            // given key is in map
            if (hasValue(splitsOffset, tmpKey)) {
                splitsOffset[tmpKey].mz.push_back(mmz);
                splitsOffset[tmpKey].intensity.push_back(iintensity);
            }
                // given key is not in map
            else {
                std::vector<double> tmpMz = {mmz};
                std::vector<int> tmpI = {iintensity};
                splitsOffset[tmpKey] = MzSpectrumPL{this->frameId, this->scanId, tmpMz, tmpI};
            }
        }
        splits.merge(splitsOffset);
    }

    std::map<int, MzSpectrumPL> retSplits;

    // check for minPeaks and minIntensity per window
    for(const auto& [bin, spectrum]: splits){

        // check minPeaks
        if(spectrum.mz.size() >= minPeaks){
            auto it = *max_element(std::begin(spectrum.intensity), std::end(spectrum.intensity));

            // check minIntensity
            if(it >= minIntensity)
                retSplits[bin] = spectrum;
        }
    }

    return retSplits;
}

std::pair<std::vector<int>, std::vector<MzSpectrumPL>> MzSpectrumPL::exportWindows(
        double windowLength,
        bool overlapping,
        int minPeaks,
        int minIntensity) const{

    auto windowMap = this->windows(windowLength, overlapping, minPeaks, minIntensity);

    std::vector<int> retBins;
    retBins.reserve(windowMap.size());
    std::vector<MzSpectrumPL> retWindows;
    retWindows.reserve(windowMap.size());

    for(const auto& [k, v]: windowMap){
        retBins.push_back(k);
        retWindows.push_back(v);
    }

    return {retBins, retWindows};
}
