#include "TimsBlock.h"

#include <algorithm>

TimsBlockPL::TimsBlockPL(std::vector<int> frames, std::vector<int> scans, std::vector<int> tofs, std::vector<int> intensities,
                         std::vector<double> retentionTimes, std::vector<double> invIonMobs, std::vector<double> mzs):
                         frame(frames), scan(scans), tof(tofs), intensity(intensities),
                         retentionTime(retentionTimes), invIonMobility(invIonMobs), mz(mzs) {}


std::vector<std::vector<int>> TimsBlockPL::getIndices() {
    return {this->frame, this->scan, this->tof, this->intensity};
}

std::vector<std::vector<double>> TimsBlockPL::getValues() {
    return {this->retentionTime, this->invIonMobility, this->mz};
}

TimsBlockPL TimsBlockPL::filterRanged(int scanMin, int scanMax, double mzMin, double mzMax, int intensityMin) {

    std::vector<int> filteredFrame;
    std::vector<int> filteredScan;
    std::vector<int> filteredTof;
    std::vector<int> filteredIntensity;
    std::vector<double> filteredRetentionTime;
    std::vector<double> filteredInvIonMobility;
    std::vector<double> filteredMz;

    for (int i = 0; i < this->frame.size(); i++) {
        if (this->scan[i] >= scanMin && this->scan[i] <= scanMax &&
            this->mz[i] >= mzMin && this->mz[i] <= mzMax &&
            this->intensity[i] >= intensityMin) {
            filteredFrame.push_back(this->frame[i]);
            filteredScan.push_back(this->scan[i]);
            filteredTof.push_back(this->tof[i]);
            filteredIntensity.push_back(this->intensity[i]);
            filteredRetentionTime.push_back(this->retentionTime[i]);
            filteredInvIonMobility.push_back(this->invIonMobility[i]);
            filteredMz.push_back(this->mz[i]);
        }
    }

    return TimsBlockPL(filteredFrame, filteredScan, filteredTof, filteredIntensity, filteredRetentionTime, filteredInvIonMobility, filteredMz);
}