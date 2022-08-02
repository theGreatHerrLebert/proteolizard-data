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