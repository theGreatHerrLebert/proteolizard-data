#ifndef CPP_SPECTRUM_H
#define CPP_SPECTRUM_H

#include "VectorizedSpectrum.h"
#include <map>
#include <vector>

/**
 * container for single timsTOF mz spectra
 */
class MzSpectrumPL {
public:
    // constructors
    MzSpectrumPL()= default;
    MzSpectrumPL(int frame, int scan, std::vector<double> m, std::vector<int> i);
    MzSpectrumPL(std::vector<MzSpectrumPL> &spectra, int resolution, bool centroid, int baselineNoiseLevel, double sigma);
    [[nodiscard]] MzSpectrumPL toResolution(int resolution) const;
    [[nodiscard]] MzVectorPL vectorize(int resolution) const;
    [[nodiscard]] MzSpectrumPL filter(double mzMin, double mzMax, int intensityMin) const;
    [[nodiscard]] std::map<int, MzSpectrumPL> windows(double windowLength, bool overlapping, int minPeaks, int minIntensity) const;
    [[nodiscard]] std::pair<std::vector<int>, std::vector<MzSpectrumPL>> exportWindows(double windowLength, bool overlapping,
                                                                          int minPeaks, int minIntensity) const;
    [[nodiscard]] MzSpectrumPL toCentroided(int baselineNoiseLevel, double sigma) const;
    [[nodiscard]] MzSpectrumPL& push(MzSpectrumPL& other);
    friend MzSpectrumPL operator+(const MzSpectrumPL &leftSpec, const MzSpectrumPL &rightSpec);
    friend MzSpectrumPL operator*(const MzSpectrumPL &Spec, const float scalar);
    friend MzSpectrumPL operator*(const float scalar, const MzSpectrumPL &Spec);
    MzSpectrumPL& operator+=(const MzSpectrumPL &rightSpec);
    int frameId{};
    // double retentionTime{};
    int scanId{};
    // double inv_ion_mobility{};
    std::vector<double> mz; // vector of mz values
    std::vector<int> intensity; // vector of intensities
};

MzSpectrumPL operator+(const MzSpectrumPL &leftSpec, const MzSpectrumPL &rightSpec);
MzSpectrumPL operator*(const MzSpectrumPL &leftSpec, const float scalar);
MzSpectrumPL operator*(const float scalar, const MzSpectrumPL &rightSpec);
#endif //CPP_SPECTRUM_H
