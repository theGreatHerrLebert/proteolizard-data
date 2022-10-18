#ifndef CPP_FRAME_H
#define CPP_FRAME_H

#include "Spectrum.h"
#include "VectorizedFrame.h"
#include "Eigen/Dense"

/**
 * create eigen sparse vector from vectorized mz spectrum
 * @param mzVector : vectorized mz spectrum to convert
 * @param numRows : dimensionality of vector
 * @return : a sparse eigen vector suited for fast vectorized operations
 */
Eigen::MatrixXd toDenseEigen(const MzVectorPL& mzVector, int numRows);

class HashBlock{
public:
    HashBlock(int c, std::vector<int> ri, std::vector<int> s, std::vector<int> b, std::vector<int> i, std::vector<int> v);

    int counter;
    std::vector<int> rowIndex;
    std::vector<int> scan;
    std::vector<int> bin;
    std::vector<int> indices;
    std::vector<int> values;
};

/**
 * same container but trying to provide cleaner OOP interface
 */
class TimsFramePL {
public:
    // constructors
    TimsFramePL()= default;
    TimsFramePL(int id, double rt, std::vector<int>  scan, std::vector<double>  mz, std::vector<int>  intensity, std::vector<int>  tof, std::vector<double>  inv_ion_mob);

    // binning to a finite resolution
    TimsFramePL filterRanged(int scanMin, int scanMax, double mzMin, double mzMax, int minIntensity=1);

    TimsFramePL toResolution(int resolution);
     friend TimsFramePL operator+(const TimsFramePL &leftFrame, const TimsFramePL &rightFrame);
    TimsFramePL fold(int resolution, int foldWidth);

    // vectorize to an integer mz value
    TimsFrameVectorizedPL vectorize(int resolution);

    // get all spectra as map from scan to spectrum
    std::map<int, MzSpectrumPL> spectra();

    std::vector<MzSpectrumPL> exportSpectra();

    std::pair<std::pair<std::vector<int>, std::vector<int>>, std::vector<MzVectorPL>> vectorizedWindows(int resolution, int minPeaksPerWindow, int minIntensity,double windowLength, bool overlapping);

    HashBlock getHashingBlocks(int resolution, int minPeaksPerWindow, int minIntensity, double windowLength, bool overlapping);

    Eigen::MatrixXd denseWindowMatrix(int resolution, int minPeaksPerWindow,
                           int minIntensity, double windowLength, bool overlapping);

    int frameId {};
    double retentionTime {};
    std::vector<int> scans {};
    std::vector<double> mzs {};
    std::vector<int> intensities {};
    std::vector<int> tofs {};
    std::vector<double> inv_ion_mobs {};
};

TimsFramePL operator+(const TimsFramePL &leftFrame, const TimsFramePL &rightFrame);

class LambdaReturn {
public:
    LambdaReturn() = default;
    LambdaReturn(std::vector<int> s, std::vector<int> b, std::vector<int> i, std::vector<int> v, std::vector<int> sw, std::vector<int> bw);

    std::vector<int> scans;
    std::vector<int> bins;
    std::vector<int> indices;
    std::vector<int> values;
    std::vector<int> scansWindow;
    std::vector<int> binsWindow;
};

#endif //CPP_FRAME_H
