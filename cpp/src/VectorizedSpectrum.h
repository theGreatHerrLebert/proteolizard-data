#ifndef CPP_VECTORIZED_SPECTRUM_H
#define CPP_VECTORIZED_SPECTRUM_H

#include <vector>

/**
 * container for single vectorized timsTOF mz spectra
 */
class MzVectorPL {
public:
    // constructors
    MzVectorPL() = default;
    MzVectorPL(int res, int frame, int scan, std::vector<int> index, std::vector<int> value);

    friend MzVectorPL operator+(const MzVectorPL &leftSpec, const MzVectorPL &rightSpec);

    int resolution{}; // number of decimals spectrum was binned to
    int frameId{};
    int scanId{};
    std::vector<int> indices; // mz indices, can be seen as vector entries
    std::vector<int> values;
};

#endif //CPP_VECTORIZED_SPECTRUM_H
