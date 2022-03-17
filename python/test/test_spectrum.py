import numpy as np

from pyproteolizard.data import MzSpectrum

# create list of mz values
mz = [389.24769132, 391.23760966, 394.22574895, 394.23824935,
      396.22520107, 401.2191587, 407.15853975, 408.23270771,
      409.21595143, 418.1376726, 560.33372185, 614.7691754,
      621.29192559, 627.79783109, 628.29483337, 628.79992595,
      634.27023685, 647.28660589, 650.81530219, 653.23075439,
      655.75950157, 659.76795231, 666.24831891, 672.30319886,
      674.82767582, 676.73846911, 683.20271353, 687.73564909,
      691.22794737, 693.27369563, 696.69691792, 704.14908821,
      707.70369977, 709.72754325, 710.18047501, 745.67061534,
      870.94955762]

# create list of intensity
intensity = [102, 114, 108, 105, 187, 107, 113, 121, 155, 109, 105, 121, 156,
             309, 130, 153, 139, 105, 135, 112, 107, 115, 100, 102, 101, 100,
             115, 109, 115, 101, 101, 105, 108, 100, 108, 105, 135]

# create spectrum from mz and intensity
spectrum = MzSpectrum(None, 1, 1, mz, intensity)


# create list of mz values
mz_2 = [389.24769132, 389.25169132, 394.22874895, 394.23124935]

# create list of intensity
intensity_2 = [1, 2, 3, 4]

# create spectrum from mz and intensity
spectrum_2 = MzSpectrum(None, 1, 1, mz_2, intensity_2)


def test_construction():
    """
    test if copy of values back and forth between C++ and python does not alter data
    """
    # get values back from C++
    mz_spec = spectrum.mz()
    intensity_spec = spectrum.intensity()

    # check if copy back and forth between python and C++ changes things
    assert set(mz) == set(mz_spec) and set(intensity) == set(intensity_spec)


def test_binning():
    """
    test if binning implemented in C++ returns same result as python implementation
    """

    assert set(spectrum_2.vectorize(2).values()) == {3, 7}


def test_grouping():
    """
    test if binning implemented in C++ returns same result as python implementation
    """
    mz_bins_1 = [int(np.round_(m, 1) * np.power(10, 1)) for m in spectrum.mz()]
    binned_mz_1 = spectrum.vectorize(1).indices()

    mz_bins_2 = [int(np.round_(m, 2) * np.power(10, 2)) for m in spectrum.mz()]
    binned_mz_2 = spectrum.vectorize(2).indices()

    print([a == b for a, b in zip(mz_bins_1, binned_mz_1)])

    assert set(mz_bins_1) == set(binned_mz_1) and set(mz_bins_2) == set(binned_mz_2)


def test_intensity_grouping():
    """

    """
    assert set(spectrum_2.vectorize(2).values()) == {3, 7}


def test_window_generation():
    """

    """

    # arrange
    mz_values = [504.326, 505.337, 509.129, 511.113, 513.845, 517.85, 521.211]
    intensity_values = [int(x) for x in np.ones_like(mz_values)]
    test_spec = MzSpectrum(None, 1, 1, mz_values, intensity_values)

    bins_true = [50, 51, 52]

    windows_true = [[504.326, 505.337, 509.129], [511.113, 513.845, 517.85], [521.211]]

    # action
    bins, windows = test_spec.windows(10, overlapping=False, min_intensity=0, min_peaks=1)
    mz_windows = [w.mz() for w in windows]
    both = np.all([np.all(x == y) for x, y in zip(windows_true, mz_windows)])

    # assert
    assert set(bins) == set(bins_true) and both


def test_window_generation_overlapping():
    """

    """

    # arrange
    mz_values = [504.326, 505.337, 509.129, 511.113, 513.845, 517.85, 521.211]
    intensity_values = [int(x) for x in np.ones_like(mz_values)]
    test_spec = MzSpectrum(None, 1, 1, mz_values, intensity_values)

    # action
    bins, windows = test_spec.windows(10, overlapping=True, min_intensity=0, min_peaks=1)
    both = [x for x in list(zip(bins, [w.mz() for w in windows])) if x[0] < 0]
    bins, windows = zip(*both)

    bins_true = [-52, -51, -50]
    windows_true = [[517.85, 521.211], [505.337, 509.129, 511.113, 513.845], [504.326]]

    both = np.all([np.all(x == y) for x, y in zip(windows_true, windows)])

    # assert
    assert set(bins) == set(bins_true) and both


def main():
    test_construction()
    test_binning()
    test_grouping()
    test_window_generation()
    test_window_generation_overlapping()


if __name__ == '__main__':
    main()
