import numpy as np
import tensorflow as tf
from pyproteolizard.data import TimsFrame

scan = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3]
mz = [500.123, 506.123, 507.127, 511.112, 500.123, 506.123, 708.58, 892.123, 500.23, 500.78]
intensity = [int(x) for x in np.ones_like(mz)]
tof = [int(x) for x in np.ones_like(mz)]
inv_ion_mob = np.ones_like(mz)

frame = TimsFrame(None, 1, scan, mz, intensity, tof, inv_ion_mob)


def test_frame_construction():
    assert set(frame.mz()) == set(mz) and set(frame.scan()) == set(scan) and set(frame.intensity()) == set(intensity)


def test_spectral_splitting():
    spec1, spec2 = frame.get_spectra()
    assert {500.123, 506.123, 507.127, 511.112} == set(spec1.mz()) \
           and {500.123, 506.123, 708.58, 892.123} == set(spec2.mz())


def test_window_block_generation():
    c, wi, si, bi, i, v = frame.get_hashing_blocks(min_peaks=1, min_intensity=1, overlapping=False)


def main():
    test_frame_construction()


if __name__ == '__main__':
    main()
