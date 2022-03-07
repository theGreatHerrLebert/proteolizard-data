import numpy as np
import libproteolizard as pl
from pyproteolizard.data import MzSpectrum, TimsFrame


def bins_to_mz(mz_bin, win_len):
    return np.abs(mz_bin) * win_len + (int(mz_bin < 0) * (0.5 * win_len))


class TimsHasher:
    def __init__(self, trials, len_trial, seed, resolution, num_dalton):
        self.__hash_ptr = pl.TimsHashGenerator(trials, len_trial, seed, resolution, num_dalton)
        self.trails = trials
        self.len_trial = len_trial
        self.seed = seed
        self.resolution = resolution
        self.num_dalton = num_dalton

    def get_matrix(self):
        return self.__hash_ptr.getMatrixCopy()

    def hash_entire_spectrum(self, spectrum: MzSpectrum):
        return self.__hash_ptr.hashMzSpectrum(spectrum.spec_ptr)

    def hash_spectrum(self, spectrum: MzSpectrum,
                      min_peaks: int = 5,
                      min_intensity: int = 150,
                      window_length: float = 10,
                      overlapping: bool = False,
                      restrict: bool = True):
        bins, hashes = self.__hash_ptr.hashMzSpectrumWindows(spectrum.spec_ptr, min_peaks, min_intensity,
                                                             window_length, overlapping, restrict)

        return np.array([bins_to_mz(b, window_length) for b in bins]), hashes

    def hash_frame(self, frame: TimsFrame,
                   min_peaks: int = 5,
                   min_intensity: int = 150,
                   window_length: float = 10,
                   overlapping: bool = False,
                   restrict: bool = True):

        scans, bins, hashes = self.__hash_ptr.hashTimsFrameWindows(frame.frame_ptr, min_peaks, min_intensity,
                                                                   window_length, overlapping, restrict)

        return scans, np.array([bins_to_mz(b, window_length) for b in bins]), hashes
