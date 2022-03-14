import tensorflow as tf
import numpy as np
import warnings
import libproteolizard as pl
from pyproteolizard.data import MzSpectrum, TimsFrame


def bins_to_mz(mz_bin, win_len):
    return np.abs(mz_bin) * win_len + (int(mz_bin < 0) * (0.5 * win_len))


class TimsHasher:
    def __init__(self, trials, len_trial, seed, resolution, num_dalton):

        if 0 < len_trial <= 32:
            self.V = tf.constant(
                np.expand_dims(np.array([np.power(2, i) for i in range(len_trial)]).astype(np.int32), 1))

        elif 32 < len_trial <= 64:
            warnings.warn(f"\nnum bits to hash set to: {len_trial}.\n" +
                          f"using int64 which might slow down computation significantly.")
            self.V = tf.constant(
                np.expand_dims(np.array([np.power(2, i) for i in range(len_trial)]).astype(np.int64), 1))

        else:
            raise ValueError(f"bit number per hash cannot be greater then 64 or smaller 1, was: {len_trial}.")

        self.trails = trials
        self.len_trial = len_trial

        self.seed = seed
        self.resolution = resolution
        self.num_dalton = num_dalton

        self.__hash_ptr = pl.TimsHashGenerator(len_trial, trials, seed, resolution, num_dalton)
        self.hash_matrix = self.get_matrix()
        self.hash_tensor = self.tf_tensor()

    def get_matrix(self):
        return self.__hash_ptr.getMatrixCopy()

    def tf_tensor(self):
        return tf.transpose(tf.constant(self.hash_matrix.astype(np.float32)))

    def calculate_keys(self, W):
        # generate signum
        S = (tf.sign(W @ self.hash_tensor) + 1) / 2

        if self.len_trial <= 32:
            # reshape into window, num_hashes, len_single_hash
            S = tf.cast(tf.reshape(S, shape=(S.shape[0], self.trails, self.len_trial)), tf.int32)

            # calculate int key from binary by base-transform
            H = tf.squeeze(S @ self.V)
            return H
        else:
            # reshape into window, num_hashes, len_single_hash
            S = tf.cast(tf.reshape(S, shape=(S.shape[0], self.trails, self.len_trial)), tf.int64)

            # calculate int key from binary by base-transform
            H = tf.squeeze(S @ self.V)
            return H

    """
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
    """