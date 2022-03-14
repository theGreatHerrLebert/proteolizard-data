import tensorflow as tf
import numpy as np
import warnings
import libproteolizard as pl
from pyproteolizard.data import MzSpectrum, TimsFrame


def bins_to_mz(mz_bin, win_len):
    return np.abs(mz_bin) * win_len + (int(mz_bin < 0) * (0.5 * win_len))


class TimsHasher:
    """
    :param trials: number of keys to be generated
    :param len_trial: bits per key, aka 0 or 1 for every signum created via random projection
    :param seed: a seed to be able to recreate random hashes for a given parameter setting
    :param resolution: decimals kept when binning the mz axis of a given mass spectrum
    :param num_dalton: length of a single window a spectrum is broken into before hashing
    """
    def __init__(self, trials, len_trial, seed=5671, resolution=1, num_dalton=10):

        # check
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

    # fetch matrix from C++ implementation
    def get_matrix(self):
        return self.__hash_ptr.getMatrixCopy()

    # create TF tensor (potentially GPU located) for hashing
    def tf_tensor(self):
        return tf.transpose(tf.constant(self.hash_matrix.astype(np.float32)))

    # create keys by random projection
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