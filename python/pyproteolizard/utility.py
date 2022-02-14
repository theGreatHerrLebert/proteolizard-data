import numpy as np


def peak_width_preserving_mz_transform(
        mz: np.array,
        M0: float = 500,
        resolution: float = 50_000):
    """
    Transform values into an index that fixes the width of the peak at half full width.
    Arguments:
        mz (np.array): An array of mass to charge ratios to transform.
        M0 (float): The m/z value at which TOFs resolution is reported.
        resolution (float): The resolution of the TOF instrument.
    """
    return (np.log(mz) - np.log(M0)) / np.log1p(1 / resolution)