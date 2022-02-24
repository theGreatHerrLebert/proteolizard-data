import libproteolizard as pl
# from data import MzSpectrum


class TimsHasher:
    def __init__(self, trials, len_trial, seed, resolution):
        self.__hash_ptr = pl.TimsHashGenerator(trials, len_trial, seed, resolution)
        self.trails = trials
        self.len_trial = len_trial
        self.seed = seed
        self.resolution = resolution

    def getMatrix(self):
        return self.__hash_ptr.getMatrixCopy()

    def hash_spectrum(self, spectrum,
                      min_peaks: int = 5,
                      min_intensity: int = 150,
                      window_length = 10,
                      overlapping: bool = False,
                      restrict: bool = True):
        return self.__hash_ptr.hashMzSpectrum(spectrum.spec_ptr, min_peaks, min_intensity, window_length, overlapping, restrict)
