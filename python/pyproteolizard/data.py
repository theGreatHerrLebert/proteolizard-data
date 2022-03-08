import numpy as np
import pandas as pd
import sqlite3
import libproteolizard as pl
import opentims_bruker_bridge as obb
import tensorflow as tf


class PyTimsDataHandle:
    def __init__(self, dp):
        """
        construct a TimsDataHandle for simple fetching of data
        :param dp: data path to bruker TimsTOF experiment
        :param bp: binary path to bruker libtimsdata.* for hidden functionality call
        """
        self.dp: str = dp
        self.bp: str = obb.get_appropriate_so_path()

        self.precursor_frames: np.array = self.__get_precursor_frame_ids()
        self.fragment_frames: np.array = self.__get_fragment_frame_ids()
        self.meta_data = self.__get_meta_data()

        try:
            self.__handle = pl.ExposedTimsDataHandle(self.dp, self.bp)

        except Exception as e:
            print(e)

    def get_frame(self, frame_id):
        """
        get frame raw data
        :param frame_id: frame id of frame which should be fetched
        :return: a TimsFrame
        """
        frame_pointer = self.__handle.getFrame(frame_id)

        return TimsFrame(frame_pointer)

    def get_slice(self, precursor_ids, fragment_ids):
        return TimsSlice(self.__handle.getSlice(precursor_ids, fragment_ids))

    def get_slice_rt_range(self, rt_min, rt_max):
        prec_ids, frag_ids = self.__get_frame_ids_by_type_rt_range(rt_min, rt_max)
        return self.get_slice(prec_ids, frag_ids)

    # TODO: move this out of the data handle, not needed for timsTOF interface
    def get_selected_precursors(self):
        """
        :return: table of peaks chosen and fragmented (DDAExperiment) of all frames in experiment
        """
        return pd.read_sql_query("SELECT * from Precursors", sqlite3.connect(self.dp + "/analysis.tdf"))

    def frames_to_rts(self, frames: np.ndarray):
        d = dict(zip(self.meta_data.Id.values, self.meta_data.Time.values))
        return [d[x] for x in frames]

    def rt_range_to_precursor_frame_ids(self, rt_min: float, rt_max: float) -> np.array:
        """
        :param rt_min: start in rt dimension in SECONDS
        :param rt_max: stop in rt dimension in SECONDS
        :return:
        """
        # extract rt frame region
        region = self.meta_data[(self.meta_data['Time'] >= rt_min) & (self.meta_data['Time'] <= rt_max)]
        return region[region['MsMsType'] == 0].Id.values

    def rt_range_to_fragment_frame_ids(self, rt_min: float, rt_max: float) -> np.array:
        """
        :param rt_min: start in rt dimension in SECONDS
        :param rt_max: stop in rt dimension in SECONDS
        :return:
        """
        # extract rt frame region
        region = self.meta_data[(self.meta_data['Time'] >= rt_min) & (self.meta_data['Time'] <= rt_max)]
        return region[region['MsMsType'] != 0].Id.values
            
    def __get_precursor_frame_ids(self):
        """
        :return: array of all precursor frame ids
        """
        return pd.read_sql_query("SELECT * from Frames WHERE MsMsType = 0",
                                 sqlite3.connect(self.dp + "/analysis.tdf")).Id.values

    def __get_fragment_frame_ids(self):
        """
        :return: array of all fragment frame ids
        """
        return pd.read_sql_query("SELECT * from Frames WHERE MsMsType != 0",
                                 sqlite3.connect(self.dp + "/analysis.tdf")).Id.values

    def __get_meta_data(self):
        """
        :return: table of complete meta data of all frames in experiment
        """
        return pd.read_sql_query("SELECT * FROM Frames", sqlite3.connect(self.dp + "/analysis.tdf"))

    def __get_frame_ids_by_type_rt_range(self, rt_start, rt_stop):
        frames = self.meta_data[(rt_start <= self.meta_data.Time) & (self.meta_data.Time <= rt_stop)]
        return frames[frames.MsMsType == 0].Id.values, frames[frames.MsMsType != 0].Id.values


class MzSpectrum:

    def __init__(self, spec_pointer, *args):

        if len(args) > 0:
            frame, scan, mz, intensity = args
            self.spec_ptr = pl.MzSpectrumPL(frame, scan, mz, intensity)

        else:
            self.spec_ptr = spec_pointer

    def frame_id(self):
        return self.spec_ptr.getFrameId()

    def scan_id(self):
        return self.spec_ptr.getScanId()

    def mz(self):
        return self.spec_ptr.getMzs()

    def intensity(self):
        return self.spec_ptr.getIntensities()

    def __add__(self, other):
        return MzSpectrum(self.spec_ptr + other.spec_ptr)

    def to_resolution(self, resolution: int = 2):
        """
        :param resolution:
        :return:
        """
        return MzVector(self.spec_ptr.vectorize(resolution))

    def windows(self, window_length=10, overlapping=True, min_peaks=3, min_intensity=50):
        bins, windows = self.spec_ptr.windows(window_length, overlapping, min_peaks, min_intensity)
        return bins, [MzSpectrum(s) for s in windows]


class MzVector:
    def __init__(self, vec_pointer):
        self.__vec = vec_pointer

    def frame_id(self):
        return self.__vec.getFrameId()

    def scan_id(self):
        return self.__vec.getScanId()

    def resolution(self):
        return self.__vec.getResolution()

    def indices(self):
        return self.__vec.getIndices()

    def values(self):
        return self.__vec.getValues()

    def __add__(self, other):
        return MzVector(self.__vec + other.__vec)


class VectorizedTimsFrame:
    def __init__(self, vec_frame_pointer):
        self.__frame = vec_frame_pointer

    def frame_id(self):
        return self.__frame.getFrameId()

    def resolution(self):
        return self.__frame.getResolution()

    def scans(self):
        return self.__frame.getScans()

    def values(self):
        return self.__frame.getValues()

    def indices(self):
        return self.__frame.getIndices()

    def filter_ranged(self, scan_min, scan_max, index_min, index_max, intensity_min):
        return VectorizedTimsFrame(self.__frame.filterRanged(scan_min, scan_max, index_min, index_max, intensity_min))

    def get_spectra(self):
        return [MzVector(s) for s in self.__frame.getSpectra()]

    def __add__(self, other):
        return VectorizedTimsFrame(self.__frame + other.frame_ptr)


class TimsFrame:
    def __init__(self, frame_pointer, *args):

        if len(args) > 0:
            frame_id, scan, mz, intensity, tof, inv_ion_mob = args
            self.frame_ptr = pl.TimsFrame(frame_id, scan, mz, intensity, tof, inv_ion_mob)

        else:
            self.frame_ptr = frame_pointer

    def frame_id(self):
        return self.frame_ptr.getFrameId()

    def scan(self):
        return self.frame_ptr.getScans()

    def mz(self):
        return self.frame_ptr.getMzs()

    def intensity(self):
        return self.frame_ptr.getIntensities()

    def tof(self):
        return self.frame_ptr.getTofs()

    def inverse_ion_mobility(self):
        return self.frame_ptr.getInverseMobilities()

    def __add__(self, other):
        return TimsFrame(self.frame_ptr + other.frame_ptr)

    def to_resolution(self, resolution=2):
        """
        :param resolution:
        :return:
        """
        return TimsFrame(self.frame_ptr.toResolution(resolution))

    def vectorized(self, resolution=2):
        """
        :param resolution:
        :return:
        """
        return VectorizedTimsFrame(self.frame_ptr.vectorize(resolution))

    def filter_ranged(self, scan_min=int(0.0), scan_max=int(1e3), mz_min=0.0, mz_max=2e3, intensity_min=50):
        return TimsFrame(self.frame_ptr.filterRanged(scan_min, scan_max, mz_min, mz_max, intensity_min))

    def fold(self, resolution=2, width=4):
        return TimsFrame(self.frame_ptr.fold(resolution, width))

    def get_spectra(self):
        spec_ptrs = self.frame_ptr.getMzSpectra()
        return [MzSpectrum(ptr) for ptr in spec_ptrs]

    def dense_windows(self, resolution=1, min_peaks=5, min_intensity=100, window_length=10, overlapping=True):
        return self.frame_ptr.getDenseWindows(resolution, min_peaks, min_intensity, window_length, overlapping)

    def get_hashing_blocks(self, resolution=1, min_peaks=3, min_intensity=50, window_length=10, overlapping=True):
        return self.frame_ptr.getHashingBlocks(resolution, min_peaks, min_intensity, window_length, overlapping)

    def hashing_block_as_dense_tensor(self, resolution=1, min_peaks=3, min_intensity=50, window_length=10,
                                      overlapping=True):

        c, wi, scan, mz_bin, i, v = self.frame_ptr.getHashingBlocks(resolution, min_peaks, min_intensity, window_length,
                                                                    overlapping)

        len_mz_vector = int(np.power(10, resolution) * window_length)

        st = tf.sparse.SparseTensor(indices=np.c_[wi, i], values=v.astype(np.float32), dense_shape=(c, len_mz_vector + 1))
        return tf.transpose(tf.sparse.to_dense(st))


class TimsSlice:
    def __init__(self, slice_ptr):
        self.__slice_ptr = slice_ptr

    def get_precursor_frames(self):
        return [TimsFrame(x) for x in self.__slice_ptr.getPrecursors()]

    def get_fragment_frames(self):
        return [TimsFrame(x) for x in self.__slice_ptr.getFragments()]

    def get_precursor_points(self):
        return Points3D(self.__slice_ptr.getPoints(True)).get_points()

    def get_fragment_points(self):
        return Points3D(self.__slice_ptr.getPoints(False)).get_points()

    def filter_ranged(self, mz_min, mz_max, scan_min=0, scan_max=1000, intensity_min=0):
        return TimsSlice(self.__slice_ptr.filterRanged(scan_min, scan_max, mz_min, mz_max, intensity_min))


class Points3D:
    def __init__(self, point_ptr):
        self.__point_ptr = point_ptr

    def get_points(self):
        ids = self.__point_ptr.getFrames()
        scans = self.__point_ptr.getScans()
        mz = self.__point_ptr.getMz()
        inv_ion_mob = self.__point_ptr.getInvIonMobility()
        intensity = self.__point_ptr.getIntensity()

        return pd.DataFrame({'frame': ids, 'scan': scans, 'invIonMob': inv_ion_mob, 'mz': mz, 'intensity': intensity})
