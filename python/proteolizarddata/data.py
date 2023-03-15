import tensorflow as tf
import numpy as np
import pandas as pd
import sqlite3
import libproteolizarddata as pl
import opentims_bruker_bridge as obb
from abc import ABC, abstractmethod
from typing import List
import json

class PyTimsDataHandle(ABC):
    def __init__(self, dp):
        """
        construct a TimsDataHandle for simple fetching of data
        :param dp: data path to bruker TimsTOF experiment
        :param bp: binary path to bruker libtimsdata.* for hidden functionality call
        """
        self.dp: str = dp
        self.bp: List[str] = obb.get_so_paths()

        self.precursor_frames: np.array = self.__get_precursor_frame_ids()
        self.fragment_frames: np.array = self.__get_fragment_frame_ids()
        self.meta_data = self.__get_meta_data()
        self.pasef_meta_data = self._get_pasef_meta_data()

        appropriate_found = False

        for so_path in self.bp:
            try:
                self.__handle = pl.ExposedTimsDataHandle(self.dp, so_path)
                appropriate_found = True
                break
            except Exception as e:
                continue
        assert appropriate_found is True, "No appropriate opentims_bruker_bridge could be found for your system"

    @property
    @abstractmethod
    def acquisition(self) -> int:
        """
        Gets acquistion as integer:
          8: ddaPasef
          9: diaPasef

        :return int: acquisition method
        """
        pass

    def get_frame(self, frame_id):
        """
        get frame raw data
        :param frame_id: frame id of frame which should be fetched
        :return: a TimsFrame
        """
        frame_pointer = self.__handle.getFrame(frame_id)

        return TimsFrame(frame_pointer)

    def get_slice(self, precursor_ids, fragment_ids):
        """

        :param precursor_ids:
        :param fragment_ids:
        :return:
        """
        return TimsSlice(self.__handle.getSlice(precursor_ids, fragment_ids))

    def get_block(self, frame_ids):
        return TimsBlock(self.__handle.getBlock(frame_ids))

    def get_slice_rt_range(self, rt_min, rt_max):
        """

        :param rt_min:
        :param rt_max:
        :return:
        """
        prec_ids, frag_ids = self.__get_frame_ids_by_type_rt_range(rt_min, rt_max)
        return self.get_slice(prec_ids, frag_ids)

    def frames_to_rts(self, frames: np.ndarray):
        """

        :param frames:
        :return:
        """
        d = dict(zip(self.meta_data.Id.values, self.meta_data.Time.values))
        return np.array([d[x] for x in frames])

    def rt_to_closest_frame_id(self, rt, precursor=True):
        """

        :param precursor:
        :param rt:
        :return:
        """
        if precursor:
            return self.precursor_frames[np.abs(self.frames_to_rts(self.precursor_frames) - rt).argmin()]
        return self.fragment_frames[np.abs(self.frames_to_rts(self.fragment_frames) - rt).argmin()]

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
        return region[region['MsMsType'] == self.acquisition].Id.values

    def get_global_mz_axis(self):
        """

        :return:
        """
        return self.__handle.getGlobalMzAxis()

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
        return pd.read_sql_query(f"SELECT * from Frames WHERE MsMsType == {self.acquisition}",
                                 sqlite3.connect(self.dp + "/analysis.tdf")).Id.values

    def __get_meta_data(self):
        """
        :return: table of complete meta data of all frames in experiment
        """
        return pd.read_sql_query("SELECT * FROM Frames", sqlite3.connect(self.dp + "/analysis.tdf"))

    @abstractmethod
    def _get_pasef_meta_data(self):
        """
        :return: table of pasef meta data
        """
        pass

    def __get_frame_ids_by_type_rt_range(self, rt_start, rt_stop):
        """

        :param rt_start:
        :param rt_stop:
        :return:
        """
        frames = self.meta_data[(rt_start <= self.meta_data.Time) & (self.meta_data.Time <= rt_stop)]
        return frames[frames.MsMsType == 0].Id.values, frames[frames.MsMsType == self.acquisition].Id.values


class PyTimsDataHandleDDA(PyTimsDataHandle):

    @property
    def acquisition(self) -> int:
        """
        gets acquisition as integer (8, DDA)
        :return (int): 8 (DDA acquisition)
        """
        return 8

    def get_selected_precursors(self):
        """
        :return: table of peaks chosen and fragmented (DDAExperiment) of all frames in experiment
        """
        return pd.read_sql_query("SELECT * from Precursors", sqlite3.connect(self.dp + "/analysis.tdf"))

    def get_selected_precursor_by_id(self, precursor_id:int) -> pd.DataFrame:
        """
        Get data of precursor by its id in precursors table.

        :param precursor_id (int): ID of precursor to get in Precursors table.
        :return: DataFrame with precursor data
        """
        # replace makes sure NULL values are np.nan
        return pd.read_sql_query(f"SELECT * from Precursors where Id={precursor_id}",
                                sqlite3.connect(self.dp + "/analysis.tdf")).replace([None],[np.nan])

    def _get_pasef_meta_data(self):
        """
        :return: table of pasef meta data
        """
        return pd.read_sql_query("SELECT * from PasefFrameMsMsInfo",
                                 sqlite3.connect(self.dp + "/analysis.tdf"))


class PyTimsDataHandleDIA(PyTimsDataHandle):

    @property
    def acquisition(self) -> int:
        """
        gets acquisition as integer (9, DIA)
        :return int: 9 (DIA acquisition)
        """
        return 9

    def _get_pasef_meta_data(self):
        """
        :return: table of pasef meta data
        """
        return pd.read_sql_query("SELECT * from DiaFrameMsMsWindows",
                                 sqlite3.connect(self.dp + "/analysis.tdf"))


class MzSpectrum:

    def __init__(self, spec_pointer, *args):

        if len(args) > 0:
            frame, scan, mz, intensity = args
            self.spec_ptr = pl.MzSpectrumPL(frame, scan, mz, intensity)

        else:
            self.spec_ptr = spec_pointer

    def __repr__(self):
        return f"MzSpectrum(frame: {self.frame_id()}, scan: {self.scan_id()}, sum intensity: {self.sum_intensity()})"

    def frame_id(self):
        """
        :return:
        """
        return self.spec_ptr.getFrameId()

    def scan_id(self):
        """
        :return:
        """
        return self.spec_ptr.getScanId()

    def mz(self):
        """
        :return:
        """
        return self.spec_ptr.getMzs()

    def intensity(self):
        """
        :return:
        """
        return self.spec_ptr.getIntensities()

    def __add__(self, other):
        """
        :param other:
        :return:
        """
        return MzSpectrum(self.spec_ptr + other.spec_ptr)

    def vectorize(self, resolution: int = 2):
        """
        :param resolution:
        :return:
        """
        return MzVector(self.spec_ptr.vectorize(resolution))

    def to_resolution(self, resolution: int = 2):
        """
        :param resolution:
        :return:
        """
        return MzSpectrum(self.spec_ptr.toResolution(resolution))

    def to_centroided(self, baseline_noise_level: int, sigma: float):
        """
        :param sigma: distance cutoff value, maximal distance for two mzs
            to be considered of same peak.
        """
        return MzSpectrum(self.spec_ptr.toCentroided(baseline_noise_level, sigma))

    def windows(self, window_length=10, overlapping=True, min_peaks=3, min_intensity=50):
        """
        :param window_length:
        :param overlapping:
        :param min_peaks:
        :param min_intensity:
        :return:
        """
        bins, windows = self.spec_ptr.windows(window_length, overlapping, min_peaks, min_intensity)
        return bins, [MzSpectrum(s) for s in windows]

    def sum_intensity(self):
        """
        :return:
        """
        return np.sum(self.intensity())

    def filter(self, mz_min=0, mz_max=2000, intensity_min=25):
        """
        to ignore one limit, it can be set to -1
        :param mz_min:
        :param mz_max:
        :param intensity_min:
        :return:
        """
        return MzSpectrum(self.spec_ptr.filter(mz_min, mz_max, intensity_min))

    def __mul__(self, scalar: float):
        """
        :param scalar: float to scale intensities by.
        """
        return MzSpectrum(self.spec_ptr * scalar)

    def __rmul__(self, scalar: float):
        """
        :param scalar: float to scale intensities by.
        """
        return MzSpectrum(scalar*self.spec_ptr)

    def to_jsons(self, only_spectrum = False):
        """
        generates json string representation of MzSpectrum
        """
        json_dict = {}
        json_dict["mz"] = self.mz().tolist()
        json_dict["intensity"] = self.intensity().tolist()

        if only_spectrum:
            return json.dumps(json_dict)

        json_dict["frame_id"] = self.frame_id()
        json_dict["scan_id"] = self.scan_id()

        return json.dumps(json_dict)

class MzVector:
    def __init__(self, vec_pointer):
        self.__vec = vec_pointer

    def frame_id(self):
        """
        :return:
        """
        return self.__vec.getFrameId()

    def scan_id(self):
        """
        :return:
        """
        return self.__vec.getScanId()

    def resolution(self):
        """
        :return:
        """
        return self.__vec.getResolution()

    def indices(self):
        """
        :return:
        """
        return self.__vec.getIndices()

    def values(self):
        """
        :return:
        """
        return self.__vec.getValues()

    def __add__(self, other):
        """
        :param other:
        :return:
        """
        return MzVector(self.__vec + other.__vec)

    def __repr__(self):
        return f"MzVector(frame: {self.frame_id()}, scan: {self.scan_id()}, sum intensity: {np.sum(self.values())})"


class VectorizedTimsFrame:
    def __init__(self, vec_frame_pointer):
        self.__frame = vec_frame_pointer

    def frame_id(self):
        """
        :return:
        """
        return self.__frame.getFrameId()

    def resolution(self):
        """
        :return:
        """
        return self.__frame.getResolution()

    def scans(self):
        """
        :return:
        """
        return self.__frame.getScans()

    def values(self):
        """
        :return:
        """
        return self.__frame.getValues()

    def indices(self):
        """
        :return:
        """
        return self.__frame.getIndices()

    def filter_ranged(self, scan_min, scan_max, index_min, index_max, intensity_min):
        """
        :param scan_min:
        :param scan_max:
        :param index_min:
        :param index_max:
        :param intensity_min:
        :return:
        """
        return VectorizedTimsFrame(self.__frame.filterRanged(scan_min, scan_max, index_min, index_max, intensity_min))

    def __repr__(self):
        return f"TimsFrameVectorized(id: {self.frame_id()}, num data-points: {self.indices().shape[0]}, sum intensity: " \
               f"{np.sum(self.values())})"

    def get_spectra(self):
        """
        :return:
        """
        return [MzVector(s) for s in self.__frame.getSpectra()]

    def __add__(self, other):
        return VectorizedTimsFrame(self.__frame + other.frame_ptr)

    def data(self):
        return pd.DataFrame({'frame': np.repeat(self.frame_id(), self.indices().shape[0]),
                             'scan': self.scans(),
                             'indices': self.indices(),
                             'values': self.values()})


class TimsFrame:
    def __init__(self, frame_pointer, *args):

        if len(args) > 0:
            frame_id, rt, scan, mz, intensity, tof, inv_ion_mob = args
            self.frame_ptr = pl.TimsFrame(frame_id, rt, scan, mz, intensity, tof, inv_ion_mob)

        else:
            self.frame_ptr = frame_pointer

    def __repr__(self):
        return f"TimsFrame(id: {self.frame_id()}, num data-points: {self.mz().shape[0]}, sum intensity: " \
               f"{np.sum(self.intensity())})"

    def frame_id(self):
        """
        :return:
        """
        return self.frame_ptr.getFrameId()

    def retention_time(self):
        """
        :return:
        """
        return self.frame_ptr.getRetentionTime()

    def scan(self):
        """
        :return:
        """
        return self.frame_ptr.getScans()

    def mz(self):
        """
        :return:
        """
        return self.frame_ptr.getMzs()

    def intensity(self):
        """
        :return:
        """
        return self.frame_ptr.getIntensities()

    def tof(self):
        """
        :return:
        """
        return self.frame_ptr.getTofs()

    def inverse_ion_mobility(self):
        """
        :return:
        """
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
        """
        reduce number of scans in a given frame by summing over consecutive scans
        :param resolution: binning in mz-dimension
        :param width: number of consecutive scans to be integrated
        :return: an integrated TimsFrame
        """
        return TimsFrame(self.frame_ptr.fold(resolution, width))

    def get_spectra(self):
        spec_ptrs = self.frame_ptr.getMzSpectra()
        return [MzSpectrum(ptr) for ptr in spec_ptrs]

    def get_dense_windows(self, resolution: int = 1, min_peaks: int = 3, min_intensity: int = 50,
                          window_length: int = 10, overlapping: bool = True):
        """
        create a collection of dense windows that represent slices of mass spectra
        :param resolution: binning of mz axis, will be 10^-resolution
        :param min_peaks: minimum number of peaks in a window
        :param min_intensity: minimum intensity of the highest peak in a window to be returned
        :param window_length: length of a window along mz-axis in Dalton
        :param overlapping: if true, overlapping windows will be created where windows are shifted by 1/2 window-length
        :return: 1D scan-indices, 1D bin-indices, 2D dense vectors of windowed mass spectra in frame
        """

        c, wi, s, b, i, v = self.frame_ptr.getHashingBlocks(resolution, min_peaks, min_intensity, window_length,
                                                            overlapping)

        len_mz_vector = int(np.power(10, resolution) * window_length)

        st = tf.sparse.SparseTensor(indices=np.c_[wi, i], values=v.astype(np.float32),
                                    dense_shape=(c, len_mz_vector + 1))

        return s.astype(np.int32), b.astype(np.int32), tf.sparse.to_dense(st)

    def get_vectorized_windows(self, resolution: int = 1, min_peaks: int = 3, min_intensity: int = 50,
                               window_length: int = 10, overlapping: bool = True):

        scans, bins, windows = self.frame_ptr.getVectorizedWindows(resolution, min_peaks, min_intensity,
                                                                   window_length, overlapping)

        return scans, bins, [MzVector(x) for x in windows]

    def data(self):
        return pd.DataFrame({'frame': np.repeat(self.frame_id(), self.mz().shape[0]),
                             'scan': self.scan(),
                             'inv_ion_mob': self.inverse_ion_mobility(),
                             'mz': self.mz(),
                             'tof': self.tof(),
                             'intensity': self.intensity()})


class TimsSlice:
    def __init__(self, slice_ptr, *args):

        if len(args) > 0:
            precursors, fragments = args
            self.__slice_ptr = pl.TimsSlicePL([f.frame_ptr for f in precursors],
                                              [f.frame_ptr for f in fragments])

        else:
            self.__slice_ptr = slice_ptr

    def get_precursor_frames(self):
        """
        :return:
        """
        return [TimsFrame(x) for x in self.__slice_ptr.getPrecursors()]

    def get_fragment_frames(self):
        """
        :return:
        """
        return [TimsFrame(x) for x in self.__slice_ptr.getFragments()]

    def get_precursor_points(self):
        """
        :return:
        """
        return Points3D(self.__slice_ptr.getPoints(True)).get_points()

    def get_fragment_points(self):
        """
        :return:
        """
        return Points3D(self.__slice_ptr.getPoints(False)).get_points()

    def filter_ranged(self, mz_min=0, mz_max=2000, scan_min=0, scan_max=1000, intensity_min=0, rt_min=0, rt_max=0):
        """
        :param mz_min:
        :param mz_max:
        :param scan_min:
        :param scan_max:
        :param intensity_min:
        :param rt_min:
        :param rt_max:
        :return:
        """

        if rt_min == 0 and rt_max == 0:

            frames = self.get_precursor_frames()
            s = frames[0].retention_time()
            e = frames[-1].retention_time()
            return TimsSlice(
                self.__slice_ptr.filterRanged(scan_min, scan_max, mz_min, mz_max, intensity_min, s, e))

        return TimsSlice(
            self.__slice_ptr.filterRanged(scan_min, scan_max, mz_min, mz_max, intensity_min, rt_min, rt_max))

    def __repr__(self):
        """
        :return:
        """
        frames = self.get_precursor_frames()

        if len(frames) == 0:
            return f"TimsSlice(None)"
        return f"TimsSlice(frame start: {frames[0].frame_id()}, frame end: {frames[len(frames)-1].frame_id()})"

    def vectorize(self, resolution=2):
        """
        :param resolution:
        :return:
        """
        return TimsSliceVectorized(self.__slice_ptr.getVectorizedSlice(resolution))


class TimsSliceVectorized:
    def __init__(self, slice_ptr):
        self.__slice_ptr = slice_ptr

    def get_precursor_frames(self):
        """
        :return:
        """
        return [VectorizedTimsFrame(x) for x in self.__slice_ptr.getPrecursors()]

    def get_fragment_frames(self):
        """
        :return:
        """
        return [VectorizedTimsFrame(x) for x in self.__slice_ptr.getFragments()]

    def get_precursor_points(self):
        """
        :return:
        """
        return Point3DVectorized(self.__slice_ptr.getPoints(True)).get_points()

    def get_fragment_points(self):
        """
        :return:
        """
        return Point3DVectorized(self.__slice_ptr.getPoints(False)).get_points()

    def filter_ranged(self, mz_min=0, mz_max=2000, scan_min=0, scan_max=1000, intensity_min=0):
        """
        :param mz_min:
        :param mz_max:
        :param scan_min:
        :param scan_max:
        :param intensity_min:
        :return:
        """
        return TimsSliceVectorized(self.__slice_ptr.filterRanged(scan_min, scan_max, mz_min, mz_max, intensity_min))

    def __repr__(self):
        """
        :return:
        """
        frames = self.get_precursor_frames()
        return f"TimsSliceVectorized(frame start: {frames[0].frame_id()}, frame end: {frames[len(frames)-1].frame_id()})"

    def get_zero_indexed_sparse_tensor(self, precursor=True, zero_index_mz=True):
        """
        translate vectorized point indices into sparse tensor representation
        :param precursor: if true, return a sparse tensor of precursor frames, otherwise fragments
        :param zero_index_mz: if true, mz indices will be zero-indexed as well
        :return: a sparse tensor of either precursor or fragment frames
        """

        if precursor:
            data = self.get_precursor_points()
        else:
            data = self.get_fragment_points()

        # need to translate frame ids to zero-indexed
        frames = data['frame'].values
        unique_frames = np.sort(np.unique(frames))
        f_dict = dict(np.c_[unique_frames, np.arange(unique_frames.shape[0])])

        # need to translate scan ids to zero-indexed
        scans = data['scan'].values
        unique_scans = np.sort(np.unique(scans))
        s_dict = dict(np.c_[unique_scans, np.arange(unique_scans.shape[0])])

        # potentially remove mz offset
        if zero_index_mz:
            indices = data.indices.values - np.min(data.indices.values)
        else:
            indices = data.indices.values

        def __fun_f(frame_id):
            return f_dict[frame_id]

        def __fun_s(scan_id):
            return s_dict[scan_id]

        f_to_index = np.vectorize(__fun_f)
        s_to_index = np.vectorize(__fun_s)

        indexed_frames = f_to_index(frames)
        indexed_scans = s_to_index(scans)

        idx = np.c_[indexed_frames, indexed_scans, indices]

        st = tf.sparse.SparseTensor(indices=idx, values=data['values'].astype(np.float32), dense_shape=np.max(idx, axis=0) + 1)

        return st, unique_frames[0], unique_frames[-1], unique_scans[0], indexed_scans[-1]


class Points3D:
    def __init__(self, point_ptr):
        self.__point_ptr = point_ptr

    def get_points(self):
        ids = self.__point_ptr.getFrames()
        scans = self.__point_ptr.getScans()
        mz = self.__point_ptr.getMz()
        inv_ion_mob = self.__point_ptr.getInvIonMobility()
        intensity = self.__point_ptr.getIntensity()

        return pd.DataFrame({'frame': ids, 'scan': scans, 'inv_ion_mob': inv_ion_mob, 'mz': mz, 'intensity': intensity})


class Point3DVectorized:
    def __init__(self, point_ptr):
        self.__point_ptr = point_ptr

    def get_points(self):
        ids = self.__point_ptr.getFrames()
        scans = self.__point_ptr.getScans()
        mz = self.__point_ptr.getMz()
        intensity = self.__point_ptr.getIntensity()
        return pd.DataFrame({'frame': ids, 'scan': scans, 'indices': mz, 'values': intensity})


class TimsBlock:
    def __init__(self, block_ptr):
        self.__block_ptr = block_ptr

    def get_indices(self):
        return self.__block_ptr.getIndices().T

    def get_values(self):
        return self.__block_ptr.getValues().T

    def filter_ranged(self, mz_min=0, mz_max=2000, scan_min=0, scan_max=1000, intensity_min=0):
        """
        :param mz_min:
        :param mz_max:
        :param scan_min:
        :param scan_max:
        :param intensity_min:
        :return:
        """
        return TimsBlock(self.__block_ptr.filterRanged(scan_min, scan_max, mz_min, mz_max, intensity_min))

    def vectorized(self, resolution: int = 2):
        """
        :param resolution:
        :return:
        """
        return TimsBlockVectorized(self.__block_ptr.getTimsBlockVectorized(resolution))


class TimsBlockVectorized:
    def __init__(self, block_ptr):
        self.__block_ptr = block_ptr

    def data(self):
        return self.__block_ptr.getData().T

    def get_zero_indexed_tensor(self):
        data = self.__block_ptr.getData().T
        rt_indices = data[:, 0]
        sc_indices = data[:, 1] - np.min(data[:, 1])
        mz_indices = data[:, 2] - np.min(data[:, 2])
        values = data[:, 3].astype(np.float32)

        return tf.sparse.SparseTensor(indices=np.c_[rt_indices, sc_indices, mz_indices], values=values,
                                      dense_shape=(rt_indices.max()+1, sc_indices.max()+1, mz_indices.max()+1))

    def filter_ranged(self, mz_min=0, mz_max=2000, scan_min=0, scan_max=1000, intensity_min=0):
        """
        :param mz_min:
        :param mz_max:
        :param scan_min:
        :param scan_max:
        :param intensity_min:
        :return:
        """
        return TimsBlockVectorized(self.__block_ptr.filterRanged(scan_min, scan_max, mz_min, mz_max, intensity_min))
