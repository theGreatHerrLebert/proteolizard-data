import numpy as np
import pandas as pd
import sqlite3
import libproteolizard as pl
import opentims_bruker_bridge as obb


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


class MzSpectrum:
    def __init__(self, spec_pointer):
        self.__spec = spec_pointer

    def frame_id(self):
        return self.__spec.getFrameId()

    def scan_id(self):
        return self.__spec.getScanId()

    def mz(self):
        return self.__spec.getMzs()

    def intensity(self):
        return self.__spec.getIntensities()

    def __add__(self, other):
        return MzSpectrum(self.__spec + other.__spec)

    def to_resolution(self, resolution: int = 2):
        """
        :param resolution:
        :return:
        """
        return self.__spec.toResolution(resolution)


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
        return VectorizedTimsFrame(self.__frame + other.__frame)


class TimsFrame:
    def __init__(self, frame_pointer):
        self.__frame = frame_pointer

    def frame_id(self):
        return self.__frame.getFrameId()

    def scan(self):
        return self.__frame.getScans()

    def mz(self):
        return self.__frame.getMzs()

    def intensity(self):
        return self.__frame.getIntensities()

    def tof(self):
        return self.__frame.getTofs()

    def invers_ion_mobility(self):
        return self.__frame.getInverseMobilities()

    def __add__(self, other):
        return TimsFrame(self.__frame + other.__frame)

    def to_resolution(self, resolution=2):
        """
        :param resolution:
        :return:
        """
        return TimsFrame(self.__frame.toResolution(resolution))

    def vectorized(self, resolution=2):
        """
        :param resolution:
        :return:
        """
        return VectorizedTimsFrame(self.__frame.vectorize(resolution))

    def filter_ranged(self, scan_min, scan_max, mz_min, mz_max):
        return TimsFrame(self.__frame.filterRanged(scan_min, scan_max, mz_min, mz_max))

    def fold(self, resolution=2, width=4):
        return TimsFrame(self.__frame.fold(resolution, width))

    def get_spectra(self):
        spec_ptrs = self.__frame.getMzSpectra()
        return [MzSpectrum(ptr) for ptr in spec_ptrs]
