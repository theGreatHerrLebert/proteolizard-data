# proteolizard-data
### An object-oriented library of C++ classes and python-wrappers to seamlessly integrate timsTOF raw-data with python

## Context
Ion-mobility enhanced tandem-MS coupled to liquid chromatography is rapidly becoming the method of 
choice for analysis of high-complexity samples generated in proteomics, lipidomics and metabolomics.
While information gain stemming from this acquisition scheme is tremendous, an additional recoding of 
ion mobility makes data processing more challenging. The extra-dimension puts a lot of pressure on the 
programmer to find solutions to raw-data processing that are fast enough to be useful.

The timsTOF mass-spectrometer from bruker daltonics is one of the latest instruments to implement
ion-mobility MS, creating raw-data files in brukers custom TDF format. 
The open-source [`opentims`](https://github.com/michalsta/opentims) project was created to provide 
easy programmatic access to such datasets, allowing the exposure of raw-data to languages including C++, 
python, R and the JVM.

## What is proteolizard-data?
`proteolizard-data` builds on top of the low-level C++ API of `opentims`, adding a thin layer of C++ classes that 
are exposed to python via bindings created with the excellent [`pybind11`](https://github.com/pybind/pybind11) 
library. The approach taken is different to that of [`timspy`](https://github.com/MatteoLacki/timspy) 
or [`alphatims`](https://github.com/MannLabs/alphatims), focussing on convenient, easy-to-use
classes that implement fast processing on dataset slice, frame or spectrum basis. It is built to 
quickly explore timsTOF datasets or prototype ideas and can be easily integrated into the well 
established, python-centric data science stack e.g. `scikit-learn`, `tensorflow` or `pytorch`.

## Navigation
* [**Build and install (py)proteolizard-data**](#build-and-install-pyproteolizard)
* [**Data handle and exposed meta data**](#data-handle-and-exposed-meta-data)
* [**Slices, frames, spectra**](#slices-frames-spectra)
* [**Binning, vectorization, pointclouds**](#binning-vectorization-pointclouds)

---
### Build and install (py)proteolizard
Right now, it's buggy and hard.

---
### Data handle and exposed meta data
To fetch data from a bruker timsTOF file, use the `PyTimsDataHandle` object. It receives a path to a bruker .d folder 
as it's only construction argument. On construction, it automatically reads experiment metadata which is generated 
during acquisition and stored by bruker software in a sqlite database.
```python
from pyproteolizard.data import PyTimsDataHandle

data_handle = PyTimsDataHandle('path/to/experiment.d')

print(data_handle.meta_data)
```

This then will give you an overview table, desribing multiple properties of each equired frame in the dataset:

|    |   Id |     Time | Polarity   |   ScanMode |   MsMsType |   TimsId |   MaxIntensity |   SummedIntensities |   NumScans |   NumPeaks |   MzCalibration |      T1 |      T2 |   TimsCalibration |   PropertyGroup |   AccumulationTime |   RampTime |
|---:|-----:|---------:|:-----------|-----------:|-----------:|---------:|---------------:|--------------------:|-----------:|-----------:|----------------:|--------:|--------:|------------------:|----------------:|-------------------:|-----------:|
|  0 |    1 | 0.597604 | +          |          8 |          0 |        0 |           7594 |            13708637 |        918 |     260611 |               1 | 25.6819 | 25.308  |                 1 |               1 |             99.953 |     99.953 |
|  1 |    2 | 0.849863 | +          |          8 |          8 |   626688 |            787 |               50776 |        918 |        609 |               1 | 25.6819 | 25.308  |                 1 |               1 |             99.953 |     99.953 |
|  2 |    3 | 0.955152 | +          |          8 |          8 |   630784 |           1066 |               54927 |        918 |        659 |               1 | 25.6819 | 25.308  |                 1 |               1 |             99.953 |     99.953 |
|  3 |    4 | 1.06146  | +          |          8 |          8 |   634880 |            721 |               39725 |        918 |        512 |               1 | 25.6819 | 25.308  |                 1 |               1 |             99.953 |     99.953 |
|  4 |    5 | 1.16897  | +          |          8 |          8 |   638976 |            471 |               41510 |        918 |        612 |               1 | 25.6819 | 25.3082 |                 1 |               1 |             99.953 |     99.953 |

Conveniently, a data handle also stores two numpy arrays, `dh.precursor_frames` and `dh.fragmet_frames`, allowing you
to collectively get all frame ids that correspond either to unfragmented or fragmented frames, respectively.

---
### Frames, spectra, slices
Mass spectrometry datasets created using liquid chromatography combined with Ion-mobility result in raw data best 
described as sparse 3D tensors (multidimensional matrices), having a retention-time, ion-mobility and mz-axis where
entries represent measured ion-intensities.

A lot of raw-data processing tasks revolve around identifying and grouping intensity values of interest. For example, 
feature detection on MS-I level means to identify peptide ions together with their charge state, monoisotopic mass and
dimensional extension for calculation of their absolute intensity.
There exist different approaches to tackle this problem, including top down (splitting marginal intensity distributions 
into individual peaks), bottom up (grouping together individual peaks into patterns) or mixed versions of those.
What they have in common is the fact that features are highly local, meaning it is sufficient to look at small pieces 
of a complete dataset to extract individual ones. 

`proteolizard-data` uses individual frames as its smallest unit of data to be loaded into memory. A frame corresponds to 
a single retention time with a collection of scans (non-normalized ion-mobility values), each scan holding a single 
mass spectrum. The frame can then easily be  split into spectra, but it also exposes a raw-data view that just contains
1D arrays of indices and intensity values. Those arrays can always be exposed individually or grouped.
```python
from pyproteolizard.data import PyTimsDataHandle, TimsFrame, MzSpectrum
data_handle = PyTimsDataHandle('path/to/experiment.d')

frame = data_handle.get_frame(1)
print(frame.data())
```

|    |   frame |   scan |       mz |   intensity |
|---:|--------:|-------:|---------:|------------:|
|  0 |       1 |     35 |  709.66  |           9 |
|  1 |       1 |     36 | 1676.65  |           9 |
|  2 |       1 |     38 | 1253.6   |          50 |
|  3 |       1 |     41 |  946.875 |           9 |
|  4 |       1 |     42 | 1191.3   |           9 |

```python
from pyproteolizard.data import PyTimsDataHandle, TimsFrame, MzSpectrum
data_handle = PyTimsDataHandle('path/to/experiment.d')

frame = data_handle.get_frame(1)
print(frame.data())
```

---
### Binning, vectorization, dense representations
DUMMY.