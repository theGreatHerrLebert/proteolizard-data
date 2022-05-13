# proteolizard-data
### An object-oriented library of C++ classes and python-wrappers to seamlessly integrate timsTOF raw-data with python

![](logo.png)

## Context
Ion-mobility enhanced tandem-MS coupled to liquid chromatography is rapidly becoming the method of 
choice for the analysis of high-complexity samples generated in proteomics, lipidomics and metabolomics.
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
* [**Frames, spectra, slices**](#frames-spectra-slices)
* [**Binning, vectorization, pointclouds**](#binning-vectorization-pointclouds)

---
### Build and install (py)proteolizard
```sh
shell> git clone --recursive https://github.com/theGreatHerrLebert/proteolizard-data
shell> cd proteolizard-data
```

```sh
shell> mkdir build && cd build
shell> cmake ../cpp -DCMAKE_BUILD_TYPE=Release
shell> make 
```

```sh
shell> cmake --install . --prefix=some/prefix/path
```
---
### Data handle and exposed meta data
To fetch data from a bruker timsTOF file, use the `PyTimsDataHandle` object. It receives a path to a bruker .d folder 
as it's only construction argument. On construction, it automatically reads experiment metadata which is generated 
during acquisition and stored by bruker software in a sqlite database.
```python
from proteolizarddata.data import PyTimsDataHandle

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
a single retention time with a collection of scans (non-normalized ion-mobility values), each scan holding one 
mass spectrum. The frame can then easily be  split into spectra, but it also exposes a raw-data view that just contains
1D arrays of indices and intensity values. Those arrays can always be retrieved individually or grouped.
```python
from proteolizarddata.data import PyTimsDataHandle, TimsFrame, MzSpectrum
data_handle = PyTimsDataHandle('path/to/experiment.d')

frame = data_handle.get_frame(1)
print(frame.data())
```

This returns a pandas table:

|    |   frame |   scan |       mz |   intensity |
|---:|--------:|-------:|---------:|------------:|
|  0 |       1 |     35 |  709.66  |           9 |
|  1 |       1 |     36 | 1676.65  |           9 |
|  2 |       1 |     38 | 1253.6   |          50 |
|  3 |       1 |     41 |  946.875 |           9 |
|  4 |       1 |     42 | 1191.3   |           9 |


But you could also query a dimension individually:
```python
print(frame.inverse_ion_mobility()[:5])
```
giving you the first five inverse ion-mobility values:

array([1.59895005, 1.59787028, 1.59571042, 1.59246983, 1.59138941])

Let us also have a look at splitting a frame into spectra:

```python
first_five_spectra = frame.get_spectra()[:5]
print(first_five_spectra)
```

this retruns a list of MzSpectrum objects:

[MzSpectrum(frame: 1, scan: 290, sum intensity: 10106),
 MzSpectrum(frame: 1, scan: 291, sum intensity: 10783),
 MzSpectrum(frame: 1, scan: 292, sum intensity: 9939),
 MzSpectrum(frame: 1, scan: 293, sum intensity: 9518),
 MzSpectrum(frame: 1, scan: 294, sum intensity: 12078)]

Both TimsFrame and MzSpectrum override pythons `+` operator, making it possible to add together spectra or frames.
Notice tough, that this is not a strictly associative and commutative operation, as resulting frame id and scan id will 
be dependent on argument order!

Finally, a data handle can also extract a collection of frames represented by an object TimsSlice, giving you the option
to load abitrary blocks of a timsTOF experiments into memory. This can either be done by a set of frame ids or retention
time start and end values, for example:
```python
# caution, retention times are stored in seconds
slice = data_handle.get_slice_rt_range(25*60, 25*60 + 30)
print(slice)
```
TimsSlice(frame start: 14143, frame end: 14418)

giving you half a minute of retention time, consisting both of precursor and fragment frames. You can either get
collections of frames via methods `get_precursor_frames()` and `get_fragment_frames()` or 1D arrays of all data-points 
with methods `get_precursor_points()` and `get_fragment_points()`. Lastly, TimsSlice, TimsFrame and MzSpectrum allow 
you to fast-filter spectra or collections of spectra using range queries. This is super useful to restrict an extraction 
to ranges of ion-mobilities, mass-over-charges or remove low intensity signals.

---
### Binning, vectorization, dense representations
DUMMY.
