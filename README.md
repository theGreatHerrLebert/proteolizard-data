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

---
### Slices, frames, spectra
Ion-mobility datasets can be interpreted as sparse 3D tensors (multidimensional matrices), having a retention-time,
ion-mobility and mz-axis where entries represent measured ion-intensities.

---
### Binning, vectorization, dense representations
DUMMY.