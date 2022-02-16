#include <vector>
#include <tuple>
#include <ostream>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <string>
#include <cmath>

#include "Frame.h"
#include "VectorizedSpectrum.h"
#include "Spectrum.h"
#include "ExposedTimsDataHandle.h"
#include "Slice.h"

namespace py = pybind11;

PYBIND11_MODULE(libproteolizard, h) {
h.doc() = "exposing fast timsTOF processing and access functions to be called from python";
// SPECTRUM
py::class_<MzSpectrum>(h, "MzSpectrum")
    .def(py::init<int &, int &, std::vector<double> &, std::vector<int> &>())
    .def("getFrameId", [](MzSpectrum &self){
        return self.frameId;
    })
    .def("getScanId", [](MzSpectrum &self){
        return self.scanId;
    })
    .def("getMzs", [](MzSpectrum &self){
        return py::array(py::cast(self.mz));
    })
    .def("getIntensities", [](MzSpectrum &self){
        return py::array(py::cast(self.intensity));
    })
    .def("toResolution", [](MzSpectrum &self, int resolution){
        return self.toResolution(resolution);
    })
    ;

// VECTORIZED SPECTRUM
py::class_<MzVector>(h, "MzVector")
    .def(py::init<int &, int &, int &, std::vector<int> &, std::vector<int> &>())
    .def("getFrameId", [](MzVector &self){
        return self.frameId;
    })
    .def("getScanId", [](MzVector &self){
        return self.scanId;
    })
    .def("getResolution", [](MzVector &self){
        return self.resolution;
    })
    .def("getIndices", [](MzVector &self){
        return py::array(py::cast(self.indices));
    })
    .def(py::self + py::self)
    .def("getValues", [](MzVector &self){
        return py::array(py::cast(self.values));
    });
// FRAME
py::class_<TimsFramePL>(h, "TimsFrame")
    .def(py::init<int &, std::vector<int> &, std::vector<double> &, std::vector<int> &, std::vector<int> &, std::vector<double> &>())
    .def("getFrameId", [](TimsFramePL &self){
        return self.frameId;
    })
    .def("getScans", [](TimsFramePL &self){
        return py::array(py::cast(self.scans));
    })
    .def("getMzs", [](TimsFramePL &self){
        return py::array(py::cast(self.mzs));
    })
    .def("getIntensities", [](TimsFramePL &self){
        return py::array(py::cast(self.intensities));
    })
    .def("getTofs", [](TimsFramePL &self){
    return py::array(py::cast(self.tofs));
    })
    .def("getInverseMobilities", [](TimsFramePL &self){
    return py::array(py::cast(self.inv_ion_mobs));
    })
    .def("vectorize", [](TimsFramePL &self, int resolution){
        return self.vectorize(resolution);
    })
    .def("toResolution", [](TimsFramePL &self, int resolution){
        return self.toResolution(resolution);
    })
    .def("fold", [](TimsFramePL &self, int resolution, int width){
        return self.fold(resolution, width);
    })
    .def("filterRanged", [](TimsFramePL &self, int scanMin, int scanMax, double mzMin, double mzMax, int intensityMin){
        return self.filterRanged(scanMin, scanMax, mzMin, mzMax, intensityMin);
    })
    .def(py::self + py::self)
    // TODO: move logic out of module
    .def("getMzSpectra", [](TimsFramePL &self){
        auto spectra = self.spectra();
        std::vector<MzSpectrum> ret;
        ret.reserve(spectra.size());

        for(auto&[key, value]: spectra)
            ret.push_back(value);

        return ret;
    })
    ;
// VECTORIZED FRAME
py::class_<TimsFrameVectorizedPL>(h, "TimsFrameVectorized")
    .def(py::init<int &, int &, std::vector<int> &, std::vector<int> &, std::vector<int> &>())
    .def("getFrameId", [](TimsFrameVectorizedPL &self){
        return self.frameId;
    })
    .def("getResolution", [](TimsFrameVectorizedPL &self){
        return self.resolution;
    })
    .def("getScans", [](TimsFrameVectorizedPL &self){
        return py::array(py::cast(self.scans));
    })
    .def("getIndices", [](TimsFrameVectorizedPL &self){
        return py::array(py::cast(self.indices));
    })
    .def(py::self + py::self)
    .def("filterRanged", [](TimsFrameVectorizedPL &self, int scanMin, int scanMax, int indexMin, int indexMax, int intensityMin){
        return self.filterRanged(scanMin, scanMax, indexMin, indexMax, intensityMin);
    })
    .def("getSpectra", [](TimsFrameVectorizedPL &self){
        auto spec = self.spectra();
        std::vector<MzVector> retVec;
        retVec.reserve(spec.size());
        for(const auto& [_, value]: spec)
            retVec.push_back(value);
        return retVec;
    })
    .def("getValues", [](TimsFrameVectorizedPL &self){
        return py::array(py::cast(self.values));
    });
// DATA ACCESS
py::class_<ExposedTimsDataHandle>(h, "ExposedTimsDataHandle")
    .def(py::init<std::string &, std::string &>())
    .def("getFrame", [](ExposedTimsDataHandle &self, int frameId){
        TimsFramePL frame =  self.getTimsFramePL(frameId);
        return frame;
    })
    .def("getSlice", [](ExposedTimsDataHandle &self, std::vector<int> precursors, std::vector<int> fragments){
        TimsSlicePL slice = self.getTimsSlicePL(precursors, fragments);
        return slice;
    })
    ;
// ----------------
py::class_<Points3D>(h, "Points3D")
    .def(py::init<std::vector<int> &, std::vector<int> &, std::vector<double> &, std::vector<double> &, std::vector<int> &>())
    .def("getFrames", [](Points3D &self){
         py::array a = py::cast(self.frame);
        return a;
    })
    .def("getScans", [](Points3D &self){
        py::array a = py::cast(self.scan);
        return a;
    })
    .def("getMz", [](Points3D &self){
        py::array a = py::cast(self.mz);
        return a;
    })
    .def("getInvIonMobility", [](Points3D &self){
        py::array a = py::cast(self.invIonMobility);
        return a;
    })
    .def("getIntensity", [](Points3D &self){
        py::array a = py::cast(self.intensity);
        return a;
    })
    ;
    // -------------- Tims-SLICE ------------
    py::class_<TimsSlicePL>(h, "TimsSlicePL")
    .def(py::init<std::vector<TimsFramePL> &, std::vector<TimsFramePL> &>())
    .def("getPrecursors", [](TimsSlicePL &self){
        return self.precursors;
    })
    .def("getFragments", [](TimsSlicePL &self, int id){
        return self.fragments;
    })
    .def("filterRanged", [](TimsSlicePL &self, int scanMin, int scanMax, double mzMin, double mzMax, int intensityMin){
        return self.filterRanged(scanMin, scanMax, mzMin, mzMax, intensityMin);
    })
    .def("getPoints", [](TimsSlicePL &self, bool precursor){
        return self.getPoints3D(precursor);
    })
    ;
}
