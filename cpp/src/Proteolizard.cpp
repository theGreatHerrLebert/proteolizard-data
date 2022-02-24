#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <string>
#include <cmath>

#include "Frame.h"
#include "Hashing.h"
#include "VectorizedSpectrum.h"
#include "Spectrum.h"
#include "ExposedTimsDataHandle.h"
#include "Slice.h"

namespace py = pybind11;

PYBIND11_MODULE(libproteolizard, h) {
h.doc() = "exposing fast timsTOF processing and access functions to be called from python";
// SPECTRUM
py::class_<MzSpectrumPL>(h, "MzSpectrumPL")
    .def(py::init<int &, int &, std::vector<double> &, std::vector<int> &>())
    .def("getFrameId", [](MzSpectrumPL &self){
        return self.frameId;
    })
    .def("getScanId", [](MzSpectrumPL &self){
        return self.scanId;
    })
    .def("getMzs", [](MzSpectrumPL &self){
        return py::array(py::cast(self.mz));
    })
    .def("getIntensities", [](MzSpectrumPL &self){
        return py::array(py::cast(self.intensity));
    })
    .def("vectorize", [](MzSpectrumPL &self, int resolution){
        return self.vectorize(resolution);
    })
    // TODO: lift logic out of bind
    .def("windows", [](MzSpectrumPL &self, double windowLength, bool overlapping, int minPeaks, int minIntensity){
        auto windows = self.windows(windowLength, overlapping, minPeaks, minIntensity);
        std::vector<int> bins;
        bins.reserve(windows.size());

        std::vector<MzSpectrumPL> win;
        win.reserve(windows.size());

        for(const auto& [key, value]: windows){
            bins.push_back(key);
            win.push_back(value);
        }

        return py::make_tuple(bins, win);
    })
    ;

// VECTORIZED SPECTRUM
py::class_<MzVectorPL>(h, "MzVectorPL")
    .def(py::init<int &, int &, int &, std::vector<int> &, std::vector<int> &>())
    .def("getFrameId", [](MzVectorPL &self){
        return self.frameId;
    })
    .def("getScanId", [](MzVectorPL &self){
        return self.scanId;
    })
    .def("getResolution", [](MzVectorPL &self){
        return self.resolution;
    })
    .def("getIndices", [](MzVectorPL &self){
        return py::array(py::cast(self.indices));
    })
    .def(py::self + py::self)
    .def("getValues", [](MzVectorPL &self){
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
        std::vector<MzSpectrumPL> ret;
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
        std::vector<MzVectorPL> retVec;
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
    .def("getFragments", [](TimsSlicePL &self){
        return self.fragments;
    })
    .def("filterRanged", [](TimsSlicePL &self, int scanMin, int scanMax, double mzMin, double mzMax, int intensityMin){
        return self.filterRanged(scanMin, scanMax, mzMin, mzMax, intensityMin);
    })
    .def("getPoints", [](TimsSlicePL &self, bool precursor){
        return self.getPoints3D(precursor);
    })
    ;
    py::class_<TimsHashGenerator>(h, "TimsHashGenerator")
    .def(py::init<int, int, int, int>())

    .def("hashMzSpectrum", [](TimsHashGenerator& self, MzSpectrumPL& spectrum, int minPeaks,
                              int minIntensity, double windowLength, bool overlapping, bool restricted){
        py::array a = py::cast(self.hashSpectrum(spectrum, minPeaks, minIntensity, windowLength, overlapping, restricted));
        return a;
    })
    .def("getMatrixCopy", &TimsHashGenerator::getMatrixCopy)
    ;
}
