#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <string>

#include "Frame.h"
#include "VectorizedSpectrum.h"
#include "Spectrum.h"
#include "ExposedTimsDataHandle.h"
#include "Slice.h"
#include "SliceVectorized.h"
#include "TimsBlock.h"
#include "TimsBlockVectorized.h"

namespace py = pybind11;

PYBIND11_MODULE(libproteolizarddata, h) {
    h.doc() = "exposing fast timsTOF processing and access functions to be called from python";

    // ---------------- CLASS SPECTRUM ------------
    py::class_<MzSpectrumPL>(h, "MzSpectrumPL")

            // -------------- CONSTRUCTOR ---------------
            .def(py::init<int &, int &, std::vector<double> &, std::vector<int> &>())

            // -------------- MEMBER ---------------
            .def("getFrameId", [](MzSpectrumPL &self) {
                return self.frameId;
            })
            .def("filter", [](MzSpectrumPL &self, double mzMin, double mzMax, int intensityMin) {
                    return self.filter(mzMin, mzMax, intensityMin);
            })
            .def("getScanId", [](MzSpectrumPL &self) {
                return self.scanId;
            })
            .def("getMzs", [](MzSpectrumPL &self) {
                return py::array(py::cast(self.mz));
            })
            .def("getIntensities", [](MzSpectrumPL &self) {
                return py::array(py::cast(self.intensity));
            })
            .def("getSize", [](MzSpectrumPL &self){
                return self.intensity.size();
            })
            .def(py::self + py::self)
            .def(py::self += py::self)
            .def(py::self * float())
            .def(float() * py::self)
            .def("vectorize", [](MzSpectrumPL &self, int resolution) {
                return self.vectorize(resolution);
            })
            .def("toResolution", [](MzSpectrumPL &self, int resolution){
                return self.toResolution(resolution);
            })
            .def("toCentroided", [](MzSpectrumPL &self, int baselineNoiseLevel, double sigma){
                return self.toCentroided(baselineNoiseLevel,sigma);
            })
            .def("push", [](MzSpectrumPL &self, MzSpectrumPL &other){
                return self.push(other);
            })
            .def("windows",
                 [](MzSpectrumPL &self, double windowLength, bool overlapping, int minPeaks, int minIntensity) {
                auto p = self.exportWindows(windowLength, overlapping, minPeaks, minIntensity);
                     return py::make_tuple(py::array(py::cast(p.first)), p.second);
                 });
    // ---------------- CLASS VECTORIZED-SPECTRUM ------------
    py::class_<MzVectorPL>(h, "MzVectorPL")

            // -------------- CONSTRUCTOR ---------------
            .def(py::init<int &, int &, int &, std::vector<int> &, std::vector<int> &>())

            // -------------- MEMBER ---------------
            .def("getFrameId", [](MzVectorPL &self) {
                return self.frameId;
            })
            .def("getScanId", [](MzVectorPL &self) {
                return self.scanId;
            })
            .def("getResolution", [](MzVectorPL &self) {
                return self.resolution;
            })
            .def("getIndices", [](MzVectorPL &self) {
                return py::array(py::cast(self.indices));
            })
            .def(py::self + py::self)
            .def("getValues", [](MzVectorPL &self) {
                return py::array(py::cast(self.values));
            });
    // ---------------- CLASS TIMS-FRAME ------------
    py::class_<TimsFramePL>(h, "TimsFrame")

            // -------------- CONSTRUCTOR ---------------
            .def(py::init<int &, double &, std::vector<int> &, std::vector<double> &, std::vector<int> &, std::vector<int> &, std::vector<double> &>())

            // -------------- MEMBER ---------------
            .def("getFrameId", [](TimsFramePL &self) {
                return self.frameId;
            })
            .def("getRetentionTime", [](TimsFramePL &self) {
                return self.retentionTime;
            })
            .def("getScans", [](TimsFramePL &self) {
                return py::array(py::cast(self.scans));
            })
            .def("getMzs", [](TimsFramePL &self) {
                return py::array(py::cast(self.mzs));
            })
            .def("getIntensities", [](TimsFramePL &self) {
                return py::array(py::cast(self.intensities));
            })
            .def("getTofs", [](TimsFramePL &self) {
                return py::array(py::cast(self.tofs));
            })
            .def("getInverseMobilities", [](TimsFramePL &self) {
                return py::array(py::cast(self.inv_ion_mobs));
            })
            .def("vectorize", [](TimsFramePL &self, int resolution) {
                return self.vectorize(resolution);
            })
            .def("toResolution", [](TimsFramePL &self, int resolution) {
                return self.toResolution(resolution);
            })
            .def("fold", [](TimsFramePL &self, int resolution, int width) {
                return self.fold(resolution, width);
            })
            .def("getDenseWindows", [](TimsFramePL &self, int resolution, int minPeaksPerWindow, int minIntensity,
                                       double windowLength, bool overlapping){

                return self.denseWindowMatrix(resolution, minPeaksPerWindow, minIntensity, windowLength, overlapping);
            })
            .def("getVectorizedWindows", [](TimsFramePL &self, int resolution, int minPeaksPerWindow, int minIntensity,
            double windowLength, bool overlapping){

                auto r = self.vectorizedWindows(resolution, minPeaksPerWindow, minIntensity, windowLength, overlapping);

                return py::make_tuple(py::array(py::cast(r.first.first)), py::array(py::cast(r.first.second)),r.second);
            })
            .def("filterRanged",
                 [](TimsFramePL &self, int scanMin, int scanMax, double mzMin, double mzMax, int intensityMin) {
                     return self.filterRanged(scanMin, scanMax, mzMin, mzMax, intensityMin);
                 })
            .def(py::self + py::self)
            .def("getHashingBlocks", [](TimsFramePL &self, int resolution,
                                     int minPeaksPerWindow,
                                     int minIntensity,
                                     double windowLength,
                                     bool overlapping){

                auto t = self.getHashingBlocks(resolution, minPeaksPerWindow, minIntensity, windowLength, overlapping);

                return py::make_tuple(t.counter, py::array(py::cast(t.rowIndex)),
                                      py::array(py::cast(t.scan)),
                                      py::array(py::cast(t.bin)),
                                      py::array(py::cast(t.indices)),
                                      py::array(py::cast(t.values)));
            })
            .def("getMzSpectra", [](TimsFramePL &self) {
                return self.exportSpectra();
            });
    // ---------------- CLASS TIMS-VECTORIZED-FRAME ------------
    py::class_<TimsFrameVectorizedPL>(h, "TimsFrameVectorized")

            // -------------- CONSTRUCTOR ---------------
            .def(py::init<int &, int &, std::vector<int> &, std::vector<int> &, std::vector<int> &>())

            // -------------- MEMBER ---------------
            .def("getFrameId", [](TimsFrameVectorizedPL &self) {
                return self.frameId;
            })
            .def("getResolution", [](TimsFrameVectorizedPL &self) {
                return self.resolution;
            })
            .def("getScans", [](TimsFrameVectorizedPL &self) {
                return py::array(py::cast(self.scans));
            })
            .def("getIndices", [](TimsFrameVectorizedPL &self) {
                return py::array(py::cast(self.indices));
            })
            .def(py::self + py::self)
            .def("filterRanged", [](TimsFrameVectorizedPL &self, int scanMin, int scanMax, int indexMin, int indexMax,
                                    int intensityMin) {
                return self.filterRanged(scanMin, scanMax, indexMin, indexMax, intensityMin);
            })
            .def("getSpectra", [](TimsFrameVectorizedPL &self) {
                auto spec = self.spectra();
                std::vector<MzVectorPL> retVec;
                retVec.reserve(spec.size());
                for (const auto&[_, value]: spec)
                    retVec.push_back(value);
                return retVec;
            })
            .def("getValues", [](TimsFrameVectorizedPL &self) {
                return py::array(py::cast(self.values));
            });
            // -------------- CLASS TIMS-BLOCK ------------
    py::class_<TimsBlockPL>(h, "TimsBlockPL")
            // -------------- CONSTRUCTOR ---------------
            .def(py::init<std::vector<int> &, std::vector<int> &, std::vector<int> &, std::vector<int> &,
            std::vector<double> &, std::vector<double> &, std::vector<double> &>())
            .def("getIndices", [](TimsBlockPL &self) {
            return py::array(py::cast(self.getIndices()));
            })
            .def("getValues", [](TimsBlockPL &self) {
            return py::array(py::cast(self.getValues()));
            })
            .def("filterRanged", [](TimsBlockPL &self, int scanMin, int scanMax,
                    double mzMin, double mzMax, int intensityMin) {
            return self.filterRanged(scanMin, scanMax, mzMin, mzMax, intensityMin);
            })
            .def("getTimsBlockVectorized", [](TimsBlockPL &self, int resolution) {
            return self.getBlockVectorized(resolution);
            })
            ;
    py::class_<TimsBlockVectorizedPL>(h, "TimsBlockVectorizedPL")
            .def(py::init<int &, int &, std::vector<int> &, std::vector<int> &, std::vector<int> &, std::vector<int> &>())
            .def("getData", [](TimsBlockVectorizedPL &self) {
                return py::array(py::cast(self.getData()));
            })
            .def("filterRanged", [](TimsBlockVectorizedPL &self, int scanMin, int scanMax, int indexMin, int indexMax,
                                    int intensityMin) {
                return self.filterRanged(scanMin, scanMax, indexMin, indexMax, intensityMin);
            })
            ;
    // ---------------- CLASS TIMS-DATA-HANDLE ------------
    py::class_<ExposedTimsDataHandle>(h, "ExposedTimsDataHandle")

            // -------------- CONSTRUCTOR ---------------
            .def(py::init<std::string &, std::string &>())

            // -------------- MEMBER ---------------
            .def("getFrame", [](ExposedTimsDataHandle &self, int frameId) {
                TimsFramePL frame = self.getTimsFramePL(frameId);
                return frame;
            })
            .def("getGlobalMzAxis", [](ExposedTimsDataHandle &self){
                std::vector<double> globalMz = self.getGlobalMzAxis();
                return py::array(py::cast(globalMz));
            })
            .def("getSlice", [](ExposedTimsDataHandle &self, std::vector<int> precursors, std::vector<int> fragments) {
                TimsSlicePL slice = self.getTimsSlicePL(precursors, fragments);
                return slice;
            })
            .def("getBlock", [](ExposedTimsDataHandle &self, std::vector<int> frameIds) {
                auto block = self.getTimsBlockPL(frameIds);
                return block;
            });
    // ---------------- CLASS TIMS-POINTS3D ------------
    py::class_<Points3D>(h, "Points3D")

            // -------------- CONSTRUCTOR ---------------
            .def(py::init<std::vector<int> &, std::vector<int> &, std::vector<double> &, std::vector<double> &, std::vector<int> &>())

                    // -------------- MEMBER ---------------
            .def("getFrames", [](Points3D &self) {
                py::array a = py::cast(self.frame);
                return a;
            })
            .def("getScans", [](Points3D &self) {
                py::array a = py::cast(self.scan);
                return a;
            })
            .def("getMz", [](Points3D &self) {
                py::array a = py::cast(self.mz);
                return a;
            })
            .def("getInvIonMobility", [](Points3D &self) {
                py::array a = py::cast(self.invIonMobility);
                return a;
            })
            .def("getIntensity", [](Points3D &self) {
                py::array a = py::cast(self.intensity);
                return a;
            });

    // ---------------- CLASS TIMS-POINTS3D ------------
    py::class_<Points3DVectorized>(h, "Points3DVectorized")

    // -------------- CONSTRUCTOR ---------------
    .def(py::init<std::vector<int> &, std::vector<int> &, std::vector<int> &, std::vector<int> &>())

    // -------------- MEMBER ---------------
    .def("getFrames", [](Points3DVectorized &self) {
    py::array a = py::cast(self.frame);
    return a;
    })
    .def("getScans", [](Points3DVectorized &self) {
    py::array a = py::cast(self.scan);
    return a;
    })
    .def("getMz", [](Points3DVectorized &self) {
    py::array a = py::cast(self.mz);
    return a;
    })
    .def("getIntensity", [](Points3DVectorized &self) {
    py::array a = py::cast(self.intensity);
    return a;
    });

    // -------------- CLASS TIMS-SLICE ------------
    py::class_<TimsSlicePL>(h, "TimsSlicePL")

            // -------------- CONSTRUCTOR ---------------
            .def(py::init<std::vector<TimsFramePL> &, std::vector<TimsFramePL> &>())

                    // -------------- MEMBER ---------------
            .def("getPrecursors", [](TimsSlicePL &self) {
                return self.precursors;
            })
            .def("getFragments", [](TimsSlicePL &self) {
                return self.fragments;
            })
            .def("filterRanged",
                 [](TimsSlicePL &self, int scanMin, int scanMax, double mzMin, double mzMax, int intensityMin, double rtMin, double rtMax) {
                     return self.filterRanged(scanMin, scanMax, mzMin, mzMax, intensityMin, rtMin, rtMax);
                 })
             .def("getVectorizedSlice", [](TimsSlicePL &self, int resolution) {
                 return self.getVectorizedSlice(resolution);
             })
            .def("getPoints", [](TimsSlicePL &self, bool precursor) {
                return self.getPoints3D(precursor);
            });

    // -------------- CLASS TIMS-SLICE-VECTORIZED ------------
    py::class_<TimsSliceVectorizedPL>(h, "TimsSliceVectorizedPL")

            // -------------- CONSTRUCTOR ---------------
           .def(py::init<std::vector<TimsFrameVectorizedPL>&, std::vector<TimsFrameVectorizedPL> &>())
           .def("getPrecursors", [](TimsSliceVectorizedPL &self) {
               return self.precursors;
           })
           .def("getFragment", [](TimsSliceVectorizedPL &self) {
               return self.fragments;
           })
           .def("filterRanged", [](TimsSliceVectorizedPL &self, int frameMin, int frameMax, int scanMin, int scanMax, int mzMin, int mzMax, int intensityMin) {
               return self.filterRanged(frameMin, frameMax, scanMin, scanMax, mzMin, mzMax, intensityMin);
           })
           .def("getPoints", [](TimsSliceVectorizedPL &self, bool precursor) {
               return self.getPoints3D(precursor);
           });
}
