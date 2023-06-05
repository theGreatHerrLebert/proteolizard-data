//
// Created by administrator on 25.05.23.
//

#include "ParquetTimsDataHandle.h"

#include <sstream>
#include <iostream>
#include <algorithm>
#include <execution>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>
#include <parquet/exception.h>

ParquetTimsDataHandle::ParquetTimsDataHandle(std::string dp):
rawDataPath(initialize(dp) + "raw/"),
metaDataPath(initialize(dp) + "meta/") {}

ParquetTimsDataHandle::~ParquetTimsDataHandle() = default;

TimsFramePL ParquetTimsDataHandle::getTimsFramePL(int blockId,
                                                  int rowGroupId,
                                                  int rowGroupIndexStart,
                                                  int rowGroupIndexStop) {
    // construct a path to the parquet file
    std::stringstream filePath;
    filePath << this->rawDataPath + "BLOCK-" << blockId << ".parquet";

    return {};
}

std::vector<TimsFramePL> createFramesFromTable(const std::shared_ptr<arrow::Table>& table) {
    // Vector to hold the frames
    std::vector<TimsFramePL> frame_vec;

    // Data vectors for each column
    std::vector<int> scans, tofs, intensities;
    std::vector<double> mzs, inv_ion_mobs;
    double retentionTime;
    int32_t previous_frameId = -2;

    for (int i = 0; i < table->num_rows(); ++i) {
        // Extract the frameId for this row
        int32_t frame_id = std::static_pointer_cast<arrow::Int32Array>(table->GetColumnByName("frameId")->chunk(0))->Value(i);

        // If the frameId changes, save the existing data and create a new frame
        if(frame_id != previous_frameId && i > 0) {
            TimsFramePL frame(previous_frameId, retentionTime, scans, mzs, intensities, tofs, inv_ion_mobs);
            frame_vec.push_back(frame);

            // Clear data vectors for new frame
            scans.clear();
            mzs.clear();
            intensities.clear();
            inv_ion_mobs.clear();
        }

        // Fetch the data for each column
        scans.push_back(std::static_pointer_cast<arrow::Int32Array>(table->GetColumnByName("scanId")->chunk(0))->Value(i));
        mzs.push_back(std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("mzId")->chunk(0))->Value(i) / 1000000.0);
        intensities.push_back(std::static_pointer_cast<arrow::Int32Array>(table->GetColumnByName("intensity")->chunk(0))->Value(i));
        inv_ion_mobs.push_back(std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("mobilityId")->chunk(0))->Value(i));

        // Update previous_frameId for next iteration
        previous_frameId = frame_id;
    }

    // Don't forget to add the last frame
    if(!scans.empty() && !mzs.empty() && !intensities.empty() && !tofs.empty() && !inv_ion_mobs.empty()){
        TimsFramePL frame(previous_frameId, retentionTime, scans, mzs, intensities, tofs, inv_ion_mobs);
        frame_vec.push_back(frame);
    }

    return frame_vec;
}

std::shared_ptr<arrow::Table> readTable(std::string path) {
    // Create Arrow file input stream
    std::shared_ptr<arrow::io::ReadableFile> infile;
    PARQUET_ASSIGN_OR_THROW(
            infile,
            arrow::io::ReadableFile::Open(path,
                                          arrow::default_memory_pool()));

    // Create a Parquet reader
    std::unique_ptr<parquet::arrow::FileReader> reader;
    PARQUET_THROW_NOT_OK(
            parquet::arrow::OpenFile(infile, arrow::default_memory_pool(), &reader));

    // Read the entire Parquet file into an Arrow table
    std::shared_ptr<arrow::Table> table;
    PARQUET_THROW_NOT_OK(reader->ReadTable(&table));

    return table;
}

TimsSlicePL ParquetTimsDataHandle::getTimsSlicePL(int frameIdStart,
                                                  int frameIdEnd,
                                                  std::vector<int> &blockIds,
                                                  std::vector<int> &msMsTypes) {

    std::vector<std::string> paths;

    for (auto i : blockIds) {
        std::stringstream ss;
        ss << this->rawDataPath + "BLOCK-" << i << ".parquet";
        std::string path_to_file = ss.str();
        paths.push_back(path_to_file);
    }

    std::vector<std::shared_ptr<arrow::Table>> tables;
    tables.resize(paths.size());

    auto p = [](std::string& path) -> std::shared_ptr<arrow::Table> {
        return readTable(path);
    };

    std::transform(std::execution::par_unseq, paths.begin(), paths.end(), tables.begin(), p);

    std::vector<std::vector<TimsFramePL>> retFrameList;
    retFrameList.resize(tables.size());

    auto f = [](std::shared_ptr<arrow::Table>& f) -> std::vector<TimsFramePL> {
        return createFramesFromTable(f);
    };

    std::transform(std::execution::par_unseq, tables.begin(), tables.end(), retFrameList.begin(), f);

    std::vector<TimsFramePL> precursors;
    std::vector<TimsFramePL> fragments;

    int counter = 0;

    for(auto &frames: retFrameList){
        for(auto &frame: frames){
            int type = msMsTypes[counter];
            if(frame.frameId >= frameIdStart && frame.frameId <= frameIdEnd){
                if(type > 0){
                    fragments.push_back(frame);
                } else {
                    precursors.push_back(frame);
                }
            }
            counter++;
        }
    }
    return {precursors, fragments};
}

std::string ParquetTimsDataHandle::initialize(const std::string &value) {
    {
        // Perform the logic here and return the string.
        if (value.back() != '/') {
            return value + '/';
        } else {
            return value;
        }
    }
}

