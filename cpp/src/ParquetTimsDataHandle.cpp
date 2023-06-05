//
// Created by administrator on 25.05.23.
//

#include "ParquetTimsDataHandle.h"

#include <sstream>
#include <arrow/api.h>
#include <arrow/dataset/api.h>
#include <arrow/io/api.h>
#include <arrow/result.h>
#include <parquet/arrow/reader.h>
#include <parquet/exception.h>

ParquetTimsDataHandle::ParquetTimsDataHandle(std::string dp):
rawDataPath(initialize(dp) + "raw/"),
metaDataPath(initialize(dp) + "meta/") {}

TimsFramePL ParquetTimsDataHandle::getTimsFramePL(int blockId,
                                                  int rowGroupId,
                                                  int rowGroupIndexStart,
                                                  int rowGroupIndexStop) {
    // construct a path to the parquet file
    std::stringstream filePath;
    ss << this->rawDataPath + "BLOCK-" << blockId << ".parquet";

    // Create a ParquetFileReader from a FileSource
    std::shared_ptr<arrow::io::ReadableFile> infile;

    PARQUET_ASSIGN_OR_THROW(
            infile,
            arrow::io::ReadableFile::Open(ss.str()));
    std::unique_ptr<parquet::arrow::FileReader> reader;
    PARQUET_THROW_NOT_OK(
            parquet::arrow::OpenFile(infile, arrow::default_memory_pool(), &reader));

    // Create a dataset from the ParquetFileReader
    auto format = std::make_shared<arrow::dataset::ParquetFileFormat>();
    auto file_source = std::make_shared<arrow::dataset::InMemoryDataset>(reader->parquet_reader());
    auto dataset = arrow::dataset::Dataset::Make({file_source}, format);

    return {};
}

static std::string ParquetTimsDataHandle::initialize(const std::string &value) {
    {
        // Perform the logic here and return the string.
        if (value.back() != '/') {
            return value + '/';
        } else {
            return value;
        }
    }
}

