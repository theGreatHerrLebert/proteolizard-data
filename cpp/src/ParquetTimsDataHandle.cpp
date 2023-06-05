//
// Created by administrator on 25.05.23.
//

#include "ParquetTimsDataHandle.h"

#include <sstream>
#include <arrow/api.h>
#include <arrow/io/api.h>
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
    filePath << this->rawDataPath + "BLOCK-" << blockId << ".parquet";

    return {};
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

