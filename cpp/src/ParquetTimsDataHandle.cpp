//
// Created by administrator on 25.05.23.
//

#include "ParquetTimsDataHandle.h"

#include <sstream>
#include <arrow/dataset/api.h>
#include <arrow/dataset/file_parquet.h>
#include <arrow/filesystem/localfs.h>
#include <arrow/io/api.h>

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

    // Use the local filesystem
    auto fs = std::make_shared<arrow::fs::LocalFileSystem>();

    // Define options
    arrow::dataset::FileSystemFactoryOptions options;

    // Create a ParquetFileFormat
    std::shared_ptr<arrow::dataset::ParquetFileFormat> format = std::make_shared<arrow::dataset::ParquetFileFormat>();

    // Create a DatasetFactory
    auto factory = arrow::dataset::FileSystemDatasetFactory::Make(fs, {filePath}, format, options);

    // Create the Dataset
    arrow::dataset::FinishOptions finish_options;

    auto result = factory->Finish(finish_options);

    /*
    if (!result.ok()) {
        // Handle error...
        return 1;
    }
     */

    std::shared_ptr<arrow::dataset::Dataset> dataset = *result;

    // Create a ScannerBuilder
    std::shared_ptr<arrow::dataset::ScannerBuilder> scanner_builder;

    auto scanner_builder_result = dataset->NewScan();

    /*
    if (!scanner_builder_result.ok()) {
        // Handle error...
        return 1;
    }
     */

    scanner_builder = *scanner_builder_result;

    // Define a predicate
    auto field = arrow::field("my_field", arrow::int32());  // Adjust this to match a field in your data

    std::shared_ptr<arrow::Scalar> literal;

    auto status = arrow::MakeScalar(10, &literal);

    auto predicate = arrow::compute::equal(arrow::compute::field_ref(field->name()), literal);

    auto status = scanner_builder->Filter(predicate);

    /*
    // Apply the predicate
    if (auto status = scanner_builder->Filter(predicate); !status.ok()) {
        // Handle error...
        return 1;
    }
     */

    // Create the Scanner
    arrow::Result<std::shared_ptr<arrow::dataset::Scanner>> scanner_result = scanner_builder->Finish();

    /*
    if (!scanner_result.ok()) {
        // Handle error...
        return 1;
    }
    */

    std::shared_ptr<arrow::dataset::Scanner> scanner = *scanner_result;

    // Scan the Dataset into a Table
    arrow::Result<std::shared_ptr<arrow::Table>> table_result = scanner->ToTable();
    /*
    if (!table_result.ok()) {
        // Handle error...
        return 1;
    }
     */
    std::shared_ptr<arrow::Table> table = *table_result;

    // Now you can use the Table...

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

