#ifndef CPP_EXPOSEDTIMSDATAHANDLE_H
#define CPP_EXPOSEDTIMSDATAHANDLE_H

#include "Frame.h"
#include "Slice.h"
#include "TimsBlock.h"
#include <memory>
#include <string>
#include <vector>

class ExposedTimsDataHandle{
public:
    // handle to read from raw data
    ExposedTimsDataHandle(const std::string& dp, const std::string& bp);
    ~ExposedTimsDataHandle();

    TimsFramePL getTimsFramePL(int frameId);
    TimsSlicePL getTimsSlicePL(std::vector<int>& precursorIds, std::vector<int>& fragmentIds);
    TimsBlockPL getTimsBlockPL(std::vector<int>& frameIds);
    std::vector<double> getGlobalMzAxis();

    // path to dataset, path to bruker binaries
    std::string datasetPath;
    std::string binaryPath;

private:
    class DataHandle;
    std::unique_ptr<DataHandle> handle_;
};

#endif //CPP_EXPOSEDTIMSDATAHANDLE_H
