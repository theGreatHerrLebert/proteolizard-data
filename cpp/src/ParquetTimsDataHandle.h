//
// Created by administrator on 25.05.23.
//

#ifndef PROTEOLIZARDDATA_PARQUETTIMSDATAHANDLE_H
#define PROTEOLIZARDDATA_PARQUETTIMSDATAHANDLE_H

#include <sstream>

#include "Frame.h"
#include "Slice.h"
#include "TimsBlock.h"

class ParquetTimsDataHandle {

public:
    ParquetTimsDataHandle(std::string dp);
    ~ParquetTimsDataHandle();

    TimsSlicePL getTimsSlicePL(int frameIdStart,
                               int frameIdEnd,
                               std::vector<int> &blockIds,
                               std::vector<int>& msMsTypes);

    TimsFramePL getTimsFramePL(int blockId,
                               int rowGroupId,
                               int rowGroupIndexStart,
                               int rowGroupIndexStop);

    // path to dataset
    const std::string rawDataPath;
    const std::string metaDataPath;

private:
    static std::string initialize(const std::string& value);
};

#endif //PROTEOLIZARDDATA_PARQUETTIMSDATAHANDLE_H