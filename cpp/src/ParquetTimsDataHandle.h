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
    ParquetTimsDataHandle(const std::string& dp);
    ~ParquetTimsDataHandle();

    TimsSlicePL getBlock(int blockId, std::vector<int>& msMsTypes);
    TimsFramePL getTimsFramePL(int blockId, int indexStart, int indexEnd);

    // path to dataset
    std::string datasetPath;
};

#endif //PROTEOLIZARDDATA_PARQUETTIMSDATAHANDLE_H