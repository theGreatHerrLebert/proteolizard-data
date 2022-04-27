find_path(ProteoLizardData_INCLUDE_DIR NAMES ProteoLizardData/ExposedTimsDataHandle.h)
mark_as_advanced(ProteoLizardData_INCLUDE_DIR)

find_library(ProteoLizardData_LIBRARY NAMES proteolizard-data)
mark_as_advanced(ProteoLizardData_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ProteoLizardData DEFAULT_MSG ProteoLizardData_INCLUDE_DIR ProteoLizardData_LIBRARY)

if(ProteoLizardData_FOUND)
    set(ProteoLizardData_INCLUDE_DIRS ${ProteoLizardData_INCLUDE_DIR})
    set(ProteoLizardData_LIBRARIES    ${ProteoLizardData_LIBRARY})
    if(NOT TARGET ProteoLizard::proteolizard-data)
        include("${CMAKE_CURRENT_LIST_DIR}/ProteoLizardDataTargets.cmake")
    endif()
endif()
