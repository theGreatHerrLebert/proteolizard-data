cmake_minimum_required(VERSION 3.16)

project(ProteoLizardData
    VERSION 0.1
    DESCRIPTION  "Integration of timsTOF raw data with Python"
    HOMEPAGE_URL "https://github.com/theGreatHerrLebert/proteolizard-data"
    LANGUAGES    CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

find_package(Eigen3 REQUIRED)
find_package(pybind11 REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Development)
find_package(SQLite3 REQUIRED)
find_package(TBB REQUIRED)
find_library(ZSTD NAMES zstd REQUIRED)

include_directories("../../include/opentims/opentims++")
include_directories(${Python3_INCLUDE_DIRS})
include_directories(${SQLITE3_INCLUDE_DIRS})
include_directories(${pybind11_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})

add_library(proteolizard SHARED
    Proteolizard.cpp
    VectorizedSpectrum.cpp
    Spectrum.cpp
    VectorizedFrame.cpp
    Frame.cpp
    Slice.cpp
    ExposedTimsDataHandle.cpp
)

target_link_libraries(proteolizard TBB::tbb)
target_link_libraries(proteolizard ${PYTHON_LIBRARIES})
link_directories("../../include/opentims/opentims++/")
target_link_libraries(proteolizard SQLite::SQLite3)
target_link_libraries(proteolizard ${ZSTD})
target_link_libraries(proteolizard Eigen3::Eigen)
