cmake_minimum_required(VERSION 3.16)
project(src)

include_directories("../../include/opentims/opentims++")
include_directories("../../include/eigen")
include_directories("../../include/pybind11/include")

find_package(TBB REQUIRED)
find_package(SQLite3 REQUIRED)
find_library(ZSTD NAMES zstd REQUIRED)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${SQLITE3_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../../../python/")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_library(proteolizard SHARED Proteolizard.cpp Hashing.h)

target_link_libraries(proteolizard TBB::tbb)
target_link_libraries(proteolizard ${PYTHON_LIBRARIES})
link_directories("../../include/opentims/opentims++/")
target_link_libraries(proteolizard SQLite::SQLite3)
target_link_libraries(proteolizard ${ZSTD})
