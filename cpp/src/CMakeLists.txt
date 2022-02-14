cmake_minimum_required(VERSION 3.16)
project(src)

include_directories("../../include/opentims/opentims++")
include_directories("../../include/eigen")

find_library(LIB_TBB NAMES tbb REQUIRED)
find_library(SQLITE NAMES sqlite3 REQUIRED)
find_library(ZSTD NAMES zstd REQUIRED)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../../../python/")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_library(proteolizard SHARED Proteolizard.cpp)

target_link_libraries(proteolizard ${LIB_TBB})
target_link_libraries(proteolizard ${PYTHON_LIBRARIES})
target_link_libraries(proteolizard ${SQLITE})
target_link_libraries(proteolizard ${ZSTD})