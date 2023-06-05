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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")

find_package(Eigen3 REQUIRED)
find_package(pybind11 REQUIRED)
find_package(SQLite3 REQUIRED)
find_package(TBB REQUIRED)
find_package(Arrow REQUIRED)
find_package(Parquet REQUIRED)
find_library(ZSTD NAMES zstd REQUIRED)

include_directories("../contrib/opentims/opentims++")
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${pybind11_INCLUDE_DIRS})
include_directories(${SQLITE3_INCLUDE_DIRS})

add_library(proteolizarddata SHARED
    src/Proteolizard.cpp
    src/VectorizedSpectrum.cpp
    src/Spectrum.cpp
    src/VectorizedFrame.cpp
    src/Frame.cpp
    src/Slice.cpp
    src/SliceVectorized.cpp
    src/TimsBlock.cpp
    src/TimsBlockVectorized.cpp
    src/ExposedTimsDataHandle.cpp
    src/ParquetTimsDataHandle.cpp
)

add_executable(ArrowParquetExample src/main.cpp
        src/VectorizedSpectrum.cpp
        src/Spectrum.cpp
        src/VectorizedFrame.cpp
        src/Frame.cpp
        src/Slice.cpp
        src/SliceVectorized.cpp
        src/TimsBlock.cpp
        src/TimsBlockVectorized.cpp
        src/ExposedTimsDataHandle.cpp
        src/ParquetTimsDataHandle.cpp
        )

set(ProteoLizardData_HEADER_FILES
    src/VectorizedSpectrum.h
    src/Spectrum.h
    src/VectorizedFrame.h
    src/Frame.h
    src/Slice.h
    src/SliceVectorized.h
    src/TimsBlock.h
    src/TimsBlockVectorized.h
    src/ExposedTimsDataHandle.h
    src/ParquetTimsDataHandle.h
)
set_target_properties(proteolizarddata PROPERTIES PUBLIC_HEADER "${ProteoLizardData_HEADER_FILES}")

target_link_libraries(proteolizarddata
    PUBLIC
        Eigen3::Eigen
        SQLite::SQLite3
        TBB::tbb
        Arrow::arrow_shared
        Parquet::parquet_shared
        ${ZSTD}
)

target_link_libraries(ArrowParquetExample
    PRIVATE
        Eigen3::Eigen
        SQLite::SQLite3
        TBB::tbb
        Arrow::arrow_shared
        Parquet::parquet_shared
        ${ZSTD}
)

include(GNUInstallDirs)
install(TARGETS proteolizarddata
    EXPORT                    ProteoLizardDataTargets
    LIBRARY DESTINATION       ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ProteoLizardData
)

install(EXPORT ProteoLizardDataTargets
    FILE        ProteoLizardDataTargets.cmake
    NAMESPACE   ProteoLizard::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ProteoLizardData
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("ProteoLizardDataConfigVersion.cmake"
    VERSION       ${ProteoLizardData_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES "cmake/ProteoLizardDataConfig.cmake" "${CMAKE_BINARY_DIR}/ProteoLizardDataConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ProteoLizardData
)
