cmake_minimum_required(VERSION 3.10)

project(sablib VERSION 0.4.0 LANGUAGES CXX)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3
find_package(Eigen3 CONFIG REQUIRED)

# Source files
set(SABLIB_SOURCES
    baseline/airpls.cpp
    baseline/aispls.cpp
    baseline/arpls.cpp
    baseline/asls.cpp
    baseline/aspls.cpp
    baseline/backcor.cpp
    baseline/beads.cpp
    baseline/cornercutting.cpp
    baseline/goldindec.cpp
    baseline/iasls.cpp
    baseline/imodpoly.cpp
    baseline/kajfoszkwiatek.cpp
    baseline/modpoly.cpp
    baseline/polynomial.cpp
    baseline/psalsa.cpp
    baseline/rubberband.cpp
    baseline/sma.cpp
    baseline/snip.cpp
    baseline/spline.cpp
    misc/functions.cpp
    smoothing/moving_average.cpp
    smoothing/moving_median.cpp
    smoothing/pspline.cpp
    smoothing/savitzky_golay.cpp
    smoothing/whittaker.cpp
    c/data.cpp
    c/baseline.cpp
    c/smoothing.cpp
)

# Header files
set(SABLIB_HEADERS
    baseline/airpls.h
    baseline/aispls.h
    baseline/arpls.h
    baseline/asls.h
    baseline/aspls.h
    baseline/backcor.h
    baseline/beads.h
    baseline/goldindec.h
    baseline/iasls.h
    baseline/imodpoly.h
    baseline/modpoly.h
    baseline/polynomial.h
    baseline/psalsa.h
    baseline/sma.h
    baseline/snip.h
    baseline/spline.h
    misc/bspline.h
    misc/convolve.h
    misc/cubic_spline.h
    misc/diff.h
    misc/expand.h
    misc/functions.h
    misc/polyfit.h
    misc/spdiags.h
    smoothing/moving_average.h
    smoothing/moving_median.h
    smoothing/pspline.h
    smoothing/savitzky_golay.h
    smoothing/whittaker.h
    c/data.h
    c/baseline.h
    c/smoothing.h
)

# Define the library
# Use a variable to decide between STATIC and SHARED (default to STATIC)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
add_library(sablib ${SABLIB_SOURCES} ${SABLIB_HEADERS})

# Generate export header for shared library support
include(GenerateExportHeader)
generate_export_header(sablib
    BASE_NAME SABLIB
    EXPORT_FILE_NAME "${CMAKE_CURRENT_BINARY_DIR}/sablib_export.h"
)

# Public include directories
target_include_directories(sablib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include/sablib>
)

# Link Eigen3
target_link_libraries(sablib PUBLIC Eigen3::Eigen)

# Installation logic
include(GNUInstallDirs)

install(TARGETS sablib
    EXPORT sablibTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install headers (preserving folder structure)
install(DIRECTORY
    baseline
    smoothing
    misc
    c
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sablib
)

# Install main header
install(FILES
    sablib.h
    csablib.h
    "${CMAKE_CURRENT_BINARY_DIR}/sablib_export.h"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sablib"
)

# Export targets for find_package
install(EXPORT sablibTargets
    FILE sablibTargets.cmake
    NAMESPACE sablib::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sablib
)

# Create sablibConfig.cmake
include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/sablibConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/sablibConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sablib
)

# Write version file
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/sablibConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

set_target_properties(sablib PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# Install config files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/sablibConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/sablibConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sablib
)
