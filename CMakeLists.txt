cmake_minimum_required(VERSION 3.10)

project(sablib VERSION 0.3.0 LANGUAGES CXX)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)

# Source files
set(SABLIB_SOURCES
    baseline/airpls.cpp
    baseline/arpls.cpp
    baseline/asls.cpp
    baseline/backcor.cpp
    baseline/beads.cpp
    baseline/goldindec.cpp
    baseline/imodpoly.cpp
    baseline/modpoly.cpp
    baseline/polynomial.cpp
    baseline/psalsa.cpp
    baseline/sma.cpp
    baseline/snip.cpp
    baseline/spline.cpp
    smoothing/moving_average.cpp
    smoothing/moving_median.cpp
    smoothing/pspline.cpp
    smoothing/savitzky_golay.cpp
    smoothing/whittaker.cpp
)

# Header files
set(SABLIB_HEADERS
    baseline/airpls.h
    baseline/arpls.h
    baseline/asls.h
    baseline/backcor.h
    baseline/beads.h
    baseline/goldindec.h
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
    misc/polyfit.h
    misc/spdiags.h
    smoothing/moving_average.h
    smoothing/moving_median.h
    smoothing/pspline.h
    smoothing/savitzky_golay.h
    smoothing/whittaker.h
)

# Define the library
add_library(sablib ${SABLIB_SOURCES} ${SABLIB_HEADERS})

# Public include directories
target_include_directories(sablib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
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
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install headers (preserving folder structure)
foreach(header ${SABLIB_HEADERS})
    get_filename_component(header_dir ${header} DIRECTORY)
    install(FILES ${header} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sablib/${header_dir}")
endforeach()

# Install main header
install(FILES sablib.h DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sablib")

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

# Install config files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/sablibConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/sablibConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sablib
)
