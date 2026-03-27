cmake_minimum_required(VERSION 3.10)

project(sablib VERSION 0.1.0 LANGUAGES CXX)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)

# Include directories
set(SABLIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Source files
set(SABLIB_SOURCES
    baseline/airpls.cpp
    baseline/arpls.cpp
    baseline/asls.cpp
    baseline/beads.cpp
    baseline/imodpoly.cpp
    baseline/modpoly.cpp
    baseline/polynomial.cpp
    baseline/sma.cpp
    baseline/spline.cpp
    smoothing/moving_average.cpp
    smoothing/moving_median.cpp
    smoothing/pspline.cpp
    smoothing/savitzky_golay.cpp
    smoothing/whittaker.cpp
)

# Header files (for IDE support)
set(SABLIB_HEADERS
    baseline/airpls.h
    baseline/arpls.h
    baseline/asls.h
    baseline/beads.h
    baseline/imodpoly.h
    baseline/modpoly.h
    baseline/polynomial.h
    baseline/sma.h
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
target_include_directories(sablib PUBLIC ${SABLIB_INCLUDE_DIR})

# Link Eigen3
target_link_libraries(sablib PUBLIC Eigen3::Eigen)

# Install rules (Optional)
install(TARGETS sablib DESTINATION lib)
install(DIRECTORY baseline/ misc/ smoothing/ DESTINATION include/sablib FILES_MATCHING PATTERN "*.h")
