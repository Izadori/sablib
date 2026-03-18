# sablib

A C++ library for signal smoothing and baseline estimation.

## Features

### Baseline Estimation

- **AsLS**: Asymmetric Least Squares.
- **airPLS**: Adaptive Iteratively Reweighted Penalized Least Squares.
- **arPLS**: Asymmetrically Reweighted Penalized Least Squares.
- **Linear/Polynomial**: Simple linear or polynomial baseline fitting to specified points.
- **SMA**: Baseline estimation using iterative Simple Moving Average.
- **Spline**: Baseline estimation using Cubic Spline interpolation.
- **BEADS**: Baseline Estimation And Denoising using Sparsity.

### Smoothing

- **Moving Average**: Simple moving average filter.
- **Weighted Moving Average**: Weighted moving average filter.
- **P-Spline**: Penalized B-Spline smoothing.
- **Savitzky-Golay**: Savitzky-Golay filter for data smoothing and differentiation.
- **Whittaker**: Whittaker smoother (Penalized Least Squares).

## Requirements

- **C++ Standard**: C++17 or later.
- **Library**: [Eigen](https://eigen.tuxfamily.org/) 3.4.0 or later.

## Build

This library is header-plus-source but intended to be built as a static library using CMake.

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

## Usage

### Integration with CMake

If you are using CMake, you can add `sablib` as a subdirectory:

```cmake
add_subdirectory(sablib)
target_link_libraries(your_project_target PUBLIC sablib)
```

Ensure Eigen3 is available and found in your `CMakeLists.txt` using `find_package(Eigen3 REQUIRED)`.

### Header Inclusion

You can include the main header to access all functionalities:

```cpp
#include "sablib/sablib.h"
```

Or include individual headers for specific algorithms:

```cpp
#include "sablib/smoothing/pspline.h"
#include "sablib/baseline/arpls.h"
```

Note: The include paths depend on how you set up your include directories. When using the CMake target `sablib`, the include paths are handled automatically, allowing you to use `#include "sablib/..."`.
