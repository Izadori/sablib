# sablib

A C++ library for signal smoothing and baseline estimation, powered by Eigen.

## Features

### Baseline Estimation

- **Linear/Polynomial**: Robust linear or polynomial baseline fitting to anchor points.
- **Spline**: Cubic Spline interpolation-based baseline estimation.
- **SMA**: Iterative Simple Moving Average for baseline estimation.
- **AsLS**: Asymmetric Least Squares.
- **airPLS**: Adaptive Iteratively Reweighted Penalized Least Squares.
- **arPLS**: Asymmetrically Reweighted Penalized Least Squares.
- **BEADS**: Baseline Estimation And Denoising using Sparsity.

### Smoothing

- **Moving Average**: Standard simple moving average filter.
- **Weighted Moving Average**: Configurable weighted moving average filter.
- **Savitzky-Golay**: Savitzky-Golay filter for simultaneous smoothing and differentiation.
- **Whittaker**: Whittaker smoother based on Penalized Least Squares.
- **P-Spline**: Penalized B-Spline smoothing.

## Requirements

- **C++ Standard**: C++17 or later.
- **Dependencies**: [Eigen](https://eigen.tuxfamily.org/) 3.4.0 or later.

## Installation & Build

`sablib` is designed to be built as a static library using CMake.

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

If Eigen3 is not found, please try `cmake .. -DCMAKE_PREFIX_PATH="path/to/Eigen3"`.

## Usage

### Integrating with CMake

You can easily integrate `sablib` into your project by adding it as a subdirectory:

```cmake
# In your CMakeLists.txt
add_subdirectory(sablib)
target_link_libraries(your_project_target PUBLIC sablib)
```

Ensure Eigen3 is configured in your project (e.g., via `find_package(Eigen3 REQUIRED)`).

### Header Inclusion

Include the main header to access all library features:

```cpp
#include "sablib/sablib.h"
```

Alternatively, you can include specific algorithm headers to minimize dependencies:

```cpp
#include "sablib/smoothing/pspline.h"
#include "sablib/baseline/arpls.h"
```

## Documentation

API documentation can be generated using Doxygen.

```bash
doxygen Doxyfile
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
