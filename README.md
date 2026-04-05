# sablib

A C++ library for signal smoothing and baseline estimation, powered by Eigen.

## Features

### Baseline Estimation

- **Linear/Polynomial**: Robust linear or polynomial baseline fitting to anchor points.
- **Spline**: Cubic Spline interpolation-based baseline estimation.
- **SMA**: Iterative Simple Moving Average for baseline estimation.
- **SNIP**: Statistics-sensitive Non-linear Iterative Peak-clipping.
- **ModPoly**: Modified Polynomial.
- **IModPoly**: Improved ModPoly.
- **Backcor**: Iterative polynomial fitting with a non-quadratic cost function.
- **Goldindec**: Goldindec algorithm (a variation of Backcor).
- **AsLS**: Asymmetric Least Squares.
- **airPLS**: Adaptive Iteratively Reweighted Penalized Least Squares.
- **arPLS**: Asymmetrically Reweighted Penalized Least Squares.
- **psalsa**: Peaked Signal's Asymmetric Least Squares Algorithm.
- **BEADS**: Baseline Estimation And Denoising using Sparsity.

### Smoothing

- **Moving Average**: Standard simple moving average filter.
- **Weighted Moving Average**: Configurable weighted moving average filter.
- **Moving Median**: Moving median.
- **Savitzky-Golay**: Savitzky-Golay filter for simultaneous smoothing and differentiation.
- **Whittaker**: Whittaker smoother based on Penalized Least Squares.
- **P-Spline**: Penalized B-Spline smoothing.

## Requirements

- **C++ Standard**: C++17 or later.
- **Dependencies**: [Eigen](https://eigen.tuxfamily.org/) 3.4.0 or later.

## Installation & Build

### Using vcpkg

`sablib` provides a vcpkg port. First, clone this repository:

```bash
git clone https://github.com/Izadori/sablib
```

Then install with the `--overlay-ports` option:

```bash
vcpkg install sablib --overlay-ports=<path/to/sablib>/vcpkg/ports
```

In your `CMakeLists.txt`:

```cmake
find_package(sablib CONFIG REQUIRED)
target_link_libraries(my_app PRIVATE sablib::sablib)
```

### Using Conan

`sablib` provides a Conan recipe. First, clone this repository:

```bash
git clone https://github.com/Izadori/sablib
```

Then install and build:

```bash
cd sablib
conan install conan/sablib --output-folder=conan/build --build=missing
conan build conan/sablib --output-folder=conan/build
```

In your `CMakeLists.txt`:

```cmake
find_package(sablib CONFIG REQUIRED)
target_link_libraries(my_app PRIVATE sablib::sablib)
```

### Building manually

`sablib` can also be built as a static library using CMake directly:

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

If Eigen3 is not found, please try `cmake .. -DCMAKE_PREFIX_PATH="path/to/Eigen3"`.

## Usage

### Integrating with CMake

#### Via vcpkg

After installing with vcpkg, use `find_package` with the vcpkg toolchain:

```cmake
find_package(sablib CONFIG REQUIRED)
target_link_libraries(my_app PRIVATE sablib::sablib)
```

#### Via Conan

After installing with Conan, use `find_package` with the Conan toolchain:

```cmake
find_package(sablib CONFIG REQUIRED)
target_link_libraries(my_app PRIVATE sablib::sablib)
```

#### Manually

You can add `sablib` to your project as follows:

```cmake
# In your CMakeLists.txt
set(SABLIB_DIR "/path/to/sablib")
include_directories(${SABLIB_DIR})
find_library(SABLIB_LIB NAMES sablib libsablib PATHS "${SABLIB_DIR}/build")
target_link_libraries(my_app PRIVATE ${SABLIB_LIB} Eigen3::Eigen)
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

Please see the [GitHub Pages](https://izadori.github.io/sablib-docs/) for documentation.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
