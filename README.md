# STL C++

Seasonal-trend decomposition for C++

[![Build Status](https://github.com/ankane/stl-cpp/actions/workflows/build.yml/badge.svg)](https://github.com/ankane/stl-cpp/actions)

## Installation

Add [the header](https://raw.githubusercontent.com/ankane/stl-cpp/v0.1.5/include/stl.hpp) to your project (supports C++17 and greater).

There is also support for CMake and FetchContent:

```cmake
include(FetchContent)

FetchContent_Declare(stl GIT_REPOSITORY https://github.com/ankane/stl-cpp.git GIT_TAG v0.1.5)
FetchContent_MakeAvailable(stl)

target_link_libraries(app PRIVATE stl::stl)
```

## Getting Started

Include the header

```cpp
#include "stl.hpp"
```

Decompose a time series

```cpp
std::vector<float> series = {
    5.0, 9.0, 2.0, 9.0, 0.0, 6.0, 3.0, 8.0, 5.0, 8.0,
    7.0, 8.0, 8.0, 0.0, 2.0, 5.0, 0.0, 5.0, 6.0, 7.0,
    3.0, 6.0, 1.0, 4.0, 4.0, 4.0, 3.0, 7.0, 5.0, 8.0
};
size_t period = 7; // period of seasonal component

auto res = stl::params().fit(series, period);
```

Get the components

```cpp
res.seasonal;
res.trend;
res.remainder;
```

## Robustness

Use robustness iterations

```cpp
auto res = stl::params().robust(true).fit(series, period);
```

Get robustness weights

```cpp
res.weights;
```

## Parameters

Set parameters

```cpp
stl::params()
    .seasonal_length(7)     // length of the seasonal smoother
    .trend_length(15)       // length of the trend smoother
    .low_pass_length(7)     // length of the low-pass filter
    .seasonal_degree(0)     // degree of locally-fitted polynomial in seasonal smoothing
    .trend_degree(1)        // degree of locally-fitted polynomial in trend smoothing
    .low_pass_degree(1)     // degree of locally-fitted polynomial in low-pass smoothing
    .seasonal_jump(1)       // skipping value for seasonal smoothing
    .trend_jump(2)          // skipping value for trend smoothing
    .low_pass_jump(1)       // skipping value for low-pass smoothing
    .inner_loops(2)         // number of loops for updating the seasonal and trend components
    .outer_loops(0)         // number of iterations of robust fitting
    .robust(false);         // if robustness iterations are to be used
```

## Strength

Get the seasonal strength

```cpp
res.seasonal_strength();
```

Get the trend strength

```cpp
res.trend_strength();
```

## Credits

This library was ported from the [Fortran implementation](https://www.netlib.org/a/stl).

## References

- [STL: A Seasonal-Trend Decomposition Procedure Based on Loess](https://www.scb.se/contentassets/ca21efb41fee47d293bbee5bf7be7fb3/stl-a-seasonal-trend-decomposition-procedure-based-on-loess.pdf)
- [Measuring strength of trend and seasonality](https://otexts.com/fpp2/seasonal-strength.html)

## History

View the [changelog](https://github.com/ankane/stl-cpp/blob/master/CHANGELOG.md)

## Contributing

Everyone is encouraged to help improve this project. Here are a few ways you can help:

- [Report bugs](https://github.com/ankane/stl-cpp/issues)
- Fix bugs and [submit pull requests](https://github.com/ankane/stl-cpp/pulls)
- Write, clarify, or fix documentation
- Suggest or add new features

To get started with development:

```sh
git clone https://github.com/ankane/stl-cpp.git
cd stl-cpp
g++ -std=c++17 -Wall -Wextra -Werror -o test/main test/main.cpp
test/main
```
