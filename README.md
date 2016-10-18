HVL_MT
===

Introduction
---

HVL_MT library can calculate values of Haaroff - van der Linde function and its derivatives by all parameters. The library can optionally make use of [MPFR library](http://www.mpfr.org/) to calculate intermediate results with arbitrary precision. If a continuous range of values is needed, HVL_MT takes advantage of multithreaded processing to speed up the computation. If multiple HVL parameters are to be evaluated for a given `x`, HVL_MT provides a set of `_prepared` functions that precalculate and reuse the terms common to all expressions.

Building
---

HVL_MT uses [CMake build system](https://cmake.org/). In order to build HVL_MT with MPFR support, MPFR library (and it's dependence [GMP library](https://gmplib.org/)) must be installed including the development files.

If you wish to use the MPFR library to perform calculations that exceed the IEEE754 `double` precision, supply addtional `-DUSE_MPFR=ON` parameter to `cmake`.

### Linux/UNIX

`cd` to the directory with HVL_MT source and run the following:

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make

**Remark:** Note that compatibility with UNIX systems other than Linux has not been tested.

### Windows

HVL_MT can be built on Windows using either the [MinGW toolchain](http://www.mingw.org/) or [Microsoft Visual Studio](https://www.visualstudio.com/). The following CMake variables have to be specified manually if you wish to use the MPFR library.

- `XGMP_INCLUDES` - Path to the directory with GMP public headers
- `MPFR_INCLUDES` - Path to the directory with MPFR public headers
- `XGMP_LIBRARIES` - Path to the directory with built GMP library
- `MPFR_LIBRARIES` - Path to the directory with built MPFR library

Assuming that the path to MinGW executables has been added to your PATH variable and you have `cd`'ed to the directory with HVL_MT source, a sample command might look like the example below. In the example GMP and MPFR are installed in `C:\gmp-bin` and `C:\mpfr-bin`, respectively.


    md build
    cd build
    cmake -G "MinGW Makefiles" .. -DXGMP_INCLUDES=c:/gmp-bin/include -MPFR_INCLUDES=c:/mpfr-bin/include -DXGMP_LIBRARIES=c:/gmp-bin/lib -DMPFR_LIBRARIES=c:/mpfr-bin/lib -DCMAKE_BUILD_TYPE=Release
    mingw32-make

**Remark 1:** In order to generate Microsoft Visual Studio solution instead, supply the appropriate value for the `-G` parameter. Please refer to CMake documentation for details.

**Remark 2:** If your MinGW toolchain uses Win32 threads instead of pthreads, supply `-DWIN32_THREADS=ON` parameter to `cmake`.

**Remark 3:**  If you intend to build HVL_MT with Microsoft Visual Studio and MPFR support, you might want to use the binary compable GMP-fork [MPIR](http://www.mpir.org/) as the mainline GMP is problematic to build with MSVS.

Licensing
---
The HVL_MT library is distributed under the terms of [The GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0.en.html).

Acknowledgments
---

HVL_MT makes use of MPFR C++ wrapper [mpreal](http://www.holoborodko.com/pavel/mpfr/) written by Pavel Holobrodko.
