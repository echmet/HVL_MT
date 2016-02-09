HVL_MT
===

Introduction
---

HVL_MT library can calculate values of Haaroff - van der Linde function including its derivatives by all its parameters. The library makes use of [MPFR library](http://www.mpfr.org/) to calculate intermediate results with arbitrary precision. If a continuous range of values is needed, HVL_MT takes advantage of multithreaded processing to speed up the computation.

Building
---

HVL_MT uses [CMake build system](https://cmake.org/). Note that HVL_MT employs arbitrary-precision library MPFR. To build HVL_MT, MPFR library (and it's dependence [GMP library](https://gmplib.org/)) must be installed including the development files.

### Linux/UNIX

cd to the directory with HVL_MT source and run the following:

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make

Remark: Note that compatibility with UNIX systems other than Linux has not been tested.

### Windows

HVL_MT can be built on Windows using the [MinGW toolchain](http://www.mingw.org/). Additionally, these CMake variables have to be specified manually

- XGMP_INCLUDES - Path to the directory with GMP public headers
- MPFR_INCLUDES - Path to the directory with MPFR public headers
- XGMP_LIBRARIES - Path to the directory with built GMP library
- MPFR_LIBRARIES - Path to the directory with build MPFR library

Assuming that the path to MinGW execuables has been added to your PATH variable and you have cd'ed to the directory with HVL_MT source, a sample command might look like this:


    md build
    cd build
    cmake -G "MinGW Makefiles" .. -DXGMP_INCLUDES=c:/gmp-bin/include -MPFR_INCLUDES=c:/mpfr-bin/include -DXGMP_LIBRARIES=c:/gmp-bin/lib -DMPFR_LIBRARIES=c:/mpfr-bin/lib -DCMAKE_BUILD_TYPE=Release
    mingw32-make

Remark: While it is possible to build HVL_MT with Microsoft Visual Studio and appropriate GMP and MPFR builds, doing so is not recommended by the developers of HVL_MT library.

Licensing
---
The HVL_MT library is released under the terms of [The GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0.en.html).

Acknowledgments
---

HVL_MT makes use of MPFR C++ wrapper [mpreal](http://www.holoborodko.com/pavel/mpfr/) written by Pavel Holobrodko.
