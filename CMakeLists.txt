cmake_minimum_required(VERSION 2.6)

project(libhvl_mt)

option(WIN32_THREADS "Use Win32 threads instead of pthreads" OFF)
option(USE_CPP17 "Take advantage of C++17 features" OFF)
option(USE_MPFR "Use MPFR high precision library" OFF)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
add_definitions("-DLIBHVL_DLLBUILD")

set(LIBHVL_INCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(hvl_mt_SRCS
    src/hvl.cpp)

if(UNIX)
  add_definitions("-Wall -Wextra -mmmx -msse -msse2 -fvisibility=hidden -DLIBHVL_PLATFORM_UNIX")

  if(USE_MPFR)
    find_package(MPFR REQUIRED)
    add_definitions("-DLIBHVL_USE_MPFR")
    include_directories(
      ${LIBHVL_INCLUDES}
      ${MPFR_INCLUDES})
    set(LINK_HIPREC
        ${MPFR_LIBRARIES})
  else()
    include_directories(
      ${LIBHVL_INCLUDES})
    set(LINK_HIPREC "")
  endif()

  if(USE_CPP17)
    add_definitions("-std=gnu++1z -DLIBHVL_USE_CPP17")
  else()
    add_definitions("-std=gnu++11")
  endif()

  set(LINK_PTHREAD "pthread")
elseif(WIN32)
  add_definitions("-DLIBHVL_PLATFORM_WIN32")

  set(hvl_mt_SRCS
      ${hvl_mt_SRCS}
      src/dllmain.cpp)

  if(USE_MPFR)
    add_definitions("-DLIBHVL_USE_MPFR")
    include_directories(
      ${LIBHVL_INCLUDES}
      ${GMP_INCLUDES}
      ${MPFR_INCLUDES})
    find_library(GMPLIB
                 gmp
                 ${GMP_LIBRARIES})
    find_library(MPFRLIB
                 mpfr
	             ${MPFR_LIBRARIES})
    set(LINK_HIPREC ${MPFRLIB} ${GMPLIB})
  else()
    include_directories(
      ${LIBHVL_INCLUDES})
    set(LINK_HIPREC "")
  endif()

  if(MSVC)
    add_definitions("-DLIBHVL_COMPILER_MSVC -DLIBHVL_THREADING_WIN32")
  elseif(MINGW)
    add_definitions("-Wall -Wextra -mmmx -msse -msse2 -fvisibility=hidden -DLIBHVL_COMPILER_MINGW")
    if (WIN32_THREADS)
      add_definitions("-DLIBHVL_THREADING_WIN32")
      set(LINK_PTHREAD "")
    else()
      add_definitions("-DLIBHVL_THREADING_PTHREAD")
      set(LINK_PTHREAD "pthread")
    endif()

    if(USE_CPP17)
      add_definitions("-std=gnu++1z -DLIBHVL_USE_CPP17")
    else()
      add_definitions("-std=gnu++11")
    endif()
  elseif(MSYS)
    add_definitions("-Wall -Wextra -mmmx -msse -msse2 -fvisibility=hidden -DLIBHVL_COMPILER_MSYS")
    if (WIN32_THREADS)
      add_definitions("-DLIBHVL_THREADING_WIN32")
      set(LINK_PTHREAD "")
    else()
      add_definitions("-DLIBHVL_THREADING_PTHREAD")
      set(LINK_PTHREAD "pthread")
    endif()

    if(USE_CPP17)
      add_definitions("-std=gnu++1z -DLIBHVL_USE_CPP17")
    else()
      add_definitions("-std=gnu++11")
    endif()
  else()
    message(FATAL_ERROR "Compiler not recognized or not supported")
  endif()
else()
  message(FATAL_ERROR "Platform not recognized or not supported")
endif()
  
add_library(libhvl_mt SHARED ${hvl_mt_SRCS})
target_link_libraries(libhvl_mt PRIVATE ${LINK_HIPREC} ${LINK_PTHREAD})
set_target_properties(libhvl_mt
                      PROPERTIES PREFIX ""
                                 VERSION 2.0
                                 SOVERSION 2.0)
