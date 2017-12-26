cmake_minimum_required(VERSION 2.8.12)

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
  set("LIBHVL_PLATFORM_UNIX" "1")
  add_definitions("-Wall -Wextra -mmmx -msse -msse2 -fvisibility=hidden")

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
  set("LIBHVL_PLATFORM_WIN32" "1")
  set(hvl_mt_SRCS
      ${hvl_mt_SRCS}
      src/dllmain.cpp)

  if(USE_MPFR)
    add_definitions("-DLIBHVL_USE_MPFR")
    include_directories(
      ${LIBHVL_INCLUDES}
      ${GMP_INCLUDE_DIR}
      ${MPFR_INCLUDE_DIR})
    set(LINK_HIPREC ${LINK_HIPREC}
        ${GMP_LIBRARY_BIN})
    set(LINK_HIPREC ${LINK_HIPREC}
        ${MPFR_LIBRARY_BIN})
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

include_directories(${INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_BINARY_DIR})
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/include/hvl.h.in" "${CMAKE_CURRENT_BINARY_DIR}/hvl.h")

add_library(hvl_mt SHARED ${hvl_mt_SRCS})
target_link_libraries(hvl_mt PRIVATE ${LINK_HIPREC} ${LINK_PTHREAD})
set_target_properties(hvl_mt
                      PROPERTIES VERSION 2.0
                                 SOVERSION 2.0
				 PUBLIC_HEADER "${CMAKE_CURRENT_BINARY_DIR}/hvl.h")
if (WIN32)
  install(TARGETS hvl_mt
          RUNTIME DESTINATION lib
          PUBLIC_HEADER DESTINATION include)
else()
  install(TARGETS hvl_mt
          LIBRARY DESTINATION lib
          PUBLIC_HEADER DESTINATION include)
endif ()
