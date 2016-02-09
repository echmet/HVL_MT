cmake_minimum_required(VERSION 2.6)

project(hvl_mt)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
add_definitions("-Wall -Wextra -mmmx -msse -msse2 -std=gnu++11")

set(hvl_mt_SRCS
    hvl.cpp
    hvl_dll.cpp)

if(UNIX)
  add_definitions("-DECHMET_MATH_HVL_PLATFORM_UNIX")
  find_package(MPFR REQUIRED)
  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MPFR_INCLUDES})
  set(LINK_HIPREC
      ${MPFR_LIBRARIES})
elseif(MINGW)
  add_definitions("-DECHMET_MATH_HVL_PLATFORM_WIN -DECHMET_MATH_HVL_MINGW32")
  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
	${XGMP_INCLUDES}
    ${MPFR_INCLUDES})
  set(hvl_mt_SRCS
      ${hvl_mt_SRCS}
      dllmain.cpp)
  find_library(GMPLIB
               gmp
	           ${XGMP_LIBRARIES})
  find_library(MPFRLIB
               mpfr
	           ${MPFR_LIBRARIES})
  set(LINK_HIPREC ${MPFRLIB} ${GMPLIB})
elseif(WIN32)
  add_definitions("-DECHMET_MATH_HVL_PLATFORM_WIN")
  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
	${XGMP_INCLUDES}
    ${MPFR_INCLUDES})
  set(hvl_mt_SRCS
      ${hvl_mt_SRCS}
      dllmain.cpp)
  find_library(GMPLIB
               gmp
	           ${XGMP_LIBRARIES})
  find_library(MPFRLIB
               mpfr
	           ${MPFR_LIBRARIES})
  set(LINK_HIPREC ${MPFRLIB} ${GMPLIB})
else()
  message(FATAL_ERROR "Platform not recognized")
endif()
  
add_library(hvl_mt SHARED ${hvl_mt_SRCS})
target_link_libraries(hvl_mt PRIVATE ${LINK_HIPREC})
set_target_properties(hvl_mt
                      PROPERTIES PREFIX ""
                                 VERSION 1.0
                                 SOVERSION 1.0)
