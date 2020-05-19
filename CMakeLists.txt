cmake_minimum_required(VERSION 3.11)

# https://cmake.org/pipermail/cmake/2008-September/023808.html
# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their
# definition
# and dump it in the cache along with proper documentation, otherwise set
# CMAKE_BUILD_TYPE
# to Debug prior to calling PROJECT()
#
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
ENDIF()
# Set the possible values of build type for cmake-gui
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
"Release" "Debug" "RelWithDebInfo" "MinSizeRel")

message(STATUS "Build type is: ${CMAKE_BUILD_TYPE}")

project(rascal_benchmark CXX)




set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(BUILD_SHARED_LIBS ON)

SET(TYPE_ARCHITECTURE "native" CACHE STRING
    "Choose the type of architecture to compile to in release mode."
)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")

    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

    # Add -march=${TYPE_ARCHITECTURE} (-march=native by default) to the release
    # build flags
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=${TYPE_ARCHITECTURE}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=${TYPE_ARCHITECTURE}")

    set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -march=${TYPE_ARCHITECTURE}")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=${TYPE_ARCHITECTURE}")
endif()



set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(rascal)


message(STATUS "rascal: ${LIBRASCAL_NAME}")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/scripts)

