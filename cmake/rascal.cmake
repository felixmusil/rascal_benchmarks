include(FetchContent)

# force new defaults for rascal
set(BUILD_BINDINGS ON CACHE INTERNAL "")
set(BUILD_EXAMPLES OFF CACHE INTERNAL "")

# set(CMAKE_BUILD_TYPE "Debug" CACHE INTERNAL "")

FetchContent_Declare(
 rascal
 GIT_REPOSITORY https://github.com/cosmo-epfl/librascal.git
 GIT_TAG        origin/feat/improve_integration_as_external
)

FetchContent_MakeAvailable(rascal)

# # This particular version of projD requires workarounds
# FetchContent_GetProperties(rascal)
# if(NOT rascal_POPULATED)
#   FetchContent_Populate(rascal)
#   add_subdirectory(${rascal_SOURCE_DIR} ${rascal_BINARY_DIR})
#   message(STATUS "rascal name: ${LIBRASCAL_NAME}")
#   # include_directories(${rascal_SOURCE_DIR}/src)
#   # include_directories(${EIGEN3_INCLUDE_DIR})
# endif()


