include(FetchContent)

# force new defaults for rascal
set(BUILD_BINDINGS OFF CACHE INTERNAL "")
set(BUILD_EXAMPLES OFF CACHE INTERNAL "")

FetchContent_Declare(
 rascal
 GIT_REPOSITORY https://github.com/cosmo-epfl/librascal.git
 GIT_TAG        origin/master
)

FetchContent_MakeAvailable(rascal)
