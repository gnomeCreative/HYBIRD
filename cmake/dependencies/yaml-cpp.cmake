############
# yaml-cpp #
############

# Include FetchContent module
include(FetchContent)

# Declare yaml-cpp dependency
FetchContent_Declare(
    yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG 39f737443b05e4135e697cb91c2b7b18095acd53  # Specific commit for stability
)

# Configure yaml-cpp build options
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "disable yaml-cpp tests" FORCE)
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "disable yaml-cpp tools" FORCE)
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "disable yaml-cpp contrib" FORCE)
set(YAML_CPP_INSTALL OFF CACHE BOOL "disable yaml-cpp install targets" FORCE)

# Make yaml-cpp available
FetchContent_MakeAvailable(yaml-cpp)
