# Minimum required CMake version
# This specifies which default policies/settings are expected
cmake_minimum_required(VERSION 3.22)
# Set the name for the project
# Can also specify VERSION, DESCRIPTION etc here if desired
project(HYBIRD LANGUAGES CXX)

# Installation target is not required
set(CMAKE_SKIP_INSTALL_RULES TRUE)

# Define headers
set(HYBIRD_INCLUDE 
    ${PROJECT_SOURCE_DIR}/src/DEM.h
    ${PROJECT_SOURCE_DIR}/src/elmt.h
    ${PROJECT_SOURCE_DIR}/src/getpot.h
    ${PROJECT_SOURCE_DIR}/src/IO.h
    ${PROJECT_SOURCE_DIR}/src/lattice.h
    ${PROJECT_SOURCE_DIR}/src/LB.h
    ${PROJECT_SOURCE_DIR}/src/macros.h
    ${PROJECT_SOURCE_DIR}/src/myvector.h
    ${PROJECT_SOURCE_DIR}/src/myvector.inl
    ${PROJECT_SOURCE_DIR}/src/node.h
    ${PROJECT_SOURCE_DIR}/src/utils.h
    )
    
# Define source files
set(HYBIRD_SRC 
    ${PROJECT_SOURCE_DIR}/src/DEM.cpp
    ${PROJECT_SOURCE_DIR}/src/elmt.cpp
    ${PROJECT_SOURCE_DIR}/src/hybird.cpp
    ${PROJECT_SOURCE_DIR}/src/IO.cpp
    ${PROJECT_SOURCE_DIR}/src/LB.cpp
    ${PROJECT_SOURCE_DIR}/src/node.cpp
    ${PROJECT_SOURCE_DIR}/src/utils.cpp
    )

# Define the target we wish to build
add_executable(hybird ${HYBIRD_INCLUDE} ${HYBIRD_SRC})

# Use C++17
target_compile_features(hybird PUBLIC cxx_std_17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#Enable OpenMP
option(HYBIRD_ENABLE_OPENMP "If Off OpenMP support will not be enabled" ON)
if(HYBIRD_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(hybird PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(hybird PRIVATE USE_OPENMP)
        # Under CLion/MinGW, CMake finds OpenMP but doesn't fill the corresponding variables
        # So do it manually
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            if(OpenMP_CXX_LIBRARIES STREQUAL "" OR OpenMP_CXX_LIB_NAMES STREQUAL "")
                target_link_options(hybird PRIVATE "-fopenmp")
            endif()
        elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
            if(MSVC_TOOLSET_VERSION GREATER_EQUAL 142)
                # Visual Studio 2019+ has experimental flag for OMP3.0 support
                set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "-openmp:experimental")
                # CMake 3.30 alt: set(OpenMP_RUNTIME_MSVC "experimental")
            endif()
        endif()
    else()    
        message(WARNING "OpenMP C++ Support was not found!")
    endif()
else()
    message(STATUS "Build configured without OpenMP support.")
endif()

# Default to release builds
get_property(GENERATOR_IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT GENERATOR_IS_MULTI_CONFIG)
    if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
        set(default_build_type "Release")
        message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
        set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING 
            "Choose the type of build, options are: Release, Debug, RelWithDebInfo, MinSizeRel or leave the value empty for the default." FORCE)
    endif()
endif()
# Ensure debug preprocessor macro is set
target_compile_definitions(hybird PRIVATE "$<$<CONFIG:Debug>:_DEBUG>")

# Update the output directories, so that they include Debug/Release
set_target_properties(hybird
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/$<CONFIG>"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/$<CONFIG>"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/$<CONFIG>"
)

# Additional compiler options
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # Extra debug info
    target_compile_options(hybird PRIVATE "$<$<CONFIG:Debug>:-ggdb3>")
    # Stricter warnings
    target_compile_options(hybird PRIVATE "$<$<COMPILE_LANGUAGE:C,CXX>:-Wall>")
    # Suppress: comparison of integer expressions of different signedness
    target_compile_options(hybird PRIVATE "$<$<COMPILE_LANGUAGE:C,CXX>:-Wno-sign-compare>")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # Stricter warnings
    target_compile_options(hybird PRIVATE "$<$<COMPILE_LANGUAGE:C,CXX>:/W4>")
endif()

# Static stdlib (required by Intel VTune)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    option(HYBIRD_STATIC_STDLIB "If enabled link against static stdlib, only necessary for Intel VTune profiler." OFF)
    # Required for intel vtune
    if (HYBIRD_STATIC_STDLIB)
        target_link_options(hybird PRIVATE -static-libgcc -static-libstdc++)
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # todo, test the correct configuration for MSVC 
    # set_property(TARGET hybird PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
