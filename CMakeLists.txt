# Minimum required CMake version
# This specifies which default policies/settings are expected
cmake_minimum_required(VERSION 3.22)

# Set CMAKE_CUDA_ARCHITECTURES defaults, now CUDA is enabled.
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  # 61: Pascal (Alessandro's Office machine)
  # 70: Volta V100 (Bessemer, Bede)
  # 80: Ampere A100 (Stanage)
  # 90: Hopper H100 (Stanage, Bede)
  # We can specify multiple target CUDA architectures
  # It merely increases compilation time and the output binary size
  # (because all CUDA is compiled for each architecture)
  # Alternatively, you can specify the minimum required directly when configuring CMake
  set(CMAKE_CUDA_ARCHITECTURES "61;70;80;90" CACHE STRING "CUDA architectures" FORCE)
endif()

# Set the name for the project
# Can also specify VERSION, DESCRIPTION etc here if desired
project(HYBIRD LANGUAGES CXX)

# Installation target is not required
set(CMAKE_SKIP_INSTALL_RULES TRUE)

# Optional CUDA support
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
    option(HYBIRD_ENABLE_CUDA "If enabled CUDA support will be used" ON)
    if(HYBIRD_ENABLE_CUDA)
        enable_language(CUDA)
    endif()
else()
    message(FATAL_ERROR "CMake did not find a CUDA compiler, it is required for this branch of HYBIRD.") 
endif()

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
    ${PROJECT_SOURCE_DIR}/src/MeasureUnits.h
    ${PROJECT_SOURCE_DIR}/src/Problem.h
    ${PROJECT_SOURCE_DIR}/src/gpu/cuda_helper.h
    ${PROJECT_SOURCE_DIR}/src/gpu/cub_temp_mem.h
    ${PROJECT_SOURCE_DIR}/src/gpu/Cylinder2.h
    ${PROJECT_SOURCE_DIR}/src/gpu/Element2.h
    ${PROJECT_SOURCE_DIR}/src/gpu/IO2.h
    ${PROJECT_SOURCE_DIR}/src/gpu/LB2.h
    ${PROJECT_SOURCE_DIR}/src/gpu/LBParams.h
    ${PROJECT_SOURCE_DIR}/src/gpu/Node2.h
    ${PROJECT_SOURCE_DIR}/src/gpu/Object2.h
    ${PROJECT_SOURCE_DIR}/src/gpu/Particle2.h
    ${PROJECT_SOURCE_DIR}/src/gpu/ValidationIO.h
    ${PROJECT_SOURCE_DIR}/src/gpu/Wall2.h
    )
    
# Define source files
set(HYBIRD_SRC 
    ${PROJECT_SOURCE_DIR}/src/DEM.cpp
    ${PROJECT_SOURCE_DIR}/src/elmt.cpp
    # src/hybird.cu includes headers that require CUDA to parse correctly
    # with some work it could be resolved to build as .cpp
    ${PROJECT_SOURCE_DIR}/src/hybird.cu
    ${PROJECT_SOURCE_DIR}/src/IO.cpp
    ${PROJECT_SOURCE_DIR}/src/LB.cpp
    ${PROJECT_SOURCE_DIR}/src/node.cpp
    ${PROJECT_SOURCE_DIR}/src/utils.cpp
    ${PROJECT_SOURCE_DIR}/src/Problem.cpp
    # @todo How can we compile these as cpp when CUDA is disabled?
    # e.g. set_source_files_properties(test.cu PROPERTIES LANGUAGE CXX)
    ${PROJECT_SOURCE_DIR}/src/gpu/IO2.cu
    ${PROJECT_SOURCE_DIR}/src/gpu/LB2.cu
    ${PROJECT_SOURCE_DIR}/src/gpu/LBParams.cu
    ${PROJECT_SOURCE_DIR}/src/gpu/Node2.cu
    )

# Define the target we wish to build
add_executable(hybird ${HYBIRD_INCLUDE} ${HYBIRD_SRC})

# Allow includes relative to src root
target_include_directories(hybird PRIVATE "${PROJECT_SOURCE_DIR}/src")

# Link with dependencies
include(${CMAKE_CURRENT_LIST_DIR}/cmake/dependencies/exprtk.cmake)
target_link_libraries(hybird PRIVATE exprtk)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/dependencies/yaml-cpp.cmake)
target_link_libraries(hybird PRIVATE yaml-cpp)

# Use C++17
target_compile_features(hybird PUBLIC cxx_std_17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable OpenMP
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
                set_target_properties(OpenMP::OpenMP_CXX PROPERTIES INTERFACE_COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:C,CXX>:-openmp:experimental>")
                # Todo how to enable OpenMP in cuda files?/Is it desired/required?
                # CMake 3.30 alt: set(OpenMP_RUNTIME_MSVC "experimental")
            endif()
        endif()
    else()    
        message(WARNING "OpenMP C++ Support was not found!")
    endif()
else()
    message(STATUS "Build configured without OpenMP support.")
endif()

if(CMAKE_CUDA_COMPILER AND HYBIRD_ENABLE_CUDA)
    # Notify the build that we are building for cuda
    target_compile_definitions(hybird PRIVATE USE_CUDA)
    target_compile_options(hybird PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
    # Suppress deprecated target warning, we know Pascal that we normally want to build for is deprecated
    target_compile_options(hybird PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Wno-deprecated-gpu-targets>)     
    # Enable relocatable device code
    set_target_properties(hybird PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    # Make CUDA include/lib available to cpp compiler
    find_package(CUDAToolkit)
    target_link_libraries(hybird PRIVATE CUDA::cudart)
else()
    # Or set .cu files to be compiled as .cpp
    foreach(HYBIRD_SRC_FILE HYBIRD_SRC)
        get_filename_component(HYBIRD_SRC_FILE_EXT HYBIRD_SRC_FILE LAST_EXT)
        if(HYBIRD_SRC_FILE_EXT STREQUAL ".cu")
            set_source_files_properties(HYBIRD_SRC_FILE PROPERTIES LANGUAGE CXX)
        endif()
    endforeach()
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
    # Don't think this is required on Windows
    # set_target_properties(hybird PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
