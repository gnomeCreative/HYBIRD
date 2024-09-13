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
    ${PROJECT_SOURCE_DIR}/src/myvector.cpp
    ${PROJECT_SOURCE_DIR}/src/node.cpp
    ${PROJECT_SOURCE_DIR}/src/utils.cpp
    )

# Define the target we wish to build
add_executable(hybird ${HYBIRD_INCLUDE} ${HYBIRD_SRC})

# Use C++17
target_compile_features(hybird PUBLIC cxx_std_17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#Enable OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(hybird PRIVATE OpenMP::OpenMP_CXX)
    # Under CLion/MinGW, CMake finds OpenMP but doesn't fill the corresponding variables
    # So do it manually
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        if(OpenMP_CXX_LIBRARIES STREQUAL "" OR OpenMP_CXX_LIB_NAMES STREQUAL "")
            target_link_options(hybird PRIVATE "-fopenmp")
        endif()
    endif()
else()    
    message(WARNING "OpenMP C++ Support was not found!")
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