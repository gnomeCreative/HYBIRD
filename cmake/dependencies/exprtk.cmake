############################################
# Mathematical Expression Toolkit (exprtk) #
############################################

# Include FetchContent module
include(FetchContent)

# Declare ExprTk dependency
FetchContent_Declare(
    exprtk
    GIT_REPOSITORY https://github.com/ArashPartow/exprtk.git
    GIT_TAG 0.0.3  # Specific commit for stability 0.0.3 9910dc988d472d17ecf309c1c9c3e38430bd9c0f
)

# Make ExprTk available
FetchContent_MakeAvailable(exprtk)

# Create an interface library target for ExprTk
add_library(exprtk INTERFACE)
target_include_directories(exprtk INTERFACE ${exprtk_SOURCE_DIR})
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(exprtk INTERFACE $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)
    #target_compile_options(exprtk INTERFACE $<$<COMPILE_LANGUAGE:CXX>:/W4>)
    #target_compile_options(exprtk INTERFACE $<$<COMPILE_LANGUAGE:CXX>:/WX>)
    target_compile_options(exprtk INTERFACE $<$<COMPILE_LANGUAGE:CXX>:/wd4334>)
    target_compile_options(exprtk INTERFACE $<$<COMPILE_LANGUAGE:CXX>:/wd4702>)

    target_compile_definitions(exprtk INTERFACE _SECURE_SCL=0            )
    target_compile_definitions(exprtk INTERFACE _CRT_SECURE_NO_WARNINGS  )
    target_compile_definitions(exprtk INTERFACE _SCL_SECURE_NO_WARNINGS  )
    #target_compile_definitions(exprtk INTERFACE _HAS_ITERATOR_DEBUGGING=0)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_SYSTEM_NAME STREQUAL "Windows")
    target_compile_options(exprtk INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wa,-mbig-obj>)
endif()
