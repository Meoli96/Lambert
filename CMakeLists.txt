cmake_minimum_required(VERSION 3.16)

project(Lambert CXX)
set(CMAKE_CXX_STANDARD 20)


set(INC_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(SRC_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
set(LMB_SOURCES
    # Source file 
    ${SRC_DIR}/main.cpp

 
)


set(LMB_HEADERS
    # Headers file 
    ${INC_DIR}/Lambert.hpp
    ${INC_DIR}/Utils.hpp
    ${INC_DIR}/ThreadPool.hpp

)


set(CMAKE_CXX_FLAGS_PROFILING "-O2 -g -pg" CACHE STRING
    "Flags used by the C++ compiler during Profiling builds."
    FORCE )
set(CMAKE_C_FLAGS_PROFILING "-O2 -g -pg" CACHE STRING
    "Flags used by the C compiler during Profiling builds."
    FORCE )
set(CMAKE_EXE_LINKER_FLAGS_PROFILING
    "-pg" CACHE STRING
    "Flags used for linking binaries during Profiling builds."
    FORCE )
set(CMAKE_SHARED_LINKER_FLAGS_PROFILING
    "-pg" CACHE STRING
    "Flags used by the shared libraries linker during Profiling builds."
    FORCE )
mark_as_advanced(
    CMAKE_CXX_FLAGS_PROFILING
    CMAKE_C_FLAGS_PROFILING
    CMAKE_EXE_LINKER_FLAGS_PROFILING
    CMAKE_SHARED_LINKER_FLAGS_PROFILING
)

# Update the documentation string of CMAKE_BUILD_TYPE for GUIs
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel Profiling."
    FORCE )
    

add_executable(
    ${PROJECT_NAME}
    ${LMB_SOURCES}
    ${LMB_HEADERS}
)


target_compile_options(
    ${PROJECT_NAME}
    PRIVATE
    -Wall
    -Wextra
    -Wpedantic
    -O0
)

target_include_directories(
    ${PROJECT_NAME}
    PRIVATE
    ${INC_DIR}
)


include_directories(/usr/local/include/IO)
find_library(IO_Astrodynamics_LIB NAMES libIO.Astrodynamics.so)
target_link_libraries(Lambert ${IO_Astrodynamics_LIB})