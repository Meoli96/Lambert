cmake_minimum_required(VERSION 3.16)

project(Lambert CXX)
set(CMAKE_CXX_STANDARD 20)
set(HOME_DIRECTORY $ENV{HOME})
set(CSPICE_PATH ${HOME_DIRECTORY}/Coding/cspice)

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



add_executable(
    ${PROJECT_NAME}
    ${LMB_SOURCES}
    ${LMB_HEADERS}
    
)


target_compile_options(
    ${PROJECT_NAME}
    PRIVATE
  
    -lm
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

include_directories(${CSPICE_PATH}/include)
set(CSPICE_LIB ${CSPICE_PATH}/lib/cspice.a)
target_link_libraries(Lambert ${CSPICE_LIB})