# __________________________________________________________________
#
# CMakeList for unit tests
# 
# __________________________________________________________________

cmake_minimum_required (VERSION 3.5)

message(" Test: field_class")

project (field_class LANGUAGES ${LANGUAGES})

set(EXE_NAME field_class)

include_directories(../../src/common)

add_executable(${EXE_NAME} ./main.cpp)

target_compile_options(${EXE_NAME} PRIVATE ${COMPILE_OPTIONS})
target_link_options(${EXE_NAME} PRIVATE ${LINK_OPTIONS})

# Target link libraries
include(${CMAKE_SOURCE_DIR}/cmake/link.cmake)


