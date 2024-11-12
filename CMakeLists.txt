# __________________________________________________________________
#
# CMakeList for Mini-PIC
# 
# __________________________________________________________________

cmake_minimum_required (VERSION 3.5)

message(" __________________________________________________________________ \n")
message(" CMakeList for Mini-PIC")
message(" __________________________________________________________________ \n")

# Disable in-source builds to prevent source tree corruption.
if( "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}" )
 message( FATAL_ERROR "FATAL: In-source builds are not allowed. You should create a separate directory for build files and delete CMakeCache.txt." )
endif()

# __________________________________________________________________
# Options

# Debug, Release, RelWithDebInfo and MinSizeRel
# set(CMAKE_BUILD_TYPE Debug)
message(" Build type: ${CMAKE_BUILD_TYPE}")

# Backend
set(BACKEND "sequential" CACHE STRING "Backend to use (sequential (default), openmp, openmp_task, openmp_target, kokkos, thrust, eventify)")

# Debug and test modes
option(DEBUG "Debug messages" OFF)
option(TEST "Simulation test for CI" OFF)
option(UNIT_TESTS "Compile Unit tests" OFF)
option(MINIPIC "Compile minipic (ON/OFF)" ON)
option(SHAMAN "Compile with Shaman" OFF)
option(PROFILER "Task profiler (ON/OFF)" OFF)
option(EVENTIFY "Build eventify if needed" OFF)


message(" Minipic compilation: ${MINIPIC}")
message(" Unit tests: ${UNIT_TESTS}")

# Check that the backend is supported
if (BACKEND STREQUAL "sequential")
  message(" Backend: sequential")
  set(EXE_NAME minipic.seq)
elseif (BACKEND STREQUAL "openmp")
  message(" Backend: openmp")
  set(EXE_NAME minipic.omp)
elseif (BACKEND STREQUAL "openmp_task")
  message(" Backend: openmp_task")
  set(EXE_NAME minipic.omp_task)
elseif (BACKEND STREQUAL "eventify")
  message(" Backend: eventify")
  set(EXE_NAME minipic.eventify)
else()
  message(FATAL_ERROR "FATAL: Backend (${BACKEND}) not supported")
endif()

message(" Executable name: ${EXE_NAME}")

# __________________________________________________________________
# Compiler type

message(" Compiler: ${CMAKE_CXX_COMPILER_ID}")

# _________________________________________________________________
# Debug mode

if (DEBUG)
 message(" Debug mode: ON")
endif()

# __________________________________________________________________
# Determine languages

set(LANGUAGES CXX)
message(" Languages: ${LANGUAGES}")

message(" __________________________________________________________________ \n")

project (minipic LANGUAGES ${LANGUAGES})

# Compiler
#set(CMAKE_CXX_COMPILER clang++)
#set(CMAKE_CXX_COMPILER icpx)
#set(CMAKE_CXX_COMPILER icpc)
#set(CMAKE_CXX_COMPILER g++)

# message("${CMAKE_CURRENT_SOURCE_DIR}")
# message("${CMAKE_SOURCE_DIR}")

# Specific behavior for test mode
if (TEST)
  set(MAIN_PATH "./main.cpp")
else()
 set(MAIN_PATH "./src/main.cpp")
endif()

set(COMPILE_OPTIONS "")
set(LINK_OPTIONS "")

# _________________________________________________________________
# Profiler task

if (PROFILER)
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_PROFILER__)
endif()

# __________________________________________________________________
# Sources and compile options

# -----> OpenMP mode <------
if (BACKEND STREQUAL "openmp")

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -fopenmp)
 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++17)
 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_OMP__ -D__MINIPIC_SIMD__)
 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} ${ARCH_OPT_FLAGS})

 include_directories(src/setups src/openmp src/common)

 # list of source files - automatic search cpp files in src/openmp and src/common directory
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/openmp/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

 find_package(OpenMP REQUIRED)

# -----> OpenMP Task mode base version <------
elseif (BACKEND STREQUAL "openmp_task")

 set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fopenmp)

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++20 -D__MINIPIC_OMP_TASK__ -D__MINIPIC_SIMD__)

 include_directories(src/setups src/openmp_task src/common)

 # list of source files - automatic search cpp files in src/openmp and src/common directory
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/openmp_task/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

 find_package(OpenMP REQUIRED)

  # list of source files - automatic search cpp files in src/openmp and src/common directory
  include_directories(src/setups src/openmp_task src/common)
  file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/openmp_task/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")
  

# -----> Eventify Task mode base version <------
elseif (BACKEND STREQUAL "eventify")
 
 set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fopenmp)
 include_directories(src/setups src/eventify src/common)

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++20 -D__MINIPIC_EVENTIFY__ -D__MINIPIC_SIMD__ ${ARCH_OPT_FLAGS})

 # list of source files - automatic search cpp files in src/openmp and src/common directory
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/eventify/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

 find_package(OpenMP REQUIRED)
 find_package (Eventify REQUIRED)

# ----> Serial mode <----
else()

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++17 -D__MINIPIC_SIMD__)

 include_directories(src/setups src/openmp src/common)

 # list of source files - automatic search cpp files in src/openmp and src/common directory
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/openmp/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

endif()

# -----> Shaman <------
if (SHAMAN)
 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__SHAMAN__)
 find_package(Shaman REQUIRED)
endif()

# -----> Debug mode <------
if (DEBUG)
 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_DEBUG__)
endif()

# __________________________________________________________________
# Create executable with options

if (MINIPIC)

  add_executable(${EXE_NAME} ${SOURCES} ${MAIN_PATH})
  target_compile_options(${EXE_NAME} PRIVATE ${COMPILE_OPTIONS})
  target_link_options(${EXE_NAME} PRIVATE ${LINK_OPTIONS})


  #Call for eventify
  if (EVENTIFY)
  target_link_libraries(${EXE_NAME} PRIVATE jsc::eventify)
  endif()

  # Target link libraries
  include(cmake/link.cmake)

endif()

# ----> Shaman <----
if (SHAMAN)

 target_link_libraries(${EXE_NAME} shaman::shaman)

 message(" Add Shaman for numerical precision analysis")

endif()

# tests

if (UNIT_TESTS)

 add_subdirectory(unit_tests/)

endif()


