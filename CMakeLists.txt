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
#set(CMAKE_BUILD_TYPE Release)
message(" Build type: ${CMAKE_BUILD_TYPE}")

# Backend
set(BACKEND "sequential" CACHE STRING "Backend to use (sequential (default), openmp, openmp_task, openmp_target, kokkos, kokkos_views,thrust, thrust_unified, eventify)")

# Atomic type for projection using Kokkos
set(KOKKOS_PROJECTION_TYPE "scatter_view" CACHE STRING "Type of projection for Kokkos (scatter_view (default), atomic_view)")

# iterator type using Thrust
set(THRUST_ITERATOR_TYPE "zip" CACHE STRING "Type of iterator for Thrust (zip (default), counting)")

if (BACKEND STREQUAL "kokkos")
  message(" Kokkos projection type: ${KOKKOS_PROJECTION_TYPE}")
endif()
if ((BACKEND STREQUAL "thrust") OR (BACKEND STREQUAL "thrust_unified"))
  message(" Thrust iterator type: ${THRUST_ITERATOR_TYPE}")
endif()

add_definitions(-DBACKEND="${BACKEND}")

# Debug and test modes
option(DEBUG "Debug messages" OFF)
option(WARNING "Warning mode (ON/OFF)" OFF)
option(TEST "Simulation test for CI" OFF)
option(UNIT_TESTS "Compile Unit tests" OFF)
option(MINIPIC "Compile minipic (ON/OFF)" ON)
option(SHAMAN "Compile with Shaman" OFF)
option(PROFILER "Task profiler (ON/OFF)" OFF)
option(EVENTIFY "Build eventify if needed" OFF)
option(NVIDIA_PROFILER "Nvidia profiler flags (ON/OFF)" OFF)


message(" Minipic compilation: ${MINIPIC}")
message(" Unit tests: ${UNIT_TESTS}")

# Check that the backend is supported
if (BACKEND STREQUAL "sequential")
  message(" Backend: sequential")
  set(EXE_NAME minipic.seq)
elseif (BACKEND STREQUAL "stdpar")
  message(" Backend: stdpar")
  set(EXE_NAME minipic.stdpar)
elseif (BACKEND STREQUAL "stdpar_cpu")
  message(" Backend: stdpar_cpu")
  set(EXE_NAME minipic.stdparcpu)
elseif (BACKEND STREQUAL "openmp")
  message(" Backend: openmp")
  set(EXE_NAME minipic.omp)
elseif (BACKEND STREQUAL "openmp_task")
  message(" Backend: openmp_task")
  set(EXE_NAME minipic.omp_task)
elseif (BACKEND STREQUAL "openmp_target")
  message(" Backend: openmp target")
  set(EXE_NAME minipic.omp_target)
elseif (BACKEND STREQUAL "openacc")
  message(" Backend: openacc")
  set(EXE_NAME minipic.acc)
elseif (BACKEND STREQUAL "kokkos")
  message(" Backend: kokkos")
  set(EXE_NAME minipic.kokkos)
elseif (BACKEND STREQUAL "kokkos_dualview_unified")
  message(" Backend: kokkos using dual views with unified memory")
  set(EXE_NAME minipic.kokkos_dualview_unified)
elseif (BACKEND STREQUAL "kokkos_unified")
  message(" Backend: kokkos with unified memory")
  set(EXE_NAME minipic.kokkos_unified)
elseif (BACKEND STREQUAL "kokkos_views")
  message(" Backend: kokkos with views")
  set(EXE_NAME minipic.kokkos_views)
elseif ((BACKEND STREQUAL "thrust") OR (BACKEND STREQUAL "thrust_unified"))
  message(" Backend: thrust")
  set(EXE_NAME minipic.thrust)
elseif (BACKEND STREQUAL "eventify")
  message(" Backend: eventify")
  set(EXE_NAME minipic.eventify)
elseif (BACKEND STREQUAL "sycl")
  message(" Backend: sycl")
  set(EXE_NAME minipic.icpx)
elseif (BACKEND STREQUAL "acpp")
  message(" Backend: acpp")
  set(EXE_NAME minipic.acpp)
else()
  message(FATAL_ERROR "FATAL: Backend (${BACKEND}) not supported")
endif()

message(" Executable name: ${EXE_NAME}")

# Specific architecture option
set(DEVICE "auto" CACHE STRING "Device architecture (auto (default), nvidia_v100)")

# Check that the architecture is supported
if ((DEVICE STREQUAL "auto") OR (DEVICE STREQUAL "nvidia_v100") OR (DEVICE STREQUAL "nvidia_a100") OR (DEVICE STREQUAL "amd_mi250") OR (DEVICE STREQUAL "nvidia_gh200") OR (DEVICE STREQUAL "nvidia_h100") OR (DEVICE STREQUAL "intel_pvc") OR (DEVICE STREQUAL "amd_mi300") OR (DEVICE STREQUAL "amd_genoa") OR (DEVICE STREQUAL "nvidia_grace") OR (DEVICE STREQUAL "nvidia_gb200"))
 message(" Device architecture: ${DEVICE}")
else()
 message(FATAL_ERROR "FATAL: Device (${DEVICE}) architecture not supported")
endif()

# ________________________________________________________________
# Determine if UVM is used

if ((BACKEND STREQUAL "stdpar") OR (BACKEND STREQUAL "thrust_unified"))
  set(MINIPIC_USES_UVM ON)
else()
  set(MINIPIC_USES_UVM OFF)
endif()

message(" Unified memory: ${MINIPIC_USES_UVM}")

# _________________________________________________________________
# Debug mode

if (DEBUG)
 message(" Debug mode: ON")
endif()

# _________________________________________________________________
# Warning mode

if (WARNING)
 message(" Warning mode: ON")
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
endif()

# __________________________________________________________________
# Determine languages

set(LANGUAGES CXX)

if ((BACKEND STREQUAL "thrust") OR (BACKEND STREQUAL "stdpar") OR (BACKEND STREQUAL "thrust_unified"))
 if ((DEVICE STREQUAL "nvidia_v100") OR (DEVICE STREQUAL "nvidia_a100") OR (DEVICE STREQUAL "nvidia_h100") OR (DEVICE STREQUAL "auto") OR (DEVICE STREQUAL "nvidia_gh200") OR (DEVICE STREQUAL "nvidia_gb200"))
 set(LANGUAGES ${LANGUAGES} CUDA)
 endif()
endif()

message(" Languages: ${LANGUAGES}")

message(" __________________________________________________________________ \n")

project (minipic LANGUAGES ${LANGUAGES})

# __________________________________________________________________
# Archiecture specific optimization flags flags

set(ARCH_OPT_FLAGS "")

if ( ((BACKEND STREQUAL "thrust") OR (BACKEND STREQUAL "thrust_unified") OR (BACKEND STREQUAL "stdpar")) AND (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")) )

  # flags for specific device architecture
  if (DEVICE STREQUAL "nvidia_v100")
    set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -gencode arch=compute_70,code=sm_70)
  elseif (DEVICE STREQUAL "nvidia_a100")
    set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -gencode arch=compute_80,code=sm_80)
  elseif (DEVICE STREQUAL "nvidia_h100")
      set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -gencode arch=compute_90,code=sm_90)
  elseif (DEVICE STREQUAL "nvidia_gh200")
      set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -gencode arch=compute_90,code=sm_90)
  elseif (DEVICE STREQUAL "nvidia_gb200")
      set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -gencode arch=compute_100,code=sm_100)
  elseif (DEVICE STREQUAL "amd_mi250")
    if (MINIPIC_USES_UVM STREQUAL ON)
      set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} --offload-arch=gfx90a:xnack+)
    else()
      set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} --offload-arch=gfx90a)
    endif()
  elseif (DEVICE STREQUAL "amd_mi300")
    if (MINIPIC_USES_UVM STREQUAL ON)
      set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} --offload-arch=gfx942:xnack+)
    else()
      set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} --offload-arch=gfx942)
    endif()
  endif()

elseif (((BACKEND STREQUAL "openacc") OR (BACKEND STREQUAL "openmp_target") OR (BACKEND STREQUAL "stdpar") OR (BACKEND STREQUAL "stdpar_cpu") OR (BACKEND STREQUAL "thrust") OR (BACKEND STREQUAL "thrust_unified")) AND ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "NVHPC") )

  # flags for specific device architecture
  if (DEVICE STREQUAL "nvidia_v100")
   set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -gpu=cc70)
  elseif (DEVICE STREQUAL "nvidia_a100")
    set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -gpu=cc80)
  elseif (DEVICE STREQUAL "nvidia_h100")
   set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -gpu=cc90)
  elseif (DEVICE STREQUAL "nvidia_gh200")
   set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -gpu=cc90)
  elseif (DEVICE STREQUAL "nvidia_gb200")
   set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -gpu=cc100)
  elseif (DEVICE STREQUAL "amd_mi250")
    if (MINIPIC_USES_UVM STREQUAL ON)
      set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} --offload-arch=gfx90a:xnack+)
    else()
      set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} --offload-arch=gfx90a)
    endif()
  elseif (DEVICE STREQUAL "amd_mi300")
  if (MINIPIC_USES_UVM STREQUAL ON)
    set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} --offload-arch=gfx942:xnack+)
  else()
    set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} --offload-arch=gfx942)
  endif()
  elseif (DEVICE STREQUAL "amd_genoa")
	  #set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -march=znver4)
 endif()

elseif (BACKEND STREQUAL "sycl")

 # flags for specific device architecture
 if (DEVICE STREQUAL "cpu_x86")
   set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -fsycl-targets=spir64_x86_64 -D__MINIPIC_ON_CPU__)
 elseif (DEVICE STREQUAL "nvidia_v100")
   set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -fsycl-targets=nvidia_gpu_sm_70 -D__MINIPIC_ON_GPU__)
 elseif (DEVICE STREQUAL "nvidia_a100")
   set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -fsycl-targets=nvidia_gpu_sm_80 -D__MINIPIC_ON_GPU__)
 elseif (DEVICE STREQUAL "nvidia_h100")
   set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -fsycl-targets=nvidia_gpu_sm_90 -D__MINIPIC_ON_GPU__)
 elseif (DEVICE STREQUAL "nvidia_gh200")
   set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -fsycl-targets=nvidia_gpu_sm_90 -D__MINIPIC_ON_GPU__)
  elseif (DEVICE STREQUAL "nvidia_gb200")
   set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -fsycl-targets=nvidia_gpu_sm_100 -D__MINIPIC_ON_GPU__)
 elseif (DEVICE STREQUAL "amd_mi250")
   set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -fsycl-targets=amd_gpu_gfx90a -D__MINIPIC_ON_GPU__)
 elseif (DEVICE STREQUAL "amd_mi300")
   set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -fsycl-targets=amdgcn-amd-amdhsa -Xsycl-target-backend=amdgcn-amd-amdhsa --offload-arch=gfx942 -D__MINIPIC_ON_GPU__)
 elseif (DEVICE STREQUAL "intel_pvc")
   set(ARCH_OPT_FLAGS ${ARCH_OPT_FLAGS} -D__MINIPIC_ON_GPU__)
 endif()

endif()

message(" Archiecture specific optimization flags: ${ARCH_OPT_FLAGS}")

# __________________________________________________________________

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

# -----> Kokkos mode <------
elseif (BACKEND STREQUAL "kokkos")

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++17 -D__MINIPIC_KOKKOS__ -D__MINIPIC_KOKKOS_COMMON__ -D__MINIPIC_KOKKOS_DUALVIEW_COMMON__ ${ARCH_OPT_FLAGS})

 if (KOKKOS_PROJECTION_TYPE STREQUAL "scatter_view")
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_KOKKOS_SCATTERVIEW__)
 else() #if (KOKKOS_PROJECTION_TYPE STREQUAL "atomic_view")
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_KOKKOS_ATOMICVIEW__)
 endif()

 include_directories(src/setups src/kokkos src/common)

 # list of source files - automatic search cpp files in src/openmp and src/common directory
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/kokkos/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

 find_package(Kokkos REQUIRED)

# -----> Kokkos mode using unified memory <------
elseif (BACKEND STREQUAL "kokkos_dualview_unified")

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++17 -D__MINIPIC_KOKKOS_DUALVIEW_UNIFIED__ -D__MINIPIC_KOKKOS_COMMON__ -D__MINIPIC_KOKKOS_DUALVIEW_COMMON__ ${ARCH_OPT_FLAGS})

 if (KOKKOS_PROJECTION_TYPE STREQUAL "scatter_view")
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_KOKKOS_SCATTERVIEW__)
 else() #if (KOKKOS_PROJECTION_TYPE STREQUAL "atomic_view")
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_KOKKOS_ATOMICVIEW__)
 endif()

 include_directories(src/setups src/kokkos src/common)

 # list of source files - automatic search cpp files in src/openmp and src/common directory
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/kokkos/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

 find_package(Kokkos REQUIRED)

# -----> Kokkos mode using unified memory <------
elseif (BACKEND STREQUAL "kokkos_unified")

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++17 -D__MINIPIC_KOKKOS_UNIFIED__ -D__MINIPIC_KOKKOS_COMMON__ ${ARCH_OPT_FLAGS})

 if (KOKKOS_PROJECTION_TYPE STREQUAL "scatter_view")
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_KOKKOS_SCATTERVIEW__)
 else() #if (KOKKOS_PROJECTION_TYPE STREQUAL "atomic_view")
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_KOKKOS_ATOMICVIEW__)
 endif()

 set(LINK_OPTIONS ${LINK_OPTIONS} ${ARCH_OPT_FLAGS})
 include_directories(src/setups src/kokkos src/common)

 # list of source files - automatic search cpp files in src/openmp and src/common directory
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/kokkos/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

 find_package(Kokkos REQUIRED)

# -----> Kokkos mode using unified memory <------
elseif (BACKEND STREQUAL "kokkos_views")

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++17 -D__MINIPIC_KOKKOS_VIEWS__ -D__MINIPIC_KOKKOS_COMMON__ ${ARCH_OPT_FLAGS})

 if (KOKKOS_PROJECTION_TYPE STREQUAL "scatter_view")
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_KOKKOS_SCATTERVIEW__)
 else() #if (KOKKOS_PROJECTION_TYPE STREQUAL "atomic_view")
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_KOKKOS_ATOMICVIEW__)
 endif()

 set(LINK_OPTIONS ${LINK_OPTIONS} ${ARCH_OPT_FLAGS})
 include_directories(src/setups src/kokkos src/common)

 # list of source files - automatic search cpp files in src/openmp and src/common directory
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/kokkos/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

 find_package(Kokkos REQUIRED)

# -----> Thrust mode <------
elseif ((BACKEND STREQUAL "thrust") OR (BACKEND STREQUAL "thrust_unified"))

 # set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} ${ARCH_OPT_FLAGS})
 # set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${ARCH_OPT_FLAGS})

  if (BACKEND STREQUAL "thrust")
    set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_THRUST__ )
  elseif (BACKEND STREQUAL "thrust_unified")
    set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_THRUST_UNIFIED__)
  endif()

  # set(CMAKE_CXX_STANDARD 20)
  # set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++20)

  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_THRUST_COMMON__ ${ARCH_OPT_FLAGS})

  if ((DEVICE STREQUAL "amd_mi250") OR (DEVICE STREQUAL "amd_mi300"))
    set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -x hip) 
  else()
    if ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      set(COMPILE_OPTIONS ${COMPILE_OPTIONS} --x cu --extended-lambda)
    else()
      set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -x cu -cuda )
      set(LINK_OPTIONS ${LINK_OPTIONS} -cuda)
    endif()
  endif()

  if (THRUST_ITERATOR_TYPE STREQUAL "zip")
    set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_THRUST_ZIP__)
  elseif (THRUST_ITERATOR_TYPE STREQUAL "counting")
    set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_THRUST_COUNTING__)
  else()
    message(FATAL_ERROR "FATAL: Thrust iterator type (${THRUST_ITERATOR_TYPE}) not supported")
  endif()

  set(LINK_OPTIONS ${LINK_OPTIONS} ${ARCH_OPT_FLAGS})

  include_directories(src/setups src/thrust src/common)

  # list of source files - automatic search cpp files in src/openmp and src/common directory
  file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/thrust/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

  if ((DEVICE STREQUAL "nvidia_v100") OR (DEVICE STREQUAL "nvidia_a100") OR (DEVICE STREQUAL "nvidia_h100") OR (DEVICE STREQUAL "auto") OR (DEVICE STREQUAL "nvidia_gh200"))
  # find_package(Thrust REQUIRED)
  find_package(CUDA REQUIRED)
  endif()

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
  
# -----> OpenMP target mode <------
elseif (BACKEND STREQUAL "openmp_target")

 # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xcompiler=-fopenmp")
 # set(CMAKE_CXX_FLAGS ${ARCH_OPT_FLAGS})

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++17 ${ARCH_OPT_FLAGS} -D__MINIPIC_OMP_TARGET__)
 set(LINK_OPTIONS ${LINK_OPTIONS} ${ARCH_OPT_FLAGS})

 # if compiler is NVHPC
 if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "NVHPC")
 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -mp=gpu)
 set(LINK_OPTIONS ${LINK_OPTIONS} -mp=gpu)
 endif()

 include_directories(src/setups src/openmp_target src/common)

 # list of source files - automatic search cpp files
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/openmp_target/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

 find_package(OpenMP REQUIRED)

# -----> OpenACC mode <------
elseif (BACKEND STREQUAL "openacc")

 # set(CMAKE_CUDA_FLAGS ${ARCH_OPT_FLAGS})
 # set(CMAKE_CXX_FLAGS ${ARCH_OPT_FLAGS})

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++17 -D__MINIPIC_OPENACC__ ${ARCH_OPT_FLAGS})
 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -acc)
 set(LINK_OPTIONS ${LINK_OPTIONS} -acc ${ARCH_OPT_FLAGS})

 include_directories(src/setups src/openmp_target src/common)

 # list of source files - automatic search cpp files
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/openmp_target/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

 find_package(OpenACC REQUIRED)

# -----> Eventify Task mode base version <------
elseif (BACKEND STREQUAL "eventify")
 
 set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fopenmp)
 include_directories(src/setups src/eventify src/common)

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++20 -D__MINIPIC_EVENTIFY__ -D__MINIPIC_SIMD__ ${ARCH_OPT_FLAGS})

 # list of source files - automatic search cpp files in src/openmp and src/common directory
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/eventify/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

 find_package(OpenMP REQUIRED)
 find_package (Eventify REQUIRED)

# ----> Sycl / OneAPI mode <----
elseif (BACKEND STREQUAL "sycl") 

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -fsycl -std=c++17 -D__MINIPIC_SYCL__ ${ARCH_OPT_FLAGS})
 set(LINK_OPTIONS ${LINK_OPTIONS} -fsycl ${ARCH_OPT_FLAGS})

  include_directories(src/setups src/sycl src/common)

  # list of source files - automatic search cpp files in src/openmp and src/common directory
  file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/sycl/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

# ----> Sycl / Adaptive C++ mode <----
elseif (BACKEND STREQUAL "acpp")

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++17 -D__MINIPIC_SYCL__ ${ARCH_OPT_FLAGS})

 include_directories(src/setups src/sycl)

 # list of source files - automatic search cpp files in src/openmp and src/common directory
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/sycl/*.cpp")


# -----> Cpp_Stdpar mode <------
elseif (BACKEND STREQUAL "stdpar")

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__MINIPIC_STDPAR__)

 include_directories(src/setups src/stdpar src/common)
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/stdpar/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")


 if ((DEVICE STREQUAL "amd_mi250") OR (DEVICE STREQUAL "amd_mi300"))
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} --hipstdpar -std=c++20 ${ARCH_OPT_FLAGS})
  set(LINK_OPTIONS ${LINK_OPTIONS} --hipstdpar ${ARCH_OPT_FLAGS})
 else()
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -x cu -cuda -std=c++20 -stdpar ${ARCH_OPT_FLAGS})
  set(LINK_OPTIONS ${LINK_OPTIONS} -stdpar -cuda -lcudart ${ARCH_OPT_FLAGS})
 endif()

 if ((DEVICE STREQUAL "nvidia_v100") OR (DEVICE STREQUAL "nvidia_a100") OR (DEVICE STREQUAL "nvidia_h100") OR (DEVICE STREQUAL "nvidia_gh200"))
 find_package(CUDA REQUIRED)
 endif()

# -----> Cpp_Stdpar mode on CPU <------
elseif (BACKEND STREQUAL "stdpar_cpu")  

 set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -std=c++20 -D__MINIPIC_STDPAR_CPU__)

 # Compiler specific flags
 if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "NVHPC")
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -stdpar -mp=multicore)
  set(LINK_OPTIONS ${LINK_OPTIONS} -stdpar -mp=multicore)
 elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -lpthread)
  set(LINK_OPTIONS ${LINK_OPTIONS} -lpthread)
 elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM")
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -fopenmp)
  set(LINK_OPTIONS ${LINK_OPTIONS} -fopenmp)
 elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -fopenmp)
  set(LINK_OPTIONS ${LINK_OPTIONS} -fopenmp)
 endif() 

 if (DEVICE STREQUAL "amd_genoa")
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} ${ARCH_OPT_FLAGS})
  set(LINK_OPTIONS ${LINK_OPTIONS} ${ARCH_OPT_FLAGS})
 endif()

 include_directories(src/setups src/stdpar src/common)

 # list of source files - automatic search cpp files in src/openmp and src/common directory
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/stdpar/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp")

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

# -----> Nvidia Profilers <------
if (NVIDIA_PROFILER)
  set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -D__NVIDIA_PROFILER__)
  set(LINK_OPTIONS ${LINK_OPTIONS} -cuda -lcudart -lnvToolsExt ${ARCH_OPT_FLAGS})
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


