# target compile options using EXE_NAME

# ---->  OpenMP mode <----
if (BACKEND STREQUAL "openmp")

  target_link_libraries("${EXE_NAME}" OpenMP::OpenMP_CXX)

# ---->  Kokkos mode <----
elseif (BACKEND STREQUAL "kokkos" OR BACKEND STREQUAL "kokkos_unified" OR BACKEND STREQUAL "kokkos_dualview_unified" OR BACKEND STREQUAL "kokkos_views")
  target_link_libraries("${EXE_NAME}" Kokkos::kokkos)

  # ---->  Thrus mode <----
elseif ((BACKEND STREQUAL "thrust") OR (BACKEND STREQUAL "thrust_unified"))

  if ((DEVICE STREQUAL "nvidia_v100") OR (DEVICE STREQUAL "nvidia_a100") OR (DEVICE STREQUAL "nvidia_h100") OR 
     (DEVICE STREQUAL "nvidia_gh200"))
    target_link_libraries("${EXE_NAME}" ${CUDA_LIBRARIES})
    #target_link_libraries("${EXE_NAME}" Thrust::thrust)
  endif()

elseif(BACKEND STREQUAL "stdpar")
  if ((DEVICE STREQUAL "nvidia_v100") OR (DEVICE STREQUAL "nvidia_a100") OR (DEVICE STREQUAL "nvidia_h100") OR 
      (DEVICE STREQUAL "nvidia_gh200"))
    target_link_libraries("${EXE_NAME}" ${CUDA_LIBRARIES})
  endif()

# ----> OpenMP Task mode <----
elseif(BACKEND STREQUAL "openmp_task")

  target_link_libraries("${EXE_NAME}" OpenMP::OpenMP_CXX)

# ----> OpenMP Target mode <----
elseif(BACKEND STREQUAL "openmp_target")
  
  target_link_libraries("${EXE_NAME}" OpenMP::OpenMP_CXX)

# ----> OpenACC mode <----
elseif(BACKEND STREQUAL "openacc")


# ----> Eventify Task mode <----
elseif(BACKEND STREQUAL "eventify")

  target_link_libraries("${EXE_NAME}" OpenMP::OpenMP_CXX)
  target_link_libraries ("${EXE_NAME}" jsc::eventify)

# ----> Sycl / OneAPI mode <----
elseif(BACKEND STREQUAL "sycl")

# Sycl / Adaptive C++ mode
elseif(BACKEND STREQUAL "acpp")

else()
   # nothing to do
endif()