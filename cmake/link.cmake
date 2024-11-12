# target compile options using EXE_NAME

# ---->  OpenMP mode <----
if (BACKEND STREQUAL "openmp")

  target_link_libraries("${EXE_NAME}" OpenMP::OpenMP_CXX)

# ----> OpenMP Task mode <----
elseif(BACKEND STREQUAL "openmp_task")

  target_link_libraries("${EXE_NAME}" OpenMP::OpenMP_CXX)

# ----> Eventify Task mode <----
elseif(BACKEND STREQUAL "eventify")

  target_link_libraries("${EXE_NAME}" OpenMP::OpenMP_CXX)
  target_link_libraries ("${EXE_NAME}" jsc::eventify)

else()
   # nothing to do
endif()