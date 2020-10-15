cmake_minimum_required(VERSION 3.14)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()

project(fr_1d
  LANGUAGES Fortran
  VERSION 2.0
  DESCRIPTION "Unstructured FV Solver for Flow Equations")

enable_language(Fortran)
  
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")  
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")  
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

#
# Compile.
#
file(GLOB_RECURSE sources *.F90)
add_executable(go ${sources})

#include_directories(${LAPACK_INCLUDE_DIRS} ${BLAS_INCLUDE_DIRS})
#target_link_libraries(go ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
