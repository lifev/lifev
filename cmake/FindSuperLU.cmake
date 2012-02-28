# Find SuperLU headers and library.
#
# This module defines the following uncached variables:
#  SuperLU_FOUND, if false, do not try to use SuperLU.
#  SuperLU_INCLUDE_DIRS, where to find the headers.
#  SuperLU_LIBRARIES, the libraries to link against to use the SuperLU library
#  SuperLU_LIBRARY_DIRS, the directory where the SuperLU library is found.

find_path(
  SuperLU_INCLUDE_DIR
  supermatrix.h
  /usr/local/include
  /usr/include
  ${SuperLU_ROOT}/include
)

if( SuperLU_INCLUDE_DIR )
  find_library(
    SuperLU_LIBRARY
    NAMES superlu
    PATHS /usr/local/lib /usr/lib ${SuperLU_ROOT}/lib
  )
  if( SuperLU_LIBRARY )
    set(SuperLU_LIBRARY_DIR "")
    get_filename_component(SuperLU_LIBRARY_DIRS ${SuperLU_LIBRARY} PATH)
    # Set uncached variables as per standard.
    set(SuperLU_FOUND ON)
    set(SuperLU_INCLUDE_DIRS ${SuperLU_INCLUDE_DIR})
    set(SuperLU_LIBRARIES ${SuperLU_LIBRARY} -lgfortran)
  endif(SuperLU_LIBRARY)
endif(SuperLU_INCLUDE_DIR)

if(SuperLU_FOUND)
  if(NOT SuperLU_FIND_QUIETLY)
    message(STATUS "FindSuperLU: Found both headers and library")
  endif(NOT SuperLU_FIND_QUIETLY)
else(SuperLU_FOUND)
  if(SuperLU_FIND_REQUIRED)
    message(FATAL_ERROR "FindSuperLU: Could not find headers and/or library")
  endif(SuperLU_FIND_REQUIRED)
endif(SuperLU_FOUND)
