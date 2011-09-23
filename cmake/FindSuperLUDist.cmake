# Find SuperLUDist headers and library.
#
# This module defines the following uncached variables:
#  SuperLUDist_FOUND, if false, do not try to use SuperLUDist.
#  SuperLUDist_INCLUDE_DIRS, where to find the headers.
#  SuperLUDist_LIBRARIES, the libraries to link against to use the SuperLUDist library
#  SuperLUDist_LIBRARY_DIRS, the directory where the SuperLUDist library is found.

find_path(
  SuperLUDist_INCLUDE_DIR
  supermatrix.h
  /usr/local/include
  /usr/include
  ${SuperLUDist_ROOT}/include
)

if( SuperLUDist_INCLUDE_DIR )
  find_library(
    SuperLUDist_LIBRARY
    NAMES superlu_dist
    PATHS /usr/local/lib /usr/lib ${SuperLUDist_ROOT}/lib
  )
  if( SuperLUDist_LIBRARY )
    set(SuperLUDist_LIBRARY_DIR "")
    get_filename_component(SuperLUDist_LIBRARY_DIRS ${SuperLUDist_LIBRARY} PATH)
    # Set uncached variables as per standard.
    set(SuperLUDist_FOUND ON)
    set(SuperLUDist_INCLUDE_DIRS ${SuperLUDist_INCLUDE_DIR})
    set(SuperLUDist_LIBRARIES ${SuperLUDist_LIBRARY} -lgfortran)
  endif(SuperLUDist_LIBRARY)
endif(SuperLUDist_INCLUDE_DIR)

if(SuperLUDist_FOUND)
  if(NOT SuperLUDist_FIND_QUIETLY)
    message(STATUS "FindSuperLUDist: Found both headers and library")
  endif(NOT SuperLUDist_FIND_QUIETLY)
else(SuperLUDist_FOUND)
  if(SuperLUDist_FIND_REQUIRED)
    message(FATAL_ERROR "FindSuperLUDist: Could not find headers and/or library")
  endif(SuperLUDist_FIND_REQUIRED)
endif(SuperLUDist_FOUND)
