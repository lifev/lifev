# Find SuperLUMT headers and library.
#
# This module defines the following uncached variables:
#  SuperLUMT_FOUND, if false, do not try to use SuperLUMT.
#  SuperLUMT_INCLUDE_DIRS, where to find the headers.
#  SuperLUMT_LIBRARIES, the libraries to link against to use the SuperLUMT library
#  SuperLUMT_LIBRARY_DIRS, the directory where the SuperLUMT library is found.

find_path(
  SuperLUMT_INCLUDE_DIR
  supermatrix.h
  /usr/local/include
  /usr/include
  ${SuperLUMT_ROOT}/include
)

if( SuperLUMT_INCLUDE_DIR )
  find_library(
    SuperLUMT_LIBRARY
    NAMES superlu_mt
    PATHS /usr/local/lib /usr/lib ${SuperLUMT_ROOT}/lib
  )
  if( SuperLUMT_LIBRARY )
    set(SuperLUMT_LIBRARY_DIR "")
    get_filename_component(SuperLUMT_LIBRARY_DIRS ${SuperLUMT_LIBRARY} PATH)
    # Set uncached variables as per standard.
    set(SuperLUMT_FOUND ON)
    set(SuperLUMT_INCLUDE_DIRS ${SuperLUMT_INCLUDE_DIR})
    set(SuperLUMT_LIBRARIES ${SuperLUMT_LIBRARY} -lpthread -lgfortran)
  endif(SuperLUMT_LIBRARY)
endif(SuperLUMT_INCLUDE_DIR)

if(SuperLUMT_FOUND)
  if(NOT SuperLUMT_FIND_QUIETLY)
    message(STATUS "FindSuperLUMT: Found both headers and library")
  endif(NOT SuperLUMT_FIND_QUIETLY)
else(SuperLUMT_FOUND)
  if(SuperLUMT_FIND_REQUIRED)
    message(FATAL_ERROR "FindSuperLUMT: Could not find headers and/or library")
  endif(SuperLUMT_FIND_REQUIRED)
endif(SuperLUMT_FOUND)
