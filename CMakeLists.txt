# @HEADER
# *******************************************************************************
# 
#     Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#     Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University
# 
#     This file is part of LifeV.
# 
#     LifeV is free software; you can redistribute it and/or modify
#     it under the terms of the GNU Lesser General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     LifeV is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
# 
#     You should have received a copy of the GNU Lesser General Public License
#     along with LifeV.  If not, see <http://www.gnu.org/licenses/>.
#
# *******************************************************************************
# @HEADER

#
# A) Define your project name and set up major project options
#

# Must set the project name as a varaible at very beginning before including anything else
# We set the project name in a separate file so CTest scripts can use it.
INCLUDE(${CMAKE_SOURCE_DIR}/ProjectName.cmake)

# CMake requires that you declare the CMake project in the top-level file and
# not in an include file :-(
PROJECT(${PROJECT_NAME} NONE)

# we do not need Fortran
SET(${PROJECT_NAME}_ENABLE_Fortran OFF CACHE BOOL "Disable Fortran compiler")

# set up a Profiling build type
SET(CMAKE_CXX_FLAGS_PROFILING "-O3 -pg -g" CACHE STRING "CXX profiling flags")
SET(CMAKE_C_FLAGS_PROFILING   "-O3 -pg -g" CACHE STRING "C profiling flags")

#
# B) Pull in the TriBITS system and execute
#

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/tribits/package_arch/TribitsProject.cmake)

# CMake requires this be in the top file and not in an include file :-(
CMAKE_MINIMUM_REQUIRED(VERSION ${TRIBITS_CMAKE_MINIMUM_REQUIRED})

# Add documentation
ADD_SUBDIRECTORY(doc)

# Add source dir to includes
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# Do all of the processing for this Tribits project
TRIBITS_PROJECT()

# Check for HDF5 APIs version
if (HDF5_IS_IN_TRILINOS AND NOT EpetraExt_FOUND)
  # H5public.h defines all we need
  find_file(LIFEV_H5PUBLIC_H H5public.h
    HINTS ${Trilinos_TPL_INCLUDE_DIRS}
    NO_DEFAULT_PATH)
  if (LIFEV_H5PUBLIC_H)
    # Found, let's open it a look for minor version
    file(STRINGS ${LIFEV_H5PUBLIC_H} LIFEV_HDF5_VERSION_TMP
      REGEX "^#define H5_VERS_MINOR")
    string(REGEX REPLACE ".*#define H5_VERS_MINOR[\t]+([0-9]+).*" "\\1" 
      LIFEV_HDF5_MINOR_VERSION ${LIFEV_HDF5_VERSION_TMP})
    if(LIFEV_HDF5_MINOR_VERSION LESS 8)
      message(WARNING ":: HDF5 version < 1.8.x! Setting to 1.6.x ...")
      set(H5_USE_16_API 1)
    endif()
  else()
    # I can't find the header ... something wrong!
    message(WARNING ":: Cannot find HDF5 version! Setting to 1.6.x ...")
    set(H5_USE_16_API 1)
  endif()
endif()

# Adding a new CMake variable for each of the enabled packages
STRING(REPLACE " " ";" ENABLED_PACKAGE_LIST ${LifeV_ENABLED_PACKAGES})

FOREACH(ENABLED_PACKAGE ${ENABLED_PACKAGE_LIST})
  STRING(TOUPPER ${ENABLED_PACKAGE} UPPER_PACKAGE)
  SET(LIFEV_HAS_${UPPER_PACKAGE} ON)
ENDFOREACH()

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/${PROJECT_NAME}_config.h.in
               ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_config.h)
# We need the next line to install LifeV_config.h.in too
INSTALL(FILES ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_config.h DESTINATION include)

INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/lifev/debug.areas
        DESTINATION share/lifev)

