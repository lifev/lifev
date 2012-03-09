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

#
# B) Pull in the TriBITS system and execute
#

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/tribits/package_arch/TribitsProject.cmake)

# CMake requires this be in the top file and not in an include file :-(
CMAKE_MINIMUM_REQUIRED(VERSION ${TRIBITS_CMAKE_MINIMUM_REQUIRED})

# Do all of the processing for this Tribits project
TRIBITS_PROJECT()

# Add documentation
ADD_SUBDIRECTORY(doc)

# Adding a new CMake variable for each of the enabled packages
STRING(REPLACE " " ";" LOCAL_PACKAGE_LIST ${${PROJECT_NAME}_ENABLED_PACKAGES})
FOREACH(LifeV_PACKAGE IN ITEMS ${LOCAL_PACKAGE_LIST})
  STRING(TOUPPER ${LifeV_PACKAGE} UPPER_PACKAGE)
  SET(LIFEV_HAS_${UPPER_PACKAGE} ON CACHE BOOL "LifeV has ${LifeV_PACKAGE}")
ENDFOREACH()

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}_config.h.in
               ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_config.h)
