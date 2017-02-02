# - Try to find the LifeV library
# Once done this will define
#
#  LIFEV_FOUND - System has LifeV
#  LIFEV_INCLUDE_DIR - The LifeV include directory
#  LIFEV_LIBRARIES - The libraries needed to use LifeV

#=============================================================================
# Copyright 2006-2009 Kitware, Inc.
# Copyright 2012 Antonio Cervone <ant.cervone@gmail.com>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# FIND_PATH(LIFEV_INCLUDE_DIRS NAMES lifev/core/LifeV.hpp
#   HINTS
#   ${LifeV_DIR}/include
#   ${LIFEV_ROOT}/include
#   ${LIFEV_INCLUDEDIR}
#   ${LIFEV_INCLUDE_DIRS}
# )

# FIND_LIBRARY(LIFEV_LIBRARIES NAMES lifevcore liblifevcore
#   HINTS
#   ${LifeV_DIR}/lib
#   ${LIFEV_ROOT}/lib
#   ${LIFEV_LIBDIR}
#   ${LIFEV_LIBRARY_DIRS}
# )

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE(LifeV
  NO_MODULE
  PATHS
  ${LifeV_DIR}
  ${LIFEV_ROOT}
)

SET(LIFEV_INCLUDE_DIRS
  ${LifeV_INCLUDE_DIRS}
  ${LifeV_TPL_INCLUDE_DIRS}
)

SET(LIFEV_LIBRARIES
  ${LifeV_LIBRARIES}
  ${LifeV_TPL_LIBRARIES}
)

# FIND_PROGRAM(LIFEV_VECTORSMALLTEST_EXECUTABLE Core_VectorSmall.exe)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(LifeV
  DEFAULT_MSG
  LIFEV_LIBRARIES
  LIFEV_INCLUDE_DIRS
)

MARK_AS_ADVANCED(
  LIFEV_INCLUDE_DIRS
  LIFEV_LIBRARIES
#  LIFEV_VECTORSMALLTEST_EXECUTABLE
)

