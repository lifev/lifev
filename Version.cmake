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
# Single file that needs to be changed on a release branch
# or on the development branch in order to configure Trilinos
# for release mode and set the version.
#

SET(LifeV_VERSION 3.10.0)
SET(LifeV_MAJOR_VERSION 3)
SET(LifeV_MINOR_VERSION 10)
SET(LifeV_MICRO_VERSION 0)
SET(LifeV_MAJOR_MINOR_VERSION 3100)
SET(LifeV_VERSION_STRING "3.10.0")
SET(LifeV_ENABLE_DEVELOPMENT_MODE_DEFAULT ON) # Change to 'OFF' for a release

# Used by testing scripts and should not be used elsewhere
SET(LifeV_REPOSITORY_BRANCH "master" CACHE INTERNAL "")
SET(LifeV_TESTING_TRACK "Nightly build 3.10" CACHE INTERNAL "")
