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
# Define the list of TPLs, their find module names, and their classification
#
# TPL_NAME:
#
#   The name of the TPL used in the CMake cache variables TPL_ENABLE_${TPL_NAME}
#
# TPL_FINDMOD:
#
#   The name of the find module under that is used to get the names of the
#   TPLs.  If ends in '/' then this gives the directory and the standard module
#   name will be used which is FindTPL${TPL_NAME}.cmake.
#
# TPL_CLASSIFICATION:
#
#   PS: Primary Stable TPL
#
#     Primary Stable TPLs are those TPLs that a LifeV developer must have
#     installed on their machine in order to be able to do LifeV
#     development.  For example, we require that you have BLAS, LAPACK, and
#     MPI installed in order to do LifeV development.  These are
#     fundamental dependencies that are needed in order to do precheckin
#     testing.
#
#   SS: Secondary Stable TPL
#
#     Secondary Stable TPLs are those TPLs that are not required in order to
#     be able to develop and test LifeV before checkins but are none the
#     less offically supported.  Support for SS TPLs is tested as part of the
#     nightly testing process.
#
#   TS: Tertiary Stable TPL
#
#     Tertiary Stable TPLs are those TPLs that are supported TPLs but can not
#     be included in the set of SS TPLs because they may conflict with other
#     SS Code.  For example, METIS is listed as a TS TPL because it conflicts
#     with ParMETIS which is declared as a SS TPL.
#
#   EX: Experimental TPL
#
#     Experimental TPLs are not offically supported.  They represent
#     experimental capabilities of LifeV packages.  Support for EX TPLs is
#     never tested as part of the main nightly testing process.  However,
#     package developers are encouraged to set up their own nightly testing
#     for their EX TPLs for their packages.
#
# The default enable for all TPLs is empty "" reguardless of the category.
# The idea is that the enabling of the TPL will be done by the package and
# other enables that the user has to set.
#
# NOTE: The TPLs must be listed in the order of increasing dependencies (if
# such dependencies exist).
#

SET( LifeV_TPLS_FINDMODS_CLASSIFICATIONS
  MPI             "${${PROJECT_NAME}_TRIBITS_DIR}/tpls/"    PS
  BLAS            "cmake/TPLs/"    PS
  LAPACK          "cmake/TPLs/"    PS
  Boost           "cmake/TPLs/"    PS
  ParMETIS        "cmake/TPLs/"    PS
  HDF5            "cmake/TPLs/"    PS
  QHull           "cmake/TPLs/"    SS
  Trilinos        "cmake/TPLs/"    PS
  muparser        "cmake/TPLs/"    PS
  )

# NOTES:
#
# (*) ParMETIS must be listed after Scotch because the
#     ParMETIS include directories must come before the
#     Scotch include directories.
#
