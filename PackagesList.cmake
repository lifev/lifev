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

INCLUDE(TribitsListHelpers)

#
# Define the LifeV package names, directories, and classification.
#
# Package classifications are:
#
#   PS: Primary Stable Package
#
#     Primary Stable Packages have at least some Primary Stable Code which is
#     expected to be fully tested before every push to the global repo.  The
#     default enable for PS packages is empty "" which allows the PS package
#     to be enabled implicitly based on other criteria.  The option
#     LifeV_ENABLE_ALL_PACKAGES=ON will cause all PS packages to be enabled
#     unless they are explicitly disabled.
#
#   SS: Secondary Stable Package
#
#     Secondary Stable Packages have no PS code or they would be classified as
#     PS packages.  A package must be classified as SS if it has a required
#     dependency on another SS package or SS TPL.  A package may also be
#     declared SS to avoid requiring it to be tested before every push to the
#     global repo.  For example, a package that does not provide any
#     significant functionally like Didasko is classified as a SS package even
#     through it could be classified as PS just based on its required package
#     and TPL dependencies.  SS packages will have their default enables set
#     to empty "".  This allows them to be enabled implicilty.  When
#     LifeV_ENABLE_ALL_PACKAGES=ON but
#     LifeV_ENABLE_SECONDARY_STABLE_CODE=OFF, the SS packages will not be
#     enabled.  However, when LifeV_ENABLE_ALL_PACKAGES=ON and
#     LifeV_ENABLE_SECONDARY_STABLE_CODE=ON, then SS packages will be
#     enabled if they are not explicitly disabled.  Packages that are SS but
#     not PS must be disabled in pre-push testing.  However, SS packages are
#     tested by the post-push CI and nightly testing processes.
#
#   EX: Experimental Package
#
#     Experimental packages are those packages that contain no PS or SS
#     code. The default enable for EX packages is always OFF which requires
#     that they be explicitly enabled in order to be turned on. EX packages
#     must be disabled in pre-push testring and are not tested as part of the
#     post-push CI or nightly testing processes.  However, package developers
#     of EX pacakges are encouraged to set up their own nightly testing for
#     thier EX packages.
#
# NOTE: These packages must be listed in strictly assending order in terms of
# package dependencies.  If you get the order wrong, then an error message
# will be printed during configuration with CMake.
#

SET( LifeV_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  Core                  lifev/core                        PS
  ETA                   lifev/eta                         PS
  BCInterface           lifev/bc_interface                PS
  OneDFSI               lifev/one_d_fsi                   PS
  LevelSet              lifev/level_set                   PS
  Darcy                 lifev/darcy                       PS
  NavierStokes          lifev/navier_stokes               PS
  NavierStokesBlocks    lifev/navier_stokes_blocks        PS
  Structure             lifev/structure                   PS
  Electrophysiology     lifev/electrophysiology           EX
  Heart                 lifev/heart                       EX
  FSI                   lifev/fsi                         PS
  FSI_blocks            lifev/fsi_blocks                  EX
  ZeroDimensional       lifev/zero_dimensional            PS
  Multiscale            lifev/multiscale                  PS
  Dummy                 lifev/dummy                       EX
)


#
# Disable certain packages on certain platforms.
#
# NOTE: This just makes the packages experimental 'EX' and therefore still
# allows the user to enable the package explicitly but the package will not
# get enabled implicitly.
#

#PACKAGE_DISABLE_ON_PLATFORMS(MOOCHO Windows)
