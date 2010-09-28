# -*- makefile -*-
#
# 
#
# Copyright (C) 2004 Christophe Prud'homme (christophe.prudhomme@epfl.ch)
#
#  Distributed under the GPL(GNU Public License):
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc

dist_pkgdata_DATA = Absorbing_1D_left.dat \
					Absorbing_1D_right.dat \
					Absorbing_NobileVergara.dat \
                    SinusoidalFlux.dat \
                    SinusoidalFluxWave.dat \
                    SinusoidalSquareFluxWave.dat \
					SinusoidalPressure.dat \
					WomersleySinusoidalFlux.dat \
					WomersleySinusoidalPressure.dat \
					UnitaryFlux.dat \
                    UnitaryPressure.dat \
                    Pressure.dat \
                    Pressure_Dynamic_STP.dat \
                    PoiseuilleFlux.dat \
                    PoiseuilleFlux_Dynamic_STP.dat \
                    PoiseuilleVelocity.dat \
                    PoiseuilleVelocity_Dynamic_STP.dat \
                    ReymondStergiopulos2009_003_R.dat \
                    ReymondStergiopulos2009_008_R.dat \
                    ReymondStergiopulos2009_010_R.dat \
                    ReymondStergiopulos2009_011_R.dat \
                    ReymondStergiopulos2009_015_R.dat \
                    ReymondStergiopulos2009_019_R.dat \
                    ReymondStergiopulos2009_022_R.dat \
                    ReymondStergiopulos2009_024_R.dat \
                    ReymondStergiopulos2009_025_R.dat \
                    ReymondStergiopulos2009_026_R.dat \
                    ReymondStergiopulos2009_031_R.dat \
                    ReymondStergiopulos2009_032_R.dat \
                    ReymondStergiopulos2009_033_R.dat \
                    ReymondStergiopulos2009_034_R.dat \
                    ReymondStergiopulos2009_036_R.dat \
                    ReymondStergiopulos2009_038_R.dat \
                    ReymondStergiopulos2009_040_R.dat \
                    ReymondStergiopulos2009_041_R.dat \
                    ReymondStergiopulos2009_045_R.dat \
                    ReymondStergiopulos2009_047_R.dat \
                    ReymondStergiopulos2009_048_R.dat \
                    ReymondStergiopulos2009_049_R.dat \
                    ReymondStergiopulos2009_051_R.dat \
                    ReymondStergiopulos2009_053_R.dat \
                    ReymondStergiopulos2009_054_R.dat \
                    ReymondStergiopulos2009_055_R.dat \
                    CarotidBifurcationFlux_Dynamic.dat \
                    HeartFlux.dat \
                    NoFlux.dat \
                    NoStress.dat \
                    TestFlux.dat \
                    TestStress.dat \
                    XY_ZeroTangentialVelocity.dat 
                    
check-local:
	@if test "$(top_srcdir)" != "$(top_builddir)"; then\
	 for i in $(dist_pkgdata_DATA); do \
           $(LN_S) -f $(srcdir)/$$i .; \
         done; \
	fi

