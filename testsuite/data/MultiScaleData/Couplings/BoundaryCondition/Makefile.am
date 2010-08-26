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
                    CarotidBifurcationFlux_Dynamic.dat \
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

