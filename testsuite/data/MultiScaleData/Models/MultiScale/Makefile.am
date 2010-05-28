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

dist_pkgdata_DATA = CarotidBifurcation_1D-3DFluid_FluxStress.dat \
                    CarotidBifurcation_1D-3DFSI_FluxStress.dat \
                    CarotidBifurcation_3DFluid.dat \
					CarotidBifurcation_3DFluid_FluxStress.dat \
                    CarotidBifurcation_3DFluid_Stress.dat \
                    CarotidBifurcation_3DFluid_Synthetic_FluxStress.dat \
                    CarotidBifurcation_3DFluid_Synthetic_Stress.dat \
                    SinusoidalWave_1D_1Cylinder.dat \
                    SinusoidalWave_1D_2Cylinders_Stress.dat \
                    SinusoidalWave_1D_2Cylinders_FluxStress.dat \
                    SinusoidalWave_1D-3DFluid_2Cylinders_FluxStress.dat \
                    SinusoidalWave_1D-3DFSI_2Cylinders_FluxStress.dat \
                    SinusoidalWave_3DFluid-1D_2Cylinders_FluxStress.dat \
                    Stationary_3DFluid_CarotidBifurcation_FluxStress.dat \
                    Stationary_3DFluid_CarotidBifurcation_Stress.dat \
                    Stationary_3DFluid_CarotidBifurcation_Synthetic_FluxStress.dat \
                    Stationary_3DFluid_CarotidBifurcation_Synthetic_Stress.dat \
                    Stationary_3DFluid_1Cylinder_tube20.dat \
                    Stationary_3DFluid_Bifurcation_FluxStress.dat \
                    Stationary_3DFluid_Bifurcation_Stress.dat \
                    Stationary_3DFluid_Hexagon_FluxStress.dat \
                    Stationary_3DFluid_Hexagon_Mixed.dat \
                    Stationary_3DFluid_Hexagon_Stress.dat \
                    TestStress_3DFSI_1Cylinder_tube20.dat \
                    TestStress_3DFSI_2Cylinders_FluxStress_tube20.dat \
                    TestStress_3DFSI_2Cylinders_Stress_tube20.dat \
                    Womersley_3DFluid_1Cylinder_tube20.dat \
                    Womersley_3DFluid_5Cylinders_FluxStress_tube20.dat \
                    Womersley_3DFluid_5Cylinders_FluxStress_tube40.dat \
                    Womersley_3DFluid_5Cylinders_FluxStress_tube60.dat \
                    Womersley_3DFluid_5Cylinders_FluxStress_tube80.dat \
                    Womersley_3DFluid_5Cylinders_Stress_tube20.dat \
                    Womersley_3DFluid_5Cylinders_Stress_tube40.dat \
                    Womersley_3DFluid_5Cylinders_Stress_tube60.dat \
                    Womersley_3DFluid_5Cylinders_Stress_tube80.dat

check-local:
	@if test "$(top_srcdir)" != "$(top_builddir)"; then\
	 for i in $(dist_pkgdata_DATA); do \
           $(LN_S) -f $(srcdir)/$$i .; \
         done; \
	fi

