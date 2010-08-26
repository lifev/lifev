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

SUBDIRS		= 	Couplings\
				Models

dist_pkgdata_DATA = MultiScale_test.dat \
					MultiScale.dat \
					Run_CarotidBifurcation_1D-3DFluid_FluxStress.dat \
					Run_CarotidBifurcation_1D-3DFSI_FluxStress.dat \
					Run_CarotidBifurcation_3DFluid.dat \
					Run_CarotidBifurcation_3DFluid_FluxStress.dat \
					Run_CarotidBifurcation_3DFluid_Stress.dat \
					Run_CarotidBifurcation_3DFluid_Synthetic_FluxStress.dat \
					Run_CarotidBifurcation_3DFluid_Synthetic_Stress.dat \
					Run_Stationary_3DFluid_1Cylinder.dat \
					Run_Stationary_3DFluid_CarotidBifurcation_FluxStress.dat \
					Run_Stationary_3DFluid_Hexagon_FluxStress.dat \
					Run_Stationary_3DFluid_Hexagon_Stress.dat \
					Run_SinusoidalWave_1D_1Cylinder.dat \
					Run_SinusoidalWave_1D_2Cylinders_Stress.dat \
					Run_SinusoidalWave_1D_2Cylinders_FluxStress.dat \
					Run_SinusoidalWave_1D_3Cylinders_FluxStress.dat \
					Run_SinusoidalWave_1D-3DFluid_2Cylinders_FluxStress.dat \
					Run_SinusoidalWave_1D-3DFluid_2Cylinders_Stress.dat \
					Run_SinusoidalWave_3DFluid-1D_2Cylinders_FluxStress.dat \
					Run_SinusoidalWave_3DFluid-1D_2Cylinders_Stress.dat \
					Run_SinusoidalWave_1D-3DFSI_2Cylinders_FluxStress.dat \
					Run_TestStress_3DFSI_1Cylinder.dat \
                    Run_TestStress_3DFSI_2Cylinders_FluxStress.dat \
                    Run_TestStress_3DFSI_2Cylinders_Stress.dat \
					Run_Womersley_3DFluid_1Cylinder.dat \
					Run_Womersley_3DFluid_5Cylinders_FluxStress_tube20.dat \
					Run_Womersley_3DFluid_5Cylinders_FluxStress_tube40.dat \
					Run_Womersley_3DFluid_5Cylinders_FluxStress_tube60.dat \
					Run_Womersley_3DFluid_5Cylinders_FluxStress_tube80.dat \
					Run_Womersley_3DFluid_5Cylinders_Stress_tube20.dat \
					Run_Womersley_3DFluid_5Cylinders_Stress_tube40.dat \
					Run_Womersley_3DFluid_5Cylinders_Stress_tube60.dat \
					Run_Womersley_3DFluid_5Cylinders_Stress_tube80.dat

check-local:
	@if test "$(top_srcdir)" != "$(top_builddir)"; then\
	 for i in $(dist_pkgdata_DATA); do \
           $(LN_S) -f $(srcdir)/$$i .; \
         done; \
	fi

