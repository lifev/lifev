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

dist_pkgdata_DATA = Fluid3D_Cylinder.dat \
					Fluid3D_Cylinder_Stokes.dat \
					Fluid3D_Carotid.dat \
					Fluid3D_Carotid_Stokes.dat \
                    MultiScale_3Dbifurcation_FluxStress.dat \
                    MultiScale_3Dbifurcation_FluxStress_Stationary.dat \
                    MultiScale_3Dbifurcation_Stress_Stationary.dat \
                    MultiScale_3DCarotidBifurcation_FluxStress.dat \
                    MultiScale_3DCarotidBifurcation_FluxStress_Stationary.dat \
                    MultiScale_3DCarotidBifurcation_Stress.dat \
                    MultiScale_3DCarotidBifurcation_Stress_Stationary.dat \
                    MultiScale_3DSyntheticCarotidBifurcation_FluxStress.dat \
                    MultiScale_3DSyntheticCarotidBifurcation_FluxStress_Stationary.dat \
                    MultiScale_3DSyntheticCarotidBifurcation_Stress.dat \
                    MultiScale_3DSyntheticCarotidBifurcation_Stress_Stationary.dat \
                    MultiScale_3Dhexagon_FluxStress_Stationary.dat \
                    MultiScale_3Dhexagon_Stress_Stationary.dat \
                    MultiScale_3Dhexagon_StressAndFluxStress_Stationary.dat \
                    MultiScale_3Dserial1.dat \
                    MultiScale_3Dserial1_Stationary.dat \
                    MultiScale_3Dserial2_FluxStress.dat \
                    MultiScale_3Dserial2_FluxStress_Stationary.dat \
                    MultiScale_3Dserial2_Stress_Stationary.dat \
                    MultiScale_3Dserial3_FluxStress.dat \
                    MultiScale_3Dserial3_FluxStress_Stationary.dat \
                    MultiScale_3Dserial3_Stress_Stationary.dat \
                    MultiScale_3Dserial5_FluxStress.dat \
                    MultiScale_3Dserial5_Stress.dat \
                    MultiScale_3Dserial6_FluxStress.dat \
                    MultiScale_3Dserial6_FluxStress_Stationary.dat \
                    MultiScale_3Dserial6_Stress.dat \
                    MultiScale_3Dserial6_Stress_Stationary.dat

check-local:
	@if test "$(top_srcdir)" != "$(top_builddir)"; then\
	 for i in $(dist_pkgdata_DATA); do \
           $(LN_S) -f $(srcdir)/$$i .; \
         done; \
	fi

