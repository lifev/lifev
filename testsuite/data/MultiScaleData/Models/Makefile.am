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

dist_pkgdata_DATA = Fluid1D.dat \
					Fluid3D_Cylinder_tube20.dat \
					Fluid3D_Cylinder_tube40.dat \
					Fluid3D_Cylinder_tube60.dat \
					Fluid3D_Cylinder_tube80.dat \
					Fluid3D_Cylinder_Stokes_tube20.dat \
					Fluid3D_Cylinder_Stokes_tube40.dat \
					Fluid3D_Cylinder_Stokes_tube60.dat \
					Fluid3D_Cylinder_Stokes_tube80.dat \
					Fluid3D_Carotid.dat \
					Fluid3D_Carotid_Stokes.dat \
                    MultiScale_CarotidBifurcation3D_FluxStress.dat \
                    MultiScale_CarotidBifurcation3D_Stress.dat \
                    MultiScale_CarotidBifurcation3D_Synthetic_FluxStress.dat \
                    MultiScale_CarotidBifurcation3D_Synthetic_Stress.dat \
                    MultiScale_Stationary3D_1Cylinder_tube20.dat \
                    MultiScale_Stationary3D_Bifurcation_FluxStress.dat \
                    MultiScale_Stationary3D_Bifurcation_Stress.dat \
                    MultiScale_Stationary3D_Hexagon_FluxStress.dat \
                    MultiScale_Stationary3D_Hexagon_Mixed.dat \
                    MultiScale_Stationary3D_Hexagon_Stress.dat \
                    MultiScale_Womersley3D_1Cylinder_tube20.dat \
                    MultiScale_Womersley3D_5Cylinders_FluxStress_tube20.dat \
                    MultiScale_Womersley3D_5Cylinders_FluxStress_tube40.dat \
                    MultiScale_Womersley3D_5Cylinders_FluxStress_tube60.dat \
                    MultiScale_Womersley3D_5Cylinders_FluxStress_tube80.dat \
                    MultiScale_Womersley3D_5Cylinders_Stress_tube20.dat \
                    MultiScale_Womersley3D_5Cylinders_Stress_tube40.dat \
                    MultiScale_Womersley3D_5Cylinders_Stress_tube60.dat \
                    MultiScale_Womersley3D_5Cylinders_Stress_tube80.dat \
                    MultiScale_3DFSIserial1.dat

check-local:
	@if test "$(top_srcdir)" != "$(top_builddir)"; then\
	 for i in $(dist_pkgdata_DATA); do \
           $(LN_S) -f $(srcdir)/$$i .; \
         done; \
	fi

