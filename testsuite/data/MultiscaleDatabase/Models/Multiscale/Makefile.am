###################################################################################################
#
#                       This file is part of the LifeV Library                       
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University   
#
#      Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
#           Date: 14-01-2010
#  License Terms: GNU LGPL
#
###################################################################################################

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc

dist_pkgdata_DATA = CarotidBifurcation_1D-3DFluid_FluxStress.dat \
                    CarotidBifurcation_1D-3DFSI_FluxStress.dat \
                    CarotidBifurcation_3DFluid.dat \
					CarotidBifurcation_3DFluid_FluxStress.dat \
                    CarotidBifurcation_3DFluid_Stress.dat \
                    CarotidBifurcation_3DFluid_Synthetic_FluxStress.dat \
                    CarotidBifurcation_3DFluid_Synthetic_Stress.dat \
                    ReymondStergiopulos2009_1D_Absorbing_Stress.dat \
                    ReymondStergiopulos2009_1D_Resistance_Stress.dat \
                    ReymondStergiopulos2009_1D_Aorta_Absorbing_Stress.dat \
                    ReymondStergiopulos2009_1D_Aorta_Resistance_Stress.dat \
                    ReymondStergiopulos2009_1D_AortaArmsLegs_Absorbing_Stress.dat \
                    ReymondStergiopulos2009_1D_AortaArmsLegs_Resistance_Stress.dat \
                    SinusoidalWave_1D_1Cylinder.dat \
                    SinusoidalWave_1D_2Cylinders_Stress.dat \
                    SinusoidalWave_1D_2Cylinders_FluxStress.dat \
                    SinusoidalWave_1D_3Cylinders_Stress.dat \
                    SinusoidalWave_1D_3Cylinders_FluxStress.dat \
                    SinusoidalWave_1D_Hexagon_Stress.dat \
                    SinusoidalWave_1D_Hexagon_FluxStress.dat \
                    SinusoidalWave_1D_Hexagon_Mixed.dat \
                    SinusoidalWave_1D-3DFluid_2Cylinders_FluxStress.dat \
                    SinusoidalWave_1D-3DFluid_2Cylinders_Stress.dat \
                    SinusoidalWave_1D-3DFSI_2Cylinders_FluxStress.dat \
                    SinusoidalWave_3DFluid-1D_2Cylinders_FluxStress.dat \
                    SinusoidalWave_3DFluid-1D_2Cylinders_Stress.dat \
                    SinusoidalSquareWave_3DFSI_1Cylinder_tube20.dat \
                    SinusoidalSquareWave_3DFSI_2Cylinders_FluxStress_tube20.dat \
                    SinusoidalSquareWave_3DFSI_2Cylinders_Stress_tube20.dat \
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
                    Womersley_3DFluid_5Cylinders_Stress_tube80.dat \
                    Template.dat

check-local:
	@if test "$(top_srcdir)" != "$(top_builddir)"; then\
	 for i in $(dist_pkgdata_DATA); do \
           $(LN_S) -f $(srcdir)/$$i .; \
         done; \
	fi

