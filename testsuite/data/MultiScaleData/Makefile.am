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

SUBDIRS		= 	Algorithm\
                Couplings\
				Models

dist_pkgdata_DATA = Run_CarotidBifurcation_1D-3DFluid_FluxStress.dat \
					Run_CarotidBifurcation_1D-3DFSI_FluxStress.dat \
					Run_CarotidBifurcation_3DFluid.dat \
					Run_CarotidBifurcation_3DFluid_FluxStress.dat \
					Run_CarotidBifurcation_3DFluid_Stress.dat \
					Run_CarotidBifurcation_3DFluid_Synthetic_FluxStress.dat \
					Run_CarotidBifurcation_3DFluid_Synthetic_Stress.dat \
					Run_ReymondStergiopulos2009_1D_Aorta.dat \
					Run_ReymondStergiopulos2009_1D_AortaArmsLegs.dat \
					Run_ReymondStergiopulos2009_1D_ArterialTree.dat \
					Run_Stationary_3DFluid_1Cylinder.dat \
					Run_Stationary_3DFluid_CarotidBifurcation_FluxStress.dat \
					Run_Stationary_3DFluid_Hexagon_FluxStress.dat \
					Run_Stationary_3DFluid_Hexagon_Stress.dat \
					Run_SinusoidalSquareWave_1D_2Cylinders_Stress.dat \
					Run_SinusoidalWave_1D_1Cylinder.dat \
					Run_SinusoidalWave_1D_2Cylinders_Stress.dat \
					Run_SinusoidalWave_1D_2Cylinders_FluxStress.dat \
					Run_SinusoidalWave_1D_3Cylinders_Stress.dat \
					Run_SinusoidalWave_1D_3Cylinders_FluxStress.dat \
					Run_SinusoidalWave_1D_Hexagon_Stress.dat \
					Run_SinusoidalWave_1D_Hexagon_FluxStress.dat \
					Run_SinusoidalWave_1D-3DFluid_2Cylinders_FluxStress.dat \
					Run_SinusoidalWave_1D-3DFluid_2Cylinders_Stress.dat \
					Run_SinusoidalWave_3DFluid-1D_2Cylinders_FluxStress.dat \
					Run_SinusoidalWave_3DFluid-1D_2Cylinders_Stress.dat \
					Run_SinusoidalWave_1D-3DFSI_2Cylinders_FluxStress.dat \
					Run_SinusoidalSquareWave_3DFSI_1Cylinder.dat \
					Run_SinusoidalSquareWave_3DFSI_2Cylinders_FluxStress.dat \
					Run_SinusoidalSquareWave_3DFSI_2Cylinders_Stress.dat \
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
					Run_Womersley_3DFluid_5Cylinders_Stress_tube80.dat \
					Test_Fluid3D.dat \
					Test_1D.dat \
					Test_FSI3D.dat

check-local:
	@if test "$(top_srcdir)" != "$(top_builddir)"; then\
	 for i in $(dist_pkgdata_DATA); do \
           $(LN_S) -f $(srcdir)/$$i .; \
         done; \
	fi

