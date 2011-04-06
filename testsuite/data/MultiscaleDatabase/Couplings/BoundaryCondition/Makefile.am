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

SUFFIXES    = .cpp .hpp .idl .c .h .f .F .o .moc

dist_pkgdata_DATA =   Absorbing_1DLeft.dat \
                      Absorbing_1DRight.dat \
                      Absorbing_NobileVergara.dat \
                      Area_One.dat \
                      FlowRate_CarotidBifurcation.dat \
                      FlowRate_Exp.dat \
                      FlowRate_Heart.dat \
                      FlowRate_Mitral.dat \
                      FlowRate_One.dat \
                      FlowRate_Poiseuille.dat \
                      FlowRate_PoiseuilleDynamicSTP.dat \
                      FlowRate_Sinusoidal.dat \
                      FlowRate_SinusoidalWave.dat \
                      FlowRate_SquareSinusoidalWave.dat \
                      FlowRate_Test.dat \
                      FlowRate_Womersley.dat \
                      FlowRate_Zero.dat \
                      Stress_One.dat \
                      Stress_PoiseuilleDynamicSTP.dat \
                      Stress_Sinusoidal.dat \
                      Stress_Test.dat \
                      Stress_Womersley.dat \
                      Stress_Zero.dat \
                      Template.dat \
                      Terminal_ReymondStergiopulos2009_003.dat \
                      Terminal_ReymondStergiopulos2009_008.dat \
                      Terminal_ReymondStergiopulos2009_010.dat \
                      Terminal_ReymondStergiopulos2009_011.dat \
                      Terminal_ReymondStergiopulos2009_015.dat \
                      Terminal_ReymondStergiopulos2009_019.dat \
                      Terminal_ReymondStergiopulos2009_022.dat \
                      Terminal_ReymondStergiopulos2009_024.dat \
                      Terminal_ReymondStergiopulos2009_025.dat \
                      Terminal_ReymondStergiopulos2009_026.dat \
                      Terminal_ReymondStergiopulos2009_031.dat \
                      Terminal_ReymondStergiopulos2009_032.dat \
                      Terminal_ReymondStergiopulos2009_033.dat \
                      Terminal_ReymondStergiopulos2009_034.dat \
                      Terminal_ReymondStergiopulos2009_036.dat \
                      Terminal_ReymondStergiopulos2009_038.dat \
                      Terminal_ReymondStergiopulos2009_040.dat \
                      Terminal_ReymondStergiopulos2009_041.dat \
                      Terminal_ReymondStergiopulos2009_045.dat \
                      Terminal_ReymondStergiopulos2009_047.dat \
                      Terminal_ReymondStergiopulos2009_048.dat \
                      Terminal_ReymondStergiopulos2009_049.dat \
                      Terminal_ReymondStergiopulos2009_051.dat \
                      Terminal_ReymondStergiopulos2009_053.dat \
                      Terminal_ReymondStergiopulos2009_054.dat \
                      Terminal_ReymondStergiopulos2009_055.dat \
                      Terminal_ReymondStergiopulos2009_057.dat \
                      Terminal_ReymondStergiopulos2009_058.dat \
                      Terminal_ReymondStergiopulos2009_064.dat \
                      Terminal_ReymondStergiopulos2009_065.dat \
                      Terminal_ReymondStergiopulos2009_071.dat \
                      Terminal_ReymondStergiopulos2009_072.dat \
                      Terminal_ReymondStergiopulos2009_074.dat \
                      Terminal_ReymondStergiopulos2009_075.dat \
                      Terminal_ReymondStergiopulos2009_076.dat \
                      Terminal_ReymondStergiopulos2009_078.dat \
                      Terminal_ReymondStergiopulos2009_080.dat \
                      Terminal_ReymondStergiopulos2009_082.dat \
                      Terminal_ReymondStergiopulos2009_084.dat \
                      Terminal_ReymondStergiopulos2009_086.dat \
                      Terminal_ReymondStergiopulos2009_088.dat \
                      Terminal_ReymondStergiopulos2009_090.dat \
                      Terminal_ReymondStergiopulos2009_091.dat \
                      Terminal_ReymondStergiopulos2009_092.dat \
                      Terminal_ReymondStergiopulos2009_093.dat \
                      Terminal_ReymondStergiopulos2009_094.dat \
                      Terminal_ReymondStergiopulos2009_096.dat \
                      Terminal_ReymondStergiopulos2009_098.dat \
                      Terminal_ReymondStergiopulos2009_099.dat \
                      Terminal_ReymondStergiopulos2009_100.dat \
               		  Terminal_ReymondStergiopulos2009_102.dat \
                      Velocity_Poiseuille.dat \
                      Velocity_PoiseuilleDynamicSTP.dat \
                      Velocity_Zero.dat \
                      Velocity_ZeroXY.dat

check-local:
	@if test "$(top_srcdir)" != "$(top_builddir)"; then\
	 for i in $(dist_pkgdata_DATA); do \
           $(LN_S) -f $(srcdir)/$$i .; \
         done; \
	fi

