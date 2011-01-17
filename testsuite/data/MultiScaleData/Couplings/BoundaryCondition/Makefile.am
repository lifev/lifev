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

