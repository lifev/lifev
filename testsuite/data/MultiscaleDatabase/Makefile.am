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

SUBDIRS		= 	Algorithms\
                Couplings\
				Models

dist_pkgdata_DATA = Run_1D_1Cylinder.dat \
					Run_1D_2Cylinders.dat \
					Run_1D_ReymondStergiopulos2009_Aorta.dat \
					Run_1D_ReymondStergiopulos2009_AortaArmsLegs.dat \
					Run_1D_ReymondStergiopulos2009_ArterialTree.dat \
					Run_1D_Rhombus.dat \
					Run_1D-Fluid3D_2Cylinders.dat \
					Run_Fluid3D_1Cylinder.dat \
					Run_Fluid3D_5Cylinders.dat \
					Run_Fluid3D_CarotidBifurcation.dat \
					Run_Fluid3D_CarotidBifurcationSynthetic.dat \
					Run_Fluid3D_Rhombus.dat \
					Run_Fluid3D-1D_2Cylinders.dat \
					Run_Fluid3D-1D_CarotidBifurcation.dat \
					Run_FSI3D_1Cylinder.dat \
					Run_FSI3D_2Cylinders.dat \
					Run_FSI3D-1D_2Cylinders.dat \
					Run_FSI3D-1D_CarotidBifurcation.dat \
					Test_Fluid3D.dat \
					Test_1D.dat \
					Test_FSI3D.dat

check-local:
	@if test "$(top_srcdir)" != "$(top_builddir)"; then\
	 for i in $(dist_pkgdata_DATA); do \
           $(LN_S) -f $(srcdir)/$$i .; \
         done; \
	fi

