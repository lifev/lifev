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

dist_pkgdata_DATA = 1D_1Cylinder.dat \
                    1D_2Cylinders_FS.dat \
                    1D_2Cylinders_S.dat \
                    1D_ReymondStergiopulos2009-Aorta_S.dat \
                    1D_ReymondStergiopulos2009-Aorta-Absorbing_S.dat \
                    1D_ReymondStergiopulos2009-AortaArmsLegs_S.dat \
                    1D_ReymondStergiopulos2009-AortaArmsLegs-Absorbing_S.dat \
                    1D_ReymondStergiopulos2009-ArterialTree_S.dat \
                    1D_ReymondStergiopulos2009-ArterialTree-Absorbing_S.dat \
                    1D_Rhombus_FS.dat \
                    1D_Rhombus_FS-S.dat \
                    1D_Rhombus_S.dat \
                    1D-Fluid3D_2Cylinders_FS.dat \
                    1D-Fluid3D_2Cylinders_S.dat \
                    1D-FSI3D_2Cylinders_FS.dat \
                    1D-FSI3D_2Cylinders_S.dat \
                    1D-FSI3D_3Cylinders_FS.dat \
                    1D-FSI3D_3Cylinders_S.dat \
                    1D-Windkessel0D_1Cylinder_FS.dat \
                    1D-Windkessel0D_1Cylinder_S.dat \
                    Fluid3D_1Cylinder-Stokes-tube20.dat \
                    Fluid3D_1Cylinder-tube20.dat \
                    Fluid3D_5Cylinders-tube20_FS.dat \
                    Fluid3D_5Cylinders-tube20_S.dat \
                    Fluid3D_5Cylinders-tube40_FS.dat \
                    Fluid3D_5Cylinders-tube40_S.dat \
                    Fluid3D_5Cylinders-tube60_FS.dat \
                    Fluid3D_5Cylinders-tube60_S.dat \
                    Fluid3D_5Cylinders-tube80_FS.dat \
                    Fluid3D_5Cylinders-tube80_S.dat \
                    Fluid3D_CarotidBifurcation.dat \
                    Fluid3D_CarotidBifurcationAnd3Cylinders_FS.dat \
                    Fluid3D_CarotidBifurcationAnd3Cylinders_S.dat \
                    Fluid3D_CarotidBifurcationSynthetic_FS.dat \
                    Fluid3D_CarotidBifurcationSynthetic_S.dat \
                    Fluid3D_Rhombus-tube20_F-FS-S.dat \
                    Fluid3D_Rhombus-tube20_FS.dat \
                    Fluid3D_Rhombus-tube20_FS-S.dat \
                    Fluid3D_Rhombus-tube20_S.dat \
                    Fluid3D-1D_2Cylinders_FS.dat \
                    Fluid3D-1D_2Cylinders_S.dat \
                    Fluid3D-1D_3Cylinders_FV.dat \
                    Fluid3D-1D_CarotidBifurcation_FS.dat \
                    FSI3D_1Cylinder-tube20.dat \
                    FSI3D_2Cylinders-tube20_FS.dat \
                    FSI3D_2Cylinders-tube20_S.dat \
                    FSI3D_3Cylinders-tube20_FS.dat \
                    FSI3D-1D_CarotidBifurcation_FS.dat \
                    FSI3D-1D_IdealLeftVentricleWithAorta_S.dat \
                    FSI3D-1D_IdealLeftVentricleWithAorta-Absorbing_S.dat \
                    FSI3D-1D_IdealLeftVentricleWithAorta-ClosedValve_S.dat \
                    FSI3D-1D_IdealLeftVentricleWithArterialTree_S.dat \
                    Template.dat \
                    Windkessel0D_Test.dat

check-local:
	@if test "$(top_srcdir)" != "$(top_builddir)"; then\
	 for i in $(dist_pkgdata_DATA); do \
           $(LN_S) -f $(srcdir)/$$i .; \
         done; \
	fi

