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

dist_pkgdata_DATA = Carotid.dat \
                    Carotid_Stokes.dat \
                    Cylinder_tube20.dat \
					Cylinder_tube40.dat \
					Cylinder_tube60.dat \
					Cylinder_tube80.dat \
					Cylinder_Stokes_tube20.dat \
					Cylinder_Stokes_tube40.dat \
					Cylinder_Stokes_tube60.dat \
					Cylinder_Stokes_tube80.dat \
					Template.dat

check-local:
	@if test "$(top_srcdir)" != "$(top_builddir)"; then\
	 for i in $(dist_pkgdata_DATA); do \
           $(LN_S) -f $(srcdir)/$$i .; \
         done; \
	fi

