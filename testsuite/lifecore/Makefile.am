# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Library                       
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University
#
#      Author(s): Umberto Villa <uvilla@emory.edu>
#           Date: 10-12-2010
#  License Terms: GNU LGPL
#
#################################################################################################

include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES		= .cpp .hpp .idl .c .h .f .F .o .moc 

check_HEADERS              = # List of all the hpp files

check_PROGRAMS 	    = test_getPot

test_getPot_SOURCES = main_getPot.cpp

EXTRA_DIST = data_getPot

link:
	test -a data_getPot || ln -s $(srcdir)/data_getPot
