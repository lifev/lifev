# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University      
#
#      Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
#           Date: 00-00-0000
#  License Terms: GNU LGPL
#
###################################################################################################

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc

include $(top_srcdir)/testsuite/Makefile.testsuite

check_HEADERS = 

check_PROGRAMS = test_LevelSetSolver
test_LevelSetSolver_SOURCES = main.cpp

EXTRA_DIST = data testsuite.at

link:
	test -a data || ln -s $(srcdir)/data

clean-results:
	-rm *.vct *.scl *.case *.geo

