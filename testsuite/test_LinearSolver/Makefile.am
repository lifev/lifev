# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications
#                Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA
#
#      Author(s): Name Surname <gwenol.grandperrin@epfl.ch>
#           Date: 09-01-2012
#  License Terms: GNU GPL
#
###################################################################################################

include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES                  = .cpp .hpp .idl .c .h .f .F .o .moc

check_HEADERS             =  # List of hpp files

check_PROGRAMS            = test_LinearSolver
test_LinearSolver_SOURCES = main.cpp

EXTRA_DIST = data

link:
	test -a data || ln -s $(srcdir)/data
	test -a SolverParamList.xml || ln -s $(srcdir)/SolverParamList.xml

clean-results:
	-rm *.vct *.scl *.case *.geo *.xmf *.h5
