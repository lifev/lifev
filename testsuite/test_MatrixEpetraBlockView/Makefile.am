# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University      
#
#      Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
#           Date: 04-01-2012
#  License Terms: GNU LGPL
#
###################################################################################################

include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES                           = .cpp .hpp .idl .c .h .f .F .o .moc

check_PROGRAMS                     = test_MatrixEpetraBlockView
test_MatrixEpetraBlockView_SOURCES = main.cpp

recheck: clean_PROGRAMS check

clean_PROGRAMS:
	-rm $(check_PROGRAMS)

clean-results:
	-rm *.vct *.scl *.case *.geo
