# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University      
#
#      Author(s): Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
#           Date: 00-00-0000
#  License Terms: GNU LGPL
#
###################################################################################################

SUFFIXES        = .cpp .hpp .idl .c .h .f .F .o .moc

include $(top_srcdir)/testsuite/Makefile.testsuite

check_HEADERS = timeAdvance.hpp linear_function.hpp nonlinear_function.hpp

check_PROGRAMS = test_timeAdvance

test_timeAdvance_SOURCES = main.cpp timeAdvance.cpp

EXTRA_DIST = data testsuite.at

link:
	test -a data || ln -s $(srcdir)/data
	test -a Mesh || ln -s $(top_srcdir)/testsuite/data/mesh/inria Mesh

clean-results:
	-rm *.vct *.scl *.case *.geo
