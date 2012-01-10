# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA          
#
#      Author(s): Samuel Quinodoz <samuel.quinodoz@epfl.ch>
#                 Umberto Villa <uvilla@emory.edu>
#           Date: 06-04-2010
#  License Terms: GNU GPL
#
###################################################################################################

include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES                   = .cpp .hpp .idl .c .h .f .F .o .moc

check_HEADERS              = test_quadrule.hpp SetOfFun.hpp test_interpolate.hpp

check_PROGRAMS             = test_currentFE test_quadRule test_interpolate

test_currentFE_SOURCES   = main.cpp # List of all the cpp files
test_quadRule_SOURCES    = test_quadrule.cpp SetOfFun.cpp
test_interpolate_SOURCES = test_interpolate.cpp

link:
	test -a data || ln -s $(srcdir)/data
	test -a Mesh || ln -s $(top_srcdir)/testsuite/data/mesh/inria/ Mesh

EXTRA_DIST = data testsuite.at

recheck: clean_PROGRAMS check

clean_PROGRAMS:
	-rm $(check_PROGRAMS)

