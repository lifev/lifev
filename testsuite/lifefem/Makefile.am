# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA          
#
#      Author(s): Name Surname <name.surname@epfl.ch>
#           Date: 00-00-0000
#  License Terms: GNU GPL
#
###################################################################################################

include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES                   = .cpp .hpp .idl .c .h .f .F .o .moc

check_HEADERS              = test_quadrule.hpp SetOfFun.hpp

check_PROGRAMS             = test_currentFE test_quadrule

test_currentFE_SOURCES  = main.cpp # List of all the cpp files
test_quadrule_SOURCES   = test_quadrule.cpp SetOfFun.cpp

EXTRA_DIST = data.txt testsuite.at

recheck: clean_PROGRAMS check

clean_PROGRAMS:
	-rm $(check_PROGRAMS)

