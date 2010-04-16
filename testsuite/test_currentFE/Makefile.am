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

check_HEADERS              = # List of all the hpp files

check_PROGRAMS             = test_currentFE
test_currentFE_SOURCES  = main.cpp # List of all the cpp files

EXTRA_DIST = testsuite.at

recheck: clean_PROGRAMS check

clean_PROGRAMS:
	-rm $(check_PROGRAMS)

