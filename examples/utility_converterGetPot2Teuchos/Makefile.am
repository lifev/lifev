# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University      
#
#      Author(s): Umberto Villa <uvilla@emory.edu>
#           Date: 04-10-2010
#  License Terms: GNU LGPL
#
###################################################################################################
include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES                   = .cpp .hpp .idl .c .h .f .F .o .moc

check_HEADERS              = converter.hpp

check_PROGRAMS             = converter_GetPot2Teuchos
converter_GetPot2Teuchos_SOURCES  = main.cpp converter.cpp

EXTRA_DIST = data

link:
	test -a data || ln -s $(srcdir)/data

recheck: clean_PROGRAMS check

clean_PROGRAMS:
	-rm $(check_PROGRAMS)

clean-results:
	-rm *.xml
