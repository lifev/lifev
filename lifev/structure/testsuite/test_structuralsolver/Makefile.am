# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University      
#
#      Author(s): Name Surname <name.surname@epfl.ch>
#           Date: 00-00-0000
#  License Terms: GNU LGPL
#
###################################################################################################

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc

include $(top_srcdir)/testsuite/Makefile.testsuite

check_HEADERS = ud_functions.hpp

check_PROGRAMS = test_structure
test_structure_SOURCES = main.cpp  ud_functions.cpp

EXTRA_DIST = data testsuite.at

link:
	test -a data || ln -s $(srcdir)/data
	test -a Mesh || ln -s $(top_srcdir)/testsuite/data/mesh/inria/ Mesh

clean-results:
	-rm velocity*.vct displacement*.vct structure.*.case structure.*.geo *.xmf *.h5 

recheck:
	-rm $(check_PROGRAMS)
	make check
