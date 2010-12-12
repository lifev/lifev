# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA          
#
#      Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
#           Date:
#  License Terms: GNU GPL
#
###################################################################################################

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc 

include $(top_srcdir)/testsuite/Makefile.testsuite

check_HEADERS = test_bdf.hpp

check_PROGRAMS = test_bdf 
test_bdf_SOURCES = main.cpp  test_bdf.cpp

link:
	test -a data || ln -s $(srcdir)/data
	test -a Mesh || ln -s $(top_srcdir)/testsuite/data/mesh/inria/ Mesh

clean-results:
	-rm *.vct *.scl *.case *.geo

EXTRA_DIST = data testsuite.at
MOSTLYCLEANFILES	= *.scl *.case *.geo
