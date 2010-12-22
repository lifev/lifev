# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University      
#
#      Author(s): Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
#           Date: 2010-07-29
#  License Terms: GNU LGPL
#
###################################################################################################

SUFFIXES        = .cpp .hpp .idl .c .h .f .F .o .moc

include $(top_srcdir)/testsuite/Makefile.testsuite

check_HEADERS = hyperbolic.hpp
check_PROGRAMS = test_hyperbolic

test_hyperbolic_SOURCES = main.cpp hyperbolic.cpp

EXTRA_DIST = data

link:
	test -a data || ln -s $(srcdir)/data
	test -a Mesh || ln -s $(top_srcdir)/testsuite/data/mesh/inria Mesh

clean-results:
	-rm *.vct *.scl *.case *.geo *.h5 *.xmf
