# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University      
#
#      Author(s): Samuel Quinodoz <samuel.quinodoz@epfl.ch>
#           Date: 00-00-0000
#  License Terms: GNU LGPL
#
###################################################################################################

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc

include $(top_srcdir)/testsuite/Makefile.testsuite

check_HEADERS = 

check_PROGRAMS = test_ADRAssembler test_ADRAssembler_2d
test_ADRAssembler_SOURCES = main.cpp
test_ADRAssembler_2d_SOURCES = main_2d.cpp

EXTRA_DIST = data testsuite.at

link:
	test -a data || ln -s $(srcdir)/data
	test -a data_2d || ln -s $(srcdir)/data_2d
	test -a Mesh || ln -s $(top_srcdir)/testsuite/data/mesh Mesh

clean-results:
	-rm *.vct *.scl *.case *.geo

