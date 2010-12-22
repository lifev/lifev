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

check_HEADERS = ensightToHdf5.hpp

check_PROGRAMS = test_filters
test_filters_SOURCES = main.cpp ensightToHdf5.cpp

link:
	test -a Mesh || ln -s $(top_srcdir)/testsuite/data/mesh/inria/ Mesh
	test -a importDir || ln -s $(top_srcdir)/testsuite/lifefilters/ importDir

EXTRA_DIST = data testsuite.at 