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

include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc 

INCLUDES	= -I$(top_srcdir)/life/lifearray

check_PROGRAMS	= test_meshentitycontainer test_vectorsmall test_ghosthandler

test_meshentitycontainer_SOURCES= test_meshentitycontainer.cpp

test_vectorsmall_SOURCES= test_vectorsmall.cpp

test_ghosthandler_SOURCES= test_ghosthandler.cpp

EXTRA_DIST  = testsuite.at

link:
	test -a data_ghosthandler || ln -s $(srcdir)/data_ghosthandler
	test -a Mesh || ln -s $(top_srcdir)/testsuite/data/mesh/inria Mesh
