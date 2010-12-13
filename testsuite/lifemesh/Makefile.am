# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University      
#
#      Author(s): 
#           Date: 00-00-0000
#  License Terms: GNU LGPL
#
###################################################################################################

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc 

check_PROGRAMS		= test_mesh

test_mesh_SOURCES	= test_mesh.cpp
 
test_mesh_LDFLAGS       =

MOSTLYCLEANFILES	= testBuilders.dat

EXTRA_DIST = data testsuite.at test_mesh.dox

include $(top_srcdir)/testsuite/Makefile.testsuite
