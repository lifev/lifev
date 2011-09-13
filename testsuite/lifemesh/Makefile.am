# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University      
#
#      Author(s): Antonio cervone <ant.cervone@gmail.com>
#           Date: 15-06-2011
#  License Terms: GNU LGPL
#
###################################################################################################

include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc 

check_PROGRAMS		= test_mesh test_baremesh

test_mesh_SOURCES	= test_mesh.cpp

test_mesh_LDFLAGS       =

test_baremesh_SOURCES	= test_baremesh.cpp

test_baremesh_LDFLAGS   =

MOSTLYCLEANFILES	= testBuilders.dat testBuilders_baremesh.dat

EXTRA_DIST = data data_baremesh testsuite.at test_mesh.dox
