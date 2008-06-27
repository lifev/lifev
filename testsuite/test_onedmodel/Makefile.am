# -*- makefile -*-
#
# Author: Vincent Martin
# Author: Tiziano Passerini
#
# Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano
#
#  Distributed under the GPL(GNU Public License):
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc 

AM_CPPFLAGS	= -DNO_BLASLAPACK

include $(top_srcdir)/testsuite/Makefile.testsuite

check_PROGRAMS = test_onedmodel test_onednet # test_twotubes 

noinst_HEADERS = interface2Vessels.hpp ud_functions.hpp

test_onedmodel_SOURCES = ud_functions.cpp main.cpp 

#	$(top_srcdir)/life/lifefem/oneDBCFunctions.cpp \
#	$(top_srcdir)/life/lifefem/oneDBCHandler.cpp \
#	$(top_srcdir)/life/lifesolver/dataOneDModel.cpp \
#	$(top_srcdir)/life/lifesolver/oneDModelHandler.cpp \
#	$(top_srcdir)/life/lifesolver/oneDModelSolver.cpp \
#	$(top_srcdir)/life/lifesolver/oneDNonLinModelParam.cpp \
#	$(top_srcdir)/life/lifesolver/vectorFunction1D.cpp 

# test_twotubes_SOURCES = interface2Vessels.cpp main_twotubes.cpp

test_onednet_SOURCES = ud_functions.cpp main_net.cpp

#	$(top_srcdir)/life/lifefem/oneDBCFunctions.cpp \
#	$(top_srcdir)/life/lifefem/oneDBCHandler.cpp \
#	$(top_srcdir)/life/lifesolver/dataOneDModel.cpp \
#	$(top_srcdir)/life/lifesolver/oneDModelHandler.cpp \
#	$(top_srcdir)/life/lifesolver/oneDModelSolver.cpp \
#	$(top_srcdir)/life/lifesolver/oneDNonLinModelParam.cpp \
#	$(top_srcdir)/life/lifesolver/vectorFunction1D.cpp 

test_onedmodel_LDADD    = $(grace_libs) $(life_libs)
# test_twotubes_LDADD     = $(grace_libs) $(life_libs)
test_onednet_LDADD     = $(grace_libs) $(life_libs)

EXTRA_DIST = data datalin2 datalin3 datanl datanl2 datanet testsuite.at 

MOSTLYCLEANFILES = *.res* *.mesh *.bb *.geo *.out *.m *.mtv *.scl *.case *.lbm *.mesh *.scl *.geo

#
