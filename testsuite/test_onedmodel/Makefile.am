# -*- makefile -*-
#
# Author: Vincent Martin
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

check_PROGRAMS = test_onedmodel
noinst_HEADERS = dataOneDModel.hpp basicOneDMesh.hpp dofOneD.hpp \
		 oneDNonLinModelParam.hpp vectorFunction1D.hpp \
		 oneDModelHandler.hpp oneDModelSolver.hpp \
		 ud_functions.hpp gracePlot.hpp 

test_onedmodel_SOURCES = dataOneDModel.cpp  basicOneDMesh.cpp \
			 oneDNonLinModelParam.cpp vectorFunction1D.cpp \
			 oneDModelHandler.cpp oneDModelSolver.cpp \
			 main.cpp gracePlot.cpp 

test_onedmodel_LDADD    = $(grace_libs) $(life_libs)

EXTRA_DIST = data datalin2 datalin3 datanl testsuite.at 

#

