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

check_PROGRAMS = test_onedmodel

noinst_HEADERS = ud_functions.hpp

test_onedmodel_SOURCES = main.cpp 

test_onedmodel_LDADD    = $(lifemc_libs) $(life_libs) 

EXTRA_DIST = data

clean_PROGRAMS:
	-rm $(check_PROGRAMS)

recheck: clean_PROGRAMS check

clean-results:
	-rm solA.m solP.m solQ.m solW1.m solW2.m *.xmf *.h5

include $(top_srcdir)/Doxyfile.am
