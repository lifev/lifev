# -*- makefile -*-
#
# Author: Christophe Prud'homme
#
# Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc 

INCLUDES		= -I$(top_srcdir)/life/lifecore -I$(top_srcdir)/life/lifearray  -I$(top_srcdir)/life/lifealg


check_HEADERS		= MatrixTest.hpp

check_LTLIBRARIES	= liblifealg_testsuite.la

liblifealg_testsuite_la_SOURCES	= MatrixTest.cpp
liblifealg_testsuite_la_LIBADD	= $(top_builddir)/life/lifealg/liblifealg.la 

check_PROGRAMS	= test_linearsolvers

test_linearsolvers_SOURCES	= test_linearsolvers.cpp
test_linearsolvers_LDADD	= liblifealg_testsuite.la


