# -*- makefile -*-
#
# Author: Martin Prosi
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

include $(top_srcdir)/testsuite/Makefile.testsuite

check_HEADERS = ud_functions.hpp ensight7-write.hpp

check_PROGRAMS = test_mass_coupling
test_mass_coupling_SOURCES = main.cpp 
EXTRA_DIST = data-lumen data-wall testsuite.at



