# -*- makefile -*-
#
# 
#
# Copyright (C) 2004 Christophe Prud'homme (christophe.prudhomme@epfl.ch)
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


noinst_HEADERS =  ud_functions.hpp 



check_PROGRAMS 	= test_parabsol
test_parabsol_SOURCES = main.cpp 

EXTRA_DIST = data testsuite.at

MOSTLYCLEANFILES=sol/u0.sol sol/u1.sol sol/u2.sol sol/u3.sol sol/u4.sol sol/u5.sol sol/u6.sol sol/u7.sol

include $(top_srcdir)/testsuite/Makefile.testsuite
