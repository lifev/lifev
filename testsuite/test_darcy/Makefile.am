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

AM_CPPFLAGS	= -DNO_BLASLAPACK

check_HEADERS = darcyHandler.hpp darcySolver.hpp dataDarcy.hpp user_fct.hpp user_diffusion.hpp 


check_PROGRAMS = test_darcy
test_darcy_SOURCES = darcyHandler.cpp darcySolver.cpp darcySolver_post.cpp dataDarcy.cpp main.cpp user_fct.cpp user_diffusion.cpp 

EXTRA_DIST = data


include $(top_srcdir)/testsuite/Makefile.testsuite
