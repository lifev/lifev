# -*- makefile -*-
#
# Author: Christophe Prud'homme
#
# Copyright (C) 2004 @COPYRIGHT@
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

INCLUDES 	= -I$(top_srcdir)/life/lifecore -I$(top_srcdir)/life/lifearray -I$(top_srcdir)/life/lifefem -I$(top_srcdir)/life/lifemesh
#INCLUDES 	= -I$(top_srcdir)/life/lifecore -I$(top_srcdir)/life/lifearray 

if HAVE_PETSC
PETSC_SRCS	= SolverPETSC.cpp
PETSC_HDRS	= SolverPETSC.hpp
endif

include_HEADERS = algebraic_facto.hpp dataAztec.hpp dataNewton.hpp inexactLU.hpp linesearch_cubic.hpp linesearch_parabolic.hpp mixed_pattern.hpp newton.hpp pattern.hpp picard.hpp saddlePointCG.hpp  values.hpp $(PETSC_HDRS)

lib_LTLIBRARIES	= liblifealg.la

liblifealg_la_SOURCES = dataAztec.cpp dataNewton.cpp pattern.cpp values.cpp $(PETSC_SRCS)

liblifealg_la_LIBADD	= ../lifefem/liblifefem.la ../lifecore/liblifecore.la $(petsc_libs) $(aztec_libs)

include $(top_srcdir)/admin/Doxyfile.am
