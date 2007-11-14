# -*- makefile -*-
#
# Author: Christophe Prud'homme
#
# Copyright (C) 2004 EPFL
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

if HAVE_PETSC
PETSC_SRCS=SolverPETSC.cpp
endif
if HAVE_UMFPACK
UMFPACK_SRCS=SolverUMFPACK.cpp
endif

lib_LTLIBRARIES	= liblifealg.la 

liblifealg_la_SOURCES	= EpetraPreconditioner.cpp SolverTrilinos.cpp  dataNewton.cpp EpetraMap.cpp $(PETSC_SRCS) $(UMFPACK_SRCS)

liblifealg_la_LIBADD	= ../lifecore/liblifecore.la ../lifearray/liblifearray.la $(petsc_libs) $(umfpack_libs) $(aztec_libs) 
include $(top_srcdir)/admin/Doxyfile.am


