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

AM_CPPFLAGS = -DNO_BLASLAPACK
INCLUDES 	= -I$(top_srcdir)/life/lifecore  -I$(top_srcdir)/life/lifearray -I$(top_srcdir)/life/lifealg \
             -I$(top_srcdir)/life/lifemesh -I$(top_srcdir)/life/lifefem -I$(top_srcdir)/life/lifefilters 

include_HEADERS	= ElasticStructureHandler.hpp NavierStokesAleHandler.hpp \
                NavierStokesAleSolverPC.hpp NavierStokesHandler.hpp \
                NavierStokesSolverPC.hpp VenantKirchhofSolver.hpp dataElasticStructure.hpp \
                dataNavierStokes.hpp dataSimplifiedStructure.hpp simplifiedStructure.hpp \
                convDiffReactHandler.hpp convDiffReactSolverPC.hpp dataConvDiffReact.hpp \
                dataDarcy.hpp darcyHandler.hpp darcySolverBase.hpp darcySolver.hpp \
                operFS.hpp steklovPoincareBase.hpp fixedPointBase.hpp exactJacobianBase.hpp


lib_LTLIBRARIES	= liblifesolver.la

liblifesolver_la_SOURCES	= darcySolver.cpp operFS.cpp steklovPoincareBase.cpp fixedPointBase.cpp exactJacobianBase.cpp 

liblifesolver_la_LIBADD    = ../lifefilters/liblifefilters.la ../lifearray/liblifearray.la ../lifecore/liblifecore.la $(petsc_libs) $(aztec_libs) $(boost_libs)

#bin_PROGRAMS=multiscale
#multiscale_SOURCES  = #multiscale.cpp
#multiscale_LDADD    = ../lifecore/liblifecore.la 

include $(top_srcdir)/admin/Doxyfile.am
