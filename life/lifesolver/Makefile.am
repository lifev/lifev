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

AM_CPPFLAGS = -DNO_BLASLAPACK

lib_LTLIBRARIES	= liblifesolver.la

liblifesolver_la_SOURCES   = FSIData.cpp ADRData.cpp VenantKirchhoffElasticData.cpp  OseenData.cpp FSIOperator.cpp FSIFixedPoint.cpp FSIMonolithic.cpp FSIMonolithicGE.cpp FSIMonolithicGI.cpp \
			                 FSIExactJacobian.cpp FSISolver.cpp HeartIonicData.cpp HeartBidomainData.cpp HeartMonodomainData.cpp HeartFunctors.cpp \
					LevelSetData.cpp VenantKirchhoffViscoelasticData.cpp MonolithicBlock.cpp MonolithicBlockComposed.cpp MonolithicBlockComposedDN.cpp \
                    MonolithicBlockComposedDND.cpp MonolithicBlockComposedDNND.cpp MonolithicBlockComposedNN.cpp MonolithicBlockMatrix.cpp \
                    MonolithicBlockMatrixRN.cpp MonolithicRobinInterface.cpp\
                    OneDimensionalData.cpp OneDimensionalFluxLinear.cpp OneDimensionalFluxNonLinear.cpp\
                    OneDimensionalSourceLinear.cpp OneDimensionalSourceNonLinear.cpp OneDimensionalPhysicsLinear.cpp\
                    OneDimensionalPhysicsNonLinear.cpp OneDimensionalSolver.cpp

liblifesolver_la_LIBADD    = ../lifefilters/liblifefilters.la ../lifearray/liblifearray.la ../lifecore/liblifecore.la ../lifefunctions/liblifefunctions.la $(petsc_libs) $(aztec_libs) $(boost_libs)

include $(top_srcdir)/Doxyfile.am
