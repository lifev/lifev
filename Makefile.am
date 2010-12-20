# -*- makefile -*-
#
# Author: Simone Deparis, Tiziano Passerini
#
# Copyright (C) 2009 EPFL POLIMI
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

lib_LTLIBRARIES	= liblifemcsolver.la

liblifemcsolver_la_SOURCES   = BCInterface1DData.cpp \
                               BCInterface3DData.cpp \
                               BCInterface3DFunctionFSI.cpp \
                               FSIMonolithic.cpp \
                               FSIMonolithicGI.cpp \
                               FSIMonolithicGE.cpp \
                               RobinInterface.cpp \
                               MonolithicBlockMatrix.cpp \
                               MonolithicBlockMatrixRN.cpp \
                               MonolithicBlock.cpp \
                               MonolithicBlockComposed.cpp \
                               MonolithicBlockComposedDN.cpp \
                               MonolithicBlockComposedDND.cpp \
                               MonolithicBlockComposedNN.cpp \
                               MonolithicBlockComposedDNND.cpp \
                               MultiscaleSolver.cpp \
                               MultiscaleAlgorithm.cpp \
                               MultiscaleAlgorithmAitken.cpp \
                               MultiscaleAlgorithmExplicit.cpp \
                               MultiscaleAlgorithmNewton.cpp \
                               MultiscaleData.cpp \
                               MultiscaleModel.cpp \
                               MultiscaleModelMultiscale.cpp \
                               MultiscaleModelFluid3D.cpp \
                               MultiscaleModelFSI3D.cpp \
                               MultiscaleModel1D.cpp \
                               MultiscaleCoupling.cpp \
                               MultiscaleCouplingBoundaryCondition.cpp \
                               MultiscaleCouplingStress.cpp \
                               MultiscaleCouplingFlowRateStress.cpp \
                               OneDimensionalData.cpp \
                               OneDimensionalPhysics.cpp \
                               OneDimensionalPhysicsLinear.cpp \
                               OneDimensionalPhysicsNonLinear.cpp \
                               OneDimensionalFluxLinear.cpp \
                               OneDimensionalFluxNonLinear.cpp \
                               OneDimensionalSolver.cpp \
                               OneDimensionalSourceLinear.cpp \
                               OneDimensionalSourceNonLinear.cpp

liblifemcsolver_la_LIBADD    = -llifesolver

include $(top_srcdir)/Doxyfile.am
