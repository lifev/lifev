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
                               Monolithic.cpp \
                               MonolithicGI.cpp \
                               MonolithicGE.cpp \
                               RobinInterface.cpp \
                               BlockMatrix.cpp \
                               BlockMatrixRN.cpp \
                               BlockInterface.cpp \
                               ComposedBlockOper.cpp \
                               ComposedDN.cpp \
                               ComposedDND.cpp \
                               ComposedNN.cpp \
                               ComposedDNND.cpp \
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
                               OneDimensionalModel_Data.cpp \
                               OneDimensionalModel_Physics.cpp \
                               OneDimensionalModel_Physics_Linear.cpp \
                               OneDimensionalModel_Physics_NonLinear.cpp \
                               OneDimensionalModel_Flux_Linear.cpp \
                               OneDimensionalModel_Flux_NonLinear.cpp \
                               OneDimensionalModel_Solver.cpp \
                               OneDimensionalModel_Source_Linear.cpp \
                               OneDimensionalModel_Source_NonLinear.cpp

liblifemcsolver_la_LIBADD    = -llifesolver

include $(top_srcdir)/Doxyfile.am
