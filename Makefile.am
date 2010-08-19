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

liblifemcsolver_la_SOURCES   = BCInterface_FSI.cpp \
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
                               MS_Solver.cpp \
                               MS_Algorithm.cpp \
                               MS_Algorithm_Aitken.cpp \
                               MS_Algorithm_Explicit.cpp \
                               MS_Algorithm_Newton.cpp \
                               MS_PhysicalData.cpp \
                               MS_PhysicalModel.cpp \
                               MS_Model_MultiScale.cpp \
                               MS_Model_Fluid3D.cpp \
                               MS_Model_FSI3D.cpp \
                               MS_Model_1D.cpp \
                               MS_PhysicalCoupling.cpp \
                               MS_Coupling_BoundaryCondition.cpp \
                               MS_Coupling_Stress.cpp \
                               MS_Coupling_FluxStress.cpp \
                               OneDimensionalModel_Data.cpp \
                               OneDimensionalModel_Physics.cpp \
                               OneDimensionalModel_Physics_Linear.cpp \
                               OneDimensionalModel_Physics_NonLinear.cpp \
                               OneDimensionalModel_Flux.cpp \
                               OneDimensionalModel_Flux_Linear.cpp \
                               OneDimensionalModel_Flux_NonLinear.cpp \
                               OneDimensionalModel_Solver.cpp \
                               OneDimensionalModel_Source.cpp \
                               OneDimensionalModel_Source_Linear.cpp \
                               OneDimensionalModel_Source_NonLinear.cpp

liblifemcsolver_la_LIBADD    = -llifesolver

include $(top_srcdir)/admin/Doxyfile.am
