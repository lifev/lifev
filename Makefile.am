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

liblifemcsolver_la_SOURCES   = Monolithic.cpp \
                               fullMonolithic.cpp \
                               dataOneDModel.cpp \
                               vectorFunction1D.cpp \
                               oneDNonLinModelParam.cpp \
                               MS_Solver.cpp \
                               MS_PhysicalData.cpp \
                               MS_PhysicalModel.cpp \
                               MS_Model_MultiScale.cpp \
                               MS_Model_Fluid3D.cpp \
                               MS_PhysicalCoupling.cpp \
                               MS_Coupling_BoundaryCondition.cpp \
                               MS_Coupling_FluxStress.cpp

liblifemcsolver_la_LIBADD    = -llifesolver

include $(top_srcdir)/admin/Doxyfile.am
