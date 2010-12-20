# -*- Mode : makefile -*-
#
#    AUTHOR: Christophe Prud'homme
#       ORG: EPFL
#    E-MAIL: christophe.prudhomme@epfl.ch
#
# DESCRIPTION:
# ============
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
# DESCRIP-END.
#

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc .cc

lib_LTLIBRARIES		= liblifefem.la

liblifefem_la_SOURCES	= bcCond.cpp bcFunction.cpp bcVector.cpp bcHandler.cpp \
          dataTime.cpp  \
          CurrentFE.cpp CurrentBoundaryFE.cpp defQuadRuleFE.cpp DOF.cpp \
          DOFInterface3Dto2D.cpp DOFInterface3Dto3D.cpp DOFInterface.cpp \
          AssemblyElemental.cpp  GeometricMap.cpp \
          localDofPattern.cpp quadRule.cpp refEle.cpp refFE.cpp refFEScalar.cpp refFEHdiv.cpp \
          refFEHybrid.cpp staticBdFE.cpp \
          bcManage.cpp quadPoint.cpp

liblifefem_la_LIBADD	= ../lifemesh/liblifemesh.la ../lifealg/liblifealg.la


include $(top_srcdir)/Doxyfile.am
