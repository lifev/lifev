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

INCLUDES 	= -I$(top_srcdir)/life/lifecore -I$(top_srcdir)/life/lifearray -I$(top_srcdir)/life/lifemesh -I$(top_srcdir)/life/lifealg

include_HEADERS		= Operator01.hpp assemb.hpp\
bcCond.hpp bcFunction.hpp bcVector.hpp bcHandler.hpp \
bdf.hpp bdfNS.hpp dataTime.hpp dataTransient.hpp \
interpolate.hpp \
postProc.hpp \
bcManage.hpp currentBdFE.hpp currentFE.hpp currentHdivFE.hpp \
dof.hpp dofInterface3Dto2D.hpp dofInterface3Dto3D.hpp dofInterfaceBase.hpp \
dofInterfaceHandler.hpp \
elemOper.hpp elemOper_ext.hpp elemOper2Fluids.hpp \
geoMap.hpp localDofPattern.hpp meshMotion.hpp quadRule.hpp \
refEle.hpp refFE.hpp refHdivFE.hpp refHybridFE.hpp \
subelements.hpp \
sobolevNorms.hpp staticBdFE.hpp regionMesh3D_ALE.hpp v2elemVec.hpp values.hpp \
refEleDG.hpp geoMapDG.hpp refFEDG.hpp \
currentFEDG.hpp currentBdDG.hpp currentBFDG.hpp currentIFDG.hpp \
dofDG.hpp dofByFace.hpp

lib_LTLIBRARIES		= liblifefem.la
liblifefem_la_SOURCES	= assemb.cpp\
 bcCond.cpp bcFunction.cpp bcVector.cpp bcHandler.cpp \
bdf.cpp bdfNS.cpp dataTime.cpp dataTransient.cpp \
currentFE.cpp currentBdFE.cpp currentHdivFE.cpp defQuadRuleFE.cpp dof.cpp \
dofInterface3Dto2D.cpp dofInterface3Dto3D.cpp dofInterfaceBase.cpp dofInterfaceHandler.cpp \
elemOper.cpp  geoMap.cpp \
localDofPattern.cpp meshMotion.cpp quadRule.cpp refEle.cpp refFE.cpp \
refHdivFE.cpp refHybridFE.cpp staticBdFE.cpp regionMesh3D_ALE.cpp v2elemVec.cpp \
refEleDG.cpp geoMapDG.cpp refFEDG.cpp \
currentFEDG.cpp currentBdDG.cpp currentBFDG.cpp currentIFDG.cpp \
dofDG.cpp dofByFace.cpp

liblifefem_la_LIBADD	= ../lifemesh/liblifemesh.la ../lifealg/liblifealg.la 


include $(top_srcdir)/admin/Doxyfile.am
