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

INCLUDES 	= -I$(top_srcdir)/life/lifecore -I$(top_srcdir)/life/lifemesh
AM_CPPFLAGS	= -DTHREEDIM

include_HEADERS		= currentBdFE.hpp currentFE.hpp currentHdivFE.hpp elemOper.hpp \
elemOper_ext.hpp geoMap.hpp localDofPattern.hpp quadRule.hpp \
refEle.hpp refFE.hpp refHdivFE.hpp refHybridFE.hpp \
sobolevNorms.hpp staticBdFE.hpp 

lib_LTLIBRARIES		= liblifefem.la
liblifefem_la_SOURCES	= currentBdFE.cpp currentFE.cpp currentHdivFE.cpp defQuadRuleFE.cpp \
elemOper.cpp geoMap.cpp localDofPattern.cpp quadRule.cpp \
refEle.cpp refFE.cpp refHdivFE.cpp refHybridFE.cpp \
staticBdFE.cpp 


liblifefem_la_LIBADD	= ../lifemesh/liblifemesh.la $(aztec_libs)
