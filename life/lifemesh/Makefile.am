# -*- Mode : makefile -*-
#
# Copyright (C) 2004 EPFL, INRIA and Politechnico di Milano
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

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc .cc

# CP temporary to make it compile
#INCLUDES 	= -I$(top_srcdir)/life/lifecore -I$(top_srcdir)/life/lifearray 
INCLUDES 	= -I$(top_srcdir)/life/lifecore -I$(top_srcdir)/life/lifearray -I$(top_srcdir)/life/lifefem

include_HEADERS		= bareItems.hpp basisElSh.hpp dataMesh.hpp geo0D.hpp \
geoElement.hpp geoND.hpp identifier.hpp markers.hpp mesh_util.hpp mesh_util_base.hpp\
markers_base.hpp meshEntity.hpp readMesh2D.hpp readMesh3D.hpp regionMesh2D.hpp \
regionMesh3D.hpp 
lib_LTLIBRARIES		= liblifemesh.la
liblifemesh_la_SOURCES	= bareItems.cpp basisElSh.cpp geo0D.cpp identifier.cpp \
markers_base.cpp readMesh3D.cpp regionMesh3D.cpp mesh_util.cpp mesh_util_base.cpp

liblifemesh_la_LIBADD	= ../lifecore/liblifecore.la


include $(top_srcdir)/admin/Doxyfile.am