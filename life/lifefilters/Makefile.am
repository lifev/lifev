# -*- makefile -*-
#
# Author: Christophe Prud'homme
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

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc 

INCLUDES	= $(all_includes)

include_HEADERS	= gmsh_wrtrs.hpp medit_wrtrs.hpp openDX_wrtrs.hpp vtk_wrtrs.hpp\
		  gmv_wrtrs.hpp mesh_util.hpp readMesh2D.hpp readMesh3D.hpp ensight7Writer.hpp\
          importer.hpp

lib_LTLIBRARIES	= liblifefilters.la

liblifefilters_la_SOURCES	= medit_wrtrs.cpp vtk_wrtrs.cpp readMesh3D.cpp\
				  openDX_wrtrs.cpp importer.cpp
liblifefilters_la_LIBADD	= ../lifefem/liblifefem.la

include $(top_srcdir)/admin/Doxyfile.am
