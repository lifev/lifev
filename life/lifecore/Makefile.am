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

include_HEADERS		= GetPot.hpp  debug.hpp\
chrono.hpp dataString.hpp dataTime.hpp fortran_wrap.hpp \
lifeV.hpp post_proc.hpp subDomainNeighbors.hpp switch.hpp \
util_string.hpp bareItems.hpp basisElSh.hpp dataTransient.hpp \
singleton.hpp policy.hpp typeInfo.hpp

lib_LTLIBRARIES		= liblifecore.la

liblifecore_la_SOURCES	= bareItems.cpp basisElSh.cpp GetPot.cpp  dataString.cpp debug.cpp \
dataTime.cpp subDomainNeighbors.cpp switch.cpp util_string.cpp dataTransient.cpp typeInfo.cpp


liblifecore_la_LIBADD	= $(aztec_libs)

include $(top_srcdir)/admin/Doxyfile.am
