# -*- makefile; indent-tabs-mode: t -*-
#
#
#
# Copyright (C) 2010 Simone Deparis (simone.deparis@epfl.ch)
#
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
#

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc

DO_NOT_COMPILE=

SUBDIRS		=  \
              $(DO_NOT_COMPILE)



$(srcdir)/package.m4: $(top_srcdir)/configure.in.in
	{									   \
	echo '# Signature of the current package.'; \
	echo 'm4_define([AT_PACKAGE_NAME],		[@PACKAGE_NAME@])'; \
	echo 'm4_define([AT_PACKAGE_TARNAME],	[@PACKAGE_TARNAME@])'; \
	echo 'm4_define([AT_PACKAGE_VERSION],	[@PACKAGE_VERSION@])'; \
	echo 'm4_define([AT_PACKAGE_STRING],	[@PACKAGE_STRING@])'; \
	echo 'm4_define([AT_PACKAGE_BUGREPORT], [@PACKAGE_BUGREPORT@])'; \
	} >$(srcdir)/package.m4


EXTRA_DIST =

EXAMPLES = $(srcdir)/examples
check-local:


include $(top_srcdir)/Doxyfile.am

MOSTLYCLEANFILES	=

DISTCLEANFILES		=



