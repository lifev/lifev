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

AUTOMAKE_OPTIONS = foreign 1.6 dejagnu

$(top_srcdir)/acinclude.m4: $(top_srcdir)/admin/acinclude.m4.in $(top_srcdir)/admin/libtool.m4.in $(top_srcdir)/libltdl/ltdl.m4
	@cd $(top_srcdir) && cat admin/acinclude.m4.in admin/libtool.m4.in libltdl/ltdl.m4 > acinclude.m4

$(top_srcdir)/configure.in: $(top_srcdir)/subdirs $(top_srcdir)/configure.files $(CONF_FILES)
	cd $(top_srcdir) && $(MAKE) -f admin/Makefile.common configure.in


EXTRA_DIST = configure.in.in configure.in.bot  configure.in.mid subdirs configure.files 

MAINTAINERCLEANFILES = configure.files subdirs


SUBDIRS = admin doc life testsuite

include $(top_srcdir)/admin/Doxyfile.am

apidox-all:
    make apidox 
	cd life && make apidox


upload-api:
	cd doc/ && scp -r api prudhomm@cmcsforge.epfl.ch:/var/lib/gforge/chroot/home/groups/lifev/htdocs/api-$(VERSION)
