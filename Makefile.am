# -*- Mode : makefile; indent-tabs-mode: t -*-
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

AUTOMAKE_OPTIONS = foreign 1.7 dist-bzip2

$(top_srcdir)/configure.in: $(top_srcdir)/subdirs $(top_srcdir)/configure.files $(CONF_FILES)
	cd $(top_srcdir) && $(MAKE) -f admin/Makefile.common configure.in


EXTRA_DIST = configure.in.in configure.in.bot  configure.in.mid subdirs configure.files Mainpage.dox BUGS

MAINTAINERCLEANFILES = configure.files subdirs


SUBDIRS = admin doc life testsuite

include $(top_srcdir)/admin/Doxyfile.am

dist-hook: changelog
# make sure that we have an up-to-date ChangeLog
	if test -w . -a $(top_srcdir) != $(top_builddir) -a ! -f ChangeLog; then\
	  cp $(srcdir)/ChangeLog .;\
	fi

changelog:
	if test -x /usr/bin/cvs2cl -a -d $(srcdir)/CVS; then\
	 (cd $(srcdir) && cvs2cl -b -t  -S --utc -P);\
    else \
     if test -w $(srcdir)/ChangeLog; then\
       touch $(srcdir)/ChangeLog;\
     fi;\
    fi

apidox-all:
	make apidox
	cd life && make apidox

all-local: 
	make apidox-all
	@echo ""
	@echo "************************************************************"
	@echo "* LifeV compilation has ended."
	@echo "*"
	@echo "* LifeV API documentation has been generated in doc/api/"
	@echo "* Open your browser at doc/api/index.html."
	@echo "*"
	@echo "* Enjoy!"
	@echo "*  C.P."
	@echo "************************************************************"

upload-api:
	-cd doc/ && rsync -avz api/ `whoami`@cmcsforge.epfl.ch:/var/lib/gforge/chroot/home/groups/lifev/htdocs/api-$(LIFEV_MAJOR_VERSION).$(LIFEV_MINOR_VERSION).$(LIFEV_MICRO_VERSION)

MOSTLYCLEANFILES	= ChangeLog
