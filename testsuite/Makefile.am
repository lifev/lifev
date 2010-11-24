# -*- makefile; indent-tabs-mode: t -*-
#
#
#
# Copyright (C) 2004 Christophe Prud'homme (christophe.prudhomme@epfl.ch)
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

SUBDIRS		= data  \
              test_TemplateTest \
              test_cylinder \
              test_restartFluid \
              test_getpot \
              test_structure \
              test_ethiersteinman \
              test_fsi \
              test_ns_flux \
              test_meshReorder \
              test_laplacian \
              test_currentFE \
              test_filters \
              test_partition \
              test_darcy \
              test_resistance \
              test_timeAdvance_orderI \
              test_timeAdvance_orderII \
	      test_ADRAssembler \
	      test_LevelSetSolver \
              $(DO_NOT_COMPILE)

TESTSUITES=$(srcdir)/test_TemplateTest/testsuite.at \
           $(srcdir)/test_cylinder/testsuite.at \
           $(srcdir)/test_structure/testsuite.at \
           $(srcdir)/test_ethiersteinman/testsuite.at\
           $(srcdir)/test_fsi/testsuite.at \
           $(srcdir)/test_filters/testsuite.at \
           $(srcdir)/test_currentFE/testsuite.at \
           $(srcdir)/test_meshReorder/testsuite.at \
           $(srcdir)/test_partition/testsuite.at \
           $(srcdir)/test_darcy/testsuite.at \
           $(srcdir)/test_resistance/testsuite.at \
	   $(srcdir)/test_ADRAssembler/testsuite.at
	   $(srcdir)/test_LevelSetSolver/testsuite.at



$(srcdir)/package.m4: $(top_srcdir)/configure.in.in
	{									   \
	echo '# Signature of the current package.'; \
	echo 'm4_define([AT_PACKAGE_NAME],		[@PACKAGE_NAME@])'; \
	echo 'm4_define([AT_PACKAGE_TARNAME],	[@PACKAGE_TARNAME@])'; \
	echo 'm4_define([AT_PACKAGE_VERSION],	[@PACKAGE_VERSION@])'; \
	echo 'm4_define([AT_PACKAGE_STRING],	[@PACKAGE_STRING@])'; \
	echo 'm4_define([AT_PACKAGE_BUGREPORT], [@PACKAGE_BUGREPORT@])'; \
	} >$(srcdir)/package.m4


EXTRA_DIST = testsuite.at testsuite package.m4

TESTSUITE = $(srcdir)/testsuite
check-local: atconfig atlocal  $(TESTSUITE) $(LN_TESTSUITE)
	$(SHELL) $(TESTSUITE)

AUTOM4TE    = autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/testsuite.at $(srcdir)/package.m4 $(TESTSUITES)
	$(AUTOTEST) -I $(srcdir) $@.at -o $@.tmp
	mv $@.tmp $@
	@if test $(top_srcdir) != $(top_builddir); then\
     ln -sf $(TESTSUITE);\
    fi

testsuite: $(srcdir)/testsuite.at $(srcdir)/package.m4 $(TESTSUITES)
	$(AUTOTEST) -I $(srcdir) $(srcdir)/$@.at -o $(srcdir)/$@.tmp
	mv $(srcdir)/$@.tmp $(srcdir)/$@
	@if test $(top_srcdir) != $(top_builddir); then\
     ln -sf $(TESTSUITE);\
    fi


include $(top_srcdir)/Doxyfile.am

MOSTLYCLEANFILES	= testsuite.log `find -type f -print | grep testsuite.dir`

DISTCLEANFILES		= atconfig



