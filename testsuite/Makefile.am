# -*- makefile -*-
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


SUBDIRS		= data  \
              lifecore lifearray lifealg \
              test_mesh test_laplacian test_naturalbc test_essentialbc \
              test_q1 test_p2 test_mixte test_fe \
              test_postproc test_bdf \
		      test_linearelasticity \
		      test_darcy test_darcytetra test_onedmodel\
		      test_ns_bdf test_ns_cyl test_ns_sstress \
		      test_fsi_newton test_fsi_picard \
		      test_masstransport test_mass_coupling test_ns_ipstab

TESTSUITES=$(srcdir)/lifecore/testsuite.at \
$(srcdir)/lifearray/testsuite.at \
$(srcdir)/lifealg/testsuite.at \
$(srcdir)/test_bdf/testsuite.at \
$(srcdir)/test_bdf/testsuite.at \
$(srcdir)/test_darcy/testsuite.at \
$(srcdir)/test_darcytetra/testsuite.at \
$(srcdir)/test_essentialbc/testsuite.at \
$(srcdir)/test_fe/testsuite.at \
$(srcdir)/test_fsi_newton/testsuite.at \
$(srcdir)/test_fsi_picard/testsuite.at \
$(srcdir)/test_linearelasticity/testsuite.at \
$(srcdir)/test_masstransport/testsuite.at \
$(srcdir)/test_mass_coupling/testsuite.at \
$(srcdir)/test_laplacian/testsuite.at \
$(srcdir)/test_mesh/testsuite.at \
$(srcdir)/test_mixte/testsuite.at \
$(srcdir)/test_naturalbc/testsuite.at \
$(srcdir)/test_ns_bdf/testsuite.at \
$(srcdir)/test_ns_cyl/testsuite.at \
$(srcdir)/test_ns_ipstab/testsuite.at \
$(srcdir)/test_ns_sstress/testsuite.at \
$(srcdir)/test_onedmodel/testsuite.at \
$(srcdir)/test_p2/testsuite.at \
$(srcdir)/test_postproc/testsuite.at \
$(srcdir)/test_q1/testsuite.at 


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
check-local: atconfig atlocal  $(TESTSUITE)
	$(SHELL) $(TESTSUITE)

AUTOM4TE    = autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/testsuite.at $(srcdir)/package.m4 $(TESTSUITES)
	$(AUTOTEST) -I $(srcdir) $@.at -o $@.tmp
	mv $@.tmp $@



include $(top_srcdir)/admin/Doxyfile.am

MOSTLYCLEANFILES	= atconfig testsuite.log `find -type f -print | grep testsuite.dir`





