###################################################################################################
#
#                       This file is part of the LifeV Applications
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University
#
#      Author(s): Name Surname <name.surname@epfl.ch>
#           Date: 00-00-0000
#  License Terms: GNU LGPL
#
###################################################################################################

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc

DO_NOT_COMPILE=

# add newly created folders to SUBDIRS to allow 'make' to enter those folders recursively
SUBDIRS		= data  \
			  lifealg \
			  lifearray \
			  lifecore \
			  lifefem \
			  lifefilters \
			  lifemesh \
			  test_MatrixEpetraStructuredFramework \
			  test_ADRAssembler \
			  test_bdf \
			  test_navierStokes \
              test_darcy \
              test_FEFunction \
              test_fsi \
              test_hyperbolic \
              test_LevelSetSolver \
              test_Parser \
              test_partition \
              test_structure \
              test_TemplateTest \
              test_timeAdvance_orderI \
              test_timeAdvance_orderII \
              test_monolithic \
              test_VectorContainer


# add newly created testsuite.at scripts to TESTSUITES to allow the automatic generation of the
# 'testsuite' script upon modification of one of them
TESTSUITES=$(srcdir)/lifealg/testsuite.at \
           $(srcdir)/lifearray/testsuite.at \
           $(srcdir)/lifecore/testsuite.at \
           $(srcdir)/lifefem/testsuite.at \
           $(srcdir)/lifefilters/testsuite.at \
           $(srcdir)/lifemesh/testsuite.at \
           $(srcdir)/test_MatrixEpetraStructuredFramework \
		   $(srcdir)/test_ADRAssembler/testsuite.at \
		   $(srcdir)/test_bdf/testsuite.at \
		   $(srcdir)/test_navierStokes/testsuite.at \
           $(srcdir)/test_darcy/testsuite.at \
           $(srcdir)/test_fsi/testsuite.at \
           $(srcdir)/test_hyperbolic/testsuite.at \
	   	   $(srcdir)/test_LevelSetSolver/testsuite.at \
           $(srcdir)/test_partition/testsuite.at \
           $(srcdir)/test_structure/testsuite.at \
           $(srcdir)/test_TemplateTest/testsuite.at \
           $(srcdir)/test_timeAdvance_orderI/testsuite.at \
           $(srcdir)/test_timeAdvance_orderII/testsuite.at \
           $(srcdir)/test_monolithic/testsuiteGCE.at \
           $(srcdir)/test_monolithic/testsuiteCE.at \
           $(srcdir)/test_Parser/testsuite.at \
           $(srcdir)/test_VectorContainer/testsuite.at


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

check-local: atconfig atlocal $(TESTSUITE) $(LN_TESTSUITE)
	$(SHELL) $(TESTSUITE)

AUTOM4TE    = autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/testsuite.at $(srcdir)/package.m4 $(TESTSUITES)
	$(AUTOTEST) -I $(srcdir) $@.at -o $@.tmp
	mv $@.tmp $@
	@if test $(top_srcdir) != $(top_builddir); then\
     ln -sf $(TESTSUITE);\
    fi


include $(top_srcdir)/Doxyfile.am

MOSTLYCLEANFILES	= testsuite.log `find -type f -print | grep testsuite.dir`

DISTCLEANFILES		= atconfig



