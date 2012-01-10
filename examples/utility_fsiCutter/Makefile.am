# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications
#                Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA
#
#      Author(s): Radu Popescu <radu.popescu@epfl.ch>
#           Date: 24-08-2010
#  License Terms: GNU GPL
#
###################################################################################################

include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES                   = .cpp .hpp .idl .c .h .f .F .o .moc

check_HEADERS              =

check_PROGRAMS             = test_fsi_cutter

test_fsi_cutter_SOURCES = main.cpp

EXTRA_DIST = data

link:
	test -a data || ln -s $(srcdir)/data
	test -a Mesh || ln -s $(top_srcdir)/testsuite/data/mesh/inria/ Mesh

recheck: clean_PROGRAMS check

clean_PROGRAMS:
	-rm $(check_PROGRAMS)

clean-results:
	-rm *.h5 *.xmf
