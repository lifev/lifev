# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications
#                Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA
#
#      Author(s): Radu Popescu <radu.popescu@epfl.ch>
#           Date: 02-07-2010
#  License Terms: GNU GPL
#
###################################################################################################

include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES                   = .cpp .hpp .idl .c .h .f .F .o .moc

check_HEADERS              = cylinder.hpp

check_PROGRAMS             = test_partition test_cylinder

test_partition_SOURCES = main_partition.cpp
test_cylinder_SOURCES = main_cylinder.cpp cylinder.cpp EqualSolutions.cpp

EXTRA_DIST = data_partition data_cylinder cylinder_ref.h5

link:
	test -a data_partition || ln -s $(srcdir)/data_partition
	test -a data_cylinder || ln -s $(srcdir)/data_cylinder
	test -a run_test.sh || ln -s $(srcdir)/run_test.sh
	test -a cylinder_ref.h5 || ln -s $(srcdir)/cylinder_ref.h5
	test -a Mesh || ln -s $(top_srcdir)/testsuite/data/mesh/inria/ Mesh

recheck: clean_PROGRAMS check

clean_PROGRAMS:
	-rm $(check_PROGRAMS)

clean-results:
	-rm *.h5 *.xmf
