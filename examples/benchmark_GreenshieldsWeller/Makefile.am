# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications
#                Copyright (C) 2001-2011 EPFL, Politecnico di Milano, INRIA
#
#      Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
#           Date: 28-10-2010
#  License Terms: GNU GPL
#
###################################################################################################

include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES                          = .cpp .hpp .idl .c .h .f .F .o .moc

check_HEADERS                     =  # List of all the hpp files

check_PROGRAMS                    = benchmark_GreenshieldsWeller

benchmark_GreenshieldsWeller_SOURCES = main.cpp

link:
	test -a data || ln -s $(top_srcdir)/examples/benchmark_GreenshieldsWeller/data data
	test -a fluid_GreenshieldsWeller_20-10-2.mesh || ln -s $(top_srcdir)/examples/benchmark_GreenshieldsWeller/fluid_GreenshieldWeller_20-10-2.mesh; ln -s $(top_srcdir)/examples/benchmark_GreenshieldsWeller/solid_GreenshieldWeller_20-10-2.mesh


recheck: clean_PROGRAMS check

clean_PROGRAMS:
	-rm $(check_PROGRAMS)

clean-results:
	-rm *.vct *.scl *.case *.geo *.xmf *.h5
