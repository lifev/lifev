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

check_PROGRAMS                    = benchmark_HPCNavierStokes
benchmark_HPCNavierStokes_SOURCES = benchmark_HPCNavierStokes.cpp

link:
	test -a Ressources || ln -s $(top_srcdir)/examples/benchmark_HPCNavierStokes/Ressources Ressources
	test -a data || ln -s $(top_srcdir)/examples/benchmark_HPCNavierStokes/Ressources/data data
	test -a 100170_cm_N2H0.008_fluid_coarse.mesh || ln -s $(top_srcdir)/examples/benchmark_HPCNavierStokes/Ressources/100170_cm_N2H0.008_fluid_coarse.mesh 100170_cm_N2H0.008_fluid_coarse.mesh

scripts:
	cp $(top_srcdir)/examples/benchmark_HPCNavierStokes/Scripts/* ./

recheck: clean_PROGRAMS check

clean_PROGRAMS:
	-rm $(check_PROGRAMS)

clean-results:
	-rm *.vct *.scl *.case *.geo *.xmf *.h5
