# -*- makefile -*-
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

include $(top_srcdir)/testsuite/Makefile.testsuite

check_HEADERS = ensightToHdf5.hpp

check_PROGRAMS = test_filters
test_filters_SOURCES = main.cpp ensightToHdf5.cpp

link:
	test -a Mesh || ln -s $(top_srcdir)/testsuite/data/mesh/inria/ Mesh
	test -a importDir || ln -s $(top_srcdir)/testsuite/lifefilters/ importDir

EXTRA_DIST = data testsuite.at \
            ethiersteinman.000.case \
            ethiersteinman.000.geo \
            ethiersteinman.001.case \
            ethiersteinman.001.geo \
            ethiersteinman.002.case \
            ethiersteinman.002.geo \
            ethiersteinman.003.case \
            ethiersteinman.003.geo \
            ethiersteinman_pressure.00000.000.scl \
            ethiersteinman_pressure.00000.001.scl \
            ethiersteinman_pressure.00000.002.scl \
            ethiersteinman_pressure.00000.003.scl \
            ethiersteinman_pressure.00001.000.scl \
            ethiersteinman_pressure.00001.001.scl \
            ethiersteinman_pressure.00001.002.scl \
            ethiersteinman_pressure.00001.003.scl \
            ethiersteinman_pressure.00002.000.scl \
            ethiersteinman_pressure.00002.001.scl \
            ethiersteinman_pressure.00002.002.scl \
            ethiersteinman_pressure.00002.003.scl \
            ethiersteinman_pressure.00003.000.scl \
            ethiersteinman_pressure.00003.001.scl \
            ethiersteinman_pressure.00003.002.scl \
            ethiersteinman_pressure.00003.003.scl \
            ethiersteinman_pressure.00004.000.scl \
            ethiersteinman_pressure.00004.001.scl \
            ethiersteinman_pressure.00004.002.scl \
            ethiersteinman_pressure.00004.003.scl \
            ethiersteinman_velocity.00000.000.vct \
            ethiersteinman_velocity.00000.001.vct \
            ethiersteinman_velocity.00000.002.vct \
            ethiersteinman_velocity.00000.003.vct \
            ethiersteinman_velocity.00001.000.vct \
            ethiersteinman_velocity.00001.001.vct \
            ethiersteinman_velocity.00001.002.vct \
            ethiersteinman_velocity.00001.003.vct \
            ethiersteinman_velocity.00002.000.vct \
            ethiersteinman_velocity.00002.001.vct \
            ethiersteinman_velocity.00002.002.vct \
            ethiersteinman_velocity.00002.003.vct \
            ethiersteinman_velocity.00003.000.vct \
            ethiersteinman_velocity.00003.001.vct \
            ethiersteinman_velocity.00003.002.vct \
            ethiersteinman_velocity.00003.003.vct \
            ethiersteinman_velocity.00004.000.vct \
            ethiersteinman_velocity.00004.001.vct \
            ethiersteinman_velocity.00004.002.vct \
            ethiersteinman_velocity.00004.003.vct