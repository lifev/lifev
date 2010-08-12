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

include $(top_srcdir)/testsuite/Makefile.testsuite

check_HEADERS = ensightToHdf5.hpp

check_PROGRAMS = test_filters
test_filters_SOURCES = main.cpp ensightToHdf5.cpp

link:
	test -a Mesh              || ln -s $(top_srcdir)/testsuite/data/mesh/inria/ Mesh

EXTRA_DIST = data testsuite.at \
  ethiersteinman.000.case                ethiersteinman_pressure.00004.000.scl\
  ethiersteinman.000.geo                 ethiersteinman_pressure.00004.001.scl\
  ethiersteinman.001.case                ethiersteinman_pressure.00004.002.scl\
  ethiersteinman.001.geo                 ethiersteinman_pressure.00004.003.scl\
  ethiersteinman.002.case                ethiersteinman_velocity.00000.000.vct\
  ethiersteinman.002.geo                 ethiersteinman_velocity.00000.001.vct\
  ethiersteinman.003.case                ethiersteinman_velocity.00000.002.vct\
  ethiersteinman.003.geo                 ethiersteinman_velocity.00000.003.vct\
  ethiersteinman_pressure.00000.000.scl  ethiersteinman_velocity.00001.000.vct\
  ethiersteinman_pressure.00000.001.scl  ethiersteinman_velocity.00001.001.vct\
  ethiersteinman_pressure.00000.002.scl  ethiersteinman_velocity.00001.002.vct\
  ethiersteinman_pressure.00000.003.scl  ethiersteinman_velocity.00001.003.vct\
  ethiersteinman_pressure.00001.000.scl  ethiersteinman_velocity.00002.000.vct\
  ethiersteinman_pressure.00001.001.scl  ethiersteinman_velocity.00002.001.vct\
  ethiersteinman_pressure.00001.002.scl  ethiersteinman_velocity.00002.002.vct\
  ethiersteinman_pressure.00001.003.scl  ethiersteinman_velocity.00002.003.vct\
  ethiersteinman_pressure.00002.000.scl  ethiersteinman_velocity.00003.000.vct\
  ethiersteinman_pressure.00002.001.scl  ethiersteinman_velocity.00003.001.vct\
  ethiersteinman_pressure.00002.002.scl  ethiersteinman_velocity.00003.002.vct\
  ethiersteinman_pressure.00002.003.scl  ethiersteinman_velocity.00003.003.vct\
  ethiersteinman_pressure.00003.000.scl  ethiersteinman_velocity.00004.000.vct\
  ethiersteinman_pressure.00003.001.scl  ethiersteinman_velocity.00004.001.vct\
  ethiersteinman_pressure.00003.002.scl  ethiersteinman_velocity.00004.002.vct\
  ethiersteinman_pressure.00003.003.scl  ethiersteinman_velocity.00004.003.vct