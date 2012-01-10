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

check_HEADERS = navierStokes.hpp

check_PROGRAMS = test_rossEthierSteinman test_kimMoin
#test_kimmoin
test_rossEthierSteinman_SOURCES = mainRossEthierSteinman.cpp ../../life/lifefunctions/RossEthierSteinmanDec.cpp
test_kimMoin_SOURCES = mainKimMoin.cpp ../../life/lifefunctions/KimMoin.cpp

EXTRA_DIST = data testsuite.at

link:
	test -a dataRossEthierSteinman || ln -s $(srcdir)/dataRossEthierSteinman
	test -a dataKimMoin || ln -s $(srcdir)/dataKimMoin
	test -a Mesh || ln -s $(top_srcdir)/testsuite/data/mesh/ Mesh

clean-results:
	-rm *.vct *.scl *.case *.geo

