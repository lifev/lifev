# -*- makefile -*-
#
# Author: Christophe Prud'homme <christophe.prudhomme@epfl.ch>
#
# Copyright (C) 2004 EPFL
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

SUFFIXES		= .cpp .hpp .idl .c .h .f .F .o .moc 

INCLUDES        = -I$(top_srcdir)/life/lifecore

LDADD           = $(top_builddir)/life/lifecore/liblifecore.la

check_PROGRAMS 	= test_casts test_copy test_mark_to_bc test_getpot\
                  test_pairs test_pf test_prova test_rotate \
                  test_seq test_size test_typedefs test_using \
                  test_valarray test_vsize test_debug



test_casts_SOURCES	= test_casts.cpp
test_copy_SOURCES	= test_copy.cpp
test_mark_to_bc_SOURCES	= test_mark_to_bc.cpp
test_pairs_SOURCES	= test_pairs.cpp
test_pf_SOURCES		= test_pf.cpp
test_prova_SOURCES	= test_prova.cpp
test_rotate_SOURCES	= test_rotate.cpp
test_seq_SOURCES	= test_seq.cpp
test_size_SOURCES	= test_size.cpp
test_typedefs_SOURCES	= test_typedefs.cpp
test_using_SOURCES	= test_using.cpp
test_valarray_SOURCES	= test_valarray.cpp
test_vsize_SOURCES	= test_vsize.cpp

test_debug_SOURCES  = test_debug.cpp
test_getpot_SOURCES	= test_getpot.cpp

POTFILES			= mesh.pot mesh2.pot mesh3.pot

EXTRA_DIST			= $(POTFILES)

check-local: $(POTFILES)
	@for potfile in $(POTFILES); do\
     ln -sf $(srcdir)/$$potfile ;  \
    done
