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
                  test_mesh test_matrix test_naturalbc test_essentialbc \
                  test_q1 test_p2 test_mixte test_fe \
                  test_postproc test_bdf \
		  test_linearelasticity \
		  test_darcy \
		  test_ns_bdf test_ns_cyl test_ns_sstress \
		  test_fsi_newton test_fsi_picard

check-qmtest: 
	qmtest -D $(srcdir) run
	qmtest -D $(srcdir) summarize

