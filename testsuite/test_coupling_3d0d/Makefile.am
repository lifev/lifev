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

AM_CPPFLAGS	= -DNO_BLASLAPACK

include $(top_srcdir)/testsuite/Makefile.testsuite

check_HEADERS  = ud_functions.hpp zeroDModelSolver.hpp meanPressure.hpp flowRate.hpp C3d0d.hpp

check_PROGRAMS  = test_coupling_NS0d test_coupling_FSI0d

test_coupling_NS0d_SOURCES = main_NS0d.cpp

test_coupling_FSI0d_SOURCES = main_FSI0d.cpp

EXTRA_DIST = data_NS0d data_FSI0d testsuite.at

MOSTLYCLEANFILES=fluid.res1 test.001.inp press.001.bb \
                 vel_x.001.bb vel_y.001.bb vel_z.001.bb fluid.res2 test.002.inp\
                 press.002.bb vel_x.002.bb vel_y.002.bb vel_z.002.bb \
                 res_DP.m res_Q.m \
                 \
                 vel_x.001.bb vel_y.001.bb vel_z.001.bb \
                 vel_x.002.bb vel_y.002.bb vel_z.002.bb \
                 vel_x.001.mesh vel_y.001.mesh vel_z.001.mesh \
                 vel_x.002.mesh vel_y.002.mesh vel_z.002.mesh \
                 dep_x.001.bb dep_y.001.bb dep_z.001.bb \
                 dep_x.002.bb dep_y.002.bb dep_z.002.bb \
                 dep_x.001.mesh dep_y.001.mesh dep_z.001.mesh \
                 dep_x.002.mesh dep_y.002.mesh dep_z.002.mesh \
                 press.001.bb press.002.bb press.001.mesh press.002.mesh\
                 defor.001.mesh defor.002.mesh \
                 iter num_iter out_iter_solid out_res_solid res \
                 res_DP.m res_Qin.m res_Qout.m
