# -*- makefile -*-
#
# Author: Christophe Prud'homme (christophe.prudhomme@epfl.ch)
#
# Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc 


include ../../Doxyfile.am



clean: clean-am clean-local

clean-local:
	-if test -d life; then rm -rf life; fi
	-if test -d testsuite; then rm -rf testsuite; fi

MOSTLYCLEANFILES    = index.hhc index.hhk files.html Mainpage_8dox.html index.hhp
