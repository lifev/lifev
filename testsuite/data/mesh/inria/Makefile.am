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

dist_pkgdata_DATA = cube_384.mesh cube_6007.mesh cylhexa.mesh cyltetra.mesh\
		    piece-of-tube-wall-tetra.mesh hexa16x16.mesh hexa2x2.mesh\
		    hexa32x32.mesh hexa2x2.mesh hexa4x4.mesh hexa8x8.mesh tube2x4x20.mesh\
		    tube_curve20.mesh vessel2x4x20.mesh vessel2x4x20_10cm.mesh\
		    vessel_curve20.mesh hexahexa10x10x10.mesh\
                    tube20.mesh vessel20.mesh \
            piece-of-tube.mesh piece-of-tube-wall.mesh hexa1x1.mesh tube.mesh vessel.mesh \
            heart13628.mesh heart1959.mesh

check-local:
	@if test "$(top_srcdir)" != "$(top_builddir)"; then\
	 for i in $(dist_pkgdata_DATA); do \
           $(LN_S) -f $(srcdir)/$$i .; \
         done; \
	fi

