# -*- Mode : makefile -*-
#
#    AUTHOR: Christophe Prud'homme
#       ORG: EPFL
#    E-MAIL: christophe.prudhomme@epfl.ch
#
# DESCRIPTION:
# ============
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
# DESCRIP-END.
#

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc .tex .pdf


PDF		= timeAdvance.pdf

timeAdvance_TEX	= timeAdvance.tex
all:$(PDF)

.tex.pdf:
	-@if test "x$(PDFLATEX)" = "xpdflatex"; then\
     echo "*** PDFLatex compiling $< ...";\
	 srcdir=$(srcdir) top_srcdir=$(top_srcdir) /bin/sh $(srcdir)/ltxcompile `basename $< .tex`;\
     test $$? -eq 0 || (cat ltxcompile.log | grep -B 5 -n Error); \
    fi

timeAdvance.pdf: $(timeAdvance_TEX)

IMAGES_SRC	= ./figures/tempo.eps ./figures/P11_order1.eps ./figures/P11_order2.eps \
./figures/P11_order3.eps  ./figures/P11_order4.eps ./figures/P21_order1.eps \
./figures/P21_order2.eps ./figures/P21_order3.eps


EXTRA_DIST = $(timeAdvance.tex) ltxcompile

MOSTLYCLEANFILES	=  timeAdvance.dvi timeAdvance.ps timeAdvance.pdf timeAdvance.aux timeAdvance.out

