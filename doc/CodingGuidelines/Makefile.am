# -*- Mode : makefile -*-
#
#    AUTHOR: Cristiano Malossi
#       ORG: EPFL
#    E-MAIL: cristiano.malossi@epfl.ch
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

PDF		= LifeV_Development_Guidelines.pdf

LifeV_Development_Guidelines_TEX	= LifeV_Development_Guidelines.tex

all:$(PDF)

.tex.pdf: 
	-@if test "x$(PDFLATEX)" = "xpdflatex"; then\
     echo "*** PDFLatex compiling $< ...";\
	 srcdir=$(srcdir) top_srcdir=$(top_srcdir) /bin/sh $(srcdir)/ltxcompile `basename $< .tex`;\
     test $$? -eq 0 || (cat ltxcompile.log | grep -B 5 -n Error); \
    fi

LifeV_Development_Guidelines.pdf: $(LifeV_Development_Guidelines_TEX)
 
IMAGES_SRC	= 

EXTRA_DIST = $(LifeV_Development_Guidelines_TEX) $(IMAGES_SRC) ltxcompile

MOSTLYCLEANFILES	= LifeV_Development_Guidelines.log LifeV_Development_Guidelines.glo LifeV_Development_Guidelines.idx LifeV_Development_Guidelines.ind LifeV_Development_Guidelines.ilg LifeV_Development_Guidelines.aux LifeV_Development_Guidelines.out \
                      LifeV_Development_Guidelines.pdf LifeV_Development_Guidelines.toc LifeV_Development_Guidelines.lof LifeV_Development_Guidelines.lot LifeV_Development_Guidelines.blg LifeV_Development_Guidelines.bbl ltxcompile.log

