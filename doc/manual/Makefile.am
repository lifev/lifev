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

SUBDIRS = cavity

PDF		= lifev-dev.pdf lifev-manual.pdf

lifev_dev_TEX	= lifev-dev.tex lifev-install.tex lifev-examples.tex lifev-dev_conventions.tex \
     lifev-dev_overview.tex lifev-dev_howto.tex

lifev_manual_TEX	= lifev-manual.tex lifev-install.tex lifev-examples.tex cavity.tex \
     navierstokes.tex fsi.tex

all:$(PDF)

.tex.pdf:
	-@if test "x$(PDFLATEX)" = "xpdflatex"; then\
     echo "*** PDFLatex compiling $< ...";\
	 srcdir=$(srcdir) top_srcdir=$(top_srcdir) /bin/sh $(srcdir)/ltxcompile `basename $< .tex`;\
     test $$? -eq 0 || (cat ltxcompile.log | grep -B 5 -n Error); \
    fi

lifev-dev.pdf: $(lifev_dev_TEX)

lifev-manual.pdf: $(lifev_manual_TEX)

IMAGES_SRC	= BCCond.fig  BCHandler.fig  BasElSha.fig  Fem_nodes.fig  mixed_pattern.fig  theFE.fig \
BCCond.eps BasElSha.eps NSCode_mesh_types.eps theFE.eps BCHandler.eps Fem_nodes.eps mixed_pattern.eps \
BCCond.pdf BasElSha.pdf NSCode_mesh_types.pdf theFE.pdf BCHandler.pdf Fem_nodes.pdf mixed_pattern.pdf \
cavityFigure.pdf_t cavityFigure.pdf ALEmapping.pdf


EXTRA_DIST = $(lifev_dev_TEX) $(lifev_manual_TEX) $(IMAGES_SRC) ltxcompile

MOSTLYCLEANFILES	= lifev-dev.log lifev-dev.glo lifev-dev.idx lifev-dev.ind lifev-dev.ilg lifev-dev.aux lifev-dev.out \
                          lifev-dev.pdf lifev-dev.toc lifev-dev.lof lifev-dev.lot lifev-dev.blg lifev-dev.bbl \
			  lifev-dev_generalities.aux lifev-dev_conventions.aux lifev-dev_overview.aux lifev-dev_howto.aux ltxcompile.log

