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

TEX_SRC	= docfe.tex jfg.tex Specifications.tex
NOTEX_SRC	= docfe jfg Specifications 
PDF		= docfe.pdf jfg.pdf Specifications.pdf

all: $(PDF)


.tex.pdf: $(NOTEX_SRC)
	-(TEXINPUTS=$$TEXINPUTS:$(srcdir):$(top_srcdir)/doc/tex//  pdflatex `basename $< .tex`)
	-(BIBINPUTS=$$BIBINPUTS:$(srcdir):$(top_srcdir)/doc/tex//  bibtex `basename $< .tex`)
	-(TEXINPUTS=$$TEXINPUTS:$(srcdir):$(top_srcdir)/doc/tex//  pdflatex `basename $< .tex`)
	-(TEXINPUTS=$$TEXINPUTS:$(srcdir):$(top_srcdir)/doc/tex//  pdflatex `basename $< .tex`)


IMAGES_SRC	= BCCond.fig  BCHandler.fig  BasElSha.fig  Fem_nodes.fig  mixed_pattern.fig  theFE.fig \
BCCond.eps BasElSha.eps NSCode_mesh_types.eps theFE.eps BCHandler.eps Fem_nodes.eps mixed_pattern.eps \
BCCond.pdf BasElSha.pdf NSCode_mesh_types.pdf theFE.pdf BCHandler.pdf Fem_nodes.pdf mixed_pattern.pdf


EXTRA_DIST = $(TEX_SRC) $(IMAGES_SRC) life.bib 

MOSTLYCLEANFILES	= docfe.log docfe.aux docfe.pdf docfe.blg docfe.bbl \
			  jfg.log jfg.aux jfg.pdf jfg.blg jfg.bbl \
			  Specifications.log Specifications.glo Specifications.idx Specifications.aux Specifications.out \
                          Specifications.pdf Specifications.toc Specifications.lof Specifications.lot Specifications.blg Specifications.bbl 