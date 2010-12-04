# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications
#                Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA
#
#      Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
#           Date: 02-07-2010
#  License Terms: GNU GPL
#
###################################################################################################

SUFFIXES    = .cpp .hpp .idl .c .h .f .F .o .moc .tex .pdf

SUBDIRS = images

PDF         = LifeV_SOTA.pdf

LifeV_SOTA_TEX = LifeV_SOTA.tex

all:$(PDF)

.tex.pdf:
	-@if test "x$(PDFLATEX)" = "xpdflatex"; then\
     echo "*** PDFLatex compiling $< ...";\
	 srcdir=$(srcdir) top_srcdir=$(top_srcdir) /bin/sh $(srcdir)/ltxcompile `basename $< .tex`;\
     test $$? -eq 0 || (cat ltxcompile.log | grep -B 5 -n Error); \
    fi

LifeV_SOTA.pdf: $(LifeV_SOTA_TEX)

IMAGES_SRC	=

EXTRA_DIST = $(LifeV_SOTA_TEX) $(IMAGES_SRC) ltxcompile

MOSTLYCLEANFILES = LifeV_SOTA.log LifeV_SOTA.aux
