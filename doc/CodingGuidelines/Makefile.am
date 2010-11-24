# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Library                       
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University   
#
#      Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
#           Date: 05-01-2010
#  License Terms: GNU LGPL
#
###################################################################################################

SUFFIXES    = .cpp .hpp .idl .c .h .f .F .o .moc .tex .pdf

PDF         = LifeV_Development_Guidelines.pdf

LifeV_Development_Guidelines_TEX = LifeV_Development_Guidelines.tex

all:$(PDF)

.tex.pdf:
	-@if test "x$(PDFLATEX)" = "xpdflatex"; then\
     echo "*** PDFLatex compiling $< ...";\
	 srcdir=$(srcdir) top_srcdir=$(top_srcdir) /bin/sh $(srcdir)/ltxcompile `basename $< .tex`;\
     test $$? -eq 0 || (cat ltxcompile.log | grep -B 5 -n Error); \
    fi

LifeV_Development_Guidelines.pdf: $(LifeV_Development_Guidelines_TEX)

IMAGES_SRC	=

LISTING_SRC = ExampleClass.hpp ExampleClass.cpp main.cpp LifeV_Makefile.am testsuite.at

EXTRA_DIST = $(LifeV_Development_Guidelines_TEX) $(IMAGES_SRC) $(LISTING_SRC) ltxcompile mcode.sty

MOSTLYCLEANFILES = LifeV_Development_Guidelines.log LifeV_Development_Guidelines.glo LifeV_Development_Guidelines.idx LifeV_Development_Guidelines.ind LifeV_Development_Guidelines.ilg LifeV_Development_Guidelines.aux LifeV_Development_Guidelines.out \
                   LifeV_Development_Guidelines.pdf LifeV_Development_Guidelines.toc LifeV_Development_Guidelines.lof LifeV_Development_Guidelines.lot LifeV_Development_Guidelines.blg LifeV_Development_Guidelines.bbl ltxcompile.log
