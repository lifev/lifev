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

IMAGES_SRC	= images/elasticityModels/cylindredomain.pdf images/fe/p03D.pdf images/fe/p13D.pdf images/fe/p23D.pdf images/fe/p1Bubble3D.pdf images/geometriesAndMeshes/aortaIRM.png  images/geometriesAndMeshes/arterz.png  images/geometriesAndMeshes/bypass.png  images/geometriesAndMeshes/carotid.png  images/geometriesAndMeshes/crosettoAorta.png  images/geometriesAndMeshes/heartfront.png  images/geometriesAndMeshes/hearttop.png  images/geometriesAndMeshes/postopAFB.png  images/geometriesAndMeshes/preopAFB1.png  images/geometriesAndMeshes/preopAFB2.png  images/geometriesAndMeshes/stent.png  images/geometriesAndMeshes/stentExpanded.png  images/geometriesAndMeshes/virtualAorta.png images/levelset/levelsetTesta.pdf images/levelset/levelsetTestb.pdf images/levelset/levelsetTestc.pdf images/monoAndBidomain/truncatedellipsoidaldomain.pdf images/monolitic/cylindredomain.pdf images/multiscale/application1.pdf images/multiscale/application2.pdf images/multiscale/application3a.pdf images/multiscale/application3b.pdf images/multiscale/bca.pdf images/multiscale/bcb.pdf images/multiscale/explicitMethoda.pdf images/multiscale/explicitMethodb.pdf images/multiscale/explicitMethodc.pdf 

EXTRA_DIST = $(LifeV_SOTA_TEX) $(IMAGES_SRC) ltxcompile mcode.sty

MOSTLYCLEANFILES = LifeV_SOTA.log LifeV_SOTA.aux
