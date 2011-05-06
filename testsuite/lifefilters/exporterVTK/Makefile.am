# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University      
#
#      Author(s): Umberto Villa <uvilla@emory.edu>
#           Date: 00-00-0000
#  License Terms: GNU LGPL
#
###################################################################################################

include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES                   = .cpp .hpp .idl .c .h .f .F .o .moc

EXTRA_DIST = data testsuite.at

MOSTLYCLEANFILES = *.vct *.scl *.case *.geo *.vtu *.pvtu *.pvd

check_HEADERS              = testExporterVTK.hpp # List of all the hpp files

check_PROGRAMS             = testExport testImport

testExport_SOURCES  = mainExport.cpp # List of all the cpp files

testImport_SOURCES  = mainImport.cpp # List of all the cpp files

link:
	test -a data || ln -s $(srcdir)/data.txt
	test -a Mesh || ln -s $(top_srcdir)/testsuite/data/mesh/inria/ Mesh

recheck: clean_PROGRAMS check

clean_PROGRAMS:
	-rm $(check_PROGRAMS)

clean-results:
	for dustfile in $(MOSTLYCLEANFILES); do \
	!( \test -e $$dustfile ) || \rm $$dustfile; \
	done
