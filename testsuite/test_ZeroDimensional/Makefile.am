###################################################################################################
#
#                       This file is part of the LifeV Library                       
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University   
#
#      Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
#           Date: 14-01-2010
#  License Terms: GNU LGPL
#
###################################################################################################

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc

include $(top_srcdir)/testsuite/Makefile.testsuite

check_HEADERS  =
check_PROGRAMS = test_ZeroDimensional

test_ZeroDimensional_SOURCES = main.cpp

EXTRA_DIST = testsuite.at

link:
	test -a RythmosAztecOOParams.xml || ln -s $(top_srcdir)/testsuite/test_ZeroDimensional/RythmosAztecOOParams.xml RythmosAztecOOParams.xml
	test -a RythmosAztecOOParamsLowsf.xml || ln -s $(top_srcdir)/testsuite/test_ZeroDimensional/RythmosAztecOOParamsLowsf.xml RythmosAztecOOParamsLowsf.xml
	test -a RythmosBelosParams.xml || ln -s $(top_srcdir)/testsuite/test_ZeroDimensional/RythmosBelosParams.xml RythmosBelosParams.xml
	test -a data.dat || ln -s $(top_srcdir)/testsuite/test_ZeroDimensional/data.dat data.dat
	test -a circuitFile.dat || ln -s $(top_srcdir)/testsuite/test_ZeroDimensional/circuitFile.dat circuitFile.dat
	test -a inputFile.dat || ln -s $(top_srcdir)/testsuite/test_ZeroDimensional/inputFile.dat inputFile.dat

recheck: clean_PROGRAMS check

clean_PROGRAMS:
	-rm $(check_PROGRAMS)

clean-results:
	-rm *.vct *.scl *.case *.geo *.mfile

include $(top_srcdir)/Doxyfile.am