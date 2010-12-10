# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Library                       
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University
#
#      Author(s): Umberto Villa <uvilla@emory.edu>
#           Date: 10-12-2010
#  License Terms: GNU LGPL
#
#################################################################################################

AM_CXXFLAGS=-O1

SUFFIXES		= .cpp .hpp .idl .c .h .f .F .o .moc 

INCLUDES        = -I$(top_srcdir)/life/lifecore $(top_srcdir)/life/lifealg

LDADD           = $(top_builddir)/life/lifecore/liblifecore.la $(top_builddir)/life/lifealg/liblifealg.la

check_PROGRAMS 	    = #add list of programs

EXTRA_DIST = 

link:
	test -a data || ln -s $(srcdir)/data