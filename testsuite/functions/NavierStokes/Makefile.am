# -*- makefile -*-
###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
#                Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University      
#
#      Author(s): Name Surname <name.surname@epfl.ch>
#           Date: 00-00-0000
#  License Terms: GNU LGPL
#
###################################################################################################

lib_LTLIBRARIES	= libtestsuite.la

libtestsuite_la_SOURCES   = RossEthierSteinmanDec.cpp RossEthierSteinmanInc.cpp Womersley.cpp \
                             ../bessel/bessik.cpp  ../bessel/bessjy.cpp \
                             ../bessel/cbessik.cpp ../bessel/cbessjy.cpp \
                             ../bessel/gamma.cpp

libtestsuite_la_LIBADD    =  $(top_builddir)/life/lifecore/liblifecore.la

include $(top_srcdir)/Doxyfile.am

