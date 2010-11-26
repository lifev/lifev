# -*- makefile -*-
#
# Author: Gilles Fourestey
# Modified: Gwenol Grandperrin
#
# Copyright (C) 2009 @COPYRIGHT
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#


# path to the compiler
CC              = /usr/bin/g++

SOURCES         = cavity.cpp
OBJECTS         = $(SOURCES:.cpp=.o)
EXECUTABLE      = cavity_example


LIFEPATH         = /u/cmcs/jtbecerr/lib/LIFEV
LIFELIBPATH      = -L$(LIFEPATH)/lib
LIFELIBS         = -llifefilters -llifesolver -llifefem \
                   -llifealg -llifearray -llifecore -llifemesh
LIFEINCLUDEPATH  = -I$(LIFEPATH)/include

TRILLIBPATH      = -L/u/cmcs/jtbecerr/lib/TRILINOS/build/lib
TRILLIBS         = -laztecoo -laztecoo -ltriutils -lml \
                   -lifpack -lamesos -lepetraext -lepetra \
                   -lteuchos -lumfpack -lamd -llapack -lblas
TRILINCLUDEPATH  = -I/u/cmcs/jtbecerr/lib/TRILINOS/build/include

MPILIBPATH       = -L/usr/lib/openmpi/lib
MPILIBS          = -lmpi_cxx -lmpi -lopen-rte
MPIINCLUDEPATH   = -I/usr/include/openmpi

METISLIBPATH	 = -L/usr/lib
METISLIBS        = -lparmetis -lmetis
METISINCLUDEPATH = -I/usr/include


# uncomment this part for optimized compilation
CXXFLAGS         = -DHAVE_CONFIG_H -DTHREEDIM  -O2
# uncomment this part for debugging
#CXXFLAGS         = -DHAVE_CONFIG_H  -DLIFEV_CHECK_ALL  -DTHREEDIM   -g -O0


LDFLAGS          = $(LIFELIBPATH) $(LIFELIBS) $(TRILLIBPATH) $(TRILLIBS) \
                   $(METISPATH) $(METISLIBS) \
	               $(MPILIBPATH) $(MPILIBS) -lhdf5 

all: $(EXECUTABLE)

$(OBJECTS): $(SOURCES)
	$(CC) -c $(CXXFLAGS) \
	$(MPIINCLUDEPATH) $(TRILINCLUDEPATH) $(LIFEINCLUDEPATH) \
	$(SOURCES) -o $@

$(EXECUTABLE): $(OBJECTS)
	echo "Linking the executable ... "
	$(CC) $< -o $@   $(LDFLAGS)

clean:
	rm -rf *o *scl *vct *vtk *case *xmf *h5 cavity_example

clean-results:
	rm -rf *o *scl *vct *vtk *case *xmf *h5

