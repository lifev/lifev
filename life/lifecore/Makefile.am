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

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc .cc


AM_CPPFLAGS	= -DTHREEDIM


include_HEADERS		= ElasticStructureHandler.hpp GetPot.hpp NavierStokesAleHandler.hpp NavierStokesAleSolverPC.hpp \
NavierStokesHandler.hpp NavierStokesHandler_miguel.hpp NavierStokesSolverPC.hpp Operator01.hpp \
RNM.hpp RNM_op.hpp RNM_opc.hpp RNM_tpl.hpp \
SimpleVect.hpp VenantKirchhofSolver.hpp algebraic_facto.hpp assemb.hpp \
bc_manage.hpp bdf.hpp bdfNS.hpp chrono.hpp \
dataAztec.hpp dataElasticStructure.hpp dataMesh.hpp dataNavierStokes.hpp \
dataNewton.hpp dataSimplifiedStructure.hpp dataString.hpp dataTime.hpp \
elemMat.hpp elemVec.hpp fortran_wrap.hpp gmsh_wrtrs.hpp \
inexactLU.hpp lifeV.hpp linesearch_cubic.hpp linesearch_parabolic.hpp \
medit_wrtrs.hpp mixed_pattern.hpp newton.hpp openDX_wrtrs.hpp \
pattern.hpp picard.hpp post_proc.hpp saddlePointCG.hpp \
simplifiedStructure.hpp subDomainNeighbors.hpp switch.hpp tab.hpp \
util_string.hpp values.hpp variables.hpp vecUnknown.hpp vtk_wrtrs.hpp


lib_LTLIBRARIES		= liblifecore.la
liblifecore_la_SOURCES	= GetPot.cpp assemb.cpp bdf.cpp bdfNS.cpp \
dataAztec.cpp dataNewton.cpp dataString.cpp dataTime.cpp \
elemMat.cpp elemVec.cpp medit_wrtrs.cpp pattern.cpp \
subDomainNeighbors.cpp switch.cpp util_string.cpp values.cpp \
vecUnknown.cpp vtk_wrtrs.cpp 

liblifecore_la_LIBADD	= $(aztec_libs)
