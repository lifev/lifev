/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include "operFS.hpp"

namespace LifeV
{

// Constructors

operFS::operFS(GetPot     &data_file):
    M_fluid       (data_file,
                   feTetraP1bubble,
                   feTetraP1,
                   quadRuleTetra64pt,
                   quadRuleTria3pt,
                   quadRuleTetra64pt,
                   quadRuleTria3pt,
                   M_BCh_u, M_BCh_mesh),
    M_solid       (data_file,
                   feTetraP1,
                   quadRuleTetra4pt,
                   quadRuleTria3pt,
                   M_BCh_d),
    M_BCh_u       (3),
    M_BCh_d       (3),
    M_BCh_mesh    (4),
    M_BCh_du      (2),
    M_BCh_dz      (3),
    M_dofFluidToStructure(feTetraP1,
                          M_solid.dDof(),
                          feTetraP1bubble,
                          M_fluid.uDof()),
    M_dofStructureToSolid(feTetraP1,
                          M_solid.dDof(),
                          feTetraP1,
                          M_solid.dDof()),
    M_dofStructureToFluidMesh(M_fluid.mesh().getRefFE(),
                            M_fluid.dofMesh(),
                            feTetraP1,
                            M_solid.dDof()),
    M_dofMeshToFluid(feTetraP1bubble,
                     M_fluid.uDof(),
                     feTetraP1bubble,
                     M_fluid.uDof()),
    M_dispStruct  ( 3*M_solid.dDof().numTotalDof() ),
    M_velo        ( 3*M_solid.dDof().numTotalDof() ),
    M_nbEval      (0)
{
    M_dofFluidToStructure.update(M_solid.mesh(),
                               1,
                               M_fluid.mesh(),
                               1,
                               0.);

    M_dofStructureToSolid.update(M_solid.mesh(),
                               1,
                               M_solid.mesh(),
                               1,
                               0.);

    M_dofStructureToFluidMesh.update(M_fluid.mesh(),
                                   1,
                                   M_solid.mesh(),
                                   1,
                                   0.0);

    M_dofMeshToFluid.update(M_fluid.mesh(),
                          1,
                          M_fluid.mesh(),
                          1,
                          0.0);

    M_solverAztec.setOptionsFromGetPot(data_file,"jacobian/aztec");
    M_method  = data_file("problem/method" ,0);

    M_fluid.initialize(u0);
    M_solid.initialize(d0,w0);
}

// Destructor
operFS::~operFS()
{
}

//

void  operFS::updateJac(Vector& sol,int iter) {
}

//


}
