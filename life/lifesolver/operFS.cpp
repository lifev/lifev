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

operFS::operFS( fluid_type& fluid,
                solid_type& solid,
                GetPot    &data_file,
                BCHandler &BCh_u,
                BCHandler &BCh_d,
                BCHandler &BCh_mesh)
    :
    M_BCh_u       (BCh_u),
    M_BCh_d       (BCh_d),
    M_BCh_mesh    (BCh_mesh),
//     M_fluid       (data_file,
//                    feTetraP1bubble,
//                    feTetraP1,
//                    quadRuleTetra64pt,
//                    quadRuleTria3pt,
//                    quadRuleTetra64pt,
//                    quadRuleTria3pt,
//                    M_BCh_u, M_BCh_mesh),
//     M_solid       (data_file,
//                    feTetraP1,
//                    quadRuleTetra4pt,
//                    quadRuleTria3pt,
//                    M_BCh_d),
    M_fluid(fluid),
    M_solid(solid),
    M_dofFluidToStructure( new DofInterface3Dto3D( feTetraP1,
                                                   M_solid->dDof(),
                                                   feTetraP1bubble,
                                                   M_fluid->uDof()) ),
    M_dofStructureToSolid( new DofInterface3Dto3D( feTetraP1,
                                                   M_solid->dDof(),
                                                   feTetraP1,
                                                   M_solid->dDof()) ),
    M_dofStructureToFluidMesh( new DofInterface3Dto3D( M_fluid->mesh().getRefFE(),
                                                       M_fluid->dofMesh(),
                                                       feTetraP1,
                                                       M_solid->dDof()) ),
    M_dofMeshToFluid( new DofInterface3Dto3D( feTetraP1bubble,
                                              M_fluid->uDof(),
                                              feTetraP1bubble,
                                              M_fluid->uDof()) ),
    M_dispStruct  ( 3*M_solid->dDof().numTotalDof() ),
    M_velo        ( 3*M_solid->dDof().numTotalDof() ),
    M_nbEval      (0)
{
    M_dofFluidToStructure->update(M_solid->mesh(),
                                  1,
                                  M_fluid->mesh(),
                                  1,
                                  0.);
    M_dofStructureToSolid->update(M_solid->mesh(),
                                  1,
                                  M_solid->mesh(),
                                  1,
                                  0.);
    M_dofStructureToFluidMesh->update(M_fluid->mesh(),
                                      1,
                                      M_solid->mesh(),
                                      1,
                                      0.0);
    M_dofMeshToFluid->update(M_fluid->mesh(),
                            1,
                            M_fluid->mesh(),
                            1,
                            0.0);
    M_solverAztec.setOptionsFromGetPot(data_file,"jacobian/aztec");
    M_method  = data_file("problem/method" ,0);
}

// Destructor
operFS::~operFS()
{
}

//
void
operFS::setup()
{
    if ( M_solid && M_fluid )
    {
        M_dofFluidToStructure->setup(feTetraP1,M_solid->dDof(),feTetraP1bubble,M_fluid->uDof());
        M_dofFluidToStructure->update(M_solid->mesh(),
                                     1,
                                     M_fluid->mesh(),
                                     1,
                                     0.);
        M_dofStructureToSolid->setup(feTetraP1,M_solid->dDof(),feTetraP1,M_solid->dDof());
        M_dofStructureToSolid->update(M_solid->mesh(),
                                     1,
                                     M_solid->mesh(),
                                     1,
                                     0.);

        M_dofStructureToFluidMesh->setup(M_fluid->mesh().getRefFE(),M_fluid->dofMesh(),feTetraP1,M_solid->dDof());
        M_dofStructureToFluidMesh->update(M_fluid->mesh(),
                                         1,
                                         M_solid->mesh(),
                                         1,
                                         0.0);

        M_dofMeshToFluid->setup(feTetraP1bubble,M_fluid->uDof(),feTetraP1bubble,M_fluid->uDof());
        M_dofMeshToFluid->update(M_fluid->mesh(),
                                1,
                                M_fluid->mesh(),
                                1,
                                0.0);
    }
}
//

void
operFS::setDataFromGetPot( GetPot const& data_file )
{
    M_solverAztec.setOptionsFromGetPot(data_file,"jacobian/aztec");
    M_method  = data_file("problem/method" ,0);
}


//

void  operFS::updateJac(Vector& sol,int iter)
{
}

//


void operFS::displacementOnInterface()
{
    UInt iBCf = M_fluid->BC_fluid().getBCbyName("Interface");
    UInt iBCd = M_solid->BC_solid().getBCbyName("Interface");

    BCBase const &BC_fluidInterface = M_fluid->BC_fluid()[iBCf];
    BCBase const &BC_solidInterface = M_solid->BC_solid()[iBCd];

    UInt nDofInterface = BC_fluidInterface.list_size();

    UInt nDimF = BC_fluidInterface.numberOfComponents();
    UInt nDimS = BC_solidInterface.numberOfComponents();

    UInt totalDofFluid = M_fluid->getDisplacement().size()/ nDimF;
    UInt totalDofSolid = M_solid->d().size()/ nDimS;

    Vector dispOnInterface(M_solid->d().size());
    dispOnInterface = ZeroVector(dispOnInterface.size());

    for (UInt iBC = 1; iBC <= nDofInterface; ++iBC)
    {
        ID IDfluid = BC_fluidInterface(iBC)->id();

        BCVectorInterface const *BCVInterface =
            static_cast <BCVectorInterface const *>
            (BC_fluidInterface.pointerToBCVector());

        ID IDsolid = BCVInterface->
            dofInterface().getInterfaceDof(IDfluid);


        for (UInt jDim = 0; jDim < nDimF; ++jDim)
        {
            dispOnInterface[IDsolid - 1 + jDim*totalDofSolid] =
                M_solid->d()              [IDsolid - 1 + jDim*totalDofSolid] -
                M_fluid->getDisplacement()[IDfluid - 1 + jDim*totalDofFluid];
        }
    }

    std::cout << "max norm disp = " << norm_inf(dispOnInterface);
    std::cout << std::endl;
    std::cout << "l2  norm disp = " << norm_2(dispOnInterface);
    std::cout << std::endl;
}


}
