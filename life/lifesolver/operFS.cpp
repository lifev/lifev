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
#include "quasiNewton.hpp"

namespace LifeV
{
// Constructors

operFS::operFS( fluid_type& fluid,
                solid_type& solid,
                GetPot    &data_file,
                bchandler_type& BCh_u,
                bchandler_type& BCh_d,
                bchandler_type& BCh_mesh):
    M_BCh_u       (BCh_u),
    M_BCh_d       (BCh_d),
    M_BCh_mesh    (BCh_mesh),
    M_fluid(fluid),
    M_solid(solid),
    M_quasiNewton(),
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
    M_dofStructureToReducedFluid( new DofInterface3Dto3D(feTetraP1,
                                                         M_fluid->pDof(),
                                                         feTetraP1,
                                                         M_solid->dDof()) ),
    M_dofReducedFluidToStructure( new DofInterface3Dto3D(feTetraP1,
                                                         solid->dDof(),
                                                         feTetraP1,
                                                         fluid->pDof()) ),
    M_dispStruct  ( 3*M_solid->dDof().numTotalDof() ),
    M_velo        ( 3*M_solid->dDof().numTotalDof() ),
    M_nbEval      (0)
{
    M_dofFluidToStructure->update(M_solid->mesh(), 1,
                                  M_fluid->mesh(), 1,
                                  0.);
    M_dofStructureToSolid->update(M_solid->mesh(), 1,
                                  M_solid->mesh(), 1,
                                  0.0);
    M_dofStructureToFluidMesh->update(M_fluid->mesh(), 1,
                                      M_solid->mesh(), 1,
                                      0.0);
    M_dofMeshToFluid->update(M_fluid->mesh(), 1,
                            M_fluid->mesh(), 1,
                            0.0);
    M_dofStructureToReducedFluid->update(M_fluid->mesh(), 1,
                                         M_solid->mesh(), 1,
                                         0.0);
    M_dofReducedFluidToStructure->update(fluid->mesh(), 1,
                                         solid->mesh(), 1,
                                         0.0);

    M_solverAztec.setOptionsFromGetPot(data_file,"jacobian/aztec");
    M_method       = data_file("problem/method" , 0);
    M_reducedFluid = data_file("problem/reducedFluid", 0);
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
        M_dispStruct.resize( 3*M_solid->dDof().numTotalDof() );
        M_velo.resize( 3*M_solid->dDof().numTotalDof() );

        M_dofFluidToStructure->setup(feTetraP1, M_solid->dDof(),
                                     feTetraP1bubble, M_fluid->uDof());
        M_dofFluidToStructure->update(M_solid->mesh(), 1,
                                      M_fluid->mesh(), 1,
                                      0.);

        M_dofStructureToSolid->setup(feTetraP1, M_solid->dDof(),
                                     feTetraP1, M_solid->dDof());
        M_dofStructureToSolid->update(M_solid->mesh(), 1,
                                      M_solid->mesh(), 1,
                                      0.);

        M_dofStructureToFluidMesh->setup(M_fluid->mesh().getRefFE(), M_fluid->dofMesh(),
                                         feTetraP1, M_solid->dDof());
        M_dofStructureToFluidMesh->update(M_fluid->mesh(), 1,
                                          M_solid->mesh(), 1,
                                          0.0);

        M_dofMeshToFluid->setup(feTetraP1bubble,M_fluid->uDof(),feTetraP1bubble,M_fluid->uDof());
        M_dofMeshToFluid->update(M_fluid->mesh(), 1,
                                 M_fluid->mesh(), 1,
                                 0.0);

        M_dofStructureToReducedFluid->setup(feTetraP1, M_fluid->pDof(),
                                            feTetraP1, M_solid->dDof());
        M_dofStructureToReducedFluid->update(M_fluid->mesh(), 1,
                                             M_solid->mesh(), 1,
                                             0.0);

        M_dofReducedFluidToStructure->setup(feTetraP1, M_solid->dDof(),
                                            feTetraP1, M_fluid->pDof());
        M_dofReducedFluidToStructure->update(M_fluid->mesh(), 1,
                                             M_solid->mesh(), 1,
                                             0.0);
    }
}
//

void
operFS::setDataFromGetPot( GetPot const& data_file )
{
    M_solverAztec.setOptionsFromGetPot(data_file,"jacobian/aztec");
    M_method  = data_file("problem/method" ,0);
    M_quasiNewton->setLinearSolver(data_file);
    M_reducedFluid = data_file("problem/reducedFluid", 0);
}


//

void
operFS::updateJac(Vector& sol,int iter)
{
}

//


Vector
operFS::displacementOnInterface()
{

    Vector dispOnInterface(M_solid->d().size());
    dispOnInterface = ZeroVector(dispOnInterface.size());

    FOR_EACH_INTERFACE_DOF( dispOnInterface[IDsolid - 1 + jDim*totalDofSolid] =
                            M_solid->d()              [IDsolid - 1 + jDim*totalDofSolid]);

    std::cout << "max norm disp = " << norm_inf(dispOnInterface);
    std::cout << std::endl;
    std::cout << "l2  norm disp = " << norm_2(dispOnInterface);
    std::cout << std::endl;

    return dispOnInterface;
}


void operFS::transferOnInterface(const Vector      &_vec1,
                                 const BCHandler   &_BC,
                                 const std::string &_BCName,
                                 Vector            &_vec2)
{
    int iBC = _BC.getBCbyName(_BCName);

    BCBase const &BCInterface = _BC[(UInt) iBC];

    UInt nDofInterface = BCInterface.list_size();

//     std::cout << "nDofInterface = " << nDofInterface << std::endl;

    UInt nDim = BCInterface.numberOfComponents();

    UInt totalDof1 = _vec1.size()/ nDim;
    UInt totalDof2 = _vec2.size()/ nDim;

    for (UInt iBC = 1; iBC <= nDofInterface; ++iBC)
    {
        ID ID1 = BCInterface(iBC)->id();

        BCVectorInterface const *BCVInterface =
            static_cast <BCVectorInterface const *>
            (BCInterface.pointerToBCVector());

        ID ID2 = BCVInterface->
            dofInterface().getInterfaceDof(ID1);

        for (UInt jDim = 0; jDim < nDim; ++jDim)
        {
            _vec2[ID2 - 1 + jDim*totalDof2] =
                _vec1[ID1 - 1 + jDim*totalDof1];
        }
    }
}


}
