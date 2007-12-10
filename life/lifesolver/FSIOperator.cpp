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

#include <life/lifesolver/FSIOperator.hpp>
//#include <life/lifesolver/reducedLinFluid.hpp>

namespace LifeV
{
// Constructors


// Destructor
FSIOperator::~FSIOperator()
{
}


void
FSIOperator::setup()
{
    const RefFE*    refFE_vel;
    const QuadRule* qR_vel;
    const QuadRule* bdQr_vel;

    const RefFE*    refFE_press;
    const QuadRule* qR_press;
    const QuadRule* bdQr_press;

    std::string uOrder = M_dataFluid->uOrder();

    int me       = M_epetraComm->MyPID();
    bool verbose = (me == 0);

    if (verbose) std::cout << "velocity order = " << uOrder << std::endl;

    if ( uOrder.compare("P2") == 0 )
    {
        if (verbose) std::cout << "P2 velocity " << std::flush;
        refFE_vel = &feTetraP2;
        qR_vel    = &quadRuleTetra15pt; // DoE 5
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else
        if ( uOrder.compare("P1") == 0 )
        {
            if (verbose) std::cout << "P1 velocity ";
            refFE_vel = &feTetraP1;
            qR_vel    = &quadRuleTetra4pt;  // DoE 2
            bdQr_vel  = &quadRuleTria3pt;   // DoE 2
        }
        else
            if ( uOrder.compare("P1Bubble") == 0 )
            {
                if (verbose) std::cout << "P1-bubble velocity " << std::flush;
                refFE_vel = &feTetraP1bubble;
                qR_vel    = &quadRuleTetra64pt;  // DoE 2
                bdQr_vel  = &quadRuleTria3pt;   // DoE 2
            }

//    Dof uDof(dataNavierStokes.mesh(), *refFE_vel);

    std::string pOrder =  M_dataFluid->pOrder();
    if ( pOrder.compare("P2") == 0 )
    {
        if (verbose) std::cout << "P2 pressure " << std::flush;
        refFE_press = &feTetraP2;
        qR_press    = &quadRuleTetra15pt; // DoE 5
        bdQr_press  = &quadRuleTria3pt;   // DoE 2
    }
    else
        if ( pOrder.compare("P1") == 0 )
        {
            if (verbose) std::cout << "P1 pressure";
            refFE_press = &feTetraP1;
            qR_press    = &quadRuleTetra4pt;  // DoE 2
            bdQr_press  = &quadRuleTria3pt;   // DoE 2
        }

    if (verbose) std::cout << std::endl;

    Dof uDof(*M_dataFluid->mesh(), *refFE_vel);
    Dof pDof(*M_dataFluid->mesh(), *refFE_press);

    const RefFE*    refFE_struct;
    const QuadRule* qR_struct;
    const QuadRule* bdQr_struct;

    std::string dOrder = M_dataSolid->order();
    if ( dOrder.compare("P2") == 0 )
    {
        if (verbose) std::cout << "P2 displacement " << std::flush;
        refFE_struct = &feTetraP2;
        qR_struct    = &quadRuleTetra15pt; // DoE 5
        bdQr_struct  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( dOrder.compare("P1") == 0 )
    {
        if (verbose) std::cout << "P1 displacement";
        refFE_struct = &feTetraP1;
        qR_struct    = &quadRuleTetra4pt;  // DoE 2
        bdQr_struct  = &quadRuleTria3pt;   // DoE 2
    }


    MPI_Barrier(MPI_COMM_WORLD);

    Dof dDof(*M_dataSolid->mesh(), *refFE_struct);


    if (verbose)
        std::cout << "fluid: building the FE space ... " << std::flush;


    if (this->isFluid())
    {

        M_fluidMeshPart.reset(new  partitionMesh< FSIOperator::mesh_type > (*M_dataFluid->mesh(), *M_epetraComm));




        M_mmFESpace.reset(new FESpace<mesh_type, EpetraMap>(*M_fluidMeshPart,
                                                            *refFE_struct,
                                                            *qR_struct,
                                                            *bdQr_struct,
                                                            3,
                                                            *M_epetraComm));


        M_meshMotion.reset(new FSIOperator::meshmotion_raw_type(*M_mmFESpace,
                                                                *M_epetraComm));

        M_uFESpace.reset( new FESpace<mesh_type, EpetraMap>(*M_fluidMeshPart,
                                                            *refFE_vel,
                                                            *qR_vel,
                                                            *bdQr_vel,
                                                            3,
                                                            *M_epetraComm));



        M_pFESpace.reset( new FESpace<mesh_type, EpetraMap>(*M_fluidMeshPart,
                                                            *refFE_press,
                                                            *qR_press,
                                                            *bdQr_press,
                                                            1,
                                                            *M_epetraComm));

        if (verbose)
            std::cout << "fluid: ok." << std::endl;

        M_fluid.reset(new FSIOperator::fluid_raw_type(dataFluid(),
                                                      *M_uFESpace,
                                                      *M_pFESpace,
                                                      *M_epetraComm));


        if (M_linearFluid)
            M_fluidLin.reset(new FSIOperator::fluidlin_raw_type(dataFluid(),
                                                                *M_uFESpace,
                                                                *M_pFESpace,
                                                                *M_epetraComm));




        vector_type u0(M_uFESpace->map() + M_pFESpace->map());

        M_bdf.reset(new BdfT<vector_type>(M_dataFluid->order_bdf()));
        M_bdf->initialize_unk(u0);


    }
    else
    {
        M_mmFESpace.reset(new FESpace<mesh_type, EpetraMap>(M_dataFluid->mesh(),
                                                            *refFE_struct,
                                                            *qR_struct,
                                                            *bdQr_struct,
                                                            3,
                                                            *M_epetraComm));


        M_meshMotion.reset(new FSIOperator::meshmotion_raw_type(*M_mmFESpace,
                                                                *M_epetraComm));

        M_uFESpace.reset( new FESpace<mesh_type, EpetraMap>(M_dataFluid->mesh(),
                                                            *refFE_vel,
                                                            *qR_vel,
                                                            *bdQr_vel,
                                                            3,
                                                            *M_epetraComm));



        M_pFESpace.reset( new FESpace<mesh_type, EpetraMap>(M_dataFluid->mesh(),
                                                            *refFE_press,
                                                            *qR_press,
                                                            *bdQr_press,
                                                            1,
                                                            *M_epetraComm));

        if (verbose)
            std::cout << "fluid: ok." << std::endl;

        M_fluid.reset(new FSIOperator::fluid_raw_type(dataFluid(),
                                                      *M_uFESpace,
                                                      *M_pFESpace,
                                                      *M_epetraComm));

        if (M_linearFluid)
            M_fluidLin.reset(new FSIOperator::fluidlin_raw_type(dataFluid(),
                                                                *M_uFESpace,
                                                                *M_pFESpace,
                                                                *M_epetraComm));

    }
//     M_fluid->setUp(dataFile);
//     M_fluid->buildSystem();


    //std::cout << M_epetraComm->MyPID() << std::endl;



    if (verbose)
        std::cout << "solid: building the FE space ... " << std::flush;

    if (this->isSolid())
    {
        M_solidMeshPart.reset( new  partitionMesh< FSIOperator::mesh_type > (*M_dataSolid->mesh(), *M_epetraComm));

        M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(*M_solidMeshPart,
                                                           *refFE_struct,
                                                           *qR_struct,
                                                           *bdQr_struct,
                                                           3,
                                                           *M_epetraComm));

//     M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(*M_solidMeshPart,
//                                                        *refFE_struct,
//                                                        *qR_struct,
//                                                        *bdQr_struct,
//                                                        3,
//                                                        *M_epetraComm));

        if (verbose)
            std::cout << "solid: ok." << std::endl;


        M_solid.reset(new FSIOperator::solid_raw_type(dataSolid(),
                                                      *M_dFESpace,
                                                      *M_epetraComm));

        if (M_linearSolid)
            M_solidLin.reset(new FSIOperator::solidlin_raw_type(dataSolid(),
                                                                *M_dFESpace,
                                                                *M_epetraComm));


    }
    else
    {
        M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(M_dataSolid->mesh(),
                                                           *refFE_struct,
                                                           *qR_struct,
                                                           *bdQr_struct,
                                                           3,
                                                           *M_epetraComm));

        if (verbose)
            std::cout << "solid: ok." << std::endl;


        M_solid.reset(new FSIOperator::solid_raw_type(dataSolid(),
                                                      *M_dFESpace,
                                                      *M_epetraComm));


        if (M_linearSolid)
            M_solidLin.reset(new FSIOperator::solidlin_raw_type(dataSolid(),
                                                                *M_dFESpace,
                                                                *M_epetraComm));

//        M_bdf.reset(new BdfT<vector_type>(M_dataFluid->order_bdf()));


    }



    MPI_Barrier(MPI_COMM_WORLD);

    if (this->isFluid())
    {
        M_dispFluidMeshOld.reset(new vector_type(*this->M_fluid->getMap().getRepeatedEpetra_Map()));
        M_veloFluidMesh.reset   (new vector_type(*this->M_fluid->getMap().getRepeatedEpetra_Map()));
        M_un.reset              (new vector_type(this->M_fluid->getMap()));
        M_rhs.reset             (new vector_type(this->M_fluid->getMap()));
    }




    M_dofFluidToStructure->setup(M_dFESpace->refFE(), dDof,
                                 M_uFESpace->refFE(), M_uFESpace->dof());
    M_dofFluidToStructure->update(*M_dataSolid->mesh(), 1,
                                  *M_uFESpace->mesh(),  1,
                                  0.);

    M_dofStructureToSolid->setup(M_dFESpace->refFE(), M_dFESpace->dof(),
                                 M_dFESpace->refFE(), M_dFESpace->dof());
    M_dofStructureToSolid->update(*M_dFESpace->mesh(), 1,
                                  *M_dFESpace->mesh(), 1,
                                  0.);

    M_dofStructureToHarmonicExtension->setup(M_uFESpace->refFE(), M_uFESpace->dof(),
                                             M_dFESpace->refFE(), M_dFESpace->dof());
    M_dofStructureToHarmonicExtension->update(*M_uFESpace->mesh(), 1,
                                              *M_dFESpace->mesh(), 1,
                                              0.0);

    M_dofHarmonicExtensionToFluid->setup(M_uFESpace->refFE(),  uDof,
                                         M_uFESpace->refFE(),  uDof);
    M_dofHarmonicExtensionToFluid->update(*M_dataFluid->mesh(), 1,
                                          *M_dataFluid->mesh(), 1,
                                          0.0);

    typedef std::map<ID, ID>::iterator Iterator;

    std::vector<int> dofInterfaceFluid;
    dofInterfaceFluid.reserve(M_dofHarmonicExtensionToFluid->locDofMap().size());

    std::vector<int> dofInterfaceSolid;
    dofInterfaceSolid.reserve(M_dofStructureToSolid->locDofMap().size());


    // now we build the sigma and lambda variables on each proc


    if (verbose)
        std::cout << "building the variables ... " << std::flush;



    std::map<ID, ID> locDofMap = M_dofStructureToHarmonicExtension->locDofMap();

    if (this->isFluid())
    {

        for (int dim = 0; dim < nDimensions; ++dim)
            for ( Iterator i = locDofMap.begin(); i != locDofMap.end(); ++i )
                {
                    dofInterfaceFluid.push_back(i->second + dim*dDof.numTotalDof()); // in solid numerotation
                }
    }


    M_fluidInterfaceMap.reset(new EpetraMap(-1,
                                            dofInterfaceFluid.size(),
                                            &dofInterfaceFluid[0],
                                            1,
                                            *M_epetraWorldComm));

    M_lambdaFluid.reset(new vector_type(*M_fluidInterfaceMap->getRepeatedEpetra_Map()));
    M_sigmaFluid.reset (new vector_type(*M_fluidInterfaceMap->getRepeatedEpetra_Map()));


    vector_type test = *M_lambdaFluid;
    M_epetraWorldComm->Barrier();


    if (this->isSolid())
    {
        std::cout << "solid" << std::endl;
        for (int dim = 0; dim < nDimensions; ++dim)
            for ( Iterator i = locDofMap.begin(); i != locDofMap.end(); ++i )
            {
                dofInterfaceSolid.push_back(i->second + dim*dDof.numTotalDof()); // in solid numerotation
            }
    }



    M_solidInterfaceMap.reset(new EpetraMap(-1,
                                            dofInterfaceSolid.size(),
                                            &dofInterfaceSolid[0],
                                            1,
                                            *M_epetraWorldComm));

    M_lambdaSolid.reset   (new vector_type(*M_solidInterfaceMap->getRepeatedEpetra_Map()));
    M_lambdaDotSolid.reset(new vector_type(*M_solidInterfaceMap->getRepeatedEpetra_Map()));
    M_sigmaSolid.reset    (new vector_type(*M_solidInterfaceMap->getRepeatedEpetra_Map()));

    M_epetraWorldComm->Barrier();

    if (verbose) std::cout << "done." << std::endl;
//    M_dofStructureToHarmonicExtension->showMe(true, std::cout);
//    M_dofHarmonicExtensionToFluid->showMe(true, std::cout);
//    M_dofStructureToReducedFluid->setup(uFESpace.refFE(), M_fluid->pressFESpace().dof(),
//                                             dFESpace.refFE(), dFESpace.dof());
//         M_dofStructureToReducedFluid->update(fluidMesh, 1,
//                                              solidMesh, 1,
//                                              0.0);

//         M_dofReducedFluidToStructure->setup(dFESpace.refFE(), dFESpace.dof(),
//                                             uFESpace.refFE(), uFESpace.dof());
//         M_dofReducedFluidToStructure->update(solidMesh, 1,
//                                              fluidMesh, 1,
//                                              0.0);
//     }

} // end setup

//

void
FSIOperator::setDataFromGetPot( GetPot const& data_file )
{
    //M_solverAztec.setOptionsFromGetPot(data_file,"jacobian/aztec");
    M_method  = data_file("problem/method" ,"steklovPoincare");

    M_dataFluid.reset(new data_fluid(data_file));
    M_dataSolid.reset(new data_solid(data_file));

}

//

void
FSIOperator::updateJacobian(vector_type& /*sol*/,int /*iter*/)
{
}

void
FSIOperator::updateSystem(fluid_source_type &fluidSource, solid_source_type &solidSource)
{
    if (this->isFluid())
    {
        double alpha = this->M_bdf->coeff_der( 0 ) / M_dataFluid->timestep();

//         vector_type beta = M_bdf->extrap();
//         vector_type rhs  = M_bdf->time_der( M_dataFluid->timestep() );

//        this->M_fluid->updateSystem( alpha, beta, rhs );

//         this->transferMeshMotionOnFluid(M_meshMotion->displacement(),
//                                         *this->M_veloFluidMesh);
        M_meshMotion->updateSystem();

        vector_type const meshDisplacement( M_meshMotion->displacement(), this->M_meshMotion->getRepeatedEpetraMap() );

        transferMeshMotionOnFluid(meshDisplacement,
                                  *this->M_dispFluidMeshOld);

        *this->M_un                = M_fluid->solution();
        *this->M_rhs               = M_fluid->matrMass()*(*this->M_un);

//         int numTotalDof = M_uFESpace->dof().numTotalDof();
//         for (int jj = 0; jj < 3; ++jj)
//             for (int ii = 1; ii <= 1050; ++ii)
//                 (*M_un)[ii + jj*numTotalDof] = M_fluid->solution()[ii + jj*numTotalDof];

//        M_fluid->postProcess();

    }

    if (this->isSolid())
    {
        this->M_solid->updateSystem();
    }

}

void
FSIOperator::shiftSolution()
{
    if (this->isFluid())
    {
        this->M_bdf->shift_right(M_fluid->solution());
    }
}

//


FSIOperator::vector_type
FSIOperator::displacementOnInterface()
{

//     vector_type dispOnInterface();
// //@     dispOnInterface = ZeroVector(dispOnInterface.size());

//  //    FOR_EACH_INTERFACE_DOF( dispOnInterface[IDsolid - 1 + jDim*totalDofSolid] =
// //                             M_solid->disp()[IDsolid - 1 + jDim*totalDofSolid]);

//     double norminf;
//     double norm2;

//     dispOnInterface.NormInf(&norminf);
//     dispOnInterface.Norm2(&norm2);

//     std::cout << "max norm disp = " << norminf;
//     std::cout << std::endl;
//     std::cout << "l2  norm disp = " << norm2;
//     std::cout << std::endl;

//     return dispOnInterface;
    assert(false);

}


void FSIOperator::sendToFluid(const vector_type& _vec)
{
   *M_lambdaFluid = _vec;

//    M_lambdaFluid->Map().Print( std::cout );

}


void FSIOperator::receiveFromSolid(vector_type& _vec)
{

    bool verbose = M_epetraComm->MyPID() == 0;

    if (verbose)
    {
        MPI_Status status;
        int length;
        int tag    = 777;

        MPI_Recv(&length, 1, MPI_INT, M_solidLeader, tag, MPI_COMM_WORLD, &status);

        std::cout << "  MPI: receiving ", length, " double from the solid ... ";

        std::vector<double> vec(length);

        MPI_Recv(&vec[0], length, MPI_DOUBLE, M_solidLeader, tag, MPI_COMM_WORLD, &status);

        std::cout << "done.";
    }
}



void
FSIOperator::transferMeshMotionOnFluid(const vector_type &_vec1,
                                       vector_type       &_vec2)
{

    std::map<ID, ID> locDofMap = M_dofFluidToStructure->locDofMap();


    int numTotalDofMesh  = M_mmFESpace->dof().numTotalDof();
    int numTotalDofFluid = M_uFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::iterator Iterator;

    for (int dim = 0; dim < nDimensions; ++dim)
        for ( Iterator it = locDofMap.begin(); it != locDofMap.end(); ++it )
            {
//                 std::cout << " doing: for " << it->second << " to " << it->second
//                           << " _vec1.Map().LID(it->second + dim*numTotalDofMesh) = " << _vec1.Map().LID(it->second + dim*numTotalDofMesh)
//                           << " _vec2.Map().LID(it->second + dim*numTotalDofFluid) = " << _vec2.Map().LID(it->second + dim*numTotalDofFluid)
//                           << std::endl;

                if (_vec1.Map().LID(it->second + dim*numTotalDofMesh) >= 0 )
                {
                    _vec2[it->second + dim*numTotalDofFluid] = _vec1[it->second + dim*numTotalDofMesh];
                }
//                 else
//                 {
//                     std::cout << " not done: from " << it->second << " to " << it->second << std::endl;
//                 }
            }

}

void
FSIOperator::interpolateVelocity(const vector_type &_vec1,
                                 vector_type       &_vec2)
{

    typedef mesh_type::VolumeShape GeoShape; // Element shape


    UInt nDofpV = M_uFESpace->refFE().nbDofPerVertex; // number of Dof per vertex
    UInt nDofpE = M_uFESpace->refFE().nbDofPerEdge;   // number of Dof per edge
    UInt nDofpF = M_uFESpace->refFE().nbDofPerFace;   // number of Dof per face
    UInt nDofpEl = M_uFESpace->refFE().nbDofPerVolume; // number of Dof per Volume

    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF = GeoShape::numFaces;    // Number of element's faces

    UInt nDofElem = M_uFESpace->refFE().nbDof; // Number of local dof per element of the M_uFESpace->mesh() (_mesh.getRefFE().nbDof)
    UInt nDofElemMesh = M_mmFESpace->refFE().nbDof;

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element
    UInt nDofElemF = nElemF * nDofpF; // number of face's Dof on a Element


    Real x, y, z;


    Vector wLoc( nDofElemMesh * nDimensions );

    ID lDof;

    // Loop on elements of the mesh
    for ( ID iElem = 1; iElem <= M_uFESpace->mesh()->numVolumes(); ++iElem )
    {


        int elemId = M_uFESpace->mesh()->volume( iElem ).localId();

        // Updating the local mesh velocity in this mesh elment
        for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
            for ( ID idof = 0; idof < nDofElemMesh; ++idof )
                {
//                     std::cout << icmp * M_mmFESpace->dof().numTotalDof()
//                         + M_mmFESpace->dof().localToGlobal( iElem, idof + 1 ) << std::endl;
//                     if (_vec1.Map().LID( icmp * M_mmFESpace->dof().numTotalDof()
//                                + M_mmFESpace->dof().localToGlobal( iElem, idof + 1)) >= 0 );
                    wLoc( icmp * nDofElemMesh + idof ) =
                        _vec1( icmp * M_mmFESpace->dof().numTotalDof()
                               + M_mmFESpace->dof().localToGlobal( iElem, idof + 1));
                }
        // Vertex based Dof
        if ( nDofpV )
        {

            // loop on element vertices
            for ( ID iVe = 1; iVe <= nElemV; ++iVe )
            {

                // Loop number of Dof per vertex
                for ( ID l = 1; l <= nDofpV; ++l )
                {
                    lDof = ( iVe - 1 ) * nDofpV + l; // Local dof in this element

                    // Nodal coordinates
                    x = M_uFESpace->refFE().xi( lDof - 1 );
                    y = M_uFESpace->refFE().eta( lDof - 1 );
                    z = M_uFESpace->refFE().zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        double __sum = 0;
                        for ( ID idof = 0; idof < nDofElemMesh; ++idof )  // Loop on local Dof on the element
                            __sum += wLoc( icmp * nDofElemMesh + idof ) * M_uFESpace->refFE().phi( idof, x, y, z );

                        // Updating interpolated mesh velocity
                        if (_vec2.Map().LID( icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( iElem, lDof )) >= 0)
                            {
                                int iDof = M_uFESpace->dof().localToGlobal( elemId, lDof  );
                                _vec2( icmp * M_uFESpace->dof().numTotalDof() + iDof ) = __sum;
                            }
                    }
                }
            }
        }

        // Edge based Dof
        if ( nDofpE )
        {

            // loop on element edges
            for ( ID iEd = 1; iEd <= nElemE; ++iEd )
            {

                // Loop number of Dof per edge
                for ( ID l = 1; l <= nDofpE; ++l )
                {
                    lDof = nDofElemV + ( iEd - 1 ) * nDofpE + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    x = M_uFESpace->refFE().xi( lDof - 1 );
                    y = M_uFESpace->refFE().eta( lDof - 1 );
                    z = M_uFESpace->refFE().zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        double __sum = 0;
                        for ( ID idof = 0; idof < nDofElem; ++idof )   // Loop on local Dof on the adjacent element
                            __sum += wLoc( icmp * nDofElem + idof ) * M_uFESpace->refFE().phi( idof, x, y, z );

                        // Updating interpolating vector
                        _vec2( icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( iElem, lDof ) ) = __sum;
                    }
                }
            }
        }

        // Face based Dof
        if ( nDofpF )
        {

            // loop on element faces
            for ( ID iFa = 1; iFa <= nElemF; ++iFa )
            {

                // Loop on number of Dof per face
                for ( ID l = 1; l <= nDofpF; ++l )
                {

                    lDof = nDofElemE + nDofElemV + ( iFa - 1 ) * nDofpF + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    x = M_uFESpace->refFE().xi( lDof - 1 );
                    y = M_uFESpace->refFE().eta( lDof - 1 );
                    z = M_uFESpace->refFE().zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        double __sum = 0;
                        for ( ID idof = 0; idof < nDofElem; ++idof )  // Loop on local Dof on the adjacent element
                            __sum += wLoc( icmp * nDofElem + idof ) * M_uFESpace->refFE().phi( idof, x, y, z );

                        // Updating interpolating vector
                        _vec2( icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( iElem, lDof + 1) ) = __sum;
                    }
                }
            }
        }

        // Element based Dof
        // Loop on number of Dof per Element
        for ( ID l = 1; l <= nDofpEl; ++l )
        {
            lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

            // Nodal coordinates
            x = M_uFESpace->refFE().xi( lDof - 1 );
            y = M_uFESpace->refFE().eta( lDof - 1 );
            z = M_uFESpace->refFE().zeta( lDof - 1 );

            // Loop on data vector components
            for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
            {

                // Interpolating data at the nodal point
                double __sum = 0;
                for ( ID idof = 0; idof < nDofElemMesh; ++idof )  // Loop on local Dof on the adjacent element
                    __sum += wLoc( icmp * nDofElemMesh + idof ) * M_uFESpace->refFE().phi( idof, x, y, z );

                // Updating interpolating vector
//                std::cout << M_uFESpace->dof().localToGlobal( elemId, lDof ) << " ";
//                std::cout << icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( elemId, lDof ) << std::endl;

                if (_vec2( icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( elemId, lDof )) >= 0.)
                    _vec2( icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( elemId, lDof ) ) = __sum;
            }
        }
    }

}





void
FSIOperator::transferFluidOnInterface(const vector_type &_vec1,
                                      vector_type       &_vec2)
{
    // e.g.: vec1=M_fluid->residual(), vec2=M_sigmaFluid
//     _vec2 = ZeroVector(_vec2.size());

    std::map<ID, ID> locDofMap = M_dofFluidToStructure->locDofMap();


    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();
    int numTotalDofFluid = M_uFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::iterator Iterator;

    for (int dim = 0; dim < nDimensions; ++dim)
        for ( Iterator it = locDofMap.begin(); it != locDofMap.end(); ++it )
            {
                 if (_vec1.Map().LID(it->second + dim*numTotalDofFluid) >= 0 )
                 {
                     _vec2[it->first + dim*numTotalDofSolid] = _vec1[it->second + dim*numTotalDofFluid];
                 }
            }
}


void
FSIOperator::transferSolidOnInterface(const vector_type &_vec1,
                                      vector_type       &_vec2)
{
    /* e.g.:
       vec1                vec2
       M_solid->disp()     M_lambdaSolid
       M_solid->vel()      M_lambdaDotSolid
       M_solid->residual() M_sigmaSolid
    */

    std::map<ID, ID> locDofMap = M_dofStructureToSolid->locDofMap();
    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::iterator Iterator;

    for (int dim = 0; dim < nDimensions; ++dim)
        for ( Iterator it = locDofMap.begin(); it != locDofMap.end(); ++it )
            {
                 if (_vec1.Map().LID(it->second + dim*numTotalDofSolid) >= 0 )
                 {
                    _vec2[it->first + dim*numTotalDofSolid] = _vec1[it->second + dim*numTotalDofSolid];
                 }
            }
}


void
FSIOperator::setFluidBC             (fluid_bchandler_type bc_fluid)
{
    if (isFluid())
        {
            M_BCh_u = bc_fluid;
        }
}

void
FSIOperator::setLinFluidBC          (fluid_bchandler_type bc_dfluid)
{
    M_BCh_du     = bc_dfluid;
}

void FSIOperator::setInvLinFluidBC (fluid_bchandler_type bc_dfluid_inv)
{
    M_BCh_du_inv = bc_dfluid_inv;
}

void FSIOperator::setHarmonicExtensionBC (fluid_bchandler_type bc_he)
{
    if ( isFluid() )
        {
            M_BCh_mesh = bc_he;
            M_meshMotion->setBC(*M_BCh_mesh);
        }
}

void FSIOperator::setSolidBC(solid_bchandler_type bc_solid)
{
    if ( isSolid() )
        {
            M_BCh_d      = bc_solid;
        }
}

void FSIOperator::setLinSolidBC(solid_bchandler_type bc_dsolid)
{
    M_BCh_dz     = bc_dsolid;
}

void FSIOperator::setInvLinSolidBC(solid_bchandler_type bc_dsolid_inv)
{
    M_BCh_dz_inv     = bc_dsolid_inv;
}

// void FSIOperator::setReducedLinFluidBC(fluid_bchandler_type bc_dredfluid)
// {
//     M_BCh_dp     = bc_dredfluid;
//     M_reducedLinFluid->setUpBC(bc_dredfluid);
// }

// void FSIOperator::setInvReducedLinFluidBC(fluid_bchandler_type bc_invdredfluid)
// {
//     M_BCh_dp_inv = bc_invdredfluid;
// }

//

void FSIOperator:: setComm     (   boost::shared_ptr<Epetra_MpiComm> comm,
                                   boost::shared_ptr<Epetra_MpiComm> worldComm)
{
    M_epetraComm       = comm;
    M_epetraWorldComm  = worldComm;
}



//

void FSIOperator::setHarmonicExtensionVelToFluid(const vector_type &vel,
                                                 UInt type)
{
    M_bcvHarmonicExtensionVelToFluid->setup(vel,
                                            M_fluid->velFESpace().dof().numTotalDof(),
                                            M_dofHarmonicExtensionToFluid,
                                            type);
}

void FSIOperator::setDerHarmonicExtensionVelToFluid(vector_type &dvel,
                                                    UInt type)
{
    M_bcvDerHarmonicExtensionVelToFluid->setup(dvel,
                                               M_fluid->velFESpace().dof().numTotalDof(),
                                               M_dofHarmonicExtensionToFluid,
                                               type);
}

void FSIOperator::setStructureDispToHarmonicExtension(vector_type &disp,
                                                      UInt type)
{
    M_bcvStructureDispToHarmonicExtension->setup(disp,
                                                 M_solid->dFESpace().dof().numTotalDof(),
                                                 M_dofStructureToHarmonicExtension,
                                                 type);
}

void FSIOperator::setStructureDispToSolid(vector_type &disp,
                                             UInt type)
{
    M_bcvStructureDispToSolid->setup(disp,
                                     M_solid->dFESpace().dof().numTotalDof(),
                                     M_dofStructureToSolid,
                                     type);
}

void FSIOperator::setDerStructureDispToSolid(vector_type &ddisp,
                                             UInt type)
{
    M_bcvDerStructureDispToSolid->setup(ddisp,
                                        M_solid->dFESpace().dof().numTotalDof(),
                                        M_dofStructureToSolid,
                                        type);
}


void FSIOperator::setFluidLoadToStructure(vector_type &load,
                                          UInt type)
{
    M_bcvFluidLoadToStructure->setup(load,
                                     M_solid->dFESpace().dof().numTotalDof(),
                                     M_dofStructureToSolid,
                                     type);
}

void FSIOperator::setDerFluidLoadToStructure(vector_type &dload,
                                             UInt type)
{
    M_bcvDerFluidLoadToStructure->setup(dload,
                                        M_fluid->velFESpace().dof().numTotalDof(),
                                        M_dofFluidToStructure,
                                        type);
}

void FSIOperator::setDerFluidLoadToFluid(vector_type &dload,
                                         UInt type)
{
    M_bcvDerFluidLoadToFluid->setup(dload,
                                    M_fluid->velFESpace().dof().numTotalDof(),
                                    M_dofHarmonicExtensionToFluid,
                                    type);
}

// void FSIOperator::setDerReducedFluidLoadToStructure(Vector &dload,
//                                                        UInt type)
// {
//     M_bcvDerReducedFluidLoadToStructure->setup(dload,
//                                                M_fluid->velFESpace().dof().numTotalDof(),
//                                                M_dofReducedFluidToStructure,
//                                                type);
// }

// void FSIOperator::setDerStructureAccToReducedFluid(Vector &acc,
//                                                        UInt type)
// {

//     M_bcvDerStructureAccToReducedFluid->setup(acc,
//                                               M_solid->dFESpace().dof().numTotalDof(),
//                                               M_dofStructureToReducedFluid,
//                                               type);
// }


}
