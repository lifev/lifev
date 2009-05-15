/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

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
#ifndef TWODIM
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
#ifndef TWODIM
    const RefFE*    refFE_vel(0);
    const QuadRule* qR_vel(0);
    const QuadRule* bdQr_vel(0);

    const RefFE*    refFE_press(0);
    const QuadRule* qR_press(0);
    const QuadRule* bdQr_press(0);

    std::string uOrder = M_dataFluid->uOrder();

    // int me       = M_epetraComm->MyPID();
    leaderPrint("velocity order = " + uOrder);

    if ( uOrder.compare("P2") == 0 )
    {
        leaderPrint( "P2 velocity ");
        refFE_vel = &feTetraP2;
        qR_vel    = &quadRuleTetra15pt; // DoE 5
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else
        if ( uOrder.compare("P1") == 0 )
        {
	    leaderPrint("P1 velocity ");
            refFE_vel = &feTetraP1;
            qR_vel    = &quadRuleTetra4pt;  // DoE 2
            bdQr_vel  = &quadRuleTria3pt;   // DoE 2
        }
        else
            if ( uOrder.compare("P1Bubble") == 0 )
            {
                leaderPrint("P1-bubble velocity ");
                refFE_vel = &feTetraP1bubble;
                qR_vel    = &quadRuleTetra64pt;  // DoE 2
                bdQr_vel  = &quadRuleTria3pt;   // DoE 2
            }
            else
            {
                std::cout << uOrder << " pressure FE not implemented yet." << std::endl;
                exit(0);
            }


//    Dof uDof(dataNavierStokes.mesh(), *refFE_vel);

    std::string pOrder =  M_dataFluid->pOrder();
    if ( pOrder.compare("P2") == 0 )
    {
        leaderPrint("P2 pressure ");
        refFE_press = &feTetraP2;
        qR_press    = qR_vel; // DoE 5
        bdQr_press  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( pOrder.compare("P1") == 0 )
    {
        leaderPrint("P1 pressure");
        refFE_press = &feTetraP1;
//        qR_press    = &quadRuleTetra64pt;  // DoE 2
        qR_press    = qR_vel;  // DoE 2
        bdQr_press  = &quadRuleTria3pt;   // DoE 2
    }
    else
    {
        std::cout << pOrder << " pressure FE not implemented yet." << std::endl;
        exit(0);
    }

    leaderPrint("\n");

    Dof uDof(*M_dataFluid->mesh(), *refFE_vel);
    Dof pDof(*M_dataFluid->mesh(), *refFE_press);

    const RefFE*    refFE_struct(0);
    const QuadRule* qR_struct(0);
    const QuadRule* bdQr_struct(0);

    std::string dOrder = M_dataSolid->order();
    if ( dOrder.compare("P2") == 0 )
    {
        leaderPrint("P2 displacement ");
        refFE_struct = &feTetraP2;
        qR_struct    = &quadRuleTetra15pt; // DoE 5
        bdQr_struct  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( dOrder.compare("P1") == 0 )
    {
        leaderPrint("P1 displacement");
        refFE_struct = &feTetraP1;
        qR_struct    = &quadRuleTetra4pt;  // DoE 2
        bdQr_struct  = &quadRuleTria3pt;   // DoE 2
    }

    leaderPrint("\n");

    MPI_Barrier(MPI_COMM_WORLD);

    Dof dDof(*M_dataSolid->mesh(), *refFE_struct);


    leaderPrint("fluid: building the FE space ... " );


    if (this->isFluid())
    {

        M_fluidMeshPart.reset(new  partitionMesh< FSIOperator::mesh_type > (*M_dataFluid->mesh(), *M_epetraComm));




        M_mmFESpace.reset(new FESpace<mesh_type, EpetraMap>(*M_fluidMeshPart,
                                                            *refFE_struct,
                                                            *qR_struct,
                                                            *bdQr_struct,
                                                            3,
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

        leaderPrint("fluid: ok.\n");


//         if (M_linearFluid)
//             M_fluidLin.reset(new FSIOperator::fluidlin_raw_type(dataFluid(),
//                                                                 *M_uFESpace,
//                                                                 *M_pFESpace,
//                                                                 *M_epetraComm));


        resetHeAndFluid();
    }
    else
    {
        M_mmFESpace.reset(new FESpace<mesh_type, EpetraMap>(M_dataFluid->mesh(),
                                                            *refFE_struct,
                                                            *qR_struct,
                                                            *bdQr_struct,
                                                            3,
                                                            *M_epetraComm));


        M_meshMotion.reset(new meshmotion_raw_type(*M_mmFESpace,
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

        leaderPrint("fluid: ok.\n");

        M_fluid.reset(new fluid_raw_type(dataFluid(),
                                                      *M_uFESpace,
                                                      *M_pFESpace,
                                                      *M_epetraComm));

//         if (M_linearFluid)
//             M_fluidLin.reset(new FSIOperator::fluidlin_raw_type(dataFluid(),
//                                                                 *M_uFESpace,
//                                                                 *M_pFESpace,
//                                                                 *M_epetraComm));

    }



    leaderPrint("solid: building the FE space ... " );

    if (this->isSolid())
    {

        solidInit(refFE_struct, bdQr_struct, qR_struct);

        leaderPrint("solid: ok.\n");

    }
    else
    {
        M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(M_dataSolid->mesh(),
                                                           *refFE_struct,
                                                           *qR_struct,
                                                           *bdQr_struct,
                                                           3,
                                                           *M_epetraComm));

        leaderPrint( "solid: ok.\n");


        M_solid.reset(new solid_raw_type(dataSolid(),
                                                      *M_dFESpace,
                                                      *M_epetraComm));


//         if (M_linearSolid)
//             M_solidLin.reset(new FSIOperator::solidlin_raw_type(dataSolid(),
//                                                                 *M_dFESpace,
//                                                                 *M_epetraComm));

    }



    MPI_Barrier(MPI_COMM_WORLD);

    M_dofFluidToStructure->setup(M_dFESpace->refFE(), dDof,
                                         M_uFESpace->refFE(), M_uFESpace->dof());
    M_dofFluidToStructure->update(*M_dataSolid->mesh(), M_fluidInterfaceFlag,
                                  *M_uFESpace->mesh(),  M_structureInterfaceFlag,
                                  M_interfaceTolerance);

    //here the solid mesh must be non partitioned in the monolithic case
    M_dofStructureToHarmonicExtension->setup(M_uFESpace->refFE(), M_uFESpace->dof(),
                                             M_dFESpace->refFE(), M_dFESpace->dof());
    M_dofStructureToHarmonicExtension->update(*M_uFESpace->mesh(), M_structureInterfaceFlag,
                                              *M_dFESpace->mesh(), M_harmonicInterfaceFlag,
                                              M_interfaceTolerance);

    M_dofStructureToSolid->setup(M_dFESpace->refFE(), M_dFESpace->dof(),
                                 M_dFESpace->refFE(), M_dFESpace->dof());
    M_dofStructureToSolid->update(*M_dFESpace->mesh(), M_structureInterfaceFlag,
                                  *M_dFESpace->mesh(), M_solidInterfaceFlag,
                                  M_interfaceTolerance);

    M_dofStructureToFluid->setup(M_uFESpace->refFE(), M_uFESpace->dof(), //modifica matteo FSI
                                 M_dFESpace->refFE(), M_dFESpace->dof());

    M_dofStructureToFluid->update(*M_uFESpace->mesh(), M_structureInterfaceFlag,
                                  //*M_dataFluid->mesh(), M_structureInterfaceFlag,
                                  *M_dataSolid->mesh(), M_fluidInterfaceFlag,
                                  M_interfaceTolerance);

    M_dofHarmonicExtensionToFluid->setup(M_uFESpace->refFE(),  uDof,
                                         M_uFESpace->refFE(),  uDof);
    M_dofHarmonicExtensionToFluid->update(*M_dataFluid->mesh(), M_harmonicInterfaceFlag,
                                          *M_dataFluid->mesh(), M_fluidInterfaceFlag,
                                          M_interfaceTolerance);

    typedef std::map<ID, ID>::const_iterator Iterator;

    std::vector<int> dofInterfaceFluid;
    dofInterfaceFluid.reserve(M_dofHarmonicExtensionToFluid->locDofMap().size());

    std::vector<int> dofInterfaceSolid;
    dofInterfaceSolid.reserve(M_dofStructureToSolid->locDofMap().size());


    // now we build the sigma and lambda variables on each proc


    leaderPrint("building the variables ... ");

    //is the interface map between HE (first) and solid (second)
    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();

    if (this->isFluid())
    {

        for (int dim = 0; dim < (int)nDimensions; ++dim)
            for ( Iterator i = locDofMap.begin(); i != locDofMap.end(); ++i )
                {
                    dofInterfaceFluid.push_back(i->second + dim*dDof.numTotalDof()); // in solid numerotation
                }
    }

    int* pointerToDofs(0);
    if (dofInterfaceFluid.size() > 0) pointerToDofs = &dofInterfaceFluid[0];

    M_fluidInterfaceMap.reset(new EpetraMap(-1,
                                            dofInterfaceFluid.size(),
                                            pointerToDofs,
                                            1,
                                            *M_epetraWorldComm));

//    vector_type test = *M_lambdaFluid;
    M_epetraWorldComm->Barrier();


    if (this->isSolid())
    {
        //std::cout << "solid" << std::endl;
        for (int dim = 0; dim < (int)nDimensions; ++dim)
            for ( Iterator i = locDofMap.begin(); i != locDofMap.end(); ++i )
                {
                    dofInterfaceSolid.push_back(i->second + dim*dDof.numTotalDof()); // in solid numerotation
                }
    }


    pointerToDofs = 0;
    if (dofInterfaceSolid.size() > 0) pointerToDofs = &dofInterfaceSolid[0];

    M_solidInterfaceMap.reset(new EpetraMap(-1,
                                            dofInterfaceSolid.size(),
                                            pointerToDofs,
                                            1,
                                            *M_epetraWorldComm));
    variablesInit( refFE_struct, bdQr_struct, qR_struct);
    M_epetraWorldComm->Barrier();

    leaderPrint(" done.\n");
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
#endif
} // end setup

//

void
FSIOperator::setDataFromGetPot( GetPot const& data_file )
{

    M_method  = data_file("problem/method" ,"steklovPoincare");
    M_algorithm  = data_file("problem/algorithm" ,"DirichletNeumann");

    M_fluidInterfaceFlag      = data_file("interface/fluid_flag",     M_fluidInterfaceFlag );
    M_solidInterfaceFlag      = data_file("interface/solid_flag",     M_fluidInterfaceFlag );
    M_structureInterfaceFlag  = data_file("interface/structure_flag", M_fluidInterfaceFlag );
    M_harmonicInterfaceFlag   = data_file("interface/harmonic_flag",  M_fluidInterfaceFlag );
    M_interfaceTolerance      = data_file("interface/tolerance",      0. );

    M_dataFluid.reset(new data_fluid(data_file));
    M_dataSolid.reset(new data_solid(data_file));
    M_monolithic  = !(M_method.compare("monolithic"));
}


void FSIOperator::leaderPrint(string const message, double const number) const
{
  if ( isLeader() )
    std::cout << message << number << std::endl;

}

void FSIOperator::leaderPrint(string const message) const
{
  if ( isLeader() )
    std::cout << message << std::flush;

}

void FSIOperator::leaderPrintMax(string const message, double const number) const
{
  double num(number);
  double globalMax;
  M_epetraWorldComm->MaxAll(&num, &globalMax, 1);

  leaderPrint( message , globalMax );

}


//

void
FSIOperator::updateJacobian(vector_type& /*sol*/,int /*iter*/)
{
}


void
FSIOperator::initializeFluid( const vector_type& velAndPressure,
                              const vector_type& displacement )
{
    this->fluid().initialize( velAndPressure );
    this->meshMotion().initialize( displacement );
    this->moveMesh( displacement);
}

void
FSIOperator::initializeSolid( const vector_type& displacement,
                          const vector_type& velocity )
{
    this->solid().initialize( displacement, velocity);
}

void
FSIOperator::moveMesh(vector_type const &dep)
{
    leaderPrint( "  Moving the mesh ... ");
    M_fluidMeshPart->mesh()->moveMesh(dep,  this->M_mmFESpace->dof().numTotalDof());
    leaderPrint(  " done.\n" );
    M_fluid->recomputeMatrix(true);
}



void
FSIOperator::setUpSystem( GetPot const& data_file )
{


    if (this->isFluid())
    {
        M_fluid->setUp(data_file);
        M_meshMotion->setUp(data_file);
//         if (M_linearFluid)
//             M_fluidLin->setUp(data_file);
    }

    if (this->isSolid())
    {
        M_solid->setUp(data_file);
//         if (M_linearSolid)
//             M_solidLin->setUp(data_file);
    }
}


void
FSIOperator::buildSystem()
{
    if (this->isFluid())
    {
        M_fluid->buildSystem();
//         if (M_linearFluid)
//             M_fluidLin->buildSystem();
    }

    if (this->isSolid())
    {
        M_solid->buildSystem();
//         if (M_linearSolid)
//             M_solidLin->buildSystem();
    }
}

void
FSIOperator::updateSystem(const vector_type& /*lambda*/)
{

    shiftSolution();

    if (this->isFluid())
    {
//        double alpha = this->M_bdf->coeff_der( 0 ) / M_dataFluid->timestep();

//         vector_type beta = M_bdf->extrap();
//         vector_type rhs  = M_bdf->time_der( M_dataFluid->timestep() );

//        this->M_fluid->updateSystem( alpha, beta, rhs );

//         this->transferMeshMotionOnFluid(M_meshMotion->displacement(),
//                                         *this->M_veloFluidMesh);

        M_meshMotion->updateSystem();

        transferMeshMotionOnFluid(M_meshMotion->disp(),
                                  *this->M_dispFluidMeshOld);

        *M_un                = M_fluid->solution();

        //*this->M_rhs               = M_fluid->matrMass()* (*this->M_un);
        //*this->M_rhs*=this->M_bdf->coeff_der( 0 ) / M_dataFluid->timestep();

        *M_rhs               = M_fluid->matrMass()*M_bdf->time_der( M_dataFluid->timestep() );

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



void
FSIOperator::transferMeshMotionOnFluid(const vector_type &_vec1,
                                       vector_type       &_vec2)
{

    //transferMeshMotionOnFluid should handle the repetition of the interface nodes.
    if (_vec1.getMaptype() == Unique)
        {
            vector_type const  vec1Repeated(_vec1, Repeated);
            transferMeshMotionOnFluid(vec1Repeated, _vec2);
            return;
        }

    if (_vec2.getMaptype() == Repeated)
        {
            vector_type  vec2Unique(_vec2, Unique);
            transferMeshMotionOnFluid(_vec1, vec2Unique);
            _vec2 = vec2Unique;
            return;
        }

    _vec2 *=0;

    interpolateVelocity(_vec1, _vec2);
    return;

    /*

    std::map<ID, ID> const& locDofMap = M_dofFluidToStructure->locDofMap();


    int numTotalDofMesh  = M_mmFESpace->dof().numTotalDof();
    int numTotalDofFluid = M_uFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::iterator Iterator;

    for (int dim = 0; dim < (int)nDimensions; ++dim)
        for ( Iterator it = locDofMap.begin(); it != locDofMap.end(); ++it )
            {
//                 std::cout << " doing: for " << it->second << " to " << it->second
//                           << " _vec1.Map().LID(it->second + dim*numTotalDofMesh) = " << _vec1.Map().LID(it->second + dim*numTotalDofMesh)
//                           << " _vec2.Map().LID(it->second + dim*numTotalDofFluid) = " << _vec2.Map().LID(it->second + dim*numTotalDofFluid)
//                           << std::endl;

                if (_vec1.BlockMap().LID(it->second + dim*numTotalDofMesh) >= 0 )
                {
                    _vec2[it->second + dim*numTotalDofFluid] = _vec1[it->second + dim*numTotalDofMesh];
                }
//                 else
//                 {
//                     std::cout << " not done: from " << it->second << " to " << it->second << std::endl;
//                 }
            }
    */
}

void
FSIOperator::interpolateVelocity(const vector_type &_vec1,
                                 vector_type       &_vec2)
{

    assert(_vec1.getMaptype() == Repeated);
    assert(_vec2.getMaptype() == Unique);

    typedef mesh_type::VolumeShape GeoShape; // Element shape


    UInt nDofpV = M_uFESpace->refFE().nbDofPerVertex; // number of Dof per vertex
    UInt nDofpE = M_uFESpace->refFE().nbDofPerEdge;   // number of Dof per edge
    UInt nDofpF = M_uFESpace->refFE().nbDofPerFace;   // number of Dof per face
    UInt nDofpEl = M_uFESpace->refFE().nbDofPerVolume; // number of Dof per Volume

    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF = GeoShape::numFaces;    // Number of element's faces

    //    UInt nDofElem = M_uFESpace->refFE().nbDof; // Number of local dof per element of the M_uFESpace->mesh() (_mesh.getRefFE().nbDof)
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


        UInt elemId = M_uFESpace->mesh()->volume( iElem ).localId();
        if (elemId != iElem)
            std::cout << " elemId = " << elemId << " iElem = " << iElem << std::endl;

        // Updating the local mesh velocity in this mesh elment
        for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
            for ( ID idof = 0; idof < nDofElemMesh; ++idof )
                {
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
                            __sum += wLoc( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi( idof, x, y, z );

                        // Updating interpolated mesh velocity
                        int iDof = icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( iElem, lDof  );
                        _vec2.checkAndSet( iDof ,__sum);

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
                        for ( ID idof = 0; idof < nDofElemMesh; ++idof )   // Loop on local Dof on the adjacent element
                            __sum += wLoc( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi( idof, x, y, z ); // Problem here with P2

                        // Updating interpolating vector
                        int iDof = icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( iElem, lDof );
                        _vec2.checkAndSet( iDof ,__sum);

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
                        for ( ID idof = 0; idof < nDofElemMesh; ++idof )  // Loop on local Dof on the adjacent element
                            __sum += wLoc( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi( idof, x, y, z ); // Problem here with P2

                        // Updating interpolating vector
                        int iDof = icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( iElem, lDof + 1);
                        _vec2.checkAndSet( iDof ,__sum);
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
                    __sum += wLoc( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi( idof, x, y, z );

                // Updating interpolating vector
//                std::cout << M_uFESpace->dof().localToGlobal( elemId, lDof ) << " ";
//                std::cout << icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( elemId, lDof ) << std::endl;

                int iDof = icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( elemId, lDof );
                _vec2.checkAndSet( iDof, __sum);
            }
        }
    }

}




// this will interpolate dofs values from fespace1 to fespace2

void
FSIOperator::interpolateInterfaceDofs(const FESpace<mesh_type, EpetraMap>& _fespace1,
                                      const vector_type&                   _vec1,
                                      const FESpace<mesh_type, EpetraMap>& _fespace2,
                                      vector_type&                         _vec2,
                                      dof_interface_type3D&                _dofInterface)
{

    assert(_vec1.getMaptype() == Repeated);
    assert(_vec2.getMaptype() == Unique);

    typedef mesh_type::VolumeShape GeoShape; // Element shape


    UInt nDofPerVert1  = _fespace1.refFE().nbDofPerVertex; // number of Dof per vertex
    UInt nDofPerEdge1  = _fespace1.refFE().nbDofPerEdge;   // number of Dof per edge
    UInt nDofPerFace1  = _fespace1.refFE().nbDofPerFace;   // number of Dof per face
    UInt nDofPerElem1  = _fespace1.refFE().nbDofPerVolume; // number of Dof per Volume

    UInt nDofPerVert2  = _fespace2.refFE().nbDofPerVertex; // number of Dof per vertex
    UInt nDofPerEdge2  = _fespace2.refFE().nbDofPerEdge;   // number of Dof per edge
    UInt nDofPerFace2  = _fespace2.refFE().nbDofPerFace;   // number of Dof per face
    UInt nDofPerElem2  = _fespace2.refFE().nbDofPerVolume; // number of Dof per Volume

    UInt nElemV        = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE        = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF        = GeoShape::numFaces;    // Number of element's faces

    UInt nBFacesVert   = GeoShape::GeoBShape::numVertices;
    UInt nBFacesEdge   = GeoShape::GeoBShape::numEdges;
    UInt nBFacesFace   = GeoShape::GeoBShape::numFaces;

    UInt nBEdges1      = _fespace1.mesh()->numBFaces();
    UInt nBEdges2      = _fespace2.mesh()->numBFaces();

    //    UInt nDofElem = M_uFESpace->refFE().nbDof; // Number of local dof per element of the M_uFESpace->mesh() (_mesh.getRefFE().nbDof)
    UInt numTotalDof1  = _fespace1.dof().numTotalDof();
    UInt numTotalDof2  = _fespace2.dof().numTotalDof();

    UInt nDofElemVert1 = nElemV * nDofPerVert1; // number of vertex's Dof on a Element
    UInt nDofElemEdge1 = nElemE * nDofPerEdge1; // number of edge's Dof on a Element
    UInt nDofElemFace1 = nElemF * nDofPerFace1; // number of face's Dof on a Element

    UInt nDofElemVert2 = nElemV * nDofPerVert2; // number of vertex's Dof on a Element
    UInt nDofElemEdge2 = nElemE * nDofPerEdge2; // number of edge's Dof on a Element
    UInt nDofElemFace2 = nElemF * nDofPerFace2; // number of face's Dof on a Element

    UInt numTotalVert1 = _fespace1.mesh()->numGlobalVertices();
    UInt numTotalEdge1 = _fespace1.mesh()->numGlobalEdges();
    UInt numTotalFace1 = _fespace1.mesh()->numGlobalFaces();
    UInt numTotalVol1  = _fespace1.mesh()->numGlobalVolumes();

    UInt numTotalVert2 = _fespace2.mesh()->numGlobalVertices();
    UInt numTotalEdge2 = _fespace2.mesh()->numGlobalEdges();
    UInt numTotalFace2 = _fespace2.mesh()->numGlobalFaces();
    UInt numTotalVol2  = _fespace2.mesh()->numGlobalVolumes();


    std::map<ID, ID> const& locDofMap = _dofInterface->locDofMap();
    std::map<ID, ID>::const_iterator iter;

    Real x, y, z;
    Real value;


    //    Vector wLoc( nDofElemMesh1 * nDimensions );


    // Loop on elements of the mesh
    if (nDofPerVert1 && nDofPerVert2)
        {
            //            std::cout << "  -> both FESpace have unknowns on their nodes" << std::endl;
            for ( ID iVert = 1; iVert <= _fespace1.mesh()->numVertices(); ++iVert )
                {
                    if (_fespace1.mesh()->pointList(iVert).marker() != M_fluidInterfaceFlag) continue;

                    ID nodeID = _fespace1.mesh()->pointList(iVert).id();
                    // Loop number of Dof per vertex
                    for ( ID l = 1; l <= nDofPerVert1; ++l )
                        {
                            // Loop on data vector components
                            for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                                {
                                    // "Interpolating" data at the nodal point
                                    iter = locDofMap.find(nodeID);
                                    double value = _vec1( icmp*numTotalDof1 + nodeID );
                                    // now what to what boundary node ( in the solid numerotation ) should we send this value ?
                                    //std::cout << "" << std::endl;
                                    int iDof = icmp*numTotalDof2 + iter->second;
                                    //                                    std::cout << " transfering " << value << " from P1 " << nodeID << " to P1 " << iDof  << " ( " << iter->second << " ) " << std::endl;
                                    // Updating interpolated mesh velocity
                                    _vec2.checkAndSet( iDof, value );
                                }
                        }
                }
        }

    // Edge based Dof
    if ( nDofPerEdge1 )
        {

            // loop on boundary edges
            for ( ID iEdge = 1; iEdge <= nBEdges1; ++iEdge )
                {

                    if (_fespace1.mesh()->edgeList(iEdge).marker() != M_fluidInterfaceFlag) continue;

                    // edge ID
                    ID edgeID = _fespace1.mesh()->edgeList(iEdge).id();
                    // dof position of the edge since unknowns are store using [ node | edges | faces | volumes ]
                    int iDofEdge = numTotalVert1 + edgeID;

                    if (nDofPerEdge2) // the second FE space has dofs on its edges
                        {
                            for ( ID l = 1; l <= nDofPerEdge1; ++l )
                                {
                                    // Loop on data vector components
                                    for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                                        {
                                            // ID of the dof in the solid numerotation
                                            iter = locDofMap.find(iDofEdge);
                                            int iDof = icmp*numTotalDof2 + iter->second;
                                            // "Interpolating" data at the nodal point
                                            double value = _vec1( icmp*numTotalDof1 + iDofEdge );
                                            // now what to what boundary node ( in the solid numerotation ) should we send this value ?
                                            //                                            std::cout << " transfering " << value << " from P2 " << edgeID << " to P2 " << iDof  << " ( " << iter->second << " ) " << std::endl;
                                            // Updating interpolated mesh velocity
                                            _vec2.checkAndSet( iDof, value );
                                        }
                                }
                        }
                    else
                        {
                            for ( ID l = 1; l <= nDofPerEdge1; ++l )
                                {
                                    // Loop on data vector components
                                    for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                                        {
                                            //
                                            // ID of the 1st node of the edge
                                            ID node1 = _fespace1.mesh()->edgeList(iEdge).point(1).id();
                                            iter = locDofMap.find(node1);
                                            // ID of the 1st dof of the edge in the solid numerotation
                                            int iDof1 = icmp*numTotalDof2 + iter->second;
                                            value = 0.5*_vec1( icmp*numTotalDof1 + iDofEdge  ) + _vec2[iDof1];
                                            //                                            std::cout << " transfering " << value << " from P2 " << iDofEdge << " to P1 " << iDof1 << " ( " << iter->second << " ) " << std::endl;
                                            _vec2.checkAndSet( iDof1, value );
                                            //
                                            // ID of the 2nd node of the edge
                                            ID node2 = _fespace1.mesh()->edgeList(iEdge).point(2).id();
                                            iter = locDofMap.find(node2);
                                            // ID of the 2nd dof of the edge in the solid numerotation
                                            int iDof2 = icmp*numTotalDof2 + iter->second;
                                            value = 0.5*_vec1( icmp*numTotalDof1 + iDofEdge ) + _vec2[iDof2];
                                            //                                            std::cout << " transfering " << value << " from P2 " << iDofEdge << " to P1 " << iDof2 << " ( " << iter->second << " ) " << std::endl;
                                            _vec2.checkAndSet( iDof2, value );
                                            // now what to what boundary node ( in the solid numerotation ) should we send this value ?
                                            //std::cout << "" << std::endl;
                                            // Updating interpolated mesh velocity
                                            //                                                     _vec2.checkAndSet( iDof2, value );
                                            //                                                     std::cout << std::endl;
                                        }
                                }
                        }



                }

        }
    else if (nDofPerEdge2)
        {
            // The first FESpace has no dofs on the edge
            // The second FESpace has dofs on the edge
            // We need to interpolate the vertex dofs on the edge dofs

            // First we need to go through the FS interface boundary edges on the second mesh

            // loop on boundary edges
            for ( ID iEdge = 1; iEdge <= nBEdges2; ++iEdge )
                {
                    if (_fespace2.mesh()->edgeList(iEdge).marker() != M_fluidInterfaceFlag) continue;
                    // Now that we have an edge on the FS interface, let's get its ID
                    ID edgeID = _fespace2.mesh()->edgeList(iEdge).id();
                    // dof position of the edge since unknowns are store using [ node | edges | faces | volumes ]
                    int iDofEdge = numTotalVert2 + edgeID;

                    for ( ID l = 1; l <= nDofPerEdge1; ++l )
                        {
                            // Loop on data vector components
                            for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                                {
                                    // interpolation of the nodal values in the middle of the segement
                                    ID node1 = _fespace2.mesh()->edgeList(iEdge).point(1).id();
                                    ID node2 = _fespace2.mesh()->edgeList(iEdge).point(2).id();
                                    value = 0.5*(_vec2(icmp*numTotalDof2 + node1) + _vec2(icmp*numTotalDof2 + node2));
                                    _vec2.checkAndSet(iDofEdge, value);
                                }
                        }
                }
        }


    //         // Face based Dof
//         if ( nDofpF )
//         {

//             // loop on element faces
//             for ( ID iFa = 1; iFa <= nElemF; ++iFa )
//             {

//                 // Loop on number of Dof per face
//                 for ( ID l = 1; l <= nDofpF; ++l )
//                 {

//                     lDof = nDofElemE + nDofElemV + ( iFa - 1 ) * nDofpF + l; // Local dof in the adjacent Element

//                     // Nodal coordinates
//                     x = _fespace1.refFE().xi( lDof - 1 );
//                     y = _fespace1.refFE().eta( lDof - 1 );
//                     z = _fespace1.refFE().zeta( lDof - 1 );

//                     // Loop on data vector components
//                     for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
//                     {

//                         // Interpolating data at the nodal point
//                         double __sum = 0;
//                         for ( ID idof = 0; idof < nDofElemMesh; ++idof )  // Loop on local Dof on the adjacent element
//                             __sum += wLoc( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi( idof, x, y, z ); // Problem here with P2

//                         // Updating interpolating vector
//                         int iDof = icmp * _fespace1.dof().numTotalDof() + _fespace1.dof().localToGlobal( iElem, lDof + 1);
//                         _vec2.checkAndSet( iDof ,__sum);
//                     }
//                 }
//             }
//         }

//         // Element based Dof
//         // Loop on number of Dof per Element
//         for ( ID l = 1; l <= nDofpEl; ++l )
//         {
//             lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

//             // Nodal coordinates
//             x = _fespace1.refFE().xi( lDof - 1 );
//             y = _fespace1.refFE().eta( lDof - 1 );
//             z = _fespace1.refFE().zeta( lDof - 1 );

//             // Loop on data vector components
//             for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
//             {

//                 // Interpolating data at the nodal point
//                 double __sum = 0;
//                 for ( ID idof = 0; idof < nDofElemMesh; ++idof )  // Loop on local Dof on the adjacent element
//                     __sum += wLoc( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi( idof, x, y, z );

//                 // Updating interpolating vector

//                 int iDof = icmp * _fespace1.dof().numTotalDof() + _fespace1.dof().localToGlobal( elemId, lDof );
//                 _vec2.checkAndSet( iDof, __sum);
//             }
//         }
}





void
FSIOperator::transferFluidOnInterface(const vector_type &_vec1,
                                      vector_type       &_vec2)
{
    // e.g.: vec1=M_fluid->residual(), vec2=M_sigmaFluid
//     _vec2 = ZeroVector(_vec2.size());

    if (_vec1.getMaptype() == Unique)
        {
            vector_type const  vec1Repeated(_vec1, Repeated);
            transferFluidOnInterface(vec1Repeated, _vec2);
            return;
        }

    if (_vec2.getMaptype() == Repeated)
        {
            vector_type  vec2Unique(_vec2, Unique);
            transferFluidOnInterface(_vec1, vec2Unique);
            _vec2 = vec2Unique;
            return;
        }

    _vec2 *= 0;

    std::map<ID, ID> const& locDofMap = M_dofFluidToStructure->locDofMap();


    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();
    int numTotalDofFluid = M_uFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::const_iterator Iterator;

    for (int dim = 0; dim < (int)nDimensions; ++dim)
        for ( Iterator it = locDofMap.begin(); it != locDofMap.end(); ++it )
            {
//                 std::cout <<  it->second + dim*numTotalDofFluid << " to "
//                           <<  it->first + dim*numTotalDofSolid << " : "
//                           << _vec1[it->second + dim*numTotalDofFluid] << std::endl;
                _vec2.checkAndSet( it->first + dim*numTotalDofSolid,
                                   _vec1[it->second + dim*numTotalDofFluid] );
            }
}

//works in serial but no yet in parallel
void
FSIOperator::transferSolidOnFluid(const vector_type &_vec1,//not working in parallel
                                      vector_type       &_vec2)
{
    //    e.g.: vec2=M_fluid->residual(), vec1=M_sigmaFluid
    //     _vec2 = ZeroVector(_vec2.size());

    if (_vec1.getMaptype() == Unique)
        {
            vector_type const  vec1Repeated(_vec1, Repeated);
            transferSolidOnInterface(vec1Repeated, _vec2);
            return;
        }

    if (_vec2.getMaptype() == Repeated)
        {
            vector_type  vec2Unique(_vec2, Unique);
            transferSolidOnInterface(_vec1, vec2Unique);
            _vec2 = vec2Unique;
            return;
        }

    _vec2 *= 0;

    std::map<ID, ID> const& locDofMap = M_dofFluidToStructure->locDofMap();


    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();
    int numTotalDofFluid = M_uFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::const_iterator Iterator;

    for (UInt dim = 0; dim < nDimensions; ++dim)
        for ( Iterator it = locDofMap.begin(); it != locDofMap.end(); ++it )
            {
                _vec2.checkAndSet( it->second + dim*numTotalDofFluid,
                                   _vec1[it->first + dim*numTotalDofSolid]);
            }

}

void
FSIOperator::transferSolidOnInterface(const vector_type &_vec1,
                                      vector_type       &_vec2)
{
    /* e.g.:
       vec1 (Unique)       vec2 (Repeated)  (On different EpetraMaps) (Changing now to Unique on both)
       M_solid->disp()     M_lambdaSolid
       M_solid->vel()      M_lambdaDotSolid
       M_solid->residual() M_sigmaSolid
    */

    if (_vec1.getMaptype() == Unique)
        {
            vector_type const  vec1Repeated(_vec1, Repeated);
            transferSolidOnInterface(vec1Repeated, _vec2);
            return;
        }

    if (_vec2.getMaptype() == Repeated)
        {
            vector_type  vec2Unique(_vec2, Unique);
            transferSolidOnInterface(_vec1, vec2Unique);
            _vec2 = vec2Unique;
            return;
        }

    _vec2 *= 0;

    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();

    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::const_iterator Iterator;

    for (UInt dim = 0; dim < nDimensions; ++dim)
        for ( Iterator it = locDofMap.begin(); it != locDofMap.end(); ++it )
            {
                _vec2.checkAndSet( it->second + dim*numTotalDofSolid,
                                   _vec1[it->second + dim*numTotalDofSolid] );
            }
}

void
FSIOperator::transferInterfaceOnSolid(const vector_type& _vec1,
                                      vector_type&       _vec2)
{
    /* e.g.:
       vec2                vec1
       M_solid->disp()     M_lambdaSolid
       M_solid->vel()      M_lambdaDotSolid
       M_solid->residual() M_sigmaSolid
    */

    if (_vec1.getMaptype() == Unique)
        {
            vector_type const  vec1Repeated(_vec1, Repeated);
            transferInterfaceOnSolid(vec1Repeated, _vec2);
            return;
        }

    if (_vec2.getMaptype() == Repeated)
        {
            vector_type  vec2Unique(_vec2, Unique);
            transferInterfaceOnSolid(_vec1, vec2Unique);
            _vec2 = vec2Unique;
            return;
        }

    _vec2 *= 0;

    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::const_iterator Iterator;

    for (int dim = 0; dim < (int)nDimensions; ++dim)
        for ( Iterator it = locDofMap.begin(); it != locDofMap.end(); ++it )
            {
                _vec2.checkAndSet( it->second + dim*numTotalDofSolid,
                                   _vec1[it->second + dim*numTotalDofSolid] );
            }
}

void
FSIOperator::setAlphafCoef( )
{
    Real pi=3.1459265358979;
    Real h=0.1, R=0.5;

    M_AlphafCoef  = 2*(this->dataSolid().rho()*h)/this->dataFluid().timestep();
    M_AlphafCoef += h*this->dataSolid().young(0)*this->dataFluid().timestep() /
                    (2*pow(R,2) *(1-pow(dataSolid().poisson(0),2)));
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

void FSIOperator::setStructureToFluid(vector_type const&velo,  UInt type)
{
    M_bcvStructureToFluid->setup(velo,
			       M_fluid->velFESpace().dof().numTotalDof(),
               		       M_dofHarmonicExtensionToFluid,
			       type);
}


void FSIOperator::setStructureToFluidParametres()
{
  this->setAlphafCoef();
  this->setAlphaf();
  if(M_Alphaf.get()==0)
    {
      this->setAlphafCoef();
      M_bcvStructureToFluid->setMixteCoef(M_AlphafCoef);
      M_bcvStructureToFluid->setBetaCoef(M_AlphafCoef);
    }
  else
    {
      M_bcvStructureToFluid->setMixteVec(this->Alphaf());
      M_bcvStructureToFluid->setBetaVec(this->Alphaf());
    }
}


void FSIOperator::setStructureDispToFluid(vector_type const&disp,  UInt type)
{
  M_bcvStructureDispToFluid->setup(disp,
				   M_fluid->velFESpace().dof().numTotalDof(),
				   M_dofStructureToFluid,
				   type);
}




void FSIOperator::setHarmonicExtensionVelToFluid(vector_type const& vel,
                                                 UInt type)
{
    M_bcvHarmonicExtensionVelToFluid->setup(vel,
                                            M_fluid->velFESpace().dof().numTotalDof(),
                                            M_dofHarmonicExtensionToFluid,
                                            type);
}

void FSIOperator::setDerHarmonicExtensionVelToFluid(vector_type const& dvel,
                                                    UInt type)
{
    M_bcvDerHarmonicExtensionVelToFluid->setup(dvel,
                                               M_fluid->velFESpace().dof().numTotalDof(),
                                               M_dofHarmonicExtensionToFluid,
                                               type);
}

void FSIOperator::setStructureDispToHarmonicExtension(vector_type const& disp,
                                                      UInt type)
{
    M_bcvStructureDispToHarmonicExtension->setup(disp,
                                                 M_solid->dFESpace().dof().numTotalDof(),
                                                 M_dofStructureToHarmonicExtension,
                                                 type);
}

void FSIOperator::setStructureDispToSolid(vector_type const& disp,
                                             UInt type)
{
    M_bcvStructureDispToSolid->setup(disp,
                                     M_solid->dFESpace().dof().numTotalDof(),
                                     M_dofStructureToSolid,
                                     type);
}


void FSIOperator::setDerStructureDispToSolid(vector_type const& ddisp,
                                             UInt type)
{
    M_bcvDerStructureDispToSolid->setup(ddisp,
                                        M_solid->dFESpace().dof().numTotalDof(),
                                        M_dofStructureToSolid,
                                        type);
}

void FSIOperator::setFluidLoadToStructure(vector_type const &load,
                                          UInt type)
{
    M_bcvFluidLoadToStructure->setup(load,
                                     M_solid->dFESpace().dof().numTotalDof(),
                                     M_dofStructureToSolid,
                                     type);
}

void FSIOperator::setSolidLoadToStructure(vector_type const &load,
                                          UInt type)
{
    M_bcvSolidLoadToStructure->setup(load,
				     M_solid->dFESpace().dof().numTotalDof(),
				     M_dofStructureToFluid,
				     type);
}

void FSIOperator::setDerFluidLoadToStructure(vector_type const& dload,
                                             UInt type)
{
    M_bcvDerFluidLoadToStructure->setup(dload,
                                        M_solid->dFESpace().dof().numTotalDof(),
                                        M_dofStructureToSolid,
                                        type);
}

void FSIOperator::setDerFluidLoadToFluid(vector_type const& dload,
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

void  FSIOperator::setLambdaFluid(const vector_type& lambda)
{
    if ( lambda.getMaptype() == Unique )
        *M_lambdaFluid         = lambda;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_lambdaFluidRepeated = *M_lambdaFluid;

}


void FSIOperator::setLambdaSolid(const vector_type& lambda)
{
    if ( lambda.getMaptype() == Unique )
        *M_lambdaSolid         = lambda;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_lambdaSolidRepeated = *M_lambdaSolid;
}

void FSIOperator::setLambdaSolidOld(const vector_type& lambda)
{
    if ( lambda.getMaptype() == Unique )
        *M_lambdaSolidOld     = lambda;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

}

void FSIOperator::setLambdaDotSolid(const vector_type& lambda)
{
    if ( lambda.getMaptype() == Unique )
        *M_lambdaDotSolid         = lambda;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_lambdaDotSolidRepeated     = *M_lambdaDotSolid;
}

void  FSIOperator::setSigmaFluid(const vector_type& sigma)
{
    if ( sigma.getMaptype() == Unique )
        *M_sigmaFluid         = sigma;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_sigmaFluidRepeated = *M_sigmaFluid;

}


void FSIOperator::setSigmaSolid(const vector_type& sigma)
{
    if ( sigma.getMaptype() == Unique )
        *M_sigmaSolid         = sigma;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_sigmaSolidRepeated = *M_sigmaSolid;
}

void FSIOperator::setMinusSigmaFluid(const vector_type& sigma)
{
    if ( sigma.getMaptype() == Unique )
      {
        *M_minusSigmaFluid        =sigma;
	*M_minusSigmaFluid       *=-1;
      }
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_minusSigmaFluidRepeated = *M_minusSigmaFluid;
}



void FSIOperator::resetHeAndFluid()
{
    M_meshMotion.reset(new meshmotion_raw_type(*M_mmFESpace,
                                               *M_epetraComm));
    M_fluid.reset(new fluid_raw_type(dataFluid(),
                                     *M_uFESpace,
                                     *M_pFESpace,
                                     *M_epetraComm));
    vector_type u0(M_fluid->getMap());
    M_bdf.reset(new BdfT<vector_type>(M_dataFluid->order_bdf()));
    M_bdf->initialize_unk(u0);
}

void FSIOperator::solidInit(const RefFE* refFE_struct, const LifeV::QuadRule* bdQr_struct, const LifeV::QuadRule* qR_struct)
{
    M_solidMeshPart.reset( new  partitionMesh< mesh_type > (*M_dataSolid->mesh(), *M_epetraComm));
    M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(*M_solidMeshPart,
                                                       *refFE_struct,
                                                       *qR_struct,
                                                       *bdQr_struct,
                                                       3,
                                                       *M_epetraComm));

    M_solid.reset(new solid_raw_type(dataSolid(),
                                     *M_dFESpace,
                                     *M_epetraComm));

    //                 if (M_linearSolid)
    //                     M_solidLin.reset(new FSIOperator::solidlin_raw_type(dataSolid(),
    //                                                                         *M_dFESpace,
    //                                                                         *M_epetraComm));
}

void FSIOperator::variablesInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct)
{
    // INITIALIZATION OF THE VARIABLES
    M_lambdaFluid.reset(new vector_type(*M_fluidInterfaceMap, Unique) );
    M_lambdaFluidRepeated.reset(new vector_type(*M_fluidInterfaceMap, Repeated) );

    if (this->isFluid())
        {
            M_dispFluidMeshOld.reset(new vector_type(M_uFESpace->map(), Repeated) );
            M_veloFluidMesh.reset   (new vector_type(M_uFESpace->map(), Repeated) );
            M_Alphaf.reset          (new vector_type(M_uFESpace->map(), Repeated));

            if (M_linearFluid)
                M_derVeloFluidMesh.reset(new vector_type(this->M_uFESpace->map(), Repeated) );

            M_un.reset (new vector_type(M_fluid->getMap()));
            M_rhs.reset(new vector_type(M_fluid->getMap()));
        }

    M_sigmaFluid.reset (new vector_type(*M_fluidInterfaceMap, Unique) );
    M_sigmaFluidRepeated.reset (new vector_type(*M_fluidInterfaceMap, Repeated) );
    M_minusSigmaFluid.reset (new vector_type(*M_fluidInterfaceMap, Unique) );
    M_minusSigmaFluidRepeated.reset (new vector_type(*M_fluidInterfaceMap, Repeated) );

    M_lambdaSolid.reset   (new vector_type(*M_solidInterfaceMap, Unique) );
    M_lambdaSolidOld.reset(new vector_type(*M_solidInterfaceMap, Unique) );
    M_lambdaDotSolid.reset(new vector_type(*M_solidInterfaceMap, Unique) );
    M_sigmaSolid.reset    (new vector_type(*M_solidInterfaceMap, Unique) );

    M_lambdaSolidRepeated.reset   (new vector_type(*M_solidInterfaceMap, Repeated) );
    M_lambdaDotSolidRepeated.reset(new vector_type(*M_solidInterfaceMap, Repeated) );
    M_sigmaSolidRepeated.reset    (new vector_type(*M_solidInterfaceMap, Repeated) );
}

void FSIOperator::couplingVariableExtrap(vector_ptrtype& lambda, vector_ptrtype& lambdaDot, bool& firstIter)
{
    *lambda      = lambdaSolid();
   if (firstIter)
    {
        firstIter = false;

        *lambda     += M_dataFluid->timestep()*lambdaDotSolid();
    }
    else
    {
        *lambda     += 1.5*M_dataFluid->timestep()*lambdaDotSolid(); // *1.5
        *lambda     -= M_dataFluid->timestep()*0.5*(*lambdaDot);
    }
        *lambdaDot   = lambdaDotSolid();


        leaderPrint("norm( disp  ) init = ", lambda->NormInf() );
        leaderPrint("norm( velo )  init = ", lambdaDot->NormInf());
}

}
#endif
