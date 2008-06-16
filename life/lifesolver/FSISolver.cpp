/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-11-18

  Copyright (C) 2004 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file FSISolver.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-11-18
 */
#include <life/lifecore/life.hpp>

#include <life/lifesolver/FSISolver.hpp>

namespace LifeV
{

Real zero_scalar( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.;
}

//! constructors

FSISolver::FSISolver( GetPot const& data_file,
                      std::string   __oper ):
//     M_dataFluid ( data_file ),
//     M_dataSolid ( data_file ),
//     M_BCh_u     (new fluid_bchandler_raw_type),
//     M_BCh_du    (new fluid_bchandler_raw_type),
//     M_BCh_d     (new fluid_bchandler_raw_type),
//     M_BCh_mesh  (new fluid_bchandler_raw_type),
    M_lambda      (),
    M_lambdaDot      (),
    M_firstIter (true),
    M_method    ( data_file("problem/method"     , "steklovPoincare") ),
    M_maxpf     ( data_file("problem/maxSubIter" , 300) ),
    M_defomega  ( data_file("problem/defOmega"   , 0.01) ),
    M_abstol    ( data_file("problem/abstol"     , 1.e-07) ),
    M_reltol    ( data_file("problem/reltol"     , 1.e-04) ),
    M_etamax    ( data_file("problem/etamax"     , 1.e-03) ),
    M_linesearch( data_file("problem/linesearch" , 0) ),
    M_epetraComm(),
    M_epetraWorldComm(),
    M_localComm (new MPI_Comm),
    M_interComm (new MPI_Comm),
    out_iter    ("iter"),
    out_res     ("res")

{
    Debug( 6220 ) << "FSISolver::FSISolver starts\n";

//     M_lambda   = ZeroVector( M_lambda.size() );
//     M_lambdaDot   = ZeroVector( M_lambdaDot.size() );

//     Debug( 6220 ) << "FSISolver::M_lambda: " << M_lambda.size() << "\n";
//     Debug( 6220 ) << "FSISolver::M_lambdaDot: " << M_lambdaDot.size() << "\n";
    int rank, numtasks;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    bool fluid = false;
    bool solid = false;

    int  fluidLeader;
    int  solidLeader;

    MPI_Group  originGroup, newGroup;
    MPI_Comm   newComm;

    MPI_Comm_group(MPI_COMM_WORLD, &originGroup);

    if (numtasks == 1)
    {
        std::cout << "Serial Fluid/Structure computation" << std::endl;
        newComm = MPI_COMM_WORLD;
        fluid = true;
        solid = true;
        fluidLeader = 0;
        solidLeader = 0;

        M_epetraWorldComm.reset(new Epetra_MpiComm(MPI_COMM_WORLD));
        M_epetraComm = M_epetraWorldComm;

    }
    else
    {
        int members[numtasks];

        solidLeader = 0;
        fluidLeader = 1-solidLeader;

        if (rank == solidLeader)
            {
                members[0] = solidLeader;
                int ierr;
                ierr = MPI_Group_incl(originGroup, 1, members, &newGroup);
                solid = true;
            }
        else
            {
                for (int ii = 0; ii <= numtasks; ++ii)
                    {
                        if ( ii < solidLeader)
                            members[ii] = ii;
                        else if ( ii > solidLeader)
                            members[ii - 1] = ii;
                    }
                int ierr;
                ierr = MPI_Group_incl(originGroup, numtasks - 1, members, &newGroup);
                fluid = true;
            }

        MPI_Comm* localComm = new MPI_Comm;
        MPI_Comm_create(MPI_COMM_WORLD, newGroup, localComm);
        M_localComm.reset(localComm);

        M_epetraComm.reset(new Epetra_MpiComm(*M_localComm.get()));
        M_epetraWorldComm.reset(new Epetra_MpiComm(MPI_COMM_WORLD));
    }

    if (fluid)
    {
        std::cout << M_epetraComm->MyPID()
                  << " ( " << rank << " ) "
                  << " out of " << M_epetraComm->NumProc()
                  << " ( " << numtasks << " ) "
                  << " is fluid." << std::endl;
        //partitionMesh< RegionMesh3D<LinearTetra> >  meshPart(M_dataFluid.mesh(), *M_epetraComm);
    }
    if (solid)
    {
        std::cout << M_epetraComm->MyPID()
                  << " ( " << rank << " ) "
                  << " out of " << M_epetraComm->NumProc()
                  << " ( " << numtasks << " ) "
                  << " is solid." << std::endl;
    }


    MPI_Barrier(MPI_COMM_WORLD);


    /*
    if (solid)
    {
        std::cout << "fluid: Building the intercommunicators ... " << std::flush;
        MPI_Comm* interComm = new MPI_Comm;
        int ierr = MPI_Intercomm_create(*M_localComm, 0, MPI_COMM_WORLD, 1, 1, interComm);
        std::cout << ierr << std::endl;
        M_interComm.reset(interComm);
    }


    if (fluid)
    {
        std::cout << "solid: Building the intercommunicators ... " << std::flush;
        MPI_Comm* interComm = new MPI_Comm;
        int ierr =  MPI_Intercomm_create(*M_localComm, 0, MPI_COMM_WORLD, 0, 1, interComm);
        std::cout << ierr << std::endl;
        M_interComm.reset(interComm);
    }

    std::cout << M_interComm.get() << std::endl;

    //M_oper->setInterComm(M_interComm.get());

    std::cout << "done." << std::endl;
    */


    MPI_Barrier(MPI_COMM_WORLD);

    Preconditioner precond  = ( Preconditioner ) data_file("problem/precond"   , DIRICHLET_NEUMANN );

    Debug( 6220 ) << "FSISolver::preconditioner: " << precond << "\n";

    if ( !__oper.empty() )
    {
        M_method = __oper;
    }

    Debug( 6220 ) << "FSISolver::setFSIOperator " << M_method << "\n";

    this->setFSIOperator( M_method );

//    M_oper = oper_fsi_ptr_mpi(new fixedPoint);

    M_oper->setFluid(fluid);
    M_oper->setSolid(solid);

    M_oper->setFluidLeader(fluidLeader);
    M_oper->setSolidLeader(solidLeader);

    Debug( 6220 ) << "FSISolver::setPreconditioner " << precond << "\n";

    std::cout << std::flush;
    M_oper->setComm(M_epetraComm, M_epetraWorldComm);

    Debug( 6220 ) << "FSISolver::setDataFromGetPot " << precond << "\n";
    std::cout << std::flush;

    M_oper->setDataFromGetPot( data_file );


    Debug( 6220 ) << "FSISolver::setPrecond " << precond << "\n";
    std::cout << std::flush;

    M_oper->setPreconditioner( precond );

    M_oper->setup();

    Debug( 6220 ) << "FSISolver:: variable setup " << precond << "\n";

    M_oper->setUpSystem(data_file);

    M_lambda.reset   (new vector_type(*M_oper->solidInterfaceMap()));
    M_lambdaDot.reset(new vector_type(*M_oper->solidInterfaceMap()));

    *M_lambda *= 0.;
    *M_lambdaDot *= 0.;



//     M_oper->
//     if (fluid)
//         {
//             M_oper->fluid().setUp(data_file);
// //            M_oper->fluidLin().setUp(data_file);
//             M_oper->meshMotion().setUp(data_file);
//         }
//     if (solid)
//         M_oper->solid().setUp(data_file);

    Debug( 6220 ) << "FSISolver:: building the fluid and solid systems " << precond << "\n";

//     if (fluid)
//     {
//         M_oper->fluid().buildSystem();
//     }

//     if (solid)
//     {
//         M_oper->solid().buildSystem();
//     }

    M_oper->buildSystem();

    MPI_Barrier(MPI_COMM_WORLD);


    Debug( 6220 ) << "FSISolver constructor ends\n";

//@     M_lambda.resize(M_oper->displacement().size());
//@     M_lambdaDot.resize(M_oper->velocity().size());
}


//


void
FSISolver::setFSIOperator( std::string const& __op )
{


    Debug( 6220 ) << "FSISolver::setFSIOperator with operator " << __op << "\n";

    M_oper = oper_fsi_ptr_mpi( FSIFactory::instance().createObject( __op ) );
//    M_oper = oper_fsi_ptr_mpi(new exactJacobian);

//    M_oper->setup();
    Debug( 6220 ) << "FSISolver::setFSIOperator done\n";
}


void
FSISolver::setFluidBC(fluid_bchandler_type bc_fluid)
{
    if (this->isFluid())
        M_oper->setFluidBC(bc_fluid);
}

void
FSISolver::setLinFluidBC(fluid_bchandler_type bc_dfluid)
{
    if (this->isFluid())
        M_oper->setLinFluidBC(bc_dfluid);
}

void
FSISolver::setInvLinFluidBC(fluid_bchandler_type bc_dfluid_inv)
{
    if (this->isFluid())
        M_oper->setInvLinFluidBC(bc_dfluid_inv);
}

void
FSISolver::setHarmonicExtensionBC(fluid_bchandler_type bc_he)
{
    if (this->isFluid())
        M_oper->setHarmonicExtensionBC(bc_he);
}

void
FSISolver::setSolidBC(solid_bchandler_type bc_solid)
{
    if (this->isSolid())
        M_oper->setSolidBC(bc_solid);
}

void
FSISolver::setLinSolidBC(solid_bchandler_type bc_dsolid)
{
    if (this->isSolid())
        M_oper->setLinSolidBC(bc_dsolid);
}

void
FSISolver::setInvLinSolidBC(solid_bchandler_type bc_dsolid_inv)
{
    if (this->isSolid())
        M_oper->setInvLinSolidBC(bc_dsolid_inv);
}

// void
// FSISolver::setReducedLinFluidBC(fluid_bchandler_type bc_dredfluid)
// {
//     M_oper->setReducedLinFluidBC(bc_dredfluid);
// }

// void
// FSISolver::setInvReducedLinFluidBC(fluid_bchandler_type bc_dredfluid_inv)
// {
//     M_oper->setInvReducedLinFluidBC(bc_dredfluid_inv);
// }


//


void
FSISolver::iterate( Real time )
{
    Debug( 6220 ) << "============================================================\n";
    Debug( 6220 ) << "Solving FSI at time " << time << " with FSIOperator: " << M_method  << "\n";
    Debug( 6220 ) << "============================================================\n";


    M_oper->setTime(time);

//     M_oper->fluid().timeAdvance( M_oper->fluid().sourceTerm(), time);
//     M_oper->solid().timeAdvance( M_oper->solid().sourceTerm(), time);

    fct_type fluidSource(zero_scalar);
    fct_type solidSource(zero_scalar);

    M_oper->updateSystem(fluidSource, solidSource);

    // displacement prediction
    MPI_Barrier(MPI_COMM_WORLD);

    if (M_firstIter)
    {
        M_firstIter = false;
        *M_lambda      = M_oper->lambdaSolid();
        *M_lambda     += timeStep()*M_oper->lambdaDotSolid();
        *M_lambdaDot   = M_oper->lambdaDotSolid();
    }
    else
    {
        *M_lambda      = M_oper->lambdaSolid();
        *M_lambda     += 1.5*timeStep()*M_oper->lambdaDotSolid(); // *1.5
        *M_lambda     -= timeStep()*0.5*(*M_lambdaDot);
        *M_lambdaDot   = M_oper->lambdaDotSolid();
    }


//     if (M_oper->isSolid())
//     {
        std::cout << "norm( disp ) init = " << M_lambda->NormInf() << std::endl;
        std::cout << "norm( velo ) init = " << M_lambdaDot->NormInf() << std::endl;
//     }

    MPI_Barrier(MPI_COMM_WORLD);

    int maxiter = M_maxpf;

    // the newton solver
    UInt status = 1;
    Debug( 6220 ) << "Calling non-linear Richardson \n";

    status = nonLinRichardson(*M_lambda,
                              *M_oper,
                              norm_inf_adaptor(),
                              M_abstol,
                              M_reltol,
                              maxiter,
                              M_etamax,
                              M_linesearch,
                              out_res,
                              time);

    if(status == 1)
    {
        std::ostringstream __ex;
        __ex << "FSISolver::iterate ( " << time << " ) Inners iterations failed to converge\n";
        throw std::logic_error( __ex.str() );
    }
    else
    {
        std::cout << "End of time "<< time << std::endl;
        std::cout << "Number of inner iterations       : "
                  << maxiter << std::endl;
        out_iter << time << " " << maxiter << " "
                 << M_oper->nbEval() << std::endl;

#warning: removed postprocessing from solver
// 	if (M_oper->isSolid()) M_oper->solid().postProcess();
//      if (M_oper->isFluid()) M_oper->fluid().postProcess();
    }

    M_oper->shiftSolution();


    Debug( 6220 ) << "FSISolver iteration at time " << time << " done\n";
    Debug( 6220 ) << "============================================================\n";
    std::cout << std::flush;
}
}
