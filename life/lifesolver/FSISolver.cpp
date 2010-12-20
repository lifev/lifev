//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/**
   @file
   @brief File containing the solver for the instances of the FSI class.
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-11-18
 */

#ifndef TWODIM

#include <life/lifecore/life.hpp>
#include <life/lifesolver/FSISolver.hpp>
//!\todo remove this header
#include <life/lifealg/nonLinRichardson.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
FSISolver::FSISolver():
        M_oper              ( ),
        M_data              ( ),
        M_fluidInterfaceMap ( ),
        M_solidInterfaceMap ( ),
        M_epetraComm        ( ),
        M_epetraWorldComm   ( ),
        M_localComm         ( new MPI_Comm ),
        M_interComm         ( new MPI_Comm )
{
#ifdef DEBUG
    Debug( 6220 ) << "FSISolver::FSISolver constructor starts\n";
#endif
}



// ===================================================
// Methods
// ===================================================
void
FSISolver::setData( const dataPtr_Type& data )
{
    M_data = data;

    int rank, numtasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    bool fluid = false;
    bool solid = false;

    int  fluidLeader(0);
    int  solidLeader(0);

    if ( ( data->method().compare("monolithicGE") && data->method().compare("monolithicGI") ) )
    {
        MPI_Group  originGroup, newGroup;
        MPI_Comm_group(MPI_COMM_WORLD, &originGroup);

        if ( numtasks == 1 )
        {
            std::cout << "Serial Fluid/Structure computation" << std::endl;
            fluid = true;
            solid = true;
            solidLeader = 0;
            fluidLeader = solidLeader;

            M_epetraWorldComm.reset( new Epetra_MpiComm(MPI_COMM_WORLD));
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
                /* int ierr = */
                MPI_Group_incl(originGroup, 1, members, &newGroup);
                solid = true;
            }
            else
            {
                for (Int ii = 0; ii <= numtasks; ++ii)
                {
                    if ( ii < solidLeader)
                        members[ii] = ii;
                    else if ( ii > solidLeader)
                        members[ii - 1] = ii;
                }

                /* int ierr = */ MPI_Group_incl(originGroup, numtasks - 1, members, &newGroup);
                fluid = true;
            }

            MPI_Comm* localComm = new MPI_Comm;
            MPI_Comm_create(MPI_COMM_WORLD, newGroup, localComm);
            M_localComm.reset(localComm);

            M_epetraComm.reset(new Epetra_MpiComm(*M_localComm.get()));
            M_epetraWorldComm.reset(new Epetra_MpiComm(MPI_COMM_WORLD));
        }
    }
    else // Monolithic or FullMonolithic
    {
        fluid = true;
        solid = true;
        solidLeader = 0;
        fluidLeader = solidLeader;

        M_epetraWorldComm.reset( new Epetra_MpiComm(MPI_COMM_WORLD));
        M_epetraComm = M_epetraWorldComm;
    }

#ifdef DEBUG
    if ( fluid )
    {
        Debug(6220) << M_epetraComm->MyPID()
        << " ( " << rank << " ) "
        << " out of " << M_epetraComm->NumProc()
        << " ( " << numtasks << " ) "
        << " is fluid." << std::endl;
    }
    if ( solid )
    {
        Debug(6220) << M_epetraComm->MyPID()
        << " ( " << rank << " ) "
        << " out of " << M_epetraComm->NumProc()
        << " ( " << numtasks << " ) "
        << " is solid." << std::endl;
    }
#endif

    M_epetraWorldComm->Barrier();

    /*
    if (solid)
    {
        std::cout << "fluid: Building the intercommunicators ... " << std::flush;
        MPI_Comm* interComm = new MPI_Comm;
        Int ierr = MPI_Intercomm_create(*M_localComm, 0, MPI_COMM_WORLD, 1, 1, interComm);
        std::cout << ierr << std::endl;
        M_interComm.reset(interComm);
    }


    if (fluid)
    {
        std::cout << "solid: Building the intercommunicators ... " << std::flush;
        MPI_Comm* interComm = new MPI_Comm;
        Int ierr =  MPI_Intercomm_create(*M_localComm, 0, MPI_COMM_WORLD, 0, 1, interComm);
        std::cout << ierr << std::endl;
        M_interComm.reset(interComm);
    }

    std::cout << M_interComm.get() << std::endl;

    //M_oper->setInterComm(M_interComm.get());

    std::cout << "done." << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    */

    this->setFSI( );

    M_oper->setFluid( fluid );
    M_oper->setSolid( solid );

    M_oper->setFluidLeader( fluidLeader );
    M_oper->setSolidLeader( solidLeader );

    M_oper->setComm( M_epetraComm, M_epetraWorldComm );

    // opening files for output on the leader only
    if (M_epetraWorldComm->MyPID() == 0)
    {
        M_out_iter.open("iter");
        M_out_res .open("res");
    }

    M_epetraWorldComm->Barrier();

#ifdef DEBUG
    Debug( 6220 ) << "FSISolver constructor ends\n";
#endif

//@     M_lambda.resize(M_oper->displacement().size());
//@     M_lambdaDot.resize(M_oper->velocity().size());

//     M_lambda   = ZeroVector( M_lambda.size() );
//     M_lambdaDot   = ZeroVector( M_lambdaDot.size() );

//     Debug( 6220 ) << "FSISolver::M_lambda: " << M_lambda.size() << "\n";
//     Debug( 6220 ) << "FSISolver::M_lambdaDot: " << M_lambdaDot.size() << "\n";

    M_oper->setData( data );
}

void
FSISolver::setup( void )
{
    M_oper->setupFluidSolid();

    M_oper->setupSystem();

    M_oper->buildSystem();
}


void
FSISolver::initialize( const std::string& /*velFName*/,
                       const std::string& /*pressName*/,
                       const std::string& /*velwName*/,
                       const std::string& /*depName*/,
                       const std::string& /*velSName*/,
                       const Real&        /*Tstart = 0.*/)
{
//             M_oper->fluid().initialize(velFName, pressName, velwName, Tstart);
//             M_oper->solid().initialize(depName, velSName, Tstart);
}


void
FSISolver::initialize(vectorPtr_Type u0, vectorPtr_Type v0)
{
    if (!u0.get())
    {
        u0.reset(new vector_Type(*M_oper->couplingVariableMap()));
        M_oper->setSolution(*u0); // couplingVariableMap()
        M_oper->initializeBDF(*u0);
    }
    else
    {
        M_oper->setSolution(*u0); // couplingVariableMap()//copy
        M_oper->initializeBDF(*u0);
    }
    if (!v0.get())
        M_oper->setSolutionDerivative(*u0); // couplingVariableMap()//copy
    //        M_oper->setSolutionDerivative(u0); // couplingVariableMap()//copy
    else
        M_oper->setSolutionDerivative(*v0);
    //M_oper->setSolutionDerivative(v0);
    //M_oper->setupBDF(*M_lambda);
}

void
FSISolver::initialize( fluidFunction_Type const& u0,
                       fluidFunction_Type const& p0,
                       solidFunction_Type const& d0,
                       solidFunction_Type const& w0,
                       fluidFunction_Type const& df0)
{
    M_oper->initialize(u0, p0 , d0, w0, df0);
}



void
FSISolver::iterate()
{
    Debug( 6220 ) << "============================================================\n";
    Debug( 6220 ) << "Solving FSI at time " << M_data->dataFluid()->dataTime()->time() << " with FSI: " << M_data->method()  << "\n";
    Debug( 6220 ) << "============================================================\n";

    // Update the system
    M_oper->updateSystem( );

    // We extract a pointer to the solution
    vectorPtr_Type lambda(new vector_Type(M_oper->solution()));
    //M_oper->solutionPtr(lambda);//copy of a shared_ptr

    // the newton solver
    UInt maxiter = M_data->maxSubIterationNumber();
    UInt status = nonLinRichardson( *lambda,
                                    *M_oper,
                                    M_data->absoluteTolerance(),
                                    M_data->relativeTolerance(),
                                    maxiter,
                                    M_data->errorTolerance(),
                                    M_data->linesearch(),
                                    M_out_res,
                                    M_data->dataFluid()->dataTime()->time() );

    // We update the solution
    M_oper->setSolution( *lambda );

    if (status == EXIT_FAILURE)
    {
        std::ostringstream __ex;
        __ex << "FSISolver::iterate ( " << M_data->dataFluid()->dataTime()->time() << " ) Inners iterations failed to converge\n";
        throw std::logic_error( __ex.str() );
    }
    else
    {
        //M_oper->displayer().leaderPrint("FSI-  Number of inner iterations:              ", maxiter, "\n" );
        if (M_epetraWorldComm->MyPID() == 0)
        {
            M_out_iter << M_data->dataFluid()->dataTime()->time() << " " << maxiter;
        }
    }

    M_oper->shiftSolution();


    Debug( 6220 ) << "FSISolver iteration at time " << M_data->dataFluid()->dataTime()->time() << " done\n";
    Debug( 6220 ) << "============================================================\n";
    std::cout << std::flush;

}

// ===================================================
// Set Functions
// ===================================================
void
FSISolver::setSourceTerms( const fluidSource_Type& fluidSource,
                           const solidSource_Type& solidSource )
{
    M_oper->fluid().setSourceTerm( fluidSource );
    M_oper->solid().setSourceTerm( solidSource );
}

void
FSISolver::setFSI( )
{
    Debug( 6220 ) << "FSISolver::setFSI with operator " << M_data->method() << "\n";
    M_oper = FSIOperPtr_Type( FSI::FSIFactory_Type::instance().createObject( M_data->method() ) );
}

void
FSISolver::setFluidBC( const fluidBchandlerPtr_Type& bc_fluid )
{
    if ( this->isFluid() )
        M_oper->setFluidBC( bc_fluid );
}

void
FSISolver::setLinFluidBC( const fluidBchandlerPtr_Type& bc_dfluid )
{
    if ( this->isFluid() )
        M_oper->setLinFluidBC( bc_dfluid );
}

void
FSISolver::setInvLinFluidBC( const fluidBchandlerPtr_Type& bc_dfluid_inv )
{
    if ( this->isFluid() )
        M_oper->setInvLinFluidBC( bc_dfluid_inv );
}

void
FSISolver::setHarmonicExtensionBC( const fluidBchandlerPtr_Type& bc_he )
{
    if ( this->isFluid() )
        M_oper->setHarmonicExtensionBC( bc_he );
}

void
FSISolver::setSolidBC( const solidBchandlerPtr_Type& bc_solid )
{
    if ( this->isSolid() )
        M_oper->setSolidBC( bc_solid );
}

void
FSISolver::setLinSolidBC( const solidBchandlerPtr_Type& bc_dsolid )
{
    if ( this->isSolid() )
        M_oper->setLinSolidBC( bc_dsolid );
}

void
FSISolver::setInvLinSolidBC( const solidBchandlerPtr_Type& bc_dsolid_inv )
{
    if ( this->isSolid() )
        M_oper->setInvLinSolidBC( bc_dsolid_inv );
}

// void
// FSISolver::setReducedLinFluidBC( const fluidBchandlerPtr_Type& bc_dredfluid )
// {
//     M_oper->setReducedLinFluidBC( bc_dredfluid );
// }

// void
// FSISolver::setInvReducedLinFluidBC( const fluidBchandlerPtr_Type& bc_dredfluid_inv )
// {
//     M_oper->setInvReducedLinFluidBC( bc_dredfluid_inv );
// }

} // Namespace LifeV
#endif /* TWODIM */
