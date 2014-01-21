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

/*!
    @file main.cpp
    @brief Application to solve different Navier-Stokes problem

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 03-12-2013
 */

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/PreconditionerLinearSolver.hpp>

// Includes the policies that I need to use
// (Here the list is quite exhaustive)
#include <lifev/navier_stokes/solver/NavierStokesSolver/NavierStokesProblem.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/NavierStokesSolver.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/AssemblyPolicyStokes.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/AssemblyPolicyGeneralizedStokes.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/AssemblyPolicyNavierStokesSemiImplicit.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/AssemblyPolicyNavierStokesNewton.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/AssemblyPolicyNavierStokesPicard.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/TimeIterationPolicyLinear.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/TimeIterationPolicyNonlinear.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/TimeIterationPolicyNonlinearIncremental.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/InitPolicyInterpolation.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/InitPolicyProjection.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/InitPolicySolver.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/ExporterPolicyHDF5.hpp>

// Includes different test cases
#include <lifev/navier_stokes/examples/TestCases/NavierStokesCavity.hpp>
#include <lifev/navier_stokes/examples/TestCases/NavierStokesEthierSteinman.hpp>

// Preconditioners for the Navier-Stokes equations
#include <lifev/navier_stokes/algorithm/PreconditionerLSC.hpp>
#include <lifev/navier_stokes/algorithm/PreconditionerPCD.hpp>
#include <lifev/navier_stokes/algorithm/PreconditionerSIMPLE.hpp>
#include <lifev/navier_stokes/algorithm/PreconditionerYosida.hpp>

using namespace LifeV;

namespace
{

typedef RegionMesh<LinearTetra>           mesh_Type;
typedef Preconditioner                    basePrec_Type;
typedef boost::shared_ptr<basePrec_Type>  basePrecPtr_Type;

typedef TimeIterationPolicyLinear< mesh_Type, AssemblyPolicyStokes< mesh_Type > > Stokes;
typedef TimeIterationPolicyLinear< mesh_Type, AssemblyPolicyGeneralizedStokes< mesh_Type > > GStokes;
typedef TimeIterationPolicyLinear< mesh_Type, AssemblyPolicyNavierStokesSemiImplicit< mesh_Type > > SemiImplicit;
typedef TimeIterationPolicyNonlinearIncremental< mesh_Type, AssemblyPolicyNavierStokesNewton< mesh_Type > > Newton;
typedef TimeIterationPolicyNonlinearIncremental< mesh_Type, AssemblyPolicyNavierStokesPicard< mesh_Type > > Picard;
typedef TimeIterationPolicyNonlinear< mesh_Type, AssemblyPolicyNavierStokesPicard< mesh_Type > > PicardOseen;
typedef InitPolicySolver< mesh_Type, Stokes > InitStokes;
typedef InitPolicySolver< mesh_Type, GStokes > InitGStokes;
typedef InitPolicyInterpolation< mesh_Type > InitInter;
typedef InitPolicyProjection<SolverPolicyLinearSolver> InitProj;
typedef ExporterPolicyNoExporter          NoExporter;
typedef ExporterPolicyHDF5< mesh_Type >   HDF5Exporter;
typedef NavierStokesSolver< mesh_Type, InitStokes, SemiImplicit, HDF5Exporter > nsSolver_Type;

// This function is just there to create a preconditioner for the problem
// Note that for now some arguments are commented.
// This is because they will be used for more advanced preconditioners
// that should be merged soon.
void setPreconditioner ( basePrecPtr_Type& precPtr,
                         const std::string& preconditionerName,
                         const std::string& precSection,
                         boost::shared_ptr<NavierStokesProblem<mesh_Type> > nsProblem,
                         const nsSolver_Type& nsSolver,
                         const GetPot& dataFile,
                         boost::shared_ptr<Epetra_Comm> Comm,
                         const bool useMinusDiv )
{
    if ( preconditionerName == "FromFile" )
    {
        std::string precName = dataFile ( "prec/prectype", "Ifpack" );
        precPtr.reset ( PRECFactory::instance().createObject ( precName ) );
        ASSERT ( precPtr.get() != 0, " Preconditioner not set" );
        precPtr->setDataFromGetPot ( dataFile, precSection );
    }
#ifdef LIFEV_HAVE_TEKO
    else if ( preconditionerName == "LSC" )
    {
        PreconditionerLSC* precLSCRawPtr ( 0 );
        precLSCRawPtr = new PreconditionerLSC;
        precLSCRawPtr->setFESpace ( nsSolver.uFESpace(), nsSolver.pFESpace() );
        precLSCRawPtr->setDataFromGetPot ( dataFile, precSection );
        precPtr.reset ( precLSCRawPtr );
    }
#endif // LIFEV_HAVE_TEKO
    else if ( preconditionerName == "PCD" )
    {
        PreconditionerPCD* precPCDRawPtr ( 0 );
        precPCDRawPtr = new PreconditionerPCD;
        precPCDRawPtr->setFESpace ( nsSolver.uFESpace(), nsSolver.pFESpace() );
        precPCDRawPtr->setBCHandler ( nsSolver.bcHandler() );
        precPCDRawPtr->setTimestep ( nsSolver.timestep() );
        precPCDRawPtr->setViscosity ( nsProblem->viscosity() );
        precPCDRawPtr->setDensity ( nsProblem->density() );
        precPCDRawPtr->setComm ( Comm );
        precPCDRawPtr->setDataFromGetPot ( dataFile, precSection );
        precPCDRawPtr->setUseMinusDivergence ( useMinusDiv );
        precPtr.reset ( precPCDRawPtr );
    }
    else if ( preconditionerName == "SIMPLE" )
    {
        PreconditionerSIMPLE* precSIMPLERawPtr ( 0 );
        precSIMPLERawPtr = new PreconditionerSIMPLE;
        precSIMPLERawPtr->setFESpace ( nsSolver.uFESpace(), nsSolver.pFESpace() );
        precSIMPLERawPtr->setDampingFactor ( 1.0 );
        precSIMPLERawPtr->setComm ( Comm );
        precSIMPLERawPtr->setDataFromGetPot ( dataFile, precSection );
        precPtr.reset ( precSIMPLERawPtr );
    }
    else if ( preconditionerName == "Yosida" )
    {
        PreconditionerYosida* precYosidaRawPtr ( 0 );
        precYosidaRawPtr = new PreconditionerYosida;
        precYosidaRawPtr->setFESpace ( nsSolver.uFESpace(), nsSolver.pFESpace() );
        precYosidaRawPtr->setTimestep ( nsSolver.timestep() );
        precYosidaRawPtr->setComm ( Comm );
        precYosidaRawPtr->setDataFromGetPot ( dataFile, precSection );
        precPtr.reset ( precYosidaRawPtr );
    }
    else
    {
        ASSERT ( false, "The preconditioner is unknown." );
    }
}

}


int
main ( int argc, char** argv )
{
    // +-----------------------------------------------+
    // |            Initialization of MPI              |
    // +-----------------------------------------------+
#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
    {
        boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_SerialComm );
#endif

        Displayer displayer ( Comm );

        displayer.leaderPrint ( " +-----------------------------------------------+\n" );
        displayer.leaderPrint ( " |            Navier-Stokes Framework            |\n" );
        displayer.leaderPrint ( " +-----------------------------------------------+\n\n" );
        displayer.leaderPrint ( " +-----------------------------------------------+\n" );
        displayer.leaderPrint ( " |           Author: Gwenol Grandperrin          |\n" );
        displayer.leaderPrint ( " |             Date: 2013-12-03                  |\n" );
        displayer.leaderPrint ( " +-----------------------------------------------+\n" );

        displayer.leaderPrint ( "\n[Initilization of MPI]\n" );
#ifdef HAVE_MPI
        displayer.leaderPrint ( "Using MPI (", Comm->NumProc(), " proc.)\n" );
#else
        displayer.leaderPrint ( "Using serial version\n" );
#endif

        LifeChrono globalChrono;
        globalChrono.start();

        // +-----------------------------------------------+
        // |               Loading the data                |
        // +-----------------------------------------------+
        displayer.leaderPrint ( "\n[Loading the data]\n" );

        // **** Stupid GetPot stuff ****
        // In principle it would be better to rely only on
        // teuchos lists but for now the preconditioners
        // still read data from GetPot.
        GetPot command_line ( argc, argv );
        const std::string datafileName = command_line.follow ( "data", 2, "-f", "--data" );
        GetPot dataFile ( datafileName );
        // *****************************
        Teuchos::ParameterList mainList = * ( Teuchos::getParametersFromXmlFile ( datafileName + ".xml" ) );

        std::string initPreconditionerName = mainList.get ( "Preconditioner for init", "none" );
        std::string preconditionerName     = mainList.get ( "Preconditioner", "none" );

        // +-----------------------------------------------+
        // |             Setup the test case               |
        // +-----------------------------------------------+
        displayer.leaderPrint ( "\n[Setup the test case]\n" );

        // Problem parameters list
        Teuchos::ParameterList problemList    = mainList.sublist ( "Navier-Stokes problem: Parameter list" );
        const std::string benchmark           = problemList.get ( "Test case name", "none" );
        const Real viscosity                  = problemList.get ( "Viscosity", 0.035 );
        const Real density                    = problemList.get ( "Density", 1.0 );
        const UInt meshRefinement             = problemList.get ( "Mesh refinement", 2 );
        std::string meshPath                  = problemList.get ( "Resources path", "./Resources" );
        meshPath.append ("/");

        // Here we create a Navier-Stokes problem object
        boost::shared_ptr< NavierStokesProblem< mesh_Type > > nsProblem;
        if ( benchmark == "Cavity" )
        {
            nsProblem.reset ( new NavierStokesCavity );
        }
        else if ( benchmark == "Ethier-Steinman" )
        {
            nsProblem.reset ( new NavierStokesEthierSteinman );
        }
        else
        {
            displayer.leaderPrint ( "[Error] This benchmark does not exist\n" );
            exit ( 1 );
        }
        nsProblem->setMesh ( meshRefinement, meshPath );
        nsProblem->setViscosity ( viscosity );
        nsProblem->setDensity ( density );

        displayer.leaderPrint ( "Test case           : ", nsProblem->name(), "\n" );

        // Mesh discretization
        displayer.leaderPrint ( "Mesh refinement     : ", meshRefinement, "\n" );
        displayer.leaderPrint ( "Mesh path           : ", meshPath, "\n" );

        // Physical quantity
        displayer.leaderPrint ( "Viscosity           : ", viscosity, "\n" );
        displayer.leaderPrint ( "Density             : ", density, "\n" );

        // +-----------------------------------------------+
        // |                Initialization                 |
        // +-----------------------------------------------+
        Teuchos::ParameterList nsSolverList   = mainList.sublist ( "Navier-Stokes solver: Parameter list" );
        bool useMinusDiv = nsSolverList.sublist ( "Time iteration: Parameter list" )
                           .sublist ( "Assembly: Parameter list" )
                           .get ( "Use minus divergence" , true );

        // We create the solver, set the problem, set the preconditioner
        // and set the parameters.
        nsSolver_Type nsSolver;
        nsSolver.setProblem ( nsProblem );
        nsSolver.setup ( nsSolverList );
        basePrecPtr_Type precPtr;
        setPreconditioner ( precPtr, initPreconditionerName,
                            "initprec", nsProblem, nsSolver, dataFile,
                            Comm, useMinusDiv );
        nsSolver.setPreconditioner ( precPtr );

        // We compute the initial condition
        nsSolver.init();

        // +-----------------------------------------------+
        // |             Solving the problem               |
        // +-----------------------------------------------+
        // We set the preconditioner
        setPreconditioner ( precPtr, preconditionerName,
                            "prec", nsProblem, nsSolver, dataFile,
                            Comm, useMinusDiv );
        nsSolver.setPreconditioner ( precPtr );

        // We solve the Navier-Stokes equations timestep after timestep
        nsSolver.solve();

        globalChrono.stop();
        displayer.leaderPrintMax ( "\nTotal simulation time: ", globalChrono.diff(), " s.\n" );
        displayer.leaderPrint ( "\n[[END_SIMULATION]]\n" );

#ifdef HAVE_MPI
    }
    MPI_Finalize();
#endif
    return ( EXIT_SUCCESS );
}
