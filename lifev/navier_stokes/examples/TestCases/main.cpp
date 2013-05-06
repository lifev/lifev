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
    @date 06-05-2013
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

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

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/PreconditionerLinearSolver.hpp>

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

#include <lifev/navier_stokes/examples/TestCases/NavierStokesCavity.hpp>
#include <lifev/navier_stokes/examples/TestCases/NavierStokesEthierSteinman.hpp>

using namespace LifeV;

namespace
{

typedef RegionMesh<LinearTetra>           mesh_Type;
typedef MatrixEpetra<Real>                matrix_Type;
typedef VectorEpetra                      vector_Type;
typedef boost::shared_ptr<VectorEpetra>   vectorPtr_Type;
typedef MapEpetra                         map_Type;
typedef FESpace< mesh_Type, map_Type >    fespace_Type;
typedef boost::shared_ptr< fespace_Type > fespacePtr_Type;
typedef Preconditioner                    basePrec_Type;
typedef boost::shared_ptr<basePrec_Type>  basePrecPtr_Type;

typedef TimeIterationPolicyLinear< AssemblyPolicyStokes > Stokes;
typedef TimeIterationPolicyLinear< AssemblyPolicyGeneralizedStokes > GStokes;
typedef TimeIterationPolicyLinear< AssemblyPolicyNavierStokesSemiImplicit > SemiImplicit;
typedef TimeIterationPolicyNonlinearIncremental< AssemblyPolicyNavierStokesNewton > Newton;
typedef TimeIterationPolicyNonlinearIncremental< AssemblyPolicyNavierStokesPicard > Picard;
typedef TimeIterationPolicyNonlinear< AssemblyPolicyNavierStokesPicard > PicardOseen;
typedef InitPolicySolver<Stokes>          InitStokes;
typedef InitPolicySolver<GStokes>         InitGStokes;
typedef InitPolicyInterpolation           InitInter;
typedef InitPolicyProjection<SolverPolicyLinearSolver> InitProj;
typedef ExporterPolicyNoExporter          NoExporter;
typedef ExporterPolicyHDF5                HDF5Exporter;
typedef NavierStokesSolver< InitStokes, SemiImplicit, HDF5Exporter > nsSolver_Type;

void setPreconditioner ( basePrecPtr_Type& precPtr,
                         const std::string& preconditionerName,
                         const std::string& precSection,
                         boost::shared_ptr<NavierStokesProblem> /*nsProblem*/,
                         const nsSolver_Type& /*nsSolver*/,
                         const GetPot& dataFile,
                         boost::shared_ptr<Epetra_Comm> /*Comm*/,
                         const bool /*useMinusDiv*/ )
{
    if ( preconditionerName == "FromFile" )
    {
        std::string precName = dataFile ( "prec/prectype", "Ifpack" );
        precPtr.reset ( PRECFactory::instance().createObject ( precName ) );
        ASSERT ( precPtr.get() != 0, " Preconditioner not set" );
        precPtr->setDataFromGetPot ( dataFile, precSection );
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
        displayer.leaderPrint ( " |             Date: 2013-05-06                  |\n" );
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

        boost::shared_ptr<NavierStokesProblem> nsProblem;
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
        nsSolver_Type nsSolver;
        nsSolver.setProblem ( nsProblem );
        nsSolver.setup ( nsSolverList );
        basePrecPtr_Type precPtr;
        setPreconditioner ( precPtr, initPreconditionerName,
                            "initprec", nsProblem, nsSolver, dataFile,
                            Comm, useMinusDiv );
        nsSolver.setPreconditioner ( precPtr );
        nsSolver.init();

        // +-----------------------------------------------+
        // |             Solving the problem               |
        // +-----------------------------------------------+
        setPreconditioner ( precPtr, preconditionerName,
                            "prec", nsProblem, nsSolver, dataFile,
                            Comm, useMinusDiv );
        nsSolver.setPreconditioner ( precPtr );
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
