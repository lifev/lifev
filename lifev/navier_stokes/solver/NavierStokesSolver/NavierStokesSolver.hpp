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
    @file NavierStokesSolver class
    @brief This class contains all the informations necessary to solve a Navier-Stokes problem

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 30-09-2011
 */

#ifndef NAVIERSTOKESSOLVER_HPP
#define NAVIERSTOKESSOLVER_HPP

#include <string>
#include <iostream>
#include <boost/shared_ptr.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#ifdef HAVE_MPI
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
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/util/LifeAssert.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>
#include <lifev/navier_stokes/solver/OseenAssembler.hpp>

#include <lifev/navier_stokes/solver/NavierStokesSolver/ExporterPolicyNoExporter.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/NavierStokesProblem.hpp>


namespace LifeV
{

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy = ExporterPolicyNoExporter >
class NavierStokesSolver : private InitPolicy, public virtual TimeIterationPolicy, private ExporterPolicy
{

public:

    typedef boost::shared_ptr< NavierStokesProblem > NSProblemPtr_Type;
    typedef MatrixEpetra<Real>                       matrix_Type;
    typedef boost::shared_ptr<matrix_Type>           matrixPtr_Type;
    typedef VectorEpetra                             vector_Type;
    typedef boost::shared_ptr<VectorEpetra>          vectorPtr_Type;
    typedef MapEpetra                                map_Type;
    typedef boost::shared_ptr<map_Type>              mapPtr_Type;
    // typedef RegionMesh<LinearTetra>                  mesh_Type;
    typedef boost::shared_ptr<mesh_Type>             meshPtr_Type;
    typedef FESpace< mesh_Type, map_Type >           fespace_Type;
    typedef boost::shared_ptr< fespace_Type >        fespacePtr_Type;
    typedef BCHandler                                bcContainer_Type;
    typedef boost::shared_ptr<bcContainer_Type>      bcContainerPtr_Type;
    //typedef LifeV::Preconditioner                  basePrec_Type;
    //typedef boost::shared_ptr<basePrec_Type>       basePrecPtr_Type;
    typedef Epetra_Comm                              comm_Type;
    typedef boost::shared_ptr<comm_Type>             commPtr_Type;
    typedef OseenAssembler< mesh_Type, matrix_Type, vector_Type > assembler_Type;
    typedef boost::shared_ptr< assembler_Type >      assemblerPtr_Type;
    typedef TimeAdvanceBDF<vector_Type>              bdf_Type;
    typedef boost::shared_ptr< bdf_Type >            bdfPtr_Type;

    //! @name  Constructors, destructor
    //@{
#ifdef HAVE_MPI
    NavierStokesSolver ( commPtr_Type comm = commPtr_Type ( new Epetra_MpiComm ( MPI_COMM_WORLD ) ) );
#else
    NavierStokesSolver ( commPtr_Type comm = commPtr_Type ( new Epetra_SerialComm ) );
#endif


    ~NavierStokesSolver();

    //@}

    //! @name  Methods
    //@{

    //! Prints the error of the finite element solution vector
    void printErrors();

    //! Setup all the parameters
    void setup ( Teuchos::ParameterList& list );

    //! Computes an initial solutions, or several solutions, if needed.
    void init();

    //! Solves the Navier-Stokes equations
    void solve();

    //@}

    //! @name  Set Methods
    //@{

    //! Setup the problem to be solved
    void setProblem ( NSProblemPtr_Type nsProblem );

    //@}

    //! @name  Get Methods
    //@{

    //! Returns the type of problem (e.g. Navier-Stokes)
    NSProblemPtr_Type problem() const;

    bcContainerPtr_Type bcHandler() const;

    //! Returns the FE space for the velocity
    fespacePtr_Type uFESpace() const;

    //! Returns the FE space for the pressure
    fespacePtr_Type pFESpace() const;

    //! Returns the initial time
    Real initialTime() const;

    //! Returns the end time
    Real endTime() const;

    //! Returns the timestep
    Real timestep() const;

    //! Returns the current time
    Real currentTime() const;

    //@}

private:

    // Informations and communications
    commPtr_Type            M_comm;
    Displayer               M_displayer;

    // Problem data
    NSProblemPtr_Type       M_nsProblem;
    meshPtr_Type            M_mesh;
    bcContainerPtr_Type     M_bcHandler;

    // Finite element discretization
    fespacePtr_Type         M_uFESpace;
    fespacePtr_Type         M_pFESpace;

    // Solution data
    mapPtr_Type             M_solutionMap;
    vectorPtr_Type          M_solution;

    // Keeping track of time
    bdfPtr_Type             M_bdf;
    Real                    M_initialTime;
    Real                    M_endTime;
    Real                    M_timestep;
    Real                    M_currentTime;

    // Parameters
    bool                    M_usePreviousSolutionAsGuess;

    // Getters for the policies
    Displayer displayer()
    {
        return M_displayer;
    }
    meshPtr_Type mesh() const
    {
        return M_mesh;
    }
    commPtr_Type comm()
    {
        return M_comm;
    }
    bdfPtr_Type bdf() const
    {
        return M_bdf;
    }
};


template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::NavierStokesSolver ( commPtr_Type comm ) :
    M_comm ( comm ), M_displayer ( comm )
{

}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::~NavierStokesSolver()
{

}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
void
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::printErrors()
{
    ASSERT ( M_nsProblem->hasExactSolution(), "The problem does not have an exact solution" );

    vector_Type velocity ( M_uFESpace->map(), Repeated );
    vector_Type pressure ( M_pFESpace->map(), Repeated );
    M_displayer.leaderPrint ( "[Computing errors]\n" );
    velocity.subset ( *M_solution );
    pressure.subset ( *M_solution, M_uFESpace->fieldDim() * M_uFESpace->dof().numTotalDof() );
    Real uRelativeError, pRelativeError, uL2Error, pL2Error;
    uL2Error = M_uFESpace->l2Error ( M_nsProblem->uexact(), velocity, M_currentTime, &uRelativeError );
    pL2Error = M_pFESpace->l20Error ( M_nsProblem->pexact(), pressure, M_currentTime, &pRelativeError );
    M_displayer.leaderPrint ( "Velocity\n" );
    M_displayer.leaderPrint ( "  L2 error      : ", uL2Error, "\n" );
    M_displayer.leaderPrint ( "  Relative error: ", uRelativeError, "\n" );
    M_displayer.leaderPrint ( "Pressure\n" );
    M_displayer.leaderPrint ( "  L2 error      : ", pL2Error, "\n");
    M_displayer.leaderPrint ( "  Relative error: ", pRelativeError, "\n" );
}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
void
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::setProblem ( NSProblemPtr_Type nsProblem )
{
    M_nsProblem = nsProblem;
}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
void
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::setup ( Teuchos::ParameterList& list )
{
    ASSERT ( M_nsProblem.get() != 0, "NavierStokesSolver::init : Error: You must set a Navier-Stokes problem first." );

    // FE parameters
    std::string uOrder = list.get ( "Velocity FE", "P2" );
    std::string pOrder = list.get ( "Pressure FE", "P1" );

    // Time parameters list
    Teuchos::ParameterList timeList = list.sublist ( "Time: Parameter list" );
    M_initialTime = timeList.get ( "Initial time", 0.0 );
    M_endTime = timeList.get ( "Final time", 1e-2 );
    M_timestep = timeList.get ( "Timestep", 1e-3 );

    // Solver parameters list
    Teuchos::ParameterList solverList = list.sublist ( "Solver: Parameter list" );
    M_usePreviousSolutionAsGuess = solverList.get ( "Use previous solution as guess", false );

    // +-----------------------------------------------+
    // |               Loading the mesh                |
    // +-----------------------------------------------+
    M_displayer.leaderPrint ( "\n[Loading the mesh]\n" );
    M_nsProblem->mesh ( M_mesh );

    // +-----------------------------------------------+
    // |            Creating the FE spaces             |
    // +-----------------------------------------------+
    M_displayer.leaderPrint ( "\n[Creating the FE spaces]\n" );
    LifeChrono feSpaceChrono;
    feSpaceChrono.start();

    M_displayer.leaderPrint ( "FE for the velocity: ", uOrder, "\n" );
    M_displayer.leaderPrint ( "FE for the pressure: ", pOrder, "\n" );

    M_displayer.leaderPrint ( "Building the velocity FE space ... " );
    M_uFESpace.reset ( new fespace_Type ( M_mesh, uOrder, 3, M_comm ) );
    M_displayer.leaderPrint ( "ok.\n" );

    M_displayer.leaderPrint ( "Building the pressure FE space ... " );
    M_pFESpace.reset ( new fespace_Type ( M_mesh, pOrder, 1, M_comm ) );
    M_displayer.leaderPrint ( "ok.\n" );

    // Creation of the total map
    M_solutionMap.reset ( new map_Type ( M_uFESpace->map() + M_pFESpace->map() ) );

    // Pressure offset in the vector
    UInt pressureOffset = M_uFESpace->fieldDim() * M_uFESpace->dof().numTotalDof();

    M_displayer.leaderPrint ( "Total Velocity Dof: ", pressureOffset, "\n" );
    M_displayer.leaderPrint ( "Total Pressure Dof: ", M_pFESpace->dof().numTotalDof(), "\n" );

    feSpaceChrono.stop();
    M_displayer.leaderPrintMax ("FE spaces building time: ", feSpaceChrono.diff(), " s.\n");

    // +-----------------------------------------------+
    // |             Boundary conditions               |
    // +-----------------------------------------------+
    M_displayer.leaderPrint ( "\n[Boundary conditions]\n" );

    LifeChrono bcSetupChrono;
    bcSetupChrono.start();

    M_bcHandler.reset ( new bcContainer_Type );

    M_displayer.leaderPrint ( "Setting BC... " );
    M_nsProblem->boundaryConditions ( M_bcHandler );
    M_displayer.leaderPrint ( "ok.\n" );

    // Update the BCHandler (internal data related to FE)
    M_bcHandler->bcUpdate ( *M_mesh, M_uFESpace->feBd(), M_uFESpace->dof() );

    bcSetupChrono.stop();
    M_displayer.leaderPrintMax ("Time to setup the BC: ", bcSetupChrono.diff(), " s.\n");

    // +-----------------------------------------------+
    // |                   Vectors                     |
    // +-----------------------------------------------+
    M_displayer.leaderPrint ( "Creation of vectors... " );
    M_solution.reset ( new vector_Type ( *M_solutionMap, Unique ) );
    M_displayer.leaderPrint ( "done\n" );

    M_bdf.reset ( new bdf_Type );
    M_bdf->setup ( TimeIterationPolicy::BDFOrder );

    // Exporter parameters list
    Teuchos::ParameterList exporterList = list.sublist ( "Exporter: Parameter list" );
    ExporterPolicy::initExporter ( exporterList,
                                   M_solution );

    // Time iteration parameters list
    Teuchos::ParameterList timeIterationList = list.sublist ( "Time iteration: Parameter list" );
    TimeIterationPolicy::initTimeIteration ( timeIterationList );
    // Init parameters list
    InitPolicy::setupInit ( timeIterationList ) ;

}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
void
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::init()
{
    // Find an initial solution
    InitPolicy::initSimulation ( M_bcHandler,
                                 M_solution );

    M_currentTime = M_initialTime;

    if ( M_nsProblem->hasExactSolution() )
    {
        printErrors();
    }

    // Export the initial solution
    ExporterPolicy::exportSolution ();
}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
void
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::solve()
{
    // Solving the problem
    M_displayer.leaderPrint ( "\n[Solving the problem]\n" );

    LifeChrono nsTimeLoopChrono;
    nsTimeLoopChrono.start();

    M_currentTime = M_initialTime + M_timestep;

    while ( M_currentTime <= M_endTime + M_timestep / 2. )
    {
        LifeChrono iterChrono;
        iterChrono.start();

        M_displayer.leaderPrint ( "\n[t = ", M_currentTime, " s.]\n" );

        if ( !M_usePreviousSolutionAsGuess )
        {
            *M_solution = 0;
        }

        TimeIterationPolicy::iterate ( M_solution,
                                       M_bcHandler,
                                       M_currentTime );

        // Updating the BDF scheme
        M_bdf->shiftRight ( *M_solution );

        ExporterPolicy::exportSolution ();

        iterChrono.stop();
        M_displayer.leaderPrintMax ( "Iteration time: ", iterChrono.diff(), " s.\n" );

        if ( M_nsProblem->hasExactSolution() )
        {
            printErrors();
        }

        M_currentTime += M_timestep;

#ifdef HAVE_MPI
        MPI_Barrier ( MPI_COMM_WORLD );
#endif
    }

    nsTimeLoopChrono.stop();
    M_displayer.leaderPrintMax ("Time of the temporal loop: ", nsTimeLoopChrono.diff(), " s.\n");

    // Ending the simulation
    ExporterPolicy::finalizeExporter();
}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
typename NavierStokesSolver< mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::NSProblemPtr_Type
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::problem() const
{
    return M_nsProblem;
}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
typename NavierStokesSolver< mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::bcContainerPtr_Type
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::bcHandler() const
{
    return M_bcHandler;
}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
typename NavierStokesSolver< mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::fespacePtr_Type
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::uFESpace() const
{
    return M_uFESpace;
}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
typename NavierStokesSolver< mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::fespacePtr_Type
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::pFESpace() const
{
    return M_pFESpace;
}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
Real
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::initialTime() const
{
    return M_initialTime;
}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
Real
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::endTime() const
{
    return M_endTime;
}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
Real
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::timestep() const
{
    return M_timestep;
}

template< class mesh_Type, class InitPolicy, class TimeIterationPolicy, class ExporterPolicy >
Real
NavierStokesSolver<mesh_Type, InitPolicy, TimeIterationPolicy, ExporterPolicy>::currentTime() const
{
    return M_currentTime;
}

} // namespace LifeV

#endif /* NAVIERSTOKESSOLVER_HPP */
