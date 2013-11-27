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
    @file InitPolicyProjection class
    @brief This class is a strategy to initialize a Navier-Stokes problem

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 11-12-2012
 */

#ifndef INITPOLICYPROJECTION_HPP
#define INITPOLICYPROJECTION_HPP

#include <iostream>
#include <boost/shared_ptr.hpp>


#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>


#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>
#include <lifev/core/fem/BCHandler.hpp>
// #include <lifev/navier_stokes/algorithm/PreconditionerPCD.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/NavierStokesProblem.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/SolverPolicyLinearSolver.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/AssemblyPolicyStokes.hpp>


namespace LifeV
{

template< class mesh_Type, class SolverPolicy = SolverPolicyLinearSolver >
struct InitPolicyProjection : public virtual SolverPolicy, public AssemblyPolicyStokes< mesh_Type >
{
    typedef VectorEpetra                             vector_Type;
    typedef boost::shared_ptr<vector_Type>           vectorPtr_Type;
    typedef MatrixEpetra<Real>                       matrix_Type;
    typedef boost::shared_ptr<matrix_Type>           matrixPtr_Type;
    typedef MeshPartitioner< mesh_Type >             meshPartitioner_Type;
    typedef MapEpetra                                map_Type;
    typedef boost::shared_ptr<map_Type>              mapPtr_Type;
    typedef FESpace< mesh_Type, map_Type >           fespace_Type;
    typedef boost::shared_ptr< fespace_Type >        fespacePtr_Type;
    typedef TimeAdvanceBDF<vector_Type>              bdf_Type;
    typedef boost::shared_ptr< bdf_Type >            bdfPtr_Type;
    typedef BCHandler                                bcContainer_Type;
    typedef boost::shared_ptr<bcContainer_Type>      bcContainerPtr_Type;
    typedef boost::shared_ptr< NavierStokesProblem<mesh_Type> > NSProblemPtr_Type;

    InitPolicyProjection() {}
    virtual ~InitPolicyProjection() {}

    void setupInit ( Teuchos::ParameterList& list );

    void initSimulation ( bcContainerPtr_Type bchandler,
                          vectorPtr_Type solution );

private:
    Teuchos::ParameterList M_list;

    virtual Displayer displayer() = 0;
    virtual fespacePtr_Type uFESpace() const = 0;
    virtual fespacePtr_Type pFESpace() const = 0;
    virtual Real timestep() const = 0;
    virtual Real initialTime() const = 0;
    virtual bdfPtr_Type bdf() const = 0;
    virtual NSProblemPtr_Type problem() const = 0;

};

template< class mesh_Type, class SolverPolicy >
void
InitPolicyProjection< mesh_Type, SolverPolicy >::
setupInit ( Teuchos::ParameterList& list )
{
    Teuchos::ParameterList assemblyList = list.sublist ( "Assembly: Parameter list" );
    M_list = list;
}

template< class mesh_Type, class SolverPolicy >
void
InitPolicyProjection< mesh_Type, SolverPolicy >::
initSimulation ( bcContainerPtr_Type bchandler,
                 vectorPtr_Type solution )
{
    ASSERT ( problem()->hasExactSolution(), "Error: You cannot use the projection method if the problem has not an analytical solution." );

    Real currentTime = initialTime() - timestep() * bdf()->order();
    UInt pressureOffset = uFESpace()->fieldDim() * uFESpace()->dof().numTotalDof();

    vectorPtr_Type rhs;
    rhs.reset ( new vector_Type ( uFESpace()->map() + pFESpace()->map(), Unique ) );

    vectorPtr_Type velocity;
    velocity.reset ( new vector_Type ( uFESpace()->map(), Unique ) );

    vectorPtr_Type pressure;
    pressure.reset ( new vector_Type ( pFESpace()->map(), Unique ) );

    // Create a first solution
    uFESpace()->interpolate ( problem()->uexact(), *velocity, currentTime );
    pFESpace()->interpolate ( problem()->pexact(), *pressure, currentTime );
    *solution = 0.0;
    solution->add ( *velocity );
    solution->add ( *pressure, pressureOffset );
    bdf()->setInitialCondition ( *solution );
    currentTime += timestep(); // we just added a first solution with interpolation

    AssemblyPolicyStokes< mesh_Type >::initAssembly ( M_list );

    displayer().leaderPrint ( "Creating the mass matrix... " );
    map_Type solutionMap ( uFESpace()->map() + pFESpace()->map() );
    matrixPtr_Type massMatrix;
    massMatrix.reset ( new matrix_Type ( solutionMap ) );
    massMatrix->zero();
    AssemblyPolicyStokes< mesh_Type >::M_assembler->addMass ( *massMatrix, 1.0 );
    massMatrix->globalAssemble();
    displayer().leaderPrint ( "done\n" );

    matrixPtr_Type systemMatrix;

    for ( ; currentTime <=  initialTime() + timestep() / 2.0; currentTime += timestep() )
    {
        *rhs      = 0.0;
        *solution = 0.0;
        *velocity = 0.0;
        *pressure = 0.0;

        uFESpace()->interpolate ( problem()->uexact(), *velocity, currentTime );
        pFESpace()->interpolate ( problem()->pexact(), *pressure, currentTime );
        solution->add ( *velocity );
        solution->add ( *pressure, pressureOffset );

        uFESpace()->interpolate ( problem()->uderexact(), *rhs, currentTime );
        rhs->globalAssemble();
        *rhs *= -1.;
        *rhs = ( *massMatrix ) * ( *rhs );

        displayer().leaderPrint ( "Updating the system... " );
        systemMatrix.reset ( new matrix_Type ( solutionMap ) );
        AssemblyPolicyStokes< mesh_Type >::assembleSystem ( systemMatrix, rhs, solution, SolverPolicy::preconditioner() );

        // We deal as in the semi-implicit way
        AssemblyPolicyStokes< mesh_Type >::M_assembler->addConvection ( *systemMatrix, 1.0, *solution );

        //        if ( SolverPolicy::preconditioner()->preconditionerType() == "PCD" )
        //        {
        //            vector_Type beta ( systemMatrix->map(), Repeated );
        //            beta += *solution;
        //            PreconditionerPCD* pcdPtr = dynamic_cast<PreconditionerPCD*> ( SolverPolicy::preconditioner().get() );
        //            pcdPtr->updateBeta ( beta );
        //        }

        displayer().leaderPrint ( "done\n" );

        LifeChrono imposingBCChrono;
        imposingBCChrono.start();
        bcManage ( *systemMatrix, *rhs, *uFESpace()->mesh(), uFESpace()->dof(), *bchandler, uFESpace()->feBd(), 1.0, currentTime );
        imposingBCChrono.stop();
        displayer().leaderPrintMax ( "Time to impose BC: ", imposingBCChrono.diff(), " s.\n" );

        systemMatrix->globalAssemble();

        displayer().leaderPrint ( "Solving the system... \n" );
        *solution = 0.0;
        SolverPolicy::solve ( systemMatrix, rhs, solution );

        // Updating bdf
        bdf()->shiftRight ( *solution );
    }
}

} // namespace LifeV

#endif /* INITPOLICYPROJECTION_HPP */
