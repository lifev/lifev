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
    @file AssemblyPolicyNavierStokesSemiImplicit class
    @brief This class contains all the informations necessary
           to assemble a Navier-Stokes problem using a
           Newton scheme

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 06-12-2012
 */

#ifndef ASSEMBLYPOLICYNAVIERSTOKESNEWTON_HPP
#define ASSEMBLYPOLICYNAVIERSTOKESNEWTON_HPP

#include <iostream>
#include <string>
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
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/util/LifeChrono.hpp>

// #include <lifev/navier_stokes/algorithm/PreconditionerPCD.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/AssemblyPolicyStokes.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/NavierStokesProblem.hpp>


namespace LifeV
{

template< class mesh_Type >
struct AssemblyPolicyNavierStokesNewton: public AssemblyPolicyStokes< mesh_Type >
{
    typedef boost::shared_ptr< NavierStokesProblem<mesh_Type> > NSProblemPtr_Type;
    typedef MatrixEpetra<Real>                       matrix_Type;
    typedef boost::shared_ptr<matrix_Type>           matrixPtr_Type;
    typedef VectorEpetra                             vector_Type;
    typedef boost::shared_ptr<VectorEpetra>          vectorPtr_Type;
    typedef MeshPartitioner< mesh_Type >             meshPartitioner_Type;
    typedef MapEpetra                                map_Type;
    typedef boost::shared_ptr<map_Type>              mapPtr_Type;
    typedef FESpace< mesh_Type, map_Type >           fespace_Type;
    typedef boost::shared_ptr< fespace_Type >        fespacePtr_Type;
    typedef TimeAdvanceBDF<vector_Type>              bdf_Type;
    typedef boost::shared_ptr< bdf_Type >            bdfPtr_Type;
    typedef Preconditioner                           preconditioner_Type;
    typedef boost::shared_ptr<preconditioner_Type>   preconditionerPtr_Type;

    enum { BDFOrder = 1 };

    void initAssembly ( Teuchos::ParameterList& list );

    void assembleSystem ( matrixPtr_Type systemMatrix,
                          vectorPtr_Type rhs,
                          vectorPtr_Type solution,
                          preconditionerPtr_Type preconditioner );

    matrixPtr_Type       M_massMatrix;

    virtual Displayer displayer() = 0;
    virtual fespacePtr_Type uFESpace() const = 0;
    virtual fespacePtr_Type pFESpace() const = 0;
    virtual NSProblemPtr_Type problem() const = 0;
    virtual Real timestep() const = 0;
    virtual bdfPtr_Type bdf() const = 0;
};

template< class mesh_Type >
void
AssemblyPolicyNavierStokesNewton< mesh_Type >::initAssembly ( Teuchos::ParameterList& list )
{
    AssemblyPolicyStokes< mesh_Type >::initAssembly ( list );

    LifeChrono assemblyChrono;
    assemblyChrono.start();

    displayer().leaderPrint ( "Creating the mass matrix... " );
    map_Type solutionMap ( uFESpace()->map() + pFESpace()->map() );
    M_massMatrix.reset ( new matrix_Type ( solutionMap ) );
    M_massMatrix->zero();
    AssemblyPolicyStokes< mesh_Type >::M_assembler->addMass ( *M_massMatrix, 1.0 );
    M_massMatrix->globalAssemble();
    displayer().leaderPrint ( "done\n" );

    assemblyChrono.stop();
    displayer().leaderPrintMax ("Matrices assembly time: ", assemblyChrono.diff(), " s.\n");
}

template< class mesh_Type >
void
AssemblyPolicyNavierStokesNewton< mesh_Type >::assembleSystem ( matrixPtr_Type systemMatrix,
                                                                vectorPtr_Type rhs,
                                                                vectorPtr_Type solution,
                                                                preconditionerPtr_Type preconditioner )
{
    AssemblyPolicyStokes< mesh_Type >::assembleSystem ( systemMatrix, rhs, solution, preconditioner );

    bdf()->updateRHSContribution ( timestep() );
    *rhs += *M_massMatrix * bdf()->rhsContributionFirstDerivative();

    double alpha = bdf()->coefficientFirstDerivative ( 0 ) / timestep();
    *systemMatrix += *M_massMatrix * alpha;

    vector_Type beta ( systemMatrix->map(), Repeated );
    beta += *solution;
    AssemblyPolicyStokes< mesh_Type >::M_assembler->addConvection ( *systemMatrix, 1.0, beta );
    AssemblyPolicyStokes< mesh_Type >::M_assembler->addNewtonConvection ( *systemMatrix, beta );
    AssemblyPolicyStokes< mesh_Type >::M_assembler->addConvectionRhs ( *rhs, 1.0, beta );

    //    if ( preconditioner->preconditionerType() == "PCD" )
    //    {
    //        PreconditionerPCD* pcdPtr = dynamic_cast<PreconditionerPCD*> ( preconditioner.get() );
    //        pcdPtr->updateBeta ( beta );
    //    }
}

} // namespace LifeV

#endif /* ASSEMBLYPOLICYNAVIERSTOKESSEMIIMPLICIT_HPP */
