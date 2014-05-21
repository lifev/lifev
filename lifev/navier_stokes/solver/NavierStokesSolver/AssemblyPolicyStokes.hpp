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
    @file AssemblyPolicyStokes class
    @brief This class contains all the informations necessary
           to assemble a Stokes problem

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 05-12-2012
 */

#ifndef ASSEMBLYPOLICYSTOKES_HPP
#define ASSEMBLYPOLICYSTOKES_HPP

#include <iostream>
#include <string>
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
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/navier_stokes/solver/OseenAssembler.hpp>

#include <lifev/navier_stokes/solver/NavierStokesSolver/NavierStokesProblem.hpp>


namespace LifeV
{

template< class mesh_Type >
struct AssemblyPolicyStokes
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
    typedef OseenAssembler< mesh_Type, matrix_Type, vector_Type > assembler_Type;
    typedef boost::shared_ptr< assembler_Type >      assemblerPtr_Type;
    typedef Preconditioner                           preconditioner_Type;
    typedef boost::shared_ptr<preconditioner_Type>   preconditionerPtr_Type;

    enum { BDFOrder = 1 };

    void initAssembly ( Teuchos::ParameterList& list );

    void assembleSystem ( matrixPtr_Type systemMatrix,
                          vectorPtr_Type rhs,
                          vectorPtr_Type solution,
                          preconditionerPtr_Type preconditioner );

    matrixPtr_Type    M_stokesMatrix;
    assemblerPtr_Type M_assembler;

    virtual Displayer displayer() = 0;
    virtual Real currentTime() const = 0;
    virtual fespacePtr_Type uFESpace() const = 0;
    virtual fespacePtr_Type pFESpace() const = 0;
    virtual NSProblemPtr_Type problem() const = 0;
};

template< class mesh_Type >
void
AssemblyPolicyStokes< mesh_Type >::initAssembly ( Teuchos::ParameterList& list )
{
    // Loading the parameters
    std::string diffusionType = list.get ( "Diffusion type", "Viscous stress" );
    bool useMinusDiv = list.get ( "Use minus divergence", true );

    // +-----------------------------------------------+
    // |              Matrices Assembly                |
    // +-----------------------------------------------+
    displayer().leaderPrint ( "\n[Matrices Assembly]\n" );

    LifeChrono assemblyChrono;
    assemblyChrono.start();

    displayer().leaderPrint ( "Setting up assembler... " );
    M_assembler.reset ( new assembler_Type );
    M_assembler->setup ( uFESpace(), pFESpace() );
    displayer().leaderPrint ( "done\n" );

    displayer().leaderPrint ( "Defining the matrices... " );
    map_Type solutionMap ( uFESpace()->map() + pFESpace()->map() );
    M_stokesMatrix.reset ( new matrix_Type ( solutionMap ) );
    M_stokesMatrix->zero();
    displayer().leaderPrint ( "done\n" );

    // Perform the assembly of the matrix
    if ( diffusionType == "Viscous stress" )
    {
        displayer().leaderPrint ( "Adding the viscous stress... " );
        M_assembler->addViscousStress ( *M_stokesMatrix, problem()->viscosity() / problem()->density() );
        displayer().leaderPrint ( "done\n" );
    }
    else if ( diffusionType == "Stiff strain" )
    {
        displayer().leaderPrint ( "Adding the stiff strain... " );
        M_assembler->addStiffStrain ( *M_stokesMatrix, problem()->viscosity() / problem()->density() );
        displayer().leaderPrint ( "done\n" );
    }
    else
    {
        displayer().leaderPrint ( "[Error] Diffusion type unknown\n" );
        exit ( 1 );
    }

    displayer().leaderPrint ( "Adding the gradient of the pressure... " );
    M_assembler->addGradPressure ( *M_stokesMatrix );
    displayer().leaderPrint ( "done\n" );

    displayer().leaderPrint ( "Adding the divergence free constraint... " );
    if ( useMinusDiv )
    {
        M_assembler->addDivergence ( *M_stokesMatrix, -1.0 );
    }
    else
    {
        M_assembler->addDivergence ( *M_stokesMatrix, 1.0 );
    }
    displayer().leaderPrint ( "done\n" );

    displayer().leaderPrint ( "Closing the matrices... " );
    M_stokesMatrix->globalAssemble();
    displayer().leaderPrint ( "done\n" );

    assemblyChrono.stop();
    displayer().leaderPrintMax ("Matrices assembly time: ", assemblyChrono.diff(), " s.\n");
}

template< class mesh_Type >
void
AssemblyPolicyStokes< mesh_Type >::assembleSystem ( matrixPtr_Type systemMatrix,
                                                    vectorPtr_Type rhs,
                                                    vectorPtr_Type /*solution*/,
                                                    preconditionerPtr_Type /*preconditioner*/ )
{
    *systemMatrix += *M_stokesMatrix;

    M_assembler->addMassRhs ( *rhs, problem()->force(), currentTime() );
}

} // namespace LifeV

#endif /* ASSEMBLYPOLICYSTOKES_HPP */
