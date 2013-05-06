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
    @file TimeIterationPolicyBlocksSolve class
    @brief This class contains all the informations necessary
           to perform a time iteration

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 06-12-2012
 */

#ifndef TIMEITERATIONPOLICYBLOCKSSOLVE_HPP
#define TIMEITERATIONPOLICYBLOCKSSOLVE_HPP

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
#include <lifev/core/array/MatrixEpetraStructuredView.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/SolverPolicyLinearSolver.hpp>


namespace LifeV
{

template< class AssemblyPolicy, class SolverPolicy = SolverPolicyLinearSolver >
struct TimeIterationPolicyBlocksSolve : private AssemblyPolicy, public virtual SolverPolicy
{
protected:
    typedef MatrixEpetra<Real>                       matrix_Type;
    typedef MatrixEpetraStructuredView<Real>         matrixBlockView_Type;
    typedef boost::shared_ptr<matrix_Type>           matrixPtr_Type;
    typedef MatrixEpetraStructured<Real>             matrixBlock_Type;
    typedef boost::shared_ptr<matrixBlock_Type>      matrixBlockPtr_Type;
    typedef VectorEpetra                             vector_Type;
    typedef boost::shared_ptr<VectorEpetra>          vectorPtr_Type;
    typedef RegionMesh<LinearTetra>                  mesh_Type;
    typedef MeshPartitioner< mesh_Type >             meshPartitioner_Type;
    typedef MapEpetra                                map_Type;
    typedef boost::shared_ptr<map_Type>              mapPtr_Type;
    typedef FESpace< mesh_Type, map_Type >           fespace_Type;
    typedef boost::shared_ptr< fespace_Type >        fespacePtr_Type;
    typedef TimeAdvanceBDF<vector_Type>              bdf_Type;
    typedef boost::shared_ptr< bdf_Type >            bdfPtr_Type;
    typedef BCHandler                                bcContainer_Type;
    typedef boost::shared_ptr<bcContainer_Type>      bcContainerPtr_Type;

    enum { BDFOrder = AssemblyPolicy::BDFOrder };

    void initTimeIteration ( Teuchos::ParameterList& list );

    void iterate ( vectorPtr_Type solution,
                   bcContainerPtr_Type bchandler,
                   const Real& currentTime );

    // Parameters
    bool                    M_computeResidual;
    matrixPtr_Type          M_systemMatrix;
    vectorPtr_Type          M_rhs;
    mapPtr_Type             M_fullSolutionMap;

    virtual Displayer displayer() = 0;
    virtual fespacePtr_Type uFESpace() const = 0;
    virtual fespacePtr_Type pFESpace() const = 0;
    virtual bcContainerPtr_Type bcHandler() const = 0;

private:
    std::string             M_blockType;

    void assembleFluidBlock ( const matrixBlockView_Type& F );
    void assembleSchurBlock ( const matrixBlockView_Type& F,
                              const matrixBlockView_Type& B,
                              const matrixBlockView_Type& Bt );
    void assembleAugmentedFluidBlock ( const matrixBlockView_Type& F,
                                       const matrixBlockView_Type& B,
                                       const matrixBlockView_Type& Bt );
    void assembleLaplacianBlock ();
};

template< class AssemblyPolicy, class SolverPolicy >
void
TimeIterationPolicyBlocksSolve<AssemblyPolicy, SolverPolicy>::
initTimeIteration ( Teuchos::ParameterList& list )
{
    // Parameters
    M_computeResidual = list.get ( "Compute exact residual", false );
    M_blockType = list.get ( "Block type", "none" );

    // Initialization
    M_fullSolutionMap.reset ( new map_Type ( uFESpace()->map() + pFESpace()->map() ) );

    // Init the assembler
    Teuchos::ParameterList assemblyList = list.sublist ( "Assembly: Parameter list" );
    AssemblyPolicy::initAssembly ( assemblyList );

    // Init the solver
    Teuchos::ParameterList solverList = list.sublist ( "Solver: Parameter list" );
    SolverPolicy::initSolver ( solverList,
                               uFESpace()->map().commPtr() );
}

template< class AssemblyPolicy, class SolverPolicy >
void
TimeIterationPolicyBlocksSolve<AssemblyPolicy, SolverPolicy>::
iterate ( vectorPtr_Type solution,
          bcContainerPtr_Type bchandler,
          const Real& currentTime )
{
    Real rhsIterNorm ( 0.0 );

    //
    // STEP 1: Updating the system
    //
    displayer().leaderPrint ( "Updating the system... " );

    matrixPtr_Type fullSystemMatrix;
    M_rhs.reset ( new vector_Type ( *M_fullSolutionMap, Unique ) );
    *M_rhs = 0.0;
    fullSystemMatrix.reset ( new matrix_Type ( *M_fullSolutionMap ) );
    AssemblyPolicy::assembleSystem ( fullSystemMatrix, M_rhs, solution, SolverPolicy::preconditioner() );
    displayer().leaderPrint ( "done\n" );

    //
    // STEP 2: Applying the boundary conditions
    //
    displayer().leaderPrint ( "Applying BC... " );
    bcManage ( *fullSystemMatrix, *M_rhs, *uFESpace()->mesh(), uFESpace()->dof(), *bchandler, uFESpace()->feBd(), 1.0, currentTime );
    fullSystemMatrix->globalAssemble();
    displayer().leaderPrint ( "done\n" );

    std::vector<UInt> blockNumRows ( 2, 0 );
    blockNumRows[0] = uFESpace()->fieldDim() * uFESpace()->dof().numTotalDof();
    blockNumRows[1] = pFESpace()->dof().numTotalDof();
    std::vector<UInt> blockNumColumns ( blockNumRows );
    matrixBlockView_Type F, Bt, B;
    F.setup ( 0, 0, blockNumRows[0], blockNumColumns[0], fullSystemMatrix.get() );
    Bt.setup ( 0, blockNumColumns[0], blockNumRows[0], blockNumColumns[1], fullSystemMatrix.get() );
    B.setup ( blockNumRows[0], 0, blockNumRows[1], blockNumColumns[0], fullSystemMatrix.get() );

    //
    // STEP 3: Computing the system to be solved
    //
    matrixPtr_Type matrixPtr;
    if ( M_blockType == "Fluid" )
    {
        assembleFluidBlock ( F );
    }
    else if ( M_blockType == "Schur" )
    {
        assembleSchurBlock ( F, B, Bt );
    }
    else if ( M_blockType == "Augmented fluid" )
    {
        assembleSchurBlock ( F, B, Bt );
    }
    else if ( M_blockType == "Laplacian" )
    {
        assembleLaplacianBlock ();
    }
    else
    {
        ASSERT ( srcBlock.matrixPtr() != 0 , "ERROR: unknown block type in TimeIterationPolicyBlocksSolve" );
    }

    // Freeing the memory
    fullSystemMatrix.reset();

    M_systemMatrix->spy ("matrix");

    // Extra information if we want to know the exact residual
    if ( M_computeResidual )
    {
        rhsIterNorm = M_rhs->norm2();
    }

    //
    // STEP 3: Solving the system
    //
    displayer().leaderPrint ( "Solving the system... \n" );
    SolverPolicy::solve ( M_systemMatrix, M_rhs, solution );

    if ( M_computeResidual )
    {
        vector_Type Ax ( solution->map() );
        vector_Type res ( *M_rhs );
        M_systemMatrix->matrixPtr()->Apply ( solution->epetraVector(), Ax.epetraVector() );
        res.epetraVector().Update ( -1, Ax.epetraVector(), 1 );
        Real residual;
        res.norm2 ( &residual );
        residual /= rhsIterNorm;
        displayer().leaderPrint ( "Scaled residual: ", residual, "\n" );
    }
}

template< class AssemblyPolicy, class SolverPolicy >
void
TimeIterationPolicyBlocksSolve<AssemblyPolicy, SolverPolicy>::
assembleFluidBlock ( const matrixBlockView_Type& F )
{
    matrixBlockPtr_Type systemMatrixTmp ( new matrixBlock_Type ( *M_fullSolutionMap ) );

    std::vector<UInt> blockNumRows ( 2, 0 );
    blockNumRows[0] = uFESpace()->fieldDim() * uFESpace()->dof().numTotalDof();
    blockNumRows[1] = pFESpace()->dof().numTotalDof();
    std::vector<UInt> blockNumColumns ( blockNumRows );

    systemMatrixTmp->setBlockStructure ( blockNumRows, blockNumColumns );

    MatrixEpetraStructuredUtility::copyBlock ( F, * (systemMatrixTmp->block (0, 0) ) );
    MatrixEpetraStructuredUtility::createIdentityBlock ( systemMatrixTmp->block ( 1, 1 ) );
    systemMatrixTmp->globalAssemble();
    M_systemMatrix = systemMatrixTmp;
}

template< class AssemblyPolicy, class SolverPolicy >
void
TimeIterationPolicyBlocksSolve<AssemblyPolicy, SolverPolicy>::
assembleSchurBlock ( const matrixBlockView_Type& F,
                     const matrixBlockView_Type& B,
                     const matrixBlockView_Type& Bt )
{
    matrixBlockPtr_Type P1c ( new matrixBlock_Type ( *M_fullSolutionMap ) );
    matrixBlockPtr_Type BBlockMat ( new matrixBlock_Type ( *M_fullSolutionMap ) );

    std::vector<UInt> blockNumRows ( 2, 0 );
    blockNumRows[0] = uFESpace()->fieldDim() * uFESpace()->dof().numTotalDof();
    blockNumRows[1] = pFESpace()->dof().numTotalDof();
    std::vector<UInt> blockNumColumns ( blockNumRows );

    BBlockMat->setBlockStructure ( blockNumRows, blockNumColumns );
    MatrixEpetraStructuredUtility::copyBlock ( B, * (BBlockMat->block ( 1, 0 ) ) );
    BBlockMat->globalAssemble();
    boost::shared_ptr<matrixBlock_Type> invDBlockMat ( new matrixBlock_Type ( *M_fullSolutionMap ) );
    invDBlockMat->setBlockStructure ( blockNumRows, blockNumColumns );
    MatrixEpetraStructuredUtility::createInvDiagBlock ( F, * (invDBlockMat->block ( 0, 0 ) ) );
    *invDBlockMat *= -1.0;
    invDBlockMat->globalAssemble();
    boost::shared_ptr<matrixBlock_Type> tmpResultMat ( new matrixBlock_Type ( *M_fullSolutionMap ) );
    BBlockMat->multiply ( false,
                          *invDBlockMat, false,
                          *tmpResultMat, true );
    BBlockMat.reset();
    invDBlockMat.reset();
    boost::shared_ptr<matrixBlock_Type> BtBlockMat ( new matrixBlock_Type ( *M_fullSolutionMap ) );
    BtBlockMat->setBlockStructure ( blockNumRows, blockNumColumns );
    MatrixEpetraStructuredUtility::copyBlock ( Bt, * (BtBlockMat->block ( 0, 1 ) ) );
    BtBlockMat->globalAssemble();
    tmpResultMat->multiply ( false,
                             *BtBlockMat, false,
                             *P1c, false );
    BtBlockMat.reset();
    tmpResultMat.reset();

    P1c->setBlockStructure ( blockNumRows, blockNumColumns );
    MatrixEpetraStructuredUtility::createIdentityBlock ( P1c->block ( 0, 0 ) );
    P1c->globalAssemble();
    M_systemMatrix = P1c;
}

template< class AssemblyPolicy, class SolverPolicy >
void
TimeIterationPolicyBlocksSolve<AssemblyPolicy, SolverPolicy>::
assembleAugmentedFluidBlock ( const matrixBlockView_Type& F,
                              const matrixBlockView_Type& B,
                              const matrixBlockView_Type& Bt )
{
    Real alphaRDF = 1e-5;
    Real invAlpha = 1 / alphaRDF;

    std::vector<UInt> blockNumRows ( 4, 0 );
    blockNumRows[0] = uFESpace()->dof().numTotalDof();
    blockNumRows[1] = blockNumRows[0];
    blockNumRows[2] = blockNumRows[0];
    blockNumRows[3] = pFESpace()->dof().numTotalDof();
    std::vector<UInt> blockNumColumns ( blockNumRows );

    displayer().leaderPrint ( "Getting the structure of A... " );
    matrixBlockView_Type A1, B1, B1t;

    A1.setup ( 0, 0, blockNumRows[0], blockNumColumns[0], M_systemMatrix.get() );
    B1.setup ( blockNumRows[0], 0, blockNumRows[3], blockNumColumns[0], M_systemMatrix.get() );
    B1t.setup ( 0, blockNumRows[0], blockNumRows[0], blockNumColumns[3], M_systemMatrix.get() );
    displayer.leaderPrint ( "done\n" );

    displayer.leaderPrint ( "Computing BtB... " );
    boost::shared_ptr<matrixBlock_Type> systemMatrixTmp ( new matrixBlock_Type ( *M_fullSolutionMap ) );
    systemMatrixTmp->setBlockStructure ( blockNumRows, blockNumColumns );

    boost::shared_ptr<matrixBlock_Type> BBlockMat ( new matrixBlock_Type ( *M_fullSolutionMap ) );
    BBlockMat->setBlockStructure ( blockNumRows, blockNumColumns );
    MatrixEpetraStructuredUtility::copyBlock ( B1, * (BBlockMat->block ( 3, 0 ) ) );
    BBlockMat->globalAssemble();

    boost::shared_ptr<matrixBlock_Type> BtBlockMat ( new matrixBlock_Type ( *M_fullSolutionMap ) );
    BtBlockMat->setBlockStructure ( blockNumRows, blockNumColumns );
    MatrixEpetraStructuredUtility::copyBlock ( B1t, * (BtBlockMat->block ( 0, 3 ) ) );
    BtBlockMat->globalAssemble();
    BtBlockMat->multiply ( false,
                           *BBlockMat, false,
                           *systemMatrixTmp, false );
    BBlockMat.reset();
    BtBlockMat.reset();

    *systemMatrixTmp *= -invAlpha;
    displayer.leaderPrint ( "done\n" );

    MatrixEpetraStructuredUtility::copyBlock ( A1, * (systemMatrixTmp->block ( 0, 0 ) ) );
    MatrixEpetraStructuredUtility::createIdentityBlock ( systemMatrixTmp->block ( 1, 1 ) );
    MatrixEpetraStructuredUtility::createIdentityBlock ( systemMatrixTmp->block ( 2, 2 ) );
    MatrixEpetraStructuredUtility::createIdentityBlock ( systemMatrixTmp->block ( 3, 3 ) );
    systemMatrixTmp->globalAssemble();

    M_systemMatrix = systemMatrixTmp;
}

template< class AssemblyPolicy, class SolverPolicy >
void
TimeIterationPolicyBlocksSolve<AssemblyPolicy, SolverPolicy>::
assembleLaplacianBlock ()
{
    displayer().leaderPrint ( "Building the laplacian operator... " );
    matrixBlockPtr_Type systemMatrixTmp ( new matrixBlock_Type ( *M_fullSolutionMap ) );
    *systemMatrixTmp *= 0.0;

    std::vector<UInt> blockNumRows ( 2, 0 );
    blockNumRows[0] = uFESpace()->fieldDim() * uFESpace()->dof().numTotalDof();
    blockNumRows[1] = pFESpace()->dof().numTotalDof();
    std::vector<UInt> blockNumColumns ( blockNumRows );
    systemMatrixTmp->setBlockStructure ( blockNumRows, blockNumColumns );
    matrixBlockView_Type B11, B22;
    systemMatrixTmp->blockView ( 0, 0, B11 );
    systemMatrixTmp->blockView ( 1, 1, B22 );

    ADRAssembler<mesh_Type, matrixBlock_Type, vector_Type> adrPressureAssembler;
    adrPressureAssembler.setup ( pFESpace(), uFESpace() );
    adrPressureAssembler.addDiffusion ( systemMatrixTmp, 1.0, B22.firstRowIndex(), B22.firstColumnIndex() );
    MatrixEpetraStructuredUtility::createIdentityBlock ( systemMatrixTmp->block ( 0, 0 ) );
    systemMatrixTmp->globalAssemble();
    M_systemMatrix = systemMatrixTmp;

    // Loop on boundary conditions
    for ( ID i = 0; i < bcHandler()->size(); ++i )
    {
        if ( bcHandler()->operator[] ( i ).type() == Natural )
        {
            for ( ID j = 0; j < bcHandler()->operator[] ( i ).list_size(); ++j )
            {
                UInt myId = bcHandler()->operator[] ( i ) [j]->id() + B22.firstRowIndex();
                M_systemMatrix->diagonalize ( myId, 1.0 );
            }
        }
    }
    displayer().leaderPrint ( "done\n" );
}

} // namespace LifeV

#endif /* TIMEITERATIONPOLICYBLOCKSSOLVE_HPP */
