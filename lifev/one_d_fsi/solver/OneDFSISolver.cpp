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
 *  @file
 *  @brief File containing a solver class for the 1D model.
 *
 *  @version 1.0
 *  @date 01-10-2006
 *  @author Vincent Martin
 *  @author Tiziano Passerini
 *  @author Lucia Mirabella
 *
 *  @version 2.0
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *  @date 01-08-2009
 *
 *  @version 2.1
 *  @date 21-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributors Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/one_d_fsi/solver/OneDFSISolver.hpp>

namespace LifeV
{

std::map< std::string, OneDFSI::physicsType_Type > OneDFSI::physicsMap;
std::map< std::string, OneDFSI::fluxTerm_Type >    OneDFSI::fluxMap;
std::map< std::string, OneDFSI::sourceTerm_Type >  OneDFSI::sourceMap;

// ===================================================
// Constructors & Destructor
// ===================================================
OneDFSISolver::OneDFSISolver():
    M_physicsPtr                   (),
    M_fluxPtr                      (),
    M_sourcePtr                    (),
    M_feSpacePtr                   (),
    M_commPtr                      (),
    M_displayer                    (),
    M_elementalMassMatrixPtr       (),
    M_elementalStiffnessMatrixPtr  (),
    M_elementalGradientMatrixPtr   (),
    M_elementalDivergenceMatrixPtr (),
    M_rhs                          (),
    M_residual                     (),
    M_fluxVector                   (),
    M_sourceVector                 (),
    M_dFdUVector                   (),
    M_dSdUVector                   (),
    M_homogeneousMassMatrixPtr     (),
    M_homogeneousGradientMatrixPtr (),
    M_dSdUMassMatrixPtr            (),
    M_dFdUStiffnessMatrixPtr       (),
    M_dFdUGradientMatrixPtr        (),
    M_dSdUDivergenceMatrixPtr      (),
    M_linearSolverPtr              (),
    M_linearViscoelasticSolverPtr  ()
{
}

// ===================================================
// Methods
// ===================================================
void
OneDFSISolver::buildConstantMatrices()
{
    std::fill( M_dFdUVector.begin(), M_dFdUVector.end(), ublas::zero_vector<Real>( M_physicsPtr->data()->numberOfNodes()) );
    std::fill( M_dSdUVector.begin(), M_dSdUVector.end(), ublas::zero_vector<Real>( M_physicsPtr->data()->numberOfNodes()) );

    for ( UInt i(0); i < 4; ++i )
    {
        M_dSdUMassMatrixPtr[i].reset( new matrix_Type( M_feSpacePtr->map() ) );
        M_dFdUStiffnessMatrixPtr[i].reset( new matrix_Type( M_feSpacePtr->map() ) );
        M_dFdUGradientMatrixPtr[i].reset( new matrix_Type( M_feSpacePtr->map() ) );
        M_dSdUDivergenceMatrixPtr[i].reset( new matrix_Type( M_feSpacePtr->map() ) );
    }

    // Elementary computation and matrix assembling
    for ( UInt iElement(0); iElement < M_physicsPtr->data()->numberOfElements(); ++iElement )
    {
        // set the elementary matrices to 0.
        M_elementalMassMatrixPtr->zero();
        M_elementalGradientMatrixPtr->zero();

        // update the current element
        M_feSpacePtr->fe().update( M_feSpacePtr->mesh()->edgeList( iElement ), UPDATE_DPHI | UPDATE_WDET );

        // update the mass and grad matrices
        mass( 1, *M_elementalMassMatrixPtr, M_feSpacePtr->fe(), 0, 0 );
        grad( 0 , -1, *M_elementalGradientMatrixPtr, M_feSpacePtr->fe(), M_feSpacePtr->fe(), 0, 0 );

        // assemble the mass and grad matrices
        assembleMatrix( *M_homogeneousMassMatrixPtr, *M_elementalMassMatrixPtr, M_feSpacePtr->fe(), M_feSpacePtr->dof() , 0, 0, 0, 0 );
        assembleMatrix( *M_homogeneousGradientMatrixPtr, *M_elementalGradientMatrixPtr, M_feSpacePtr->fe(), M_feSpacePtr->dof() , 0, 0, 0, 0  );
    }

    // Dirichlet boundary conditions set in the mass matrix
    M_homogeneousMassMatrixPtr->globalAssemble();
    M_homogeneousGradientMatrixPtr->globalAssemble();

    // In the classical case the linear system use a mass matrix (with Dirichlet BC)
    //matrixPtr_Type systemMatrix( new matrix_Type( M_feSpacePtr->map() ) );
    //systemMatrix->insertValueDiagonal( M_physicsPtr->data()->mesh()->meanH() );
    matrix_Type systemMatrix( *M_homogeneousMassMatrixPtr );
    applyDirichletBCToMatrix( systemMatrix );
    M_linearSolverPtr->setMatrix( systemMatrix );
}

void
OneDFSISolver::setupSolution( solution_Type& solution, const MapEpetra& map, const bool& onlyMainQuantities )
{
    solution["Q"].reset( new vector_Type( map ) );
    solution["P"].reset( new vector_Type( map ) );
    solution["AoverA0minus1"].reset( new vector_Type( map ) );

    if( onlyMainQuantities )
        return;

    solution["A"].reset( new vector_Type( map ) );
    solution["W1"].reset( new vector_Type( map ) );
    solution["W2"].reset( new vector_Type( map ) );

    // Flux correction with viscoelastic term
    if ( M_physicsPtr->data()->viscoelasticWall() )
    {
        solution["Q_visc"].reset( new vector_Type( map ) ); // viscoelastic contribution to the flux
        solution["P_visc"].reset( new vector_Type( map ) ); // viscoelastic contribution to the pressure
    }

    // correction flux with inertial term
    if ( M_physicsPtr->data()->inertialWall() )
        solution["Q_inert"].reset( new vector_Type( map ) );

    // correction flux with longitudinal term
    if ( M_physicsPtr->data()->longitudinalWall() )
        solution["Q_long"].reset( new vector_Type( map ) );

    // Initialize solution to zero
    for ( solutionConstIterator_Type i = solution.begin(); i != solution.end(); ++i )
        *i->second = 0.;
}

void
OneDFSISolver::initialize( solution_Type& solution )
{
    for ( UInt iNode(0); iNode < M_physicsPtr->data()->numberOfNodes() ; ++iNode )
    {
        (*solution["A"])[iNode] = M_physicsPtr->data()->area0( iNode );
        (*solution["Q"])[iNode] = 0;
    }

    // Compute W1 and W2 from A and Q
    computeW1W2( solution );

    // Compute A/A0 from A
    computeAreaRatio( solution );

    // Compute initial pressure (taking into account the viscoelastic wall)
    M_physicsPtr->setArea_tn( *solution["A"] );
    computePressure( solution, M_physicsPtr->data()->dataTime()->timeStep() );
}

void
OneDFSISolver::computeW1W2( solution_Type& solution )
{
    for ( UInt iNode(0); iNode < M_physicsPtr->data()->numberOfNodes() ; ++iNode )
    {
        M_physicsPtr->fromUToW( ( *solution["W1"]) [iNode], (*solution["W2"]) [iNode],
                                ( *solution["A"] ) [iNode], (*solution["Q"] ) [iNode], iNode );
    }
}

void
OneDFSISolver::computePressure( solution_Type& solution, const Real& timeStep )
{
    for ( UInt iNode(0); iNode < M_physicsPtr->data()->numberOfNodes() ; ++iNode )
    {
        ( *solution["P"] ) [iNode] = M_physicsPtr->elasticPressure( ( *solution["A"] ) [iNode], iNode )
                                   + M_physicsPtr->externalPressure();

        if ( M_physicsPtr->data()->viscoelasticWall() )
        {
            ( *solution["P_visc"] ) [iNode] = M_physicsPtr->viscoelasticPressure( ( *solution["A"]) [iNode], timeStep, iNode );
            ( *solution["P"] )      [iNode] += ( *solution["P_visc"] ) [iNode];
        }
    }
}

void
OneDFSISolver::computeAreaRatio( solution_Type& solution )
{
    for ( UInt iNode(0); iNode < M_physicsPtr->data()->numberOfNodes() ; ++iNode )
        ( *solution["AoverA0minus1"] ) [iNode] = ( *solution["A"] ) [iNode] / M_physicsPtr->data()->area0( iNode ) - 1;
}

void
OneDFSISolver::computeArea( solution_Type& solution )
{
    for ( UInt iNode(0); iNode < M_physicsPtr->data()->numberOfNodes() ; ++iNode )
        ( *solution["A"] ) [iNode] = ( (*solution["AoverA0minus1"] ) [iNode] + 1 ) * M_physicsPtr->data()->area0( iNode );
}

void
OneDFSISolver::updateRHS( const solution_Type& solution, const Real& timeStep )
{
    updatedFdU( solution ); // Update the vector containing the values of the flux at the nodes and its jacobian
    updatedSdU( solution ); // Update the vector containing the values of the source term at the nodes and its jacobian
    updateMatrices();       // Update the matrices for the non-linear terms

    // Taylor-Galerkin scheme: (explicit, U = [U1,U2]^T, with U1=A, U2=Q )
    // (Un+1, phi) =          (               Un,     phi     )-> massFactor^{-1} * Un+1 = mass * U
    //             + dt     * (           Fh(Un),     dphi/dz )->            grad * F(U)
    //             - dt^2/2 * (diffFh(Un) Sh(Un),     dphi/dz )-> gradDiffFlux(U) * S(U)
    //             + dt^2/2 * (diffSh(Un) dFh/dz(Un), phi     )->   divDiffSrc(U) * F(U)
    //             - dt^2/2 * (diffFh(Un) dFh/dz(Un), dphi/dz )->stiffDiffFlux(U) * F(U)
    //             - dt     * (           Sh(Un),     phi     )->            mass * S(U)
    //             + dt^2/2 * (diffSh(Un) Sh(Un),     phi     )->  massDiffSrc(U) * S(U)

    Real dt2over2 = timeStep * timeStep * 0.5;

    // Initialize residual to 0
    *M_residual[0] *= 0;
    *M_residual[1] *= 0;

    for ( UInt i(0); i < 2; ++i )
    {
        // rhs = rhs + dt * grad * F(Un)
        *M_residual[i] += ( *M_homogeneousGradientMatrixPtr ) * ( timeStep * *M_fluxVector[i] );

        // rhs = rhs - dt * mass * S(Un)
        *M_residual[i] += ( *M_homogeneousMassMatrixPtr ) * ( - timeStep * *M_sourceVector[i] );

        for ( UInt j(0); j < 2; ++j )
        {
            // rhs = rhs - dt^2/2 * gradDiffFlux * S(Un)
            M_dFdUGradientMatrixPtr[2*i + j]->globalAssemble();
            *M_residual[i] += *M_dFdUGradientMatrixPtr[2*i + j] * ( -dt2over2 * *M_sourceVector[j] );

            // rhs = rhs + dt^2/2 * divDiffSrc * F(Un)
            M_dSdUDivergenceMatrixPtr[2*i + j]->globalAssemble();
            *M_residual[i] += *M_dSdUDivergenceMatrixPtr[2*i + j] * ( dt2over2 * *M_fluxVector[j] );

            // rhs = rhs - dt^2/2 * stiffDiffFlux * F(Un)
            M_dFdUStiffnessMatrixPtr[2*i + j]->globalAssemble();
            *M_residual[i] += *M_dFdUStiffnessMatrixPtr[2*i + j]*( -dt2over2 * *M_fluxVector[j] );

            // rhs = rhs + dt^2/2 * massDiffSrc * S(Un)
            M_dSdUMassMatrixPtr[2*i + j]->globalAssemble();
            *M_residual[i] += *M_dSdUMassMatrixPtr[2*i + j] * ( dt2over2 * *M_sourceVector[j] );
        }
    }

    // rhs = mass * Un + residual
    *M_rhs[0]  = *M_residual[0] + ( *M_homogeneousMassMatrixPtr ) * *solution.find("A")->second;
    *M_rhs[1]  = *M_residual[1] + ( *M_homogeneousMassMatrixPtr ) * *solution.find("Q")->second;

    if ( M_physicsPtr->data()->viscoelasticWall() )
        *M_rhs[1]  -= ( *M_homogeneousMassMatrixPtr ) * *solution.find("Q_visc")->second;
}

void
OneDFSISolver::iterate( OneDFSIBCHandler& bcHandler, solution_Type& solution, const Real& time, const Real& timeStep )
{
    // Apply BC to RHS
    bcHandler.applyBC( time, timeStep, solution, M_fluxPtr, M_rhs );

    // Compute A^n+1
    vector_Type area( *M_rhs[0] );
    M_linearSolverPtr->solveSystem( *M_rhs[0], area, M_homogeneousMassMatrixPtr );

    // Compute Q^n+1
    vector_Type flowRate( *M_rhs[1] );
    M_linearSolverPtr->solveSystem( *M_rhs[1], flowRate, M_homogeneousMassMatrixPtr );

    // Correct flux with inertial, viscoelastic and longitudinal terms
    if ( M_physicsPtr->data()->inertialWall() )
    {
        *solution["Q_inert"] = inertialFlowRateCorrection( flowRate );
        flowRate += *solution["Q_inert"];
    }

    if ( M_physicsPtr->data()->viscoelasticWall() )
    {
        *solution["Q_visc"] = viscoelasticFlowRateCorrection( area, flowRate, *solution.find("Q_visc")->second, timeStep, bcHandler );
        flowRate += *solution["Q_visc"];
    }

    if ( M_physicsPtr->data()->longitudinalWall() )
    {
        *solution["Q_long"] = longitudinalFlowRateCorrection();
        flowRate += *solution["Q_long"];
    }

    // Update the solution container
    *solution["A"] = area;
    *solution["Q"] = flowRate;

    // Compute W1 and W2 from A and Q
    computeW1W2( solution );

    // Compute A/A0 from A
    computeAreaRatio( solution );

    // Update the pressure (taking into account the viscoelastic wall)
    computePressure( solution, timeStep );
}

OneDFSISolver::vector_Type
OneDFSISolver::viscoelasticFlowRateCorrection( const vector_Type& newArea, const vector_Type& newElasticFlowRate,
                                               const vector_Type& oldViscoelasticFlowRate,
                                               const Real& timeStep, OneDFSIBCHandler& bcHandler,
                                               const bool& updateSystemMatrix )
{
    // Matrix
    matrix_Type massMatrix( M_feSpacePtr->map() );
    matrix_Type stiffnessMatrix( M_feSpacePtr->map() );

    Real massCoefficient;
    Real stiffnessCoefficient;

    // RHS
    vector_Type rhs( M_feSpacePtr->map() );
    rhs = 0;

    // Elementary computation and matrix assembling
    for ( UInt iElement(0); iElement < M_physicsPtr->data()->numberOfElements() ; ++iElement )
    {
        // Update the current element
        M_feSpacePtr->fe().update( M_feSpacePtr->mesh()->edgeList(iElement), UPDATE_DPHI | UPDATE_WDET );

        // Compute mass coefficient
        massCoefficient = 1 / ( 0.5 * ( newArea[ iElement ] + newArea[ iElement + 1 ] ) );
//      massCoefficient = 1 / ( 0.5 * ( M_physicsPtr->data()->area0( iElement ) + M_physicsPtr->data()->area0( iElement + 1 ) ) ); // For debug purposes

        // Compute stiffness coefficient
        stiffnessCoefficient  = timeStep * 0.5 * ( M_physicsPtr->data()->viscoelasticCoefficient( iElement ) + M_physicsPtr->data()->viscoelasticCoefficient( iElement + 1 ) )
                              / M_physicsPtr->data()->densityRho() * massCoefficient * std::sqrt( massCoefficient );

        // Set the elementary matrices to 0.
        M_elementalMassMatrixPtr->zero();
        M_elementalStiffnessMatrixPtr->zero();

        // Assemble the elemental matrix
        mass(  massCoefficient,      *M_elementalMassMatrixPtr,      M_feSpacePtr->fe(), 0, 0 );
        stiff( stiffnessCoefficient, *M_elementalStiffnessMatrixPtr, M_feSpacePtr->fe(), 0, 0 );

        // Assemble the stiffness matrix
        assembleMatrix( massMatrix,      *M_elementalMassMatrixPtr,      M_feSpacePtr->fe(), M_feSpacePtr->dof(), 0, 0, 0, 0 );
        assembleMatrix( stiffnessMatrix, *M_elementalStiffnessMatrixPtr, M_feSpacePtr->fe(), M_feSpacePtr->dof(), 0, 0, 0, 0 );

        // Natural BC
        if ( iElement == 0 )
            rhs( 0 )            -= stiffnessCoefficient * M_fluxPtr->physics()->data()->computeSpatialDerivativeAtNode( newElasticFlowRate, 0, 1 );
        if ( iElement == M_physicsPtr->data()->numberOfElements() - 1 )
            rhs( iElement + 1 ) += stiffnessCoefficient * M_fluxPtr->physics()->data()->computeSpatialDerivativeAtNode( newElasticFlowRate, iElement + 1, 1 );
    }

    // Matrices global assemble
    massMatrix.globalAssemble();
    stiffnessMatrix.globalAssemble();

    // RHS
    rhs += massMatrix * (oldViscoelasticFlowRate) + stiffnessMatrix * (-newElasticFlowRate);

    // System matrix
    massMatrix += stiffnessMatrix;

    // Apply BC to Matrix and RHS
    bcHandler.applyViscoelasticBC( M_fluxPtr, massMatrix, rhs );

    // Compute flow rate correction at t^n+1
    vector_Type newViscoelasticFlowRate( rhs );

    if ( updateSystemMatrix )
        M_linearViscoelasticSolverPtr->setMatrix( massMatrix );
    M_linearViscoelasticSolverPtr->solveSystem( rhs, newViscoelasticFlowRate, M_homogeneousMassMatrixPtr );

    return newViscoelasticFlowRate;
}

Real
OneDFSISolver::computeCFL( const solution_Type& solution, const Real& timeStep ) const
{
    Real lambdaMax( 0. );

    container2D_Type eigenvalues;
    container2D_Type leftEigenvector1;
    container2D_Type leftEigenvector2;

    for ( UInt iNode(0); iNode < M_physicsPtr->data()->numberOfNodes() ; ++iNode )
    {
        // compute the eigenvalues at node
        M_fluxPtr->eigenValuesEigenVectors( ( *solution.find("A")->second ) ( iNode ),
                                         ( *solution.find("Q")->second ) ( iNode ),
                                           eigenvalues, leftEigenvector1, leftEigenvector2, iNode );

        lambdaMax = std::max<Real>( std::max<Real>( std::fabs(eigenvalues[0]), std::fabs(eigenvalues[1]) ), lambdaMax );
    }

    Real minH = MeshUtility::MeshStatistics::computeSize(*M_feSpacePtr->mesh()).minH;
    return lambdaMax * timeStep / minH;
}

void
OneDFSISolver::resetOutput( const solution_Type& solution )
{
    std::ofstream outfile;
    for ( solutionConstIterator_Type i = solution.begin(); i != solution.end(); ++i )
    {
        std::string file = M_physicsPtr->data()->postprocessingDirectory() + "/" + M_physicsPtr->data()->postprocessingFile() + "_" + i->first + ".mfile";
        outfile.open( file.c_str(), std::ios::trunc );
        outfile.close();
    }
}

void
OneDFSISolver::postProcess( const solution_Type& solution, const Real& time )
{
    std::ofstream outfile;
    for ( solutionConstIterator_Type i = solution.begin(); i != solution.end(); ++i )
    {
        std::string file = M_physicsPtr->data()->postprocessingDirectory() + "/" + M_physicsPtr->data()->postprocessingFile() + "_" + i->first + ".mfile";
        outfile.open( file.c_str(), std::ios::app );
        outfile.setf( ios::scientific, ios::floatfield );

        outfile << time << " ";
        for ( UInt iNode(0); iNode < static_cast< UInt > ( (*i->second).size() ); ++iNode )
            outfile << (*i->second)(iNode) << " ";

        outfile << std::endl;
        outfile.close();
    }
}

// ===================================================
// Set Methods
// ===================================================
void
OneDFSISolver::setProblem( const physicsPtr_Type& physicsPtr,
                           const fluxPtr_Type&    fluxPtr,
                           const sourcePtr_Type&  sourcePtr )
{
    M_physicsPtr = physicsPtr;
    M_fluxPtr    = fluxPtr;
    M_sourcePtr  = sourcePtr;
}

void
OneDFSISolver::setCommunicator( const commPtr_Type& commPtr )
{
    M_commPtr = commPtr;
    M_displayer.setCommunicator( commPtr );
}

void
OneDFSISolver::setFESpace( const feSpacePtr_Type& feSpacePtr )
{
    M_feSpacePtr = feSpacePtr;

    //Elementary Matrices
    M_elementalMassMatrixPtr.reset ( new MatrixElemental( M_feSpacePtr->fe().nbFEDof(), 1, 1 ) );
    M_elementalStiffnessMatrixPtr.reset( new MatrixElemental( M_feSpacePtr->fe().nbFEDof(), 1, 1 ) );
    M_elementalGradientMatrixPtr.reset ( new MatrixElemental( M_feSpacePtr->fe().nbFEDof(), 1, 1 ) );
    M_elementalDivergenceMatrixPtr.reset  ( new MatrixElemental( M_feSpacePtr->fe().nbFEDof(), 1, 1 ) );

    //Vectors
    for ( UInt i(0) ; i < 2 ; ++i )
    {
        M_rhs[i].reset( new vector_Type( M_feSpacePtr->map() ) );
        M_residual[i].reset( new vector_Type( M_feSpacePtr->map() ) );
        M_fluxVector[i].reset( new vector_Type( M_feSpacePtr->map() ) );
        M_sourceVector[i].reset( new vector_Type( M_feSpacePtr->map() ) );
    }

    //Matrix
    M_homogeneousMassMatrixPtr.reset( new matrix_Type( M_feSpacePtr->map() ) );
    M_homogeneousGradientMatrixPtr.reset( new matrix_Type( M_feSpacePtr->map() ) );
}

void
OneDFSISolver::setLinearSolver( const linearSolverPtr_Type& linearSolverPtr )
{
    M_linearSolverPtr = linearSolverPtr;
}

void
OneDFSISolver::setLinearViscoelasticSolver( const linearSolverPtr_Type& linearViscoelasticSolverPtr )
{
    M_linearViscoelasticSolverPtr = linearViscoelasticSolverPtr;
}

// ===================================================
// Get Methods
// ===================================================
UInt
OneDFSISolver::boundaryDOF( const bcSide_Type& bcSide ) const
{
    switch ( bcSide )
    {
    case OneDFSI::left:

        return 0;

        break;

    case OneDFSI::right:

        return M_physicsPtr->data()->numberOfNodes() - 1;

        break;

    default:

        std::cout << "Warning: bcSide \"" << bcSide << "\" not available!" << std::endl;

        return 0;
    }
}

Real
OneDFSISolver::boundaryValue( const solution_Type& solution, const bcType_Type& bcType, const bcSide_Type& bcSide ) const
{
    UInt boundaryDof( boundaryDOF( bcSide ) );

    switch ( bcType )
    {
    case OneDFSI::A:

        return (*solution.find("A")->second)( boundaryDof );

    case OneDFSI::Q:

        // Flow rate is positive with respect to the outgoing normal
        return (*solution.find("Q")->second)( boundaryDof ) * ( ( bcSide == OneDFSI::left ) ? -1. : 1. );

    case OneDFSI::W1:

        return (*solution.find("W1")->second)( boundaryDof );

    case OneDFSI::W2:

        return (*solution.find("W2")->second)( boundaryDof );

    case OneDFSI::P:

        return (*solution.find("P")->second)( boundaryDof );

    case OneDFSI::S:

        return -(*solution.find("P")->second)( boundaryDof );

    case OneDFSI::T:
    {
        Real P     = ( *solution.find("P")->second )( boundaryDof );
        Real rho   = M_physicsPtr->data()->densityRho();
        Real alpha = M_physicsPtr->data()->alpha( boundaryDof );
        Real Q     = ( *solution.find("Q")->second )( boundaryDof );
        Real A     = ( *solution.find("A")->second )( boundaryDof );

        // Note that the kinetic contribution should account for the real velocity profile
        // through the alpha coefficient (i.e., the Coriolis coefficient)
        return -P - 0.5 * rho * alpha * Q * Q / ( A * A );
    }
    default:

        std::cout << "Warning: bcType \"" << bcType << "\"not available!" << std::endl;
        return 0.;

    }
}

void
OneDFSISolver::boundaryEigenValuesEigenVectors( const bcSide_Type& bcSide,
                                                const solution_Type& solution,
                                                container2D_Type& eigenvalues,
                                                container2D_Type& leftEigenvector1,
                                                container2D_Type& leftEigenvector2 )
{
    UInt boundaryDof( boundaryDOF( bcSide ) );

    M_fluxPtr->eigenValuesEigenVectors( (*solution.find("A")->second)( boundaryDof ),
                                     (*solution.find("Q")->second)( boundaryDof ),
                                     eigenvalues, leftEigenvector1, leftEigenvector2,
                                     boundaryDof );
}

// ===================================================
// Private Methods
// ===================================================
void
OneDFSISolver::updateFlux( const solution_Type& solution )
{
    Real Ai, Qi;

    for ( UInt iNode(0); iNode < M_physicsPtr->data()->numberOfNodes() ; ++iNode )
    {
        Ai = (*solution.find("A")->second)( iNode );
        Qi = (*solution.find("Q")->second)( iNode );

        (*M_fluxVector[0])( iNode ) = M_fluxPtr->flux( Ai, Qi, 0, iNode );
        (*M_fluxVector[1])( iNode ) = M_fluxPtr->flux( Ai, Qi, 1, iNode );
    }
}

void
OneDFSISolver::updatedFdU( const solution_Type& solution )
{
    // first update the Flux vector
    updateFlux( solution );

    // then update the derivative of the Flux vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;

    for ( UInt iElement(0); iElement < M_physicsPtr->data()->numberOfElements(); ++iElement )
    {
        // for P1Seg and appropriate mesh only!
        Aii   = (*solution.find("A")->second)( iElement );
        Qii   = (*solution.find("Q")->second)( iElement );

        Aiip1 = (*solution.find("A")->second)( iElement + 1 );
        Qiip1 = (*solution.find("Q")->second)( iElement + 1 );

        for ( UInt ii=0; ii<2; ++ii )
        {
            for ( UInt jj=0; jj<2; ++jj )
            {
                tmp  = M_fluxPtr->dFdU(   Aii,   Qii, ii, jj, iElement );     // left node of current element
                tmp += M_fluxPtr->dFdU( Aiip1, Qiip1, ii, jj, iElement + 1 ); // right node of current element

                M_dFdUVector[ 2*ii + jj ]( iElement ) = 0.5 * tmp;
            }
        }
    }
}

void
OneDFSISolver::updateSource( const solution_Type& solution )
{
    Real Ai, Qi;

    for ( UInt iNode(0); iNode < M_physicsPtr->data()->numberOfNodes() ; ++iNode )
    {
        Ai = (*solution.find("A")->second)( iNode );
        Qi = (*solution.find("Q")->second)( iNode );

        (*M_sourceVector[0])( iNode ) = M_sourcePtr->source( Ai, Qi, 0, iNode );
        (*M_sourceVector[1])( iNode ) = M_sourcePtr->source( Ai, Qi, 1, iNode );
    }
}

void
OneDFSISolver::updatedSdU( const solution_Type& solution )
{
    // first update the Source vector
    updateSource( solution );

    // then update the derivative of the Source vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;

    for ( UInt iElement(0); iElement < M_physicsPtr->data()->numberOfElements(); ++iElement )
    {
        // for P1Seg and appropriate mesh only!
        Aii   = (*solution.find("A")->second)( iElement);
        Qii   = (*solution.find("Q")->second)( iElement);
        Aiip1 = (*solution.find("A")->second)( iElement + 1 );
        Qiip1 = (*solution.find("Q")->second)( iElement + 1 );

        for ( UInt ii=0; ii<2; ++ii )
        {
            for ( UInt jj=0; jj<2; ++jj )
            {
                tmp =  M_sourcePtr->dSdU(   Aii,   Qii, ii, jj, iElement );     // left node of current element
                tmp += M_sourcePtr->dSdU( Aiip1, Qiip1, ii, jj, iElement + 1 ); // right node of current element

                M_dSdUVector[ 2*ii + jj ]( iElement ) = 0.5 * tmp;
            }
        }
    }
}

void
OneDFSISolver::updateMatrices()
{
    // Matrices initialization
    for ( UInt i(0); i < 4; ++i )
    {
        M_dSdUMassMatrixPtr[i].reset( new matrix_Type( M_feSpacePtr->map() ) );
        M_dFdUStiffnessMatrixPtr[i].reset( new matrix_Type( M_feSpacePtr->map() ) );
        M_dFdUGradientMatrixPtr[i].reset( new matrix_Type( M_feSpacePtr->map() ) );
        M_dSdUDivergenceMatrixPtr[i].reset( new matrix_Type( M_feSpacePtr->map() ) );
    }

    // Elementary computation and matrix assembling
    for ( UInt iElement(0); iElement < M_physicsPtr->data()->numberOfElements(); ++iElement )
    {
        // Update the current element
        M_feSpacePtr->fe().update( M_feSpacePtr->mesh()->edgeList( iElement), UPDATE_DPHI | UPDATE_WDET );

        for ( UInt ii(0); ii < 2; ++ii )
        {
            for ( UInt jj(0); jj < 2; ++jj )
            {
                // Update the elemental matrices
                updateElementalMatrices( M_dFdUVector[ 2*ii + jj ]( iElement ), M_dSdUVector[ 2*ii + jj ]( iElement ) );

                // Assemble the global matrices
                matrixAssemble( ii, jj );
            }
        }
    }
}

void
OneDFSISolver::updateElementalMatrices( const Real& dFdU, const Real& dSdU )
{
    // Set the elementary matrices to 0.
    M_elementalMassMatrixPtr->zero();
    M_elementalStiffnessMatrixPtr->zero();
    M_elementalGradientMatrixPtr->zero();
    M_elementalDivergenceMatrixPtr->zero();

    // Update the mass matrix
    mass( dSdU, *M_elementalMassMatrixPtr, M_feSpacePtr->fe(), 0, 0 );
//    std::cout << "Elem Stiff matrix :" << std::endl;
//    M_elementalMassMatrixPtr->showMe( std::cout );

    // Update the stiffness matrix
    stiff( dFdU, *M_elementalStiffnessMatrixPtr, M_feSpacePtr->fe(), 0 ,0 );
    //AssemblyElemental::stiffness( *M_elementalStiffnessMatrixPtr, M_feSpacePtr->fe(), M_coeffStiff, 0 );
//    std::cout << "Elem Stiff matrix :" << std::endl;
//    M_elementalStiffnessMatrixPtr->showMe( std::cout );

    /*! Update the gradient matrix
      gradient operator:
      grad_{ij} = \int_{fe} coeff \phi_j \frac{d \phi_i}{d x}

      BEWARE :
      \param 0: the first argument "0" corresponds to the first
      and only coordinate (1D!), and HERE it starts from 0... (Damm'!)

      \param - M_coeffGrad: the sign "-" in the second argument
      is added to correspond to the described operator.
      (There is a minus in the elemOper implementation).
    */
    grad( 0, -dFdU, *M_elementalGradientMatrixPtr, M_feSpacePtr->fe(), M_feSpacePtr->fe(), 0, 0 );
//    std::cout << "Elem Grad matrix :" << std::endl;
//    M_elementalGradientMatrixPtr->showMe( std::cout );

    /*! update the divergence matrix
      divergence operator: (transpose of the gradient)
      div_{ij} = \int_{fe} coeff \frac{d \phi_j}{d x} \phi_i

      \note formally this M_elementalDivergenceMatrixPtr is not necessary
      as it is the transpose of the M_elementalGradientMatrixPtr.
      But for the sake of clarity, I prefer to keep it. (low cost!)

      BEWARE : same remarks as grad (see above).
    */
    div( 0, -dSdU, *M_elementalDivergenceMatrixPtr, M_feSpacePtr->fe(), M_feSpacePtr->fe(), 0, 0 );
//    std::cout << "Elem Div matrix :" << std::endl;
//    M_elementalDivergenceMatrixPtr->showMe( std::cout );
}

void
OneDFSISolver::matrixAssemble( const UInt& ii, const UInt& jj )
{
    // Assemble the mass matrix
    assembleMatrix( *M_dSdUMassMatrixPtr[ 2*ii + jj ], *M_elementalMassMatrixPtr, M_feSpacePtr->fe(), M_feSpacePtr->dof(), 0, 0, 0, 0 );

    // Assemble the stiffness matrix
    assembleMatrix( *M_dFdUStiffnessMatrixPtr[ 2*ii + jj ], *M_elementalStiffnessMatrixPtr, M_feSpacePtr->fe(), M_feSpacePtr->dof(), 0, 0, 0, 0 );

    // Assemble the gradient matrix
    assembleMatrix( *M_dFdUGradientMatrixPtr[ 2*ii + jj ], *M_elementalGradientMatrixPtr, M_feSpacePtr->fe(), M_feSpacePtr->dof(), 0, 0, 0, 0 );

    // Assemble the divergence matrix
    assembleMatrix( *M_dSdUDivergenceMatrixPtr[ 2*ii + jj ], *M_elementalDivergenceMatrixPtr, M_feSpacePtr->fe(), M_feSpacePtr->dof(), 0, 0, 0, 0 );
}

void
OneDFSISolver::applyDirichletBCToMatrix( matrix_Type& matrix )
{
    // Dirichlet BC
    matrix.globalAssemble();
    matrix.diagonalize( 0, 1, 0 );
    matrix.diagonalize( M_physicsPtr->data()->numberOfNodes() - 1, 1, 0 );

    //matrix.spy("SystemMatrix");
}

OneDFSISolver::vector_Type
OneDFSISolver::inertialFlowRateCorrection( const vector_Type& flux )
{
    matrix_Type matrixLHS(M_feSpacePtr->map());
    matrix_Type stiffRHS (M_feSpacePtr->map());

    MatrixElemental elmatMassLHS  (M_feSpacePtr->fe().nbFEDof(), 1, 1);
    MatrixElemental elmatStiffLHS (M_feSpacePtr->fe().nbFEDof(), 1, 1);
    MatrixElemental elmatStiffRHS (M_feSpacePtr->fe().nbFEDof(), 1, 1);

    vector_Type rhs(M_feSpacePtr->map());

    Real coeffMass;
    Real coeffStiff;

    Real m, meanA0;

    matrixLHS *= 0.;
    stiffRHS  *= 0.;

    // Elementary computation and matrix assembling
    for ( UInt iElement(0); iElement < M_physicsPtr->data()->numberOfElements(); ++iElement )
    {
        // set the elementary matrices to 0.
        elmatMassLHS. zero();
        elmatStiffLHS.zero();
        elmatStiffRHS.zero();

        coeffMass  = (*M_rhs[0])( iElement ) + (*M_rhs[0])( iElement + 1 );
        coeffMass /= 2;
        coeffMass  = 1./ coeffMass;

        meanA0  = M_physicsPtr->data()->area0(iElement) + M_physicsPtr->data()->area0(iElement + 1);
        meanA0 /= 2;

        m = M_physicsPtr->data()->densityWall()*M_physicsPtr->data()->thickness(iElement + 1)/
            ( 2*std::sqrt(4*std::atan(1.))*std::sqrt(meanA0) );

        coeffStiff = m/M_physicsPtr->data()->densityRho();

        // Update the current element
        M_feSpacePtr->fe().update( M_feSpacePtr->mesh()->edgeList(iElement), UPDATE_DPHI | UPDATE_WDET );

        mass (   coeffMass,  elmatMassLHS,  M_feSpacePtr->fe(), 0, 0 );
        stiff(   coeffStiff, elmatStiffLHS, M_feSpacePtr->fe(), 0, 0 );
        stiff( - coeffStiff, elmatStiffRHS, M_feSpacePtr->fe(), 0, 0 );

        // assemble the mass and grad matrices
        assembleMatrix( matrixLHS, elmatMassLHS, M_feSpacePtr->fe(), M_feSpacePtr->dof(), 0, 0, 0, 0 );
        assembleMatrix( matrixLHS, elmatStiffLHS, M_feSpacePtr->fe(), M_feSpacePtr->dof(), 0, 0, 0, 0 );
        assembleMatrix( stiffRHS, elmatStiffRHS, M_feSpacePtr->fe(), M_feSpacePtr->dof(), 0, 0, 0, 0 );

#ifdef HAVE_LIFEV_DEBUG
        debugStream( 6310 ) << "\n\tm = "           << m << "\n\t_coeffMass = "  << coeffMass << "\n\t_coeffStiff = " << coeffStiff << "\n";
#endif
    } // end loop on elements

    // update rhs
    //_stiffRHS.Axpy( 1., flux , 0., _rhs );

    rhs = stiffRHS*flux;

    UInt firstDof = 0;
    UInt lastDof  = rhs.size()-1;

    // symmetric treatment (cholesky can be used)
    // first row modified (Dirichlet)
    rhs( firstDof ) = 0.;
    // second row modified (for symmetry)
    //  _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    // last row modified (Dirichlet)
    rhs( lastDof ) = 0.;
    // penultimate row modified (for symmetry)
    //  _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    applyDirichletBCToMatrix( matrixLHS );

    //_tridiagsolver.Factor( _matrixLHS );

    // cholesky or lapack lu solve
    // solve the system: rhs1 = massFactor^{-1} * rhs1
    //_tridiagsolver.Solve( _matrixLHS, _rhs );

    vector_Type sol( rhs);

    M_linearSolverPtr->setMatrix( matrixLHS );
    //@int numIter = M_linearSolverPtr->solveSystem( _rhs, _sol, _matrixLHS, true);

    //std::cout <<" iterations number :  " << numIter << std::endl;

    return sol;
}

OneDFSISolver::vector_Type
OneDFSISolver::longitudinalFlowRateCorrection()
{
    matrix_Type massLHS(M_feSpacePtr->map());
    matrix_Type massRHS(M_feSpacePtr->map());

    //TriDiagCholesky< Real, matrix_Type, Vector > _tridiagsolver(M_physicsPtr->data()->numberOfNodes());

    MatrixElemental elmatMassLHS (M_feSpacePtr->fe().nbFEDof(),1,1);
    MatrixElemental elmatMassRHS (M_feSpacePtr->fe().nbFEDof(),1,1);

    vector_Type rhs(M_feSpacePtr->map());
    // let g = sqrt(A) - sqrt(A0)
    // now f = _d3g_dz3
    vector_Type g(M_feSpacePtr->map());
    vector_Type f(M_feSpacePtr->map());

    //          _g = *M_rhs[0];
    for ( UInt iNode(0); iNode < M_physicsPtr->data()->numberOfNodes(); ++iNode )
        g(iNode) = std::sqrt((*M_rhs[0])(iNode)) - std::sqrt(M_physicsPtr->data()->area0(iNode));

    UInt iNode;

    Real coeffMassLHS;
    Real coeffMassRHS;

    Real a;
    //    std::ostringstream output;
    massLHS *= 0.;
    massRHS *= 0.;

    Real h( M_physicsPtr->data()->length() / static_cast<Real>(M_physicsPtr->data()->numberOfElements() - 1) );

    // Elementary computation and matrix assembling
    // Loop on elements
    for ( UInt iElement(0); iElement < M_physicsPtr->data()->numberOfElements() ; ++iElement )
    {
        iNode = iElement;

        // set the elementary matrices to 0.
        elmatMassLHS.zero();
        elmatMassRHS.zero();

        // coeff (1/A) (average values over the element)
        coeffMassLHS = (*M_rhs[0])( iNode ) + (*M_rhs[0])( iNode+1 );
        coeffMassLHS /= 2;
        coeffMassLHS = 1./ coeffMassLHS;

        a = M_physicsPtr->data()->inertialModulus() / std::sqrt(4*std::atan(1.));
        coeffMassRHS = M_physicsPtr->data()->dataTime()->timeStep() * a / M_physicsPtr->data()->densityRho();

        // backward differentiation when near to the left boundary
        // d3 A (xi)/ dz3 = (1/h^3) * ( -A(xi) + A(xi+3) - 3A(xi+2) + 3A(xi+1) )

        // forward differentiation when near to the right boundary
        // d3 A (xi)/ dz3 = (1/h^3) * ( A(xi) - A(xi-3) + 3A(xi-2) - 3A(xi-1) )

        // central differentiation otherwise
        // d3 A (xi)/ dz3 = (1/h^3) * ( -A(xi-2) + 2A(xi-1) - 2A(xi+1) + A(xi+2) )

#ifdef HAVE_LIFEV_DEBUG
        debugStream( 6310 ) << "\ninode = " << iNode << "\n";
#endif
        if (iNode<2)
        { // backward differentiation
            f( iNode ) = - g( iNode ) + g( iNode+3 ) - 3 * g( iNode+2 ) + 3 * g( iNode+1 );
#ifdef HAVE_LIFEV_DEBUG
            debugStream( 6310 ) << "\n\tbackward differentiation = " << coeffMassLHS << "\n";
#endif
        }
        else if (iNode>M_feSpacePtr->mesh()->numEdges()-2)
        { // forward differentiation
            f( iNode ) = g( iNode ) - g( iNode-3 ) + 3 * g( iNode-2 ) - 3 * g( iNode-1 );
#ifdef HAVE_LIFEV_DEBUG
            debugStream( 6310 ) << "\n\forward differentiation = " << coeffMassLHS << "\n";
#endif
        }
        else
        { // central differentiation
            f( iNode ) = - g( iNode - 2 ) + 2 * g( iNode-1 ) - 2 * g( iNode+1 ) + g( iNode+2 );
#ifdef HAVE_LIFEV_DEBUG
            debugStream( 6310 ) << "\n\tcentral differentiation = " << coeffMassLHS << "\n";
#endif
        }

        f(iNode) *= 1 / ( 2 * OneDFSI::pow30(h, 3) );

        // Update the current element
        M_feSpacePtr->fe().update( M_feSpacePtr->mesh()->edgeList(iElement), UPDATE_DPHI | UPDATE_WDET );

        mass( coeffMassLHS, elmatMassLHS, M_feSpacePtr->fe(),0, 0 );
        mass( coeffMassRHS, elmatMassRHS, M_feSpacePtr->fe(),0, 0 );

        // assemble the mass and grad matrices
        //assemb_mat( massLHS, elmatMassLHS, M_feSpacePtr->fe(), M_feSpacePtr->dof() , 0, 0 );
        assembleMatrix( massLHS, elmatMassLHS, M_feSpacePtr->fe(), M_feSpacePtr->dof() , 0, 0, 0, 0 );
        //assemb_mat( massRHS, elmatMassRHS, M_feSpacePtr->fe(), M_feSpacePtr->dof() , 0, 0 );
        assembleMatrix( massRHS, elmatMassRHS, M_feSpacePtr->fe(), M_feSpacePtr->dof() , 0, 0, 0, 0 );

#ifdef HAVE_LIFEV_DEBUG
        debugStream( 6310 ) << "\n\t_coeffMassLHS = " << coeffMassLHS << "\n";
        debugStream( 6310 ) << "\n\t_coeffMassRHS = " << coeffMassRHS << "\n";
#endif
    } // end loop on elements

    // update rhs
    //    _massLHS.Axpy( 1., (*solution["P"]) , 0., _rhs );
    rhs = massRHS*f;

    UInt firstDof = 0;
    UInt lastDof  = rhs.size() - 1;

    // symmetric treatment (cholesky can be used)
    // first row modified (Dirichlet)
    rhs( firstDof ) = 0.;
    // second row modified (for symmetry)
    //    _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    // last row modified (Dirichlet)
    rhs( lastDof ) = 0.;
    // penultimate row modified (for symmetry)
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    applyDirichletBCToMatrix( massLHS);

    //@_tridiagsolver.Factor( _massLHS );

    // cholesky or lapack lu solve
    // solve the system: rhs1 = massFactor^{-1} * rhs1
    //@_tridiagsolver.Solve( _massLHS, _rhs );

    return rhs;
}

}
