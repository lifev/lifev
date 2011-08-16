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
 *  @mantainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <life/lifesolver/OneDimensionalSolver.hpp>

namespace LifeV
{

std::map< std::string, OneDimensional::physicsType_Type > OneDimensional::physicsMap;
std::map< std::string, OneDimensional::fluxTerm_Type >    OneDimensional::fluxMap;
std::map< std::string, OneDimensional::sourceTerm_Type >  OneDimensional::sourceMap;

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalSolver::OneDimensionalSolver():
    M_physics                   (),
    M_flux                      (),
    M_source                    (),
    M_feSpace                   (),
    M_comm                      (),
    M_displayer                 (),
    M_elementalMassMatrix       (),
    M_elementalStiffnessMatrix  (),
    M_elementalGradientMatrix   (),
    M_elementalDivergenceMatrix (),
    M_rhs                       (),
    M_residual                  (),
    M_fluxVector                (),
    M_sourceVector              (),
    M_dFdUVector                (),
    M_dSdUVector                (),
    M_homogeneousMassMatrix     (),
    M_homogeneousGradientMatrix (),
    M_dSdUMassMatrix            (),
    M_dFdUStiffnessMatrix       (),
    M_dFdUGradientMatrix        (),
    M_dSdUDivergenceMatrix      (),
    M_linearSolver              (),
    M_linearViscoelasticSolver  ()
{
}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalSolver::buildConstantMatrices()
{
    std::fill( M_dFdUVector.begin(), M_dFdUVector.end(), ublas::zero_vector<Real>( M_physics->data()->numberOfNodes()) );
    std::fill( M_dSdUVector.begin(), M_dSdUVector.end(), ublas::zero_vector<Real>( M_physics->data()->numberOfNodes()) );

    for ( UInt i(0); i < 4; ++i )
    {
        M_dSdUMassMatrix[i].reset( new matrix_Type( M_feSpace->map() ) );
        M_dFdUStiffnessMatrix[i].reset( new matrix_Type( M_feSpace->map() ) );
        M_dFdUGradientMatrix[i].reset( new matrix_Type( M_feSpace->map() ) );
        M_dSdUDivergenceMatrix[i].reset( new matrix_Type( M_feSpace->map() ) );
    }

    // Elementary computation and matrix assembling
    for ( UInt iElement(0); iElement < M_physics->data()->numberOfElements(); ++iElement )
    {
        // set the elementary matrices to 0.
        M_elementalMassMatrix->zero();
        M_elementalGradientMatrix->zero();

        // update the current element
        M_feSpace->fe().update( M_feSpace->mesh()->edgeList( iElement ), UPDATE_DPHI | UPDATE_WDET );

        // update the mass and grad matrices
        mass( 1, *M_elementalMassMatrix, M_feSpace->fe(), 0, 0 );
        grad( 0 , -1, *M_elementalGradientMatrix, M_feSpace->fe(), M_feSpace->fe(), 0, 0 );

        // assemble the mass and grad matrices
        assembleMatrix( *M_homogeneousMassMatrix, *M_elementalMassMatrix, M_feSpace->fe(), M_feSpace->dof() , 0, 0, 0, 0 );
        assembleMatrix( *M_homogeneousGradientMatrix, *M_elementalGradientMatrix, M_feSpace->fe(), M_feSpace->dof() , 0, 0, 0, 0  );
    }

    // Dirichlet boundary conditions set in the mass matrix
    M_homogeneousMassMatrix->globalAssemble();
    M_homogeneousGradientMatrix->globalAssemble();

    // In the classical case the linear system use a mass matrix (with Dirichlet BC)
    //matrixPtr_Type systemMatrix( new matrix_Type( M_feSpace->map() ) );
    //systemMatrix->insertValueDiagonal( M_physics->data()->mesh()->meanH() );
    matrix_Type systemMatrix( *M_homogeneousMassMatrix );
    applyDirichletBCToMatrix( systemMatrix );
    M_linearSolver->setMatrix( systemMatrix );
}

void
OneDimensionalSolver::setupSolution( solution_Type& solution, const MapEpetra& map, const bool& onlyMainQuantities )
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
    if ( M_physics->data()->viscoelasticWall() )
    {
        solution["Q_visc"].reset( new vector_Type( map ) ); // viscoelastic contribution to the flux
        solution["P_visc"].reset( new vector_Type( map ) ); // viscoelastic contribution to the pressure
    }

    // correction flux with inertial term
    if ( M_physics->data()->inertialWall() )
        solution["Q_inert"].reset( new vector_Type( map ) );

    // correction flux with longitudinal term
    if ( M_physics->data()->longitudinalWall() )
        solution["Q_long"].reset( new vector_Type( map ) );

    // Initialize solution to zero
    for ( solutionConstIterator_Type i = solution.begin(); i != solution.end(); ++i )
        *i->second = 0.;
}

void
OneDimensionalSolver::initialize( solution_Type& solution )
{
    for ( UInt iNode(0); iNode < M_physics->data()->numberOfNodes() ; ++iNode )
    {
        (*solution["A"])[iNode] = M_physics->data()->area0( iNode );
        (*solution["Q"])[iNode] = 0;
    }

    // Compute W1 and W2 from A and Q
    computeW1W2( solution );

    // Compute A/A0 from A
    computeAreaRatio( solution );

    // Compute initial pressure (taking into account the viscoelastic wall)
    M_physics->setArea_tn( *solution["A"] );
    computePressure( solution, M_physics->data()->dataTime()->timeStep() );
}

void
OneDimensionalSolver::computeW1W2( solution_Type& solution )
{
    for ( UInt iNode(0); iNode < M_physics->data()->numberOfNodes() ; ++iNode )
    {
        M_physics->fromUToW( ( *solution["W1"]) [iNode], (*solution["W2"]) [iNode],
                             ( *solution["A"] ) [iNode], (*solution["Q"] ) [iNode], iNode );
    }
}

void
OneDimensionalSolver::computePressure( solution_Type& solution, const Real& timeStep )
{
    for ( UInt iNode(0); iNode < M_physics->data()->numberOfNodes() ; ++iNode )
    {
        ( *solution["P"] ) [iNode] = M_physics->elasticPressure( ( *solution["A"] ) [iNode], iNode )
                                   + M_physics->externalPressure();

        if ( M_physics->data()->viscoelasticWall() )
        {
            ( *solution["P_visc"] ) [iNode] = M_physics->viscoelasticPressure( ( *solution["A"]) [iNode], timeStep, iNode );
            ( *solution["P"] )      [iNode] += ( *solution["P_visc"] ) [iNode];
        }
    }
}

void
OneDimensionalSolver::computeAreaRatio( solution_Type& solution )
{
    for ( UInt iNode(0); iNode < M_physics->data()->numberOfNodes() ; ++iNode )
        ( *solution["AoverA0minus1"] ) [iNode] = ( *solution["A"] ) [iNode] / M_physics->data()->area0( iNode ) - 1;
}

void
OneDimensionalSolver::computeArea( solution_Type& solution )
{
    for ( UInt iNode(0); iNode < M_physics->data()->numberOfNodes() ; ++iNode )
        ( *solution["A"] ) [iNode] = ( (*solution["AoverA0minus1"] ) [iNode] + 1 ) * M_physics->data()->area0( iNode );
}

void
OneDimensionalSolver::updateRHS( const solution_Type& solution, const Real& timeStep )
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
        *M_residual[i] += ( *M_homogeneousGradientMatrix ) * ( timeStep * *M_fluxVector[i] );

        // rhs = rhs - dt * mass * S(Un)
        *M_residual[i] += ( *M_homogeneousMassMatrix ) * ( - timeStep * *M_sourceVector[i] );

        for ( UInt j(0); j < 2; ++j )
        {
            // rhs = rhs - dt^2/2 * gradDiffFlux * S(Un)
            M_dFdUGradientMatrix[2*i + j]->globalAssemble();
            *M_residual[i] += *M_dFdUGradientMatrix[2*i + j] * ( -dt2over2 * *M_sourceVector[j] );

            // rhs = rhs + dt^2/2 * divDiffSrc * F(Un)
            M_dSdUDivergenceMatrix[2*i + j]->globalAssemble();
            *M_residual[i] += *M_dSdUDivergenceMatrix[2*i + j] * ( dt2over2 * *M_fluxVector[j] );

            // rhs = rhs - dt^2/2 * stiffDiffFlux * F(Un)
            M_dFdUStiffnessMatrix[2*i + j]->globalAssemble();
            *M_residual[i] += *M_dFdUStiffnessMatrix[2*i + j]*( -dt2over2 * *M_fluxVector[j] );

            // rhs = rhs + dt^2/2 * massDiffSrc * S(Un)
            M_dSdUMassMatrix[2*i + j]->globalAssemble();
            *M_residual[i] += *M_dSdUMassMatrix[2*i + j] * ( dt2over2 * *M_sourceVector[j] );
        }
    }

    // rhs = mass * Un + residual
    *M_rhs[0]  = *M_residual[0] + ( *M_homogeneousMassMatrix ) * *solution.find("A")->second;
    *M_rhs[1]  = *M_residual[1] + ( *M_homogeneousMassMatrix ) * *solution.find("Q")->second;
}

void
OneDimensionalSolver::iterate( OneDimensionalBCHandler& bcHandler, solution_Type& solution, const Real& time, const Real& timeStep )
{
    // Apply BC to RHS
    bcHandler.applyBC( time, timeStep, solution, M_flux, M_rhs );

    // Compute A^n+1
    vector_Type area( *M_rhs[0] );
    M_linearSolver->solveSystem( *M_rhs[0], area, M_homogeneousMassMatrix );

    // Compute Q^n+1
    vector_Type flowRate( *M_rhs[1] );
    M_linearSolver->solveSystem( *M_rhs[1], flowRate, M_homogeneousMassMatrix );

    // Correct flux with inertial, viscoelastic and longitudinal terms
    if ( M_physics->data()->inertialWall() )
    {
        *solution["Q_inert"] = inertialFlowRateCorrection( flowRate );
        flowRate += *solution["Q_inert"];
    }

    if ( M_physics->data()->viscoelasticWall() )
    {
        *solution["Q_visc"] = viscoelasticFlowRateCorrection( area, flowRate, timeStep, bcHandler );
        flowRate += *solution["Q_visc"];
    }

    if ( M_physics->data()->longitudinalWall() )
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

OneDimensionalSolver::vector_Type
OneDimensionalSolver::viscoelasticFlowRateCorrection( const vector_Type& area, const vector_Type& flowRate, const Real& timeStep, OneDimensionalBCHandler& bcHandler, const bool& updateSystemMatrix )
{
    // Matrix
    matrix_Type systemMatrix( M_feSpace->map() );
    matrix_Type stiffnessMatrix( M_feSpace->map() );

    Real massCoefficient;
    Real stiffnessCoefficient;

    // RHS
    vector_Type rhs( M_feSpace->map() );
    rhs = 0;

    // Elementary computation and matrix assembling
    for ( UInt iElement(0); iElement < M_physics->data()->numberOfElements() ; ++iElement )
    {
        // Update the current element
        M_feSpace->fe().update( M_feSpace->mesh()->edgeList(iElement), UPDATE_DPHI | UPDATE_WDET );

        // Compute mass coefficient
        massCoefficient = 1 / ( 0.5 * ( area[ iElement ] + area[ iElement + 1 ] ) );
//      massCoefficient = 1 / ( 0.5 * ( M_physics->data()->area0( iElement ) + M_physics->data()->area0( iElement + 1 ) ) );

        // Compute stiffness coefficient
        stiffnessCoefficient  = timeStep * 0.5 * ( M_physics->data()->viscoelasticCoefficient( iElement ) + M_physics->data()->viscoelasticCoefficient( iElement + 1 ) )
                              / M_physics->data()->densityRho() * massCoefficient * std::sqrt( massCoefficient );

        // Set the elementary matrices to 0.
        M_elementalMassMatrix->zero();
        M_elementalStiffnessMatrix->zero();

        // Assemble the elemental matrix
        mass(  massCoefficient,      *M_elementalMassMatrix,      M_feSpace->fe(), 0, 0 );
        stiff( stiffnessCoefficient, *M_elementalStiffnessMatrix, M_feSpace->fe(), 0, 0 );

        // Assemble the stiffness matrix
        assembleMatrix( systemMatrix,    *M_elementalMassMatrix,      M_feSpace->fe(), M_feSpace->dof(), 0, 0, 0, 0 );
        assembleMatrix( stiffnessMatrix, *M_elementalStiffnessMatrix, M_feSpace->fe(), M_feSpace->dof(), 0, 0, 0, 0 );

        // Natural BC
        if ( iElement == 0 )
            rhs( 0 )            -= stiffnessCoefficient * M_flux->physics()->data()->computeSpatialDerivativeAtNode( flowRate, 0, 1 );
        if ( iElement == M_physics->data()->numberOfElements() - 1 )
            rhs( iElement + 1 ) += stiffnessCoefficient * M_flux->physics()->data()->computeSpatialDerivativeAtNode( flowRate, iElement + 1, 1 );
    }

    // System Matrix = MassMatrix + stiffnessCoefficient * StiffnessMatrix
    stiffnessMatrix.globalAssemble();
    systemMatrix.globalAssemble();
    systemMatrix += stiffnessMatrix;

    // RHS
    rhs += stiffnessMatrix * (-flowRate);

    // Apply BC to Matrix and RHS
    bcHandler.applyViscoelasticBC( M_flux, systemMatrix, rhs );

    // Compute flow rate correction at t^n+1
    vector_Type flowRateCorrection( rhs );

    if ( updateSystemMatrix )
        M_linearViscoelasticSolver->setMatrix( systemMatrix );
    M_linearViscoelasticSolver->solveSystem( rhs, flowRateCorrection, M_homogeneousMassMatrix );

    return flowRateCorrection;
}

Real
OneDimensionalSolver::computeCFL( const solution_Type& solution, const Real& timeStep ) const
{
    Real lambdaMax( 0. );

    container2D_Type eigenvalues;
    container2D_Type leftEigenvector1;
    container2D_Type leftEigenvector2;

    for ( UInt iNode(0); iNode < M_physics->data()->numberOfNodes() ; ++iNode )
    {
        // compute the eigenvalues at node
        M_flux->eigenValuesEigenVectors( ( *solution.find("A")->second ) ( iNode ),
                                         ( *solution.find("Q")->second ) ( iNode ),
                                           eigenvalues, leftEigenvector1, leftEigenvector2, iNode );

        lambdaMax = std::max<Real>( std::max<Real>( std::fabs(eigenvalues[0]), std::fabs(eigenvalues[1]) ), lambdaMax );
    }

    return lambdaMax * timeStep / M_feSpace->mesh()->minH();
}

void
OneDimensionalSolver::resetOutput( const solution_Type& solution )
{
    std::ofstream outfile;
    for ( solutionConstIterator_Type i = solution.begin(); i != solution.end(); ++i )
    {
        std::string file = M_physics->data()->postprocessingDirectory() + "/" + M_physics->data()->postprocessingFile() + "_" + i->first + ".m";
        outfile.open( file.c_str(), std::ios::trunc );
        outfile.close();
    }
}

void
OneDimensionalSolver::postProcess( const solution_Type& solution, const Real& time )
{
    std::ofstream outfile;
    for ( solutionConstIterator_Type i = solution.begin(); i != solution.end(); ++i )
    {
        std::string file = M_physics->data()->postprocessingDirectory() + "/" + M_physics->data()->postprocessingFile() + "_" + i->first + ".m";
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
OneDimensionalSolver::setProblem( const physicsPtr_Type& physics,
                                  const fluxPtr_Type&    flux,
                                  const sourcePtr_Type&  source )
{
    M_physics = physics;
    M_flux    = flux;
    M_source  = source;
}

void
OneDimensionalSolver::setCommunicator( const commPtr_Type& comm )
{
    M_comm = comm;
    M_displayer.setCommunicator( comm );
}

void
OneDimensionalSolver::setFESpace( const feSpacePtr_Type& feSpace )
{
    M_feSpace = feSpace;

    //Elementary Matrices
    M_elementalMassMatrix.reset ( new MatrixElemental( M_feSpace->fe().nbFEDof(), 1, 1 ) );
    M_elementalStiffnessMatrix.reset( new MatrixElemental( M_feSpace->fe().nbFEDof(), 1, 1 ) );
    M_elementalGradientMatrix.reset ( new MatrixElemental( M_feSpace->fe().nbFEDof(), 1, 1 ) );
    M_elementalDivergenceMatrix.reset  ( new MatrixElemental( M_feSpace->fe().nbFEDof(), 1, 1 ) );

    //Vectors
    for ( UInt i(0) ; i < 2 ; ++i )
    {
        M_rhs[i].reset( new vector_Type( M_feSpace->map() ) );
        M_residual[i].reset( new vector_Type( M_feSpace->map() ) );
        M_fluxVector[i].reset( new vector_Type( M_feSpace->map() ) );
        M_sourceVector[i].reset( new vector_Type( M_feSpace->map() ) );
    }

    //Matrix
    M_homogeneousMassMatrix.reset( new matrix_Type( M_feSpace->map() ) );
    M_homogeneousGradientMatrix.reset( new matrix_Type( M_feSpace->map() ) );
}

void
OneDimensionalSolver::setLinearSolver( const linearSolverPtr_Type& linearSolver )
{
    M_linearSolver = linearSolver;
}

void
OneDimensionalSolver::setLinearViscoelasticSolver( const linearSolverPtr_Type& linearViscoelasticSolver )
{
    M_linearViscoelasticSolver = linearViscoelasticSolver;
}

// ===================================================
// Get Methods
// ===================================================
UInt
OneDimensionalSolver::boundaryDOF( const bcSide_Type& bcSide ) const
{
    switch ( bcSide )
    {
    case OneDimensional::left:

        return 0;

        break;

    case OneDimensional::right:

        return M_physics->data()->numberOfNodes() - 1;

        break;

    default:

        std::cout << "Warning: bcSide \"" << bcSide << "\" not available!" << std::endl;

        return 0;
    }
}

Real
OneDimensionalSolver::boundaryValue( const solution_Type& solution, const bcType_Type& bcType, const bcSide_Type& bcSide ) const
{
    UInt boundaryDof( boundaryDOF( bcSide ) );

    switch ( bcType )
    {
    case OneDimensional::A:

        return (*solution.find("A")->second)( boundaryDof );

    case OneDimensional::Q:

        // Flow rate is positive with respect to the outgoing normal
        return (*solution.find("Q")->second)( boundaryDof ) * ( ( bcSide == OneDimensional::left ) ? -1. : 1. );

    case OneDimensional::W1:

        return (*solution.find("W1")->second)( boundaryDof );

    case OneDimensional::W2:

        return (*solution.find("W2")->second)( boundaryDof );

    case OneDimensional::P:

        return (*solution.find("P")->second)( boundaryDof );

    case OneDimensional::S:

        return -(*solution.find("P")->second)( boundaryDof );

    default:

        std::cout << "Warning: bcType \"" << bcType << "\"not available!" << std::endl;
        return 0.;

    }
}

void
OneDimensionalSolver::boundaryEigenValuesEigenVectors( const bcSide_Type& bcSide,
                                                       const solution_Type& solution,
                                                             container2D_Type& eigenvalues,
                                                             container2D_Type& leftEigenvector1,
                                                             container2D_Type& leftEigenvector2 )
{
    UInt boundaryDof( boundaryDOF( bcSide ) );

    M_flux->eigenValuesEigenVectors( (*solution.find("A")->second)( boundaryDof ),
                                     (*solution.find("Q")->second)( boundaryDof ),
                                     eigenvalues, leftEigenvector1, leftEigenvector2,
                                     boundaryDof );
}

// ===================================================
// Private Methods
// ===================================================
void
OneDimensionalSolver::updateFlux( const solution_Type& solution )
{
    Real Ai, Qi;

    for ( UInt iNode(0); iNode < M_physics->data()->numberOfNodes() ; ++iNode )
    {
        Ai = (*solution.find("A")->second)( iNode );
        Qi = (*solution.find("Q")->second)( iNode );

        (*M_fluxVector[0])( iNode ) = M_flux->flux( Ai, Qi, 0, iNode );
        (*M_fluxVector[1])( iNode ) = M_flux->flux( Ai, Qi, 1, iNode );
    }
}

void
OneDimensionalSolver::updatedFdU( const solution_Type& solution )
{
    // first update the Flux vector
    updateFlux( solution );

    // then update the derivative of the Flux vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;

    for ( UInt iElement(0); iElement < M_physics->data()->numberOfElements(); ++iElement )
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
                tmp  = M_flux->dFdU(   Aii,   Qii, ii, jj, iElement );     // left node of current element
                tmp += M_flux->dFdU( Aiip1, Qiip1, ii, jj, iElement + 1 ); // right node of current element

                M_dFdUVector[ 2*ii + jj ]( iElement ) = 0.5 * tmp;
            }
        }
    }
}

void
OneDimensionalSolver::updateSource( const solution_Type& solution )
{
    Real Ai, Qi;

    for ( UInt iNode(0); iNode < M_physics->data()->numberOfNodes() ; ++iNode )
    {
        Ai = (*solution.find("A")->second)( iNode );
        Qi = (*solution.find("Q")->second)( iNode );

        (*M_sourceVector[0])( iNode ) = M_source->source( Ai, Qi, 0, iNode );
        (*M_sourceVector[1])( iNode ) = M_source->source( Ai, Qi, 1, iNode );
    }
}

void
OneDimensionalSolver::updatedSdU( const solution_Type& solution )
{
    // first update the Source vector
    updateSource( solution );

    // then update the derivative of the Source vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;

    for ( UInt iElement(0); iElement < M_physics->data()->numberOfElements(); ++iElement )
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
                tmp =  M_source->dSdU(   Aii,   Qii, ii, jj, iElement );     // left node of current element
                tmp += M_source->dSdU( Aiip1, Qiip1, ii, jj, iElement + 1 ); // right node of current element

                M_dSdUVector[ 2*ii + jj ]( iElement ) = 0.5 * tmp;
            }
        }
    }
}

void
OneDimensionalSolver::updateMatrices()
{
    // Matrices initialization
    for ( UInt i(0); i < 4; ++i )
    {
        M_dSdUMassMatrix[i].reset( new matrix_Type( M_feSpace->map() ) );
        M_dFdUStiffnessMatrix[i].reset( new matrix_Type( M_feSpace->map() ) );
        M_dFdUGradientMatrix[i].reset( new matrix_Type( M_feSpace->map() ) );
        M_dSdUDivergenceMatrix[i].reset( new matrix_Type( M_feSpace->map() ) );
    }

    // Elementary computation and matrix assembling
    for ( UInt iElement(0); iElement < M_physics->data()->numberOfElements(); ++iElement )
    {
        // Update the current element
        M_feSpace->fe().update( M_feSpace->mesh()->edgeList( iElement), UPDATE_DPHI | UPDATE_WDET );

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
OneDimensionalSolver::updateElementalMatrices( const Real& dFdU, const Real& dSdU )
{
    // Set the elementary matrices to 0.
    M_elementalMassMatrix->zero();
    M_elementalStiffnessMatrix->zero();
    M_elementalGradientMatrix->zero();
    M_elementalDivergenceMatrix->zero();

    // Update the mass matrix
    mass( dSdU, *M_elementalMassMatrix, M_feSpace->fe(), 0, 0 );
//    std::cout << "Elem Stiff matrix :" << std::endl;
//    M_elementalMassMatrix->showMe( std::cout );

    // Update the stiffness matrix
    stiff( dFdU, *M_elementalStiffnessMatrix, M_feSpace->fe(), 0 ,0 );
    //AssemblyElemental::stiffness( *M_elementalStiffnessMatrix, M_feSpace->fe(), M_coeffStiff, 0 );
//    std::cout << "Elem Stiff matrix :" << std::endl;
//    M_elementalStiffnessMatrix->showMe( std::cout );

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
    grad( 0, -dFdU, *M_elementalGradientMatrix, M_feSpace->fe(), M_feSpace->fe(), 0, 0 );
//    std::cout << "Elem Grad matrix :" << std::endl;
//    M_elementalGradientMatrix->showMe( std::cout );

    /*! update the divergence matrix
      divergence operator: (transpose of the gradient)
      div_{ij} = \int_{fe} coeff \frac{d \phi_j}{d x} \phi_i

      \note formally this M_elementalDivergenceMatrix is not necessary
      as it is the transpose of the M_elementalGradientMatrix.
      But for the sake of clarity, I prefer to keep it. (low cost!)

      BEWARE : same remarks as grad (see above).
    */
    div( 0, -dSdU, *M_elementalDivergenceMatrix, M_feSpace->fe(), M_feSpace->fe(), 0, 0 );
//    std::cout << "Elem Div matrix :" << std::endl;
//    M_elementalDivergenceMatrix->showMe( std::cout );
}

void
OneDimensionalSolver::matrixAssemble( const UInt& ii, const UInt& jj )
{
    // Assemble the mass matrix
    assembleMatrix( *M_dSdUMassMatrix[ 2*ii + jj ], *M_elementalMassMatrix, M_feSpace->fe(), M_feSpace->dof(), 0, 0, 0, 0 );

    // Assemble the stiffness matrix
    assembleMatrix( *M_dFdUStiffnessMatrix[ 2*ii + jj ], *M_elementalStiffnessMatrix, M_feSpace->fe(), M_feSpace->dof(), 0, 0, 0, 0 );

    // Assemble the gradient matrix
    assembleMatrix( *M_dFdUGradientMatrix[ 2*ii + jj ], *M_elementalGradientMatrix, M_feSpace->fe(), M_feSpace->dof(), 0, 0, 0, 0 );

    // Assemble the divergence matrix
    assembleMatrix( *M_dSdUDivergenceMatrix[ 2*ii + jj ], *M_elementalDivergenceMatrix, M_feSpace->fe(), M_feSpace->dof(), 0, 0, 0, 0 );
}

void
OneDimensionalSolver::applyDirichletBCToMatrix( matrix_Type& matrix )
{
    // Dirichlet BC
    matrix.globalAssemble();
    matrix.diagonalize( 0, 1, 0 );
    matrix.diagonalize( M_physics->data()->numberOfNodes() - 1, 1, 0 );

    //matrix.spy("SystemMatrix");
}

OneDimensionalSolver::vector_Type
OneDimensionalSolver::inertialFlowRateCorrection( const vector_Type& flux )
{
    matrix_Type matrixLHS(M_feSpace->map());
    matrix_Type stiffRHS (M_feSpace->map());

    MatrixElemental elmatMassLHS  (M_feSpace->fe().nbFEDof(), 1, 1);
    MatrixElemental elmatStiffLHS (M_feSpace->fe().nbFEDof(), 1, 1);
    MatrixElemental elmatStiffRHS (M_feSpace->fe().nbFEDof(), 1, 1);

    vector_Type rhs(M_feSpace->map());

    Real coeffMass;
    Real coeffStiff;

    Real m, meanA0;

    matrixLHS *= 0.;
    stiffRHS  *= 0.;

    // Elementary computation and matrix assembling
    for ( UInt iElement(0); iElement < M_physics->data()->numberOfElements(); ++iElement )
    {
        // set the elementary matrices to 0.
        elmatMassLHS. zero();
        elmatStiffLHS.zero();
        elmatStiffRHS.zero();

        coeffMass  = (*M_rhs[0])( iElement ) + (*M_rhs[0])( iElement + 1 );
        coeffMass /= 2;
        coeffMass  = 1./ coeffMass;

        meanA0  = M_physics->data()->area0(iElement) + M_physics->data()->area0(iElement + 1);
        meanA0 /= 2;

        m = M_physics->data()->densityWall()*M_physics->data()->thickness(iElement + 1)/
            ( 2*std::sqrt(4*std::atan(1))*std::sqrt(meanA0) );

        coeffStiff = m/M_physics->data()->densityRho();

        // Update the current element
        M_feSpace->fe().update( M_feSpace->mesh()->edgeList(iElement), UPDATE_DPHI | UPDATE_WDET );

        mass (   coeffMass,  elmatMassLHS,  M_feSpace->fe(), 0, 0 );
        stiff(   coeffStiff, elmatStiffLHS, M_feSpace->fe(), 0, 0 );
        stiff( - coeffStiff, elmatStiffRHS, M_feSpace->fe(), 0, 0 );

        // assemble the mass and grad matrices
        assembleMatrix( matrixLHS, elmatMassLHS, M_feSpace->fe(), M_feSpace->dof(), 0, 0, 0, 0 );
        assembleMatrix( matrixLHS, elmatStiffLHS, M_feSpace->fe(), M_feSpace->dof(), 0, 0, 0, 0 );
        assembleMatrix( stiffRHS, elmatStiffRHS, M_feSpace->fe(), M_feSpace->dof(), 0, 0, 0, 0 );

#ifdef HAVE_LIFEV_DEBUG
        Debug( 6310 ) << "\n\tm = "           << m << "\n\t_coeffMass = "  << coeffMass << "\n\t_coeffStiff = " << coeffStiff << "\n";
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

    M_linearSolver->setMatrix( matrixLHS );
    //@int numIter = M_linearSolver->solveSystem( _rhs, _sol, _matrixLHS, true);

    //std::cout <<" iterations number :  " << numIter << std::endl;

    return sol;
}

OneDimensionalSolver::vector_Type
OneDimensionalSolver::longitudinalFlowRateCorrection()
{
    matrix_Type massLHS(M_feSpace->map());
    matrix_Type massRHS(M_feSpace->map());

    //TriDiagCholesky< Real, matrix_Type, Vector > _tridiagsolver(M_physics->data()->numberOfNodes());

    MatrixElemental elmatMassLHS (M_feSpace->fe().nbFEDof(),1,1);
    MatrixElemental elmatMassRHS (M_feSpace->fe().nbFEDof(),1,1);

    vector_Type rhs(M_feSpace->map());
    // let g = sqrt(A) - sqrt(A0)
    // now f = _d3g_dz3
    vector_Type g(M_feSpace->map());
    vector_Type f(M_feSpace->map());

    //          _g = *M_rhs[0];
    for ( UInt iNode(0); iNode < M_physics->data()->numberOfNodes(); ++iNode )
        g(iNode) = std::sqrt((*M_rhs[0])(iNode)) - std::sqrt(M_physics->data()->area0(iNode));

    UInt iNode;

    Real coeffMassLHS;
    Real coeffMassRHS;

    Real a;
    //    std::ostringstream output;
    massLHS *= 0.;
    massRHS *= 0.;

    Real h( M_physics->data()->length() / static_cast<Real>(M_physics->data()->numberOfElements() - 1) );

    // Elementary computation and matrix assembling
    // Loop on elements
    for ( UInt iElement(0); iElement < M_physics->data()->numberOfElements() ; ++iElement )
    {
        iNode = iElement;

        // set the elementary matrices to 0.
        elmatMassLHS.zero();
        elmatMassRHS.zero();

        // coeff (1/A) (average values over the element)
        coeffMassLHS = (*M_rhs[0])( iNode ) + (*M_rhs[0])( iNode+1 );
        coeffMassLHS /= 2;
        coeffMassLHS = 1./ coeffMassLHS;

        a = M_physics->data()->inertialModulus() / std::sqrt(4*std::atan(1));
        coeffMassRHS = M_physics->data()->dataTime()->timeStep() * a / M_physics->data()->densityRho();

        // backward differentiation when near to the left boundary
        // d3 A (xi)/ dz3 = (1/h^3) * ( -A(xi) + A(xi+3) - 3A(xi+2) + 3A(xi+1) )

        // forward differentiation when near to the right boundary
        // d3 A (xi)/ dz3 = (1/h^3) * ( A(xi) - A(xi-3) + 3A(xi-2) - 3A(xi-1) )

        // central differentiation otherwise
        // d3 A (xi)/ dz3 = (1/h^3) * ( -A(xi-2) + 2A(xi-1) - 2A(xi+1) + A(xi+2) )

#ifdef HAVE_LIFEV_DEBUG
        Debug( 6310 ) << "\ninode = " << iNode << "\n";
#endif
        if (iNode<2)
        { // backward differentiation
            f( iNode ) = - g( iNode ) + g( iNode+3 ) - 3 * g( iNode+2 ) + 3 * g( iNode+1 );
#ifdef HAVE_LIFEV_DEBUG
            Debug( 6310 ) << "\n\tbackward differentiation = " << coeffMassLHS << "\n";
#endif
        }
        else if (iNode>M_feSpace->mesh()->numEdges()-2)
        { // forward differentiation
            f( iNode ) = g( iNode ) - g( iNode-3 ) + 3 * g( iNode-2 ) - 3 * g( iNode-1 );
#ifdef HAVE_LIFEV_DEBUG
            Debug( 6310 ) << "\n\forward differentiation = " << coeffMassLHS << "\n";
#endif
        }
        else
        { // central differentiation
            f( iNode ) = - g( iNode - 2 ) + 2 * g( iNode-1 ) - 2 * g( iNode+1 ) + g( iNode+2 );
#ifdef HAVE_LIFEV_DEBUG
            Debug( 6310 ) << "\n\tcentral differentiation = " << coeffMassLHS << "\n";
#endif
        }

        f(iNode) *= 1 / ( 2 * OneDimensional::pow30(h, 3) );

        // Update the current element
        M_feSpace->fe().update( M_feSpace->mesh()->edgeList(iElement), UPDATE_DPHI | UPDATE_WDET );

        mass( coeffMassLHS, elmatMassLHS, M_feSpace->fe(),0, 0 );
        mass( coeffMassRHS, elmatMassRHS, M_feSpace->fe(),0, 0 );

        // assemble the mass and grad matrices
        //assemb_mat( massLHS, elmatMassLHS, M_feSpace->fe(), M_feSpace->dof() , 0, 0 );
        assembleMatrix( massLHS, elmatMassLHS, M_feSpace->fe(), M_feSpace->dof() , 0, 0, 0, 0 );
        //assemb_mat( massRHS, elmatMassRHS, M_feSpace->fe(), M_feSpace->dof() , 0, 0 );
        assembleMatrix( massRHS, elmatMassRHS, M_feSpace->fe(), M_feSpace->dof() , 0, 0, 0, 0 );

#ifdef HAVE_LIFEV_DEBUG
        Debug( 6310 ) << "\n\t_coeffMassLHS = " << coeffMassLHS << "\n";
        Debug( 6310 ) << "\n\t_coeffMassRHS = " << coeffMassRHS << "\n";
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
