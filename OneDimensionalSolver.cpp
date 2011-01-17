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

#include <lifemc/lifesolver/OneDimensionalSolver.hpp>

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
    M_leftNodeId                (),
    M_leftInternalNodeId        (),
    M_rightNodeId               (),
    M_rightInternalNodeId       (),
    M_elementalMassMatrix       (),
    M_elementalStiffnessMatrix  (),
    M_elementalGradientMatrix   (),
    M_elementalDivergenceMatrix (),
    M_rhs                       (),
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
    M_bcDirLeft                 (),
    M_bcDirRight                ()
{
}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalSolver::buildConstantMatrices()
{
    std::fill( M_dFdUVector.begin(), M_dFdUVector.end(), ublas::zero_vector<Real>( M_feSpace->dim()) );
    std::fill( M_dSdUVector.begin(), M_dSdUVector.end(), ublas::zero_vector<Real>( M_feSpace->dim()) );

    for ( UInt i(0); i < 4; ++i )
    {
        M_dSdUMassMatrix[i].reset( new matrix_Type( M_feSpace->map() ) );
        M_dFdUStiffnessMatrix[i].reset( new matrix_Type( M_feSpace->map() ) );
        M_dFdUGradientMatrix[i].reset( new matrix_Type( M_feSpace->map() ) );
        M_dSdUDivergenceMatrix[i].reset( new matrix_Type( M_feSpace->map() ) );
    }

    // Elementary computation and matrix assembling
    for ( UInt iEdge(1); iEdge <= M_feSpace->mesh()->numEdges(); ++iEdge )
    {
        // set the elementary matrices to 0.
        M_elementalMassMatrix->zero();
        M_elementalGradientMatrix->zero();

        // update the current element
        M_feSpace->fe().update( M_feSpace->mesh()->edgeList( iEdge ), UPDATE_DPHI | UPDATE_WDET );

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
    //matrixPtr_Type matrFull( new matrix_Type( M_feSpace->map() ) );
    //matrFull->insertValueDiagonal( M_physics->data()->mesh()->meanH() );
    //matrFull->GlobalAssemble();
    matrixPtr_Type matrFull( new matrix_Type( *M_homogeneousMassMatrix ));
    updateBCDirichletMatrix( *matrFull );
    M_linearSolver->setMatrix( *matrFull );
}

void
OneDimensionalSolver::setupSolution( solution_Type& solution )
{
    solution["A"].reset( new vector_Type( M_feSpace->map() ) );
    solution["A/A0-1"].reset( new vector_Type( M_feSpace->map() ) );
    solution["Q"].reset( new vector_Type( M_feSpace->map() ) );
    solution["W1"].reset( new vector_Type( M_feSpace->map() ) );
    solution["W2"].reset( new vector_Type( M_feSpace->map() ) );
    solution["P"].reset( new vector_Type( M_feSpace->map() ) );

    // Flux correction with viscoelastic term
    if ( M_physics->data()->viscoelasticWall() )
    {
        solution["Q_visc"].reset( new vector_Type( M_feSpace->map() ) ); // viscoelastic contribution to the flux
        solution["P_visc"].reset( new vector_Type( M_feSpace->map() ) ); // viscoelastic contribution to the pressure
    }

    // flux second derivative
    if ( M_physics->data()->fluxSecondDer() )
        solution["d2Q_dx2"].reset( new vector_Type( M_feSpace->map() ) );

    // correction flux with inertial term
    if ( M_physics->data()->inertialWall() )
        solution["Q_inert"].reset( new vector_Type( M_feSpace->map() ) );

    // correction flux with longitudinal term
    if ( M_physics->data()->longitudinalWall() )
        solution["Q_long"].reset( new vector_Type( M_feSpace->map() ) );

    // Initialize solution to zero
    for ( solutionConstIterator_Type i = solution.begin(); i != solution.end(); ++i )
        *i->second = 0.;
}

void
OneDimensionalSolver::initialize( solution_Type& solution )
{
    for ( UInt inode(M_leftNodeId); inode <= M_rightNodeId ; ++inode )
    {
        (*solution["A"])[inode] = M_physics->data()->area0(inode - 1);
        (*solution["Q"])[inode] = 0;
    }

    // Compute W1 and W2 from A and Q
    computeW1W2( solution );

    // Compute A/A0 from A
    computeAreaRatio( solution );

    // Compute initial pressure (taking into account the viscoelastic wall)
    computePressure( solution, M_physics->data()->dataTime()->timeStep() );
}

void
OneDimensionalSolver::computeW1W2( solution_Type& solution )
{
    for (UInt iElement(0); iElement < M_feSpace->dim() ; ++iElement )
    {
        M_physics->fromUToW( ( *solution["W1"]) [iElement + 1], (*solution["W2"]) [iElement + 1],
                               ( *solution["A"] ) [iElement + 1], (*solution["Q"] ) [iElement + 1], iElement);
    }
}

void
OneDimensionalSolver::computePressure( solution_Type& solution, const Real& timeStep )
{
    for ( UInt i(0); i < M_feSpace->dim() ; ++i )
    {
        ( *solution["P"] ) [i+1] = M_physics->elasticPressure( ( *solution["A"] ) [i+1], i );

        if ( M_physics->data()->viscoelasticWall() )
        {
            // Viscoelastic pressure is not working right now!
            // We need to pass from outside U_prevtime and U_2prevtime
            ( *solution["P_visc"] ) [i+1] = M_physics->viscoelasticPressure( ( *solution["A"]) [i+1],
                                                                            ( *M_UPreviousTime["A"]) [i+1],
                                                                            ( *M_U2PreviousTime["A"]) [i+1], timeStep, i );
            ( *solution["P"] ) [i+1] += ( *solution["P_visc"] ) [i+1];
        }
    }
}

void
OneDimensionalSolver::computeAreaRatio( solution_Type& solution )
{
    for ( UInt inode(M_leftNodeId); inode <= M_rightNodeId ; ++inode )
        ( *solution["A/A0-1"] ) [inode] = ( *solution["A"] ) [inode] / M_physics->data()->area0(inode - 1) - 1;
}

void
OneDimensionalSolver::computeArea( solution_Type& solution )
{
    for ( UInt inode(M_leftNodeId); inode <= M_rightNodeId ; ++inode )
        ( *solution["A"] ) [inode] = ( (*solution["A/A0-1"] ) [inode] + 1 ) * M_physics->data()->area0(inode - 1);
}

void
OneDimensionalSolver::updateRHS( const solution_Type& solution, const Real& timeStep )
{
    updatedFdU( solution );   // Update the vector containing the values of the flux at the nodes and its jacobian
    updatedSdU( solution ); // Update the vector containing the values of the source term at the nodes and its jacobian
    updateMatrices();            // Update the matrices for the non-linear terms

    // Taylor-Galerkin scheme: (explicit, U = [U1,U2]^T, with U1=A, U2=Q )
    // (Un+1, phi) =          (               Un,     phi     )-> massFactor^{-1} * Un+1 = mass * U
    //             + dt     * (           Fh(Un),     dphi/dz )->            grad * F(U)
    //             - dt^2/2 * (diffFh(Un) Sh(Un),     dphi/dz )-> gradDiffFlux(U) * S(U)
    //             + dt^2/2 * (diffSh(Un) dFh/dz(Un), phi     )->   divDiffSrc(U) * F(U)
    //             - dt^2/2 * (diffFh(Un) dFh/dz(Un), dphi/dz )->stiffDiffFlux(U) * F(U)
    //             - dt     * (           Sh(Un),     phi     )->            mass * S(U)
    //             + dt^2/2 * (diffSh(Un) Sh(Un),     phi     )->  massDiffSrc(U) * S(U)

    Real dt2over2 = timeStep * timeStep * 0.5;

    // rhs = mass * Un
    *M_rhs[0]  = ( *M_homogeneousMassMatrix ) * *solution.find("A")->second;
    *M_rhs[1]  = ( *M_homogeneousMassMatrix ) * *solution.find("Q")->second;
    for ( UInt i(0); i < 2; ++i )
    {
        // rhs = rhs + dt * grad * F(Un)
        *M_rhs[i] += ( *M_homogeneousGradientMatrix ) * ( timeStep * *M_fluxVector[i] );

        // rhs = rhs - dt * mass * S(Un)
        *M_rhs[i] += ( *M_homogeneousMassMatrix ) * ( - timeStep * *M_sourceVector[i] );

        for ( UInt j(0); j < 2; ++j )
        {
            // rhs = rhs - dt^2/2 * gradDiffFlux * S(Un)
            M_dFdUGradientMatrix[2*i + j]->globalAssemble();
            *M_rhs[i] += *M_dFdUGradientMatrix[2*i + j] * ( -dt2over2 * *M_sourceVector[j] );

            // rhs = rhs + dt^2/2 * divDiffSrc * F(Un)
            M_dSdUDivergenceMatrix[2*i + j]->globalAssemble();
            *M_rhs[i] += *M_dSdUDivergenceMatrix[2*i + j] * ( dt2over2 * *M_fluxVector[j] );

            // rhs = rhs - dt^2/2 * stiffDiffFlux * F(Un)
            M_dFdUStiffnessMatrix[2*i + j]->globalAssemble();
            *M_rhs[i] += *M_dFdUStiffnessMatrix[2*i + j]*( -dt2over2 * *M_fluxVector[j] );

            // rhs = rhs + dt^2/2 * massDiffSrc * S(Un)
            M_dSdUMassMatrix[2*i + j]->globalAssemble();
            *M_rhs[i] += *M_dSdUMassMatrix[2*i + j] * ( dt2over2 * *M_sourceVector[j] );
        }
    }
}

void
OneDimensionalSolver::iterate( OneDimensionalBCHandler& bcH, solution_Type& solution, const Real& time, const Real& timeStep )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 6310 ) << "[timeAdvance] \tcompute BC\n";
#endif
    // Compute BC
    bcH.applyBC( time, timeStep, solution, M_flux, M_bcDirLeft, M_bcDirRight );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 6310 ) << "[timeAdvance] \tcompute BC dirichlet vector\n";
#endif
    // Apply BC to RHS
    updateBCDirichletVector();

    // Compute A^n+1
    vector_Type area( *M_rhs[0] );
    M_linearSolver->solveSystem( *M_rhs[0], area, M_homogeneousMassMatrix );

    // Compute Q^n+1
    vector_Type flowrate( *M_rhs[1] );
    M_linearSolver->solveSystem( *M_rhs[1], flowrate, M_homogeneousMassMatrix );

    // Correct flux with inertial, viscoelastic and longitudinal terms
    if ( M_physics->data()->inertialWall() )
    {
        *solution["Q_inert"] = inertialFluxCorrection( flowrate );
        flowrate += *solution["Q_inert"];
    }

    if ( M_physics->data()->viscoelasticWall() )
    {
        *solution["Q_visc"] = viscoelasticFluxCorrection( flowrate, timeStep );
        flowrate += *solution["Q_visc"];
    }

    if ( M_physics->data()->longitudinalWall() )
    {
        *solution["Q_long"] = longitudinalFluxCorrection();
        flowrate += *solution["Q_long"];
    }

    // compute L2 projection of d2Q_dx2
    //    if( M_physics->data()->fluxSecondDer() )
    //      M_d2_U2_dx2 = _compute_d2Q_dx2( flowrate );

    // Update the solution container
    *solution["A"] = area;
    *solution["Q"] = flowrate;

    // Compute W1 and W2 from A and Q
    computeW1W2( solution );

    // Compute A/A0 from A
    computeAreaRatio( solution );

    // Update the pressure (taking into account the viscoelastic wall)
    computePressure( solution, timeStep );
}

Real
OneDimensionalSolver::computeCFL( const solution_Type& solution, const Real& timeStep ) const
{
    Real lambdaMax = 0.;

    container2D_Type eigenvalues;
    container2D_Type leftEigenvector1;
    container2D_Type leftEigenvector2;

    for ( UInt inode(0); inode < M_feSpace->dim() ; ++inode )
    {
        // compute the eigenvalues at node
        M_flux->eigenValuesEigenVectors( ( *solution.find("A")->second ) ( inode + 1 ),
                                           ( *solution.find("Q")->second ) ( inode + 1 ),
                                             eigenvalues, leftEigenvector1, leftEigenvector2, inode );

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
OneDimensionalSolver::postProcess( const solution_Type& solution )
{
    std::ofstream outfile;
    for ( solutionConstIterator_Type i = solution.begin(); i != solution.end(); ++i )
    {
        std::string file = M_physics->data()->postprocessingDirectory() + "/" + M_physics->data()->postprocessingFile() + "_" + i->first + ".m";
        outfile.open( file.c_str(), std::ios::app );

        for ( UInt ii(0); ii < M_feSpace->dim(); ++ii )
            outfile << (*i->second)(ii + 1) << " ";

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

    //Id of left and right bc nodes
    M_leftNodeId          = 1;
    M_leftInternalNodeId  = M_leftNodeId + 1;
    M_rightNodeId         = M_feSpace->dim();
    M_rightInternalNodeId = M_rightNodeId - 1;

    //Elementary Matrices
    M_elementalMassMatrix.reset ( new MatrixElemental( M_feSpace->fe().nbFEDof(), 1, 1 ) );
    M_elementalStiffnessMatrix.reset( new MatrixElemental( M_feSpace->fe().nbFEDof(), 1, 1 ) );
    M_elementalGradientMatrix.reset ( new MatrixElemental( M_feSpace->fe().nbFEDof(), 1, 1 ) );
    M_elementalDivergenceMatrix.reset  ( new MatrixElemental( M_feSpace->fe().nbFEDof(), 1, 1 ) );

    //Vectors
    for ( UInt i(0) ; i < 2 ; ++i )
    {
        M_rhs[i].reset( new vector_Type( M_feSpace->map() ) );
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
OneDimensionalSolver::setBCValuesLeft( const Real& bcL1, const Real& bcL2 )
{
    M_bcDirLeft[0] = bcL1;
    M_bcDirLeft[1] = bcL2;
}

void
OneDimensionalSolver::setBCValuesRight( const Real& bcR1, const Real& bcR2 )
{
    M_bcDirRight[0] = bcR1;
    M_bcDirRight[1] = bcR2;
}

// ===================================================
// Get Methods
// ===================================================
OneDimensionalSolver::container2D_Type
OneDimensionalSolver::bcValuesLeft( const solution_Type& solution ) const
{
    container2D_Type temp;

    temp[0] = (*solution.find("A")->second)( leftNodeId() );
    temp[1] = (*solution.find("Q")->second)( leftNodeId() );

    return temp;
}

OneDimensionalSolver::container2D_Type
OneDimensionalSolver::bcValuesInternalLeft( const solution_Type& solution ) const
{
    container2D_Type temp;

    temp[0] = (*solution.find("A")->second)( leftInternalNodeId() );
    temp[1] = (*solution.find("Q")->second)( leftInternalNodeId() );

    return temp;
}

OneDimensionalSolver::container2D_Type
OneDimensionalSolver::bcValuesRight( const solution_Type& solution ) const
{
    container2D_Type temp;

    temp[0] = (*solution.find("A")->second)( rightNodeId() );
    temp[1] = (*solution.find("Q")->second)( rightNodeId() );

    return temp;
}

OneDimensionalSolver::container2D_Type
OneDimensionalSolver::bcValuesInternalRight( const solution_Type& solution ) const
{
    container2D_Type temp;

    temp[0] = (*solution.find("A")->second)( rightInternalNodeId() );
    temp[1] = (*solution.find("Q")->second)( rightInternalNodeId() );

    return temp;
}

Real
OneDimensionalSolver::boundaryValue( const solution_Type& solution, const bcType_Type& bcType, const bcSide_Type& bcSide ) const
{
    UInt boundaryDof;

    switch ( bcSide )
    {
    case OneDimensional::left:
        boundaryDof = 1;
        break;
    case OneDimensional::right:

        boundaryDof = M_flux->physics()->data()->numberOfElements() + 1;

        boundaryDof = M_flux->physics()->data()->numberOfElements() + 1;

        break;
    default:
        std::cout << "Warning: bcSide \"" << bcSide << "\" not available!" << std::endl;
        return 0;
        break;
    }

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
    UInt boundaryDof;

    switch ( bcSide )
    {
    case OneDimensional::left:

        boundaryDof = 0;

        break;

    case OneDimensional::right:

        boundaryDof = M_flux->physics()->data()->numberOfElements();
        boundaryDof = M_flux->physics()->data()->numberOfElements();

        break;

    default:

        std::cout << "Warning: bcSide \"" << bcSide << "\" not available!" << std::endl;

        return;
    }

    M_flux->eigenValuesEigenVectors( (*solution.find("A")->second)( boundaryDof + 1 ),
                                     (*solution.find("Q")->second)( boundaryDof + 1 ),
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

    for ( UInt iNode(0); iNode < M_feSpace->dim() ; ++iNode )
    {
        Ai = (*solution.find("A")->second)( iNode + 1 );
        Qi = (*solution.find("Q")->second)( iNode + 1 );

        (*M_fluxVector[0])( iNode + 1 ) = M_flux->flux( Ai, Qi, 1, iNode );
        (*M_fluxVector[1])( iNode + 1 ) = M_flux->flux( Ai, Qi, 2, iNode );
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

    for ( UInt iElement(0); iElement < M_feSpace->dim() - 1; ++iElement )
    {
        // for P1Seg and appropriate mesh only!
        Aii   = (*solution.find("A")->second)( iElement + 1 );
        Qii   = (*solution.find("Q")->second)( iElement + 1 );

        Aiip1 = (*solution.find("A")->second)( iElement + 2 );
        Qiip1 = (*solution.find("Q")->second)( iElement + 2 );

        for ( UInt ii=1; ii<3; ++ii )
        {
            for ( UInt jj=1; jj<3; ++jj )
            {
                tmp  = M_flux->dFdU(   Aii,   Qii, ii, jj, iElement );     // left node of current element
                tmp += M_flux->dFdU( Aiip1, Qiip1, ii, jj, iElement + 1 ); // right node of current element

                M_dFdUVector[ 2*(ii - 1) + jj - 1 ]( iElement ) = 0.5 * tmp;
            }
        }
    }
}

void
OneDimensionalSolver::updateSource( const solution_Type& solution )
{
    Real Ai, Qi;

    for ( UInt iNode(0); iNode < M_feSpace->dim() ; ++iNode )
    {
        Ai = (*solution.find("A")->second)( iNode + 1 );
        Qi = (*solution.find("Q")->second)( iNode + 1 );

        (*M_sourceVector[0])( iNode + 1 ) = M_source->source( Ai, Qi, 1, iNode );
        (*M_sourceVector[1])( iNode + 1 ) = M_source->source( Ai, Qi, 2, iNode );
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

    for ( UInt iElement(0); iElement < M_feSpace->dim() - 1 ; ++iElement )
    {
        // for P1Seg and appropriate mesh only!
        Aii   = (*solution.find("A")->second)( iElement + 1);
        Qii   = (*solution.find("Q")->second)( iElement + 1);
        Aiip1 = (*solution.find("A")->second)( iElement + 2 );
        Qiip1 = (*solution.find("Q")->second)( iElement + 2 );

        for ( UInt ii=1; ii<3; ++ii )
        {
            for ( UInt jj=1; jj<3; ++jj )
            {
                tmp =  M_source->dSdU(   Aii,   Qii, ii, jj, iElement );     // left node of current element
                tmp += M_source->dSdU( Aiip1, Qiip1, ii, jj, iElement + 1 ); // right node of current element

                M_dSdUVector[ 2*(ii - 1) + jj - 1 ]( iElement ) = 0.5 * tmp;
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
    for (UInt iEdge = 1; iEdge <= M_feSpace->mesh()->numEdges(); ++iEdge)
    {
        // Update the current element
        M_feSpace->fe().update( M_feSpace->mesh()->edgeList(iEdge), UPDATE_DPHI | UPDATE_WDET );

        for ( UInt ii(1); ii <= 2; ++ii )
        {
            for ( UInt jj(1); jj <= 2; ++jj )
            {
                // Update the elemental matrices
                updateElementalMatrices( M_dFdUVector[ 2*(ii-1) + jj-1 ]( iEdge - 1 ), M_dSdUVector[ 2*(ii-1) + jj-1 ]( iEdge - 1 ) );

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
    assembleMatrix( *M_dSdUMassMatrix[ 2*(ii-1) + jj-1 ], *M_elementalMassMatrix, M_feSpace->fe(), M_feSpace->dof(), 0, 0, 0, 0 );

    // Assemble the stiffness matrix
    assembleMatrix( *M_dFdUStiffnessMatrix[ 2*(ii-1) + jj-1 ], *M_elementalStiffnessMatrix, M_feSpace->fe(), M_feSpace->dof(), 0, 0, 0, 0 );

    // Assemble the gradient matrix
    assembleMatrix( *M_dFdUGradientMatrix[ 2*(ii-1) + jj-1 ], *M_elementalGradientMatrix, M_feSpace->fe(), M_feSpace->dof(), 0, 0, 0, 0 );

    // Assemble the divergence matrix
    assembleMatrix( *M_dSdUDivergenceMatrix[ 2*(ii-1) + jj-1 ], *M_elementalDivergenceMatrix, M_feSpace->fe(), M_feSpace->dof(), 0, 0, 0, 0 );
}

void
OneDimensionalSolver::updateBCDirichletVector()
{
#ifdef HAVE_LIFEV_DEBUG
    Debug( 6310 ) << "[updateBCDirichletVector] \t firstDof = " << M_leftNodeId
    << " lastDof = " << M_rightNodeId << ";\n";

    Debug( 6310 ) << "[updateBCDirichletVector] \t bcDirLeft[0] = " << M_bcDirLeft[0]
    << " bcDirLeft[1] = " << M_bcDirLeft[1] << ";\n";

    Debug( 6310 ) << "[updateBCDirichletVector] \t bcDirRight[0] = " << M_bcDirRight[0]
    << " bcDirRight[1] = " << M_bcDirRight[1] << ";\n";
#endif

    // first row modified
    (*M_rhs[0])( M_leftNodeId ) = M_bcDirLeft[0];
    (*M_rhs[1])( M_leftNodeId ) = M_bcDirLeft[1];

    // last row modified
    (*M_rhs[0])( M_rightNodeId ) = M_bcDirRight[0];
    (*M_rhs[1])( M_rightNodeId ) = M_bcDirRight[1];
}

void
OneDimensionalSolver::updateBCDirichletMatrix( matrix_Type& mat )
{
    mat.diagonalize(M_leftNodeId - 1, 1, 0);
    mat.diagonalize(M_rightNodeId - 1, 1, 0);
}

OneDimensionalSolver::vector_Type
OneDimensionalSolver::inertialFluxCorrection( const vector_Type& flux )
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
    for (UInt iEdge = 1; iEdge <= M_feSpace->mesh()->numEdges(); ++iEdge) // Loop on elements
    {
        // set the elementary matrices to 0.
        elmatMassLHS. zero();
        elmatStiffLHS.zero();
        elmatStiffRHS.zero();

        coeffMass  = (*M_rhs[0])( iEdge - 1 ) + (*M_rhs[0])( iEdge );
        coeffMass /= 2;
        coeffMass  = 1./ coeffMass;

        meanA0  = M_physics->data()->area0(iEdge - 1) + M_physics->data()->area0(iEdge);
        meanA0 /= 2;

        m = M_physics->data()->densityWall()*M_physics->data()->thickness(iEdge)/
            ( 2*std::sqrt(4*std::atan(1))*std::sqrt(meanA0) );

        coeffStiff = m/M_physics->data()->densityRho();

        // Update the current element
        M_feSpace->fe().update( M_feSpace->mesh()->edgeList(iEdge), UPDATE_DPHI | UPDATE_WDET );

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

    UInt firstDof = 1;
    UInt lastDof  = rhs.size();

    // symmetric treatment (cholesky can be used)
    // first row modified (Dirichlet)
    rhs( firstDof ) = 0.;
    // second row modified (for symmetry)
    //  _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    // last row modified (Dirichlet)
    rhs( lastDof ) = 0.;
    // penultimate row modified (for symmetry)
    //  _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    updateBCDirichletMatrix( matrixLHS );

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
OneDimensionalSolver::viscoelasticFluxCorrection( const vector_Type& flux, const Real& timeStep )
{
    matrix_Type matrixLHS(M_feSpace->map());
    matrix_Type stiffRHS (M_feSpace->map());

    //TriDiagCholesky< Real, matrix_Type, Vector > _tridiagsolver(M_feSpace->dim());

    MatrixElemental elmatMassLHS  (M_feSpace->fe().nbFEDof(),1,1);
    MatrixElemental elmatStiffLHS (M_feSpace->fe().nbFEDof(),1,1);
    MatrixElemental elmatStiffRHS (M_feSpace->fe().nbFEDof(),1,1);

    vector_Type rhs(M_feSpace->map());

    Real coeffMass;
    Real coeffStiff;

    Real gamma, meanA0(1.);

    matrixLHS *= 0.;
    stiffRHS  *= 0.;

    // Elementary computation and matrix assembling
    // Loop on elements
    for (UInt iEdge = 1; iEdge <= M_feSpace->mesh()->numEdges(); ++iEdge)
    {
        // set the elementary matrices to 0.
        elmatMassLHS.zero();
        elmatStiffLHS.zero();
        elmatStiffRHS.zero();

        // this comes from the exact derivation of generalized string model
        // + Voigt viscoelasticity
        coeffMass = (*M_rhs[0])( iEdge - 1 ) + (*M_rhs[0])( iEdge );
        coeffMass *= 0.5;
        coeffMass = 1./std::sqrt( coeffMass);

        if (M_physics->data()->linearizeStringModel())
        {
            // this is to recover the linearized version (_coeffMass = 1/A)
            coeffMass *= coeffMass;

            meanA0  = M_physics->data()->area0( iEdge - 1 ) + M_physics->data()->area0( iEdge );
            meanA0 *= 0.5;

            if (M_physics->data()->linearizeEquations())
            {
                // when using linearized equations, A \simeq A0
                coeffMass = 1./ meanA0;
            }
        }
        gamma = M_physics->data()->viscoelasticModulus() / ( 2 * std::sqrt(4*std::atan(1)) );
        gamma *= 1 / std::sqrt( meanA0 );
        coeffStiff = timeStep * gamma / M_physics->data()->densityRho();

        // Update the current element
        M_feSpace->fe().update( M_feSpace->mesh()->edgeList(iEdge), UPDATE_DPHI | UPDATE_WDET );

        mass (   coeffMass,  elmatMassLHS,  M_feSpace->fe(),0, 0 );
        stiff(   coeffStiff, elmatStiffLHS, M_feSpace->fe(),0, 0 );
        stiff( - coeffStiff, elmatStiffRHS, M_feSpace->fe(),0, 0 );

        // assemble the mass and grad matrices
        assembleMatrix( matrixLHS, elmatMassLHS, M_feSpace->fe(), M_feSpace->dof() , 0, 0, 0, 0 );
        assembleMatrix( matrixLHS, elmatStiffLHS, M_feSpace->fe(), M_feSpace->dof() , 0, 0, 0, 0 );
        assembleMatrix( stiffRHS, elmatStiffRHS, M_feSpace->fe(), M_feSpace->dof() , 0, 0, 0, 0 );

#ifdef HAVE_LIFEV_DEBUG
        Debug( 6310 ) << "\n\tgamma = " << gamma << "\n\t_coeffMass = " << coeffMass << "\n\t_coeffStiff = " << coeffStiff << "\n";
#endif
    } // end loop on elements

    // update rhs
    // rhs = _stiffRHS * rhs
    rhs = stiffRHS * flux;

    // NOTE I should add to rhs the value of boundary integral ("neumann-like")
    // BUT it's useless since I'm going to impose dirichlet conditions on all boundaries

    UInt firstDof = 1;
    UInt lastDof  = rhs.size();

    // symmetric treatment (cholesky can be used)
    // first row modified (Dirichlet)
    rhs( firstDof ) = 0.;
    // second row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    // _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    // last row modified (Dirichlet)
    rhs( lastDof ) = 0.;
    // penultimate row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    updateBCDirichletMatrix(matrixLHS);

    vector_Type sol(rhs);

    M_linearSolver->setMatrix(matrixLHS);
    //int numIter = M_linearSolver->solveSystem( _rhs, _sol, _matrixLHS, true);

    return sol;
}

OneDimensionalSolver::vector_Type
OneDimensionalSolver::longitudinalFluxCorrection()
{
    matrix_Type massLHS(M_feSpace->map());
    matrix_Type massRHS(M_feSpace->map());

    //TriDiagCholesky< Real, matrix_Type, Vector > _tridiagsolver(M_feSpace->dim());

    MatrixElemental elmatMassLHS (M_feSpace->fe().nbFEDof(),1,1);
    MatrixElemental elmatMassRHS (M_feSpace->fe().nbFEDof(),1,1);

    vector_Type rhs(M_feSpace->map());
    // let g = sqrt(A) - sqrt(A0)
    // now f = _d3g_dz3
    vector_Type g(M_feSpace->map());
    vector_Type f(M_feSpace->map());

    //          _g = *M_rhs[0];
    for ( UInt i=0; i<M_feSpace->dim(); ++i )
        g(i) = std::sqrt((*M_rhs[0])(i)) - std::sqrt(M_physics->data()->area0(i));

    UInt inode;

    Real coeffMassLHS;
    Real coeffMassRHS;

    Real a;
    //    std::ostringstream output;
    massLHS *= 0.;
    massRHS *= 0.;

    Real h( M_physics->data()->length() / static_cast<Real>(M_physics->data()->numberOfElements() - 1) );

    // Elementary computation and matrix assembling
    // Loop on elements
    for (UInt iEdge = 1; iEdge <= M_feSpace->mesh()->numEdges(); ++iEdge)
    {
        inode = iEdge - 1;

        // set the elementary matrices to 0.
        elmatMassLHS.zero();
        elmatMassRHS.zero();

        // coeff (1/A) (average values over the element)
        coeffMassLHS = (*M_rhs[0])( inode ) + (*M_rhs[0])( inode+1 );
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
        Debug( 6310 ) << "\ninode = " << inode << "\n";
#endif
        if (inode<2)
        { // backward differentiation
            f( inode ) = - g( inode ) + g( inode+3 ) - 3 * g( inode+2 ) + 3 * g( inode+1 );
#ifdef HAVE_LIFEV_DEBUG
            Debug( 6310 ) << "\n\tbackward differentiation = " << coeffMassLHS << "\n";
#endif
        }
        else if (inode>M_feSpace->mesh()->numEdges()-2)
        { // forward differentiation
            f( inode ) = g( inode ) - g( inode-3 ) + 3 * g( inode-2 ) - 3 * g( inode-1 );
#ifdef HAVE_LIFEV_DEBUG
            Debug( 6310 ) << "\n\forward differentiation = " << coeffMassLHS << "\n";
#endif
        }
        else
        { // central differentiation
            f( inode ) = - g( inode - 2 ) + 2 * g( inode-1 ) - 2 * g( inode+1 ) + g( inode+2 );
#ifdef HAVE_LIFEV_DEBUG
            Debug( 6310 ) << "\n\tcentral differentiation = " << coeffMassLHS << "\n";
#endif
        }

        f(inode) *= 1 / ( 2 * std::pow(h,3) );

        // Update the current element
        M_feSpace->fe().update( M_feSpace->mesh()->edgeList(iEdge), UPDATE_DPHI | UPDATE_WDET );

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

    UInt firstDof = 1;
    UInt lastDof  = rhs.size();

    // symmetric treatment (cholesky can be used)
    // first row modified (Dirichlet)
    rhs( firstDof ) = 0.;
    // second row modified (for symmetry)
    //    _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    // last row modified (Dirichlet)
    rhs( lastDof ) = 0.;
    // penultimate row modified (for symmetry)
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    updateBCDirichletMatrix( massLHS);

    //@_tridiagsolver.Factor( _massLHS );

    // cholesky or lapack lu solve
    // solve the system: rhs1 = massFactor^{-1} * rhs1
    //@_tridiagsolver.Solve( _massLHS, _rhs );

    return rhs;
}

/*
template< class Params, class Flux, class Source >
ScalVec
OneDimensionalSolver::_compute_d2Q_dx2( const ScalVec& flux )
{
    matrix_Type _massLHS (M_feSpace->dim());
    matrix_Type _stiffRHS(M_feSpace->dim());

    TriDiagCholesky< Real, matrix_Type, Vector > _tridiagsolver(M_feSpace->dim());

<<<<<<< HEAD
    MatrixElemental _elmatMassLHS (M_FESpace -> fe().nbFEDof(),1,1);
    MatrixElemental _elmatStiffRHS (M_FESpace -> fe().nbFEDof(),1,1);
=======
    ElemMat _elmatMassLHS (M_feSpace->fe().nbFEDof(),1,1);
    ElemMat _elmatStiffRHS (M_feSpace->fe().nbFEDof(),1,1);
>>>>>>> Renaming some members + some optimizations.

    ScalVec _rhs(M_feSpace->dim());

    _massLHS  *= 0.;
    _stiffRHS *= 0.;

    // Elementary computation and matrix assembling
    // Loop on elements
    for(UInt iEdge = 1; iEdge <= M_feSpace->mesh()->numEdges(); ++iEdge){

        // set the elementary matrices to 0.
        _elmatMassLHS *= 0.;
        _elmatStiffRHS *= 0.;

        // update the current element
        M_feSpace->fe().update( M_feSpace->mesh()->edgeList( iEdge ), UPDATE_DPHI | UPDATE_WDET );
        //std::cout << M_feSpace->fe().currentId() << std::endl;

        mass( 1., _elmatMassLHS, M_feSpace->fe(),0, 0 );
        stiff( -1., _elmatStiffRHS, M_feSpace->fe(),0, 0 );

        // assemble the mass and grad matrices
        assemb_mat( _massLHS, _elmatMassLHS, M_feSpace->fe(), M_feSpace->dof() , 0, 0 );
        assemb_mat( _stiffRHS, _elmatStiffRHS, M_feSpace->fe(), M_feSpace->dof() , 0, 0 );

    } // end loop on elements

    // update rhs
    _rhs = _stiffRHS*flux;

    UInt firstDof = 0;
    UInt lastDof  = _rhs.size()-1;

    // symmetric treatment (cholesky can be used)
    // first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    // second row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    // _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    // last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    // penultimate row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    updateBCDirichletMatrix(_massLHS);

//     Debug( 6310 ) << "\n[_compute_d2Q_dx2] _massLHS.Diag():\n";
//     for(int i=0; i<_massLHS.OrderMatrix(); ++i)
//         Debug( 6310 ) << "\t" << _massLHS.Diag()(i);
//     Debug( 6310 ) << "\n[_compute_d2Q_dx2] _massLHS.UpDiag():\n";
//     for(int i=0; i<_massLHS.OrderMatrix()-1; ++i)
//         Debug( 6310 ) << "\t" << _massLHS.UpDiag()(i);
//     Debug( 6310 ) << "\n[_compute_d2Q_dx2] _massLHS.LowDiag():\n";
//     for(int i=0; i<_massLHS.OrderMatrix()-1; ++i)
//         Debug( 6310 ) << "\t" << _massLHS.LowDiag()(i);
//     Debug( 6310 ) << "\n";

//     _tridiagsolver.Factor( _massLHS );

//     Debug( 6310 ) << "\n[_compute_d2Q_dx2] solving with rhs:\n";
//     for(int i=0; i<_massLHS.OrderMatrix(); ++i)
//         Debug( 6310 ) << "\t" << _rhs(i);
//     Debug( 6310 ) << "\n";

    // cholesky or lapack lu solve
    // solve the system: rhs1 = massFactor^{-1} * rhs1
    _tridiagsolver.Solve( _massLHS, _rhs );

    return _rhs;
}
*/

}
