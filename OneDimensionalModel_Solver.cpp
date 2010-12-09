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
    @file
    @brief File containing a solver class for the 1D model.

    @version 1.0
    @date 01-10-2006
    @author Vincent Martin
    @author Tiziano Passerini
    @author Lucia Mirabella

    @version 2.0
    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @date 01-08-2009

    @version 2.1
    @date 21-04-2010
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>

    @contributors Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>

    @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/OneDimensionalModel_Solver.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_Solver::OneDimensionalModel_Solver():
    M_physics               (),
    M_flux                  (),
    M_source                (),
    M_FESpace               (),
    M_comm                  (),
    M_leftNodeId            (),
    M_leftInternalNodeId    (),
    M_rightNodeId           (),
    M_rightInternalNodeId   (),
    M_coeffMass             (),
    M_coeffStiff            (),
    M_coeffGrad             (),
    M_coeffDiv              (),
    M_elmatMass             (),
    M_elmatStiff            (),
    M_elmatGrad             (),
    M_elmatDiv              (),
    M_rhs                   (),
    M_fluxVector            (),
    M_diffFlux              ( 4 ),
    M_sourceVector          (),
    M_diffSrc               ( 4 ),
    M_massMatrix            (),
    M_massMatrixDiffSrc     ( 4  ),
    M_stiffMatrixDiffFlux   ( 4  ),
    M_gradMatrix            (),
    M_gradMatrixDiffFlux    ( 4  ),
    M_divMatrixDiffSrc      ( 4 ),
    M_linearSolver          (),
    M_bcDirLeft             (),
    M_bcDirRight            ()
{
}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalModel_Solver::buildConstantMatrices()
{
    std::fill( M_diffFlux.begin(), M_diffFlux.end(), ublas::zero_vector<Real>( M_FESpace -> dim()) );
    std::fill( M_diffSrc.begin(),  M_diffSrc.end(),  ublas::zero_vector<Real>( M_FESpace -> dim()) );

    for ( UInt i = 0; i < 4; ++i )
    {
        M_massMatrixDiffSrc[i].reset( new matrix_Type( M_FESpace -> map() ));
        M_stiffMatrixDiffFlux[i].reset( new matrix_Type( M_FESpace -> map() ));
        M_gradMatrixDiffFlux[i].reset( new matrix_Type( M_FESpace -> map() ));
        M_divMatrixDiffSrc[i].reset( new matrix_Type( M_FESpace -> map() ));
    }

    // set the coeff to 1.
    M_coeffMass = 1.;
    M_coeffGrad = 1.;

    // Elementary computation and matrix assembling
    for ( UInt iedge(1); iedge <= M_FESpace -> mesh() -> numEdges(); ++iedge )
    {
        // set the elementary matrices to 0.
        M_elmatMass -> zero();
        M_elmatGrad -> zero();

        // update the current element
        M_FESpace -> fe().updateFirstDerivQuadPt( M_FESpace -> mesh() -> edgeList(iedge) );

        // update the mass and grad matrices
        mass(      M_coeffMass, *M_elmatMass, M_FESpace -> fe(),0, 0 );
        grad( 0 , -M_coeffGrad, *M_elmatGrad, M_FESpace -> fe(), M_FESpace -> fe(), 0, 0 );

        // assemble the mass and grad matrices
        assemb_mat( *M_massMatrix, *M_elmatMass, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );
        assemb_mat( *M_gradMatrix, *M_elmatGrad, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );
    }

    // Dirichlet boundary conditions set in the mass matrix
    M_massMatrix -> GlobalAssemble();
    M_gradMatrix -> GlobalAssemble();

    // In the classical case the linear system use a mass matrix (with Dirichlet BC)
    //matrixPtr_Type matrFull( new matrix_Type( M_FESpace -> map() ) );
    //matrFull -> insertValueDiagonal( M_physics -> Data() -> mesh() -> meanH() );
    //matrFull -> GlobalAssemble();
    matrixPtr_Type matrFull( new matrix_Type( *M_massMatrix ));
    updateBCDirichletMatrix( *matrFull );
    M_linearSolver -> setMatrix( *matrFull );
}

void
OneDimensionalModel_Solver::setupSolution( solution_Type& solution )
{
    solution["A"].reset( new vector_Type( M_FESpace -> map() ) );
    solution["A/A0-1"].reset( new vector_Type( M_FESpace -> map() ) );
    solution["Q"].reset( new vector_Type( M_FESpace -> map() ) );
    solution["W1"].reset( new vector_Type( M_FESpace -> map() ) );
    solution["W2"].reset( new vector_Type( M_FESpace -> map() ) );
    solution["P"].reset( new vector_Type( M_FESpace -> map() ) );

    // Flux correction with viscoelastic term
    if ( M_physics -> Data() -> viscoelasticWall() )
    {
        solution["Q_visc"].reset( new vector_Type( M_FESpace -> map() ) ); // viscoelastic contribution to the flux
        solution["P_visc"].reset( new vector_Type( M_FESpace -> map() ) ); // viscoelastic contribution to the pressure
    }

    // flux second derivative
    if ( M_physics -> Data() -> fluxSecondDer() )
        solution["d2Q_dx2"].reset( new vector_Type( M_FESpace -> map() ) );

    // correction flux with inertial term
    if ( M_physics -> Data() -> inertialWall() )
        solution["Q_inert"].reset( new vector_Type( M_FESpace -> map() ) );

    // correction flux with longitudinal term
    if ( M_physics -> Data() -> longitudinalWall() )
        solution["Q_long"].reset( new vector_Type( M_FESpace -> map() ) );

    // Initialize solution to zero
    for ( solutionConstIterator_Type i = solution.begin(); i != solution.end(); ++i )
        *i -> second = 0.;
}

void
OneDimensionalModel_Solver::initialize( solution_Type& solution )
{
    for ( UInt inode(M_leftNodeId); inode <= M_rightNodeId ; ++inode )
    {
        (*solution["A"])[inode] = M_physics -> Data() -> Area0(inode - 1);
        (*solution["Q"])[inode] = 0;
    }

    // Compute W1 and W2 from A and Q
    computeW1W2( solution );

    // Compute A/A0 from A
    computeAreaRatio( solution );

    // Compute initial pressure (taking into account the viscoelastic wall)
    computePressure( solution, M_physics -> Data() -> dataTime() -> getTimeStep() );
}

void
OneDimensionalModel_Solver::computeW1W2( solution_Type& solution )
{
    for (UInt ielem = 0; ielem < M_FESpace -> dim() ; ++ielem )
    {
        M_physics -> W_from_U( ( *solution["W1"]) [ielem + 1], (*solution["W2"]) [ielem + 1],
                               ( *solution["A"] ) [ielem + 1], (*solution["Q"] ) [ielem + 1], ielem);
    }
}

void
OneDimensionalModel_Solver::computePressure( solution_Type& solution, const Real& TimeStep )
{
    for ( UInt i(0); i < M_FESpace -> dim() ; ++i )
    {
        ( *solution["P"] ) [i+1] = M_physics -> elasticPressure( ( *solution["A"] ) [i+1], i );

        if ( M_physics -> Data() -> viscoelasticWall() )
        {
            // Viscoelastic pressure is not working right now!
            // We need to pass from outside U_prevtime and U_2prevtime
            ( *solution["P_visc"] ) [i+1] = M_physics -> viscoelasticPressure( ( *solution["A"]) [i+1],
                                                                            ( *M_UPreviousTime["A"]) [i+1],
                                                                            ( *M_U2PreviousTime["A"]) [i+1], TimeStep, i );
            ( *solution["P"] ) [i+1] += ( *solution["P_visc"] ) [i+1];
        }
    }
}

void
OneDimensionalModel_Solver::computeAreaRatio( solution_Type& solution )
{
    for ( UInt inode(M_leftNodeId); inode <= M_rightNodeId ; ++inode )
        ( *solution["A/A0-1"] ) [inode] = ( *solution["A"] ) [inode] / M_physics -> Data() -> Area0(inode - 1) - 1;
}

void
OneDimensionalModel_Solver::computeArea( solution_Type& solution )
{
    for ( UInt inode(M_leftNodeId); inode <= M_rightNodeId ; ++inode )
        ( *solution["A"] ) [inode] = ( (*solution["A/A0-1"] ) [inode] + 1 ) * M_physics -> Data() -> Area0(inode - 1);
}

/*
void
OneDimensionalModel_Solver::initialize( solution_Type& solution )
{
    ASSERT_PRE( (M_leftNodeId <= M_rightNodeId) && (M_rightNodeId <= M_FESpace -> dim()),
                "[initialize] outside tube boundaries" );

    Real value1, value2;

    Vector exponent( M_FESpace -> dim() );
    for ( UInt inode(M_leftNodeId); inode <= M_rightNodeId ; ++inode )
    {
        exponent[inode - 1] *= 0.;

        // first half of a gaussian signal, centered in firstnode;
        // second half of a gaussian signal, centered in lastnode;
        // width represents the total duration of the gaussian signal
        // (rise + decay)
        if( (inode < M_leftNodeId) || (inode > M_rightNodeId) )
        {
            exponent[inode - 1]  = - std::pow( Real(int( inode - M_leftNodeId )), 2 );
            exponent[inode - 1] /= 2*std::pow( M_physics -> Data() -> Length(), 2 );
        }
    }

    switch( M_physics -> Data() -> initialVariable() )
    {
        case OneD_InitializeArea:
            Debug( 6310 ) << "[initialize] 0- OneDInitArea\n";

            value1 = M_physics -> Data() -> initialValue();
            value2 = 0;

            *solution["Q"] = value2;

            for ( UInt inode(M_leftNodeId); inode <= M_rightNodeId ; ++inode )
            {
                (*solution["A"])[inode] = value1 * ( 1 + M_physics -> Data() -> multiplier() * std::exp( exponent[inode - 1] ) );

                M_physics -> W_from_U( (*solution["W1"])[inode], (*solution["W2"])[inode],
                                     (*solution["A"])[inode],  (*solution["Q"])[inode], inode - 1 );
            }
            break;

        case OneD_InitializeFlux:
            // HYPOTHESIS: when initializing flux, area is equal to Area0
            Debug( 6310 ) << "[initialize] 0- OneDInitFlux\n";

            value1 = M_physics -> Data() -> Area0(0); // this if Area0 is constant
            value2 = M_physics -> Data() -> initialValue();

            for ( UInt inode(M_leftNodeId); inode <= M_rightNodeId ; ++inode )
            {
                (*solution["A"])[inode] = M_physics -> Data() -> Area0(inode - 1);
                (*solution["Q"])[inode] = value2*( 1 + M_physics -> Data() -> multiplier()*std::exp( exponent[inode - 1] ) );

                M_physics -> W_from_U( (*solution["W1"])[inode], (*solution["W2"])[inode],
                                     (*solution["A"])[inode],  (*solution["Q"])[inode], inode - 1 );
            }
            break;

        case OneD_InitializeRiemann1:
            Debug( 6310 ) << "[initialize] 0- OneDInitRiemann1\n";

            value1 = M_physics -> Data() -> initialValue();
            value2 = -value1;

            std::cout << "[initialize] WARNING! Initializing W2 = - W1" << " (assuming Q = 0)" << std::endl;

            for ( UInt inode(M_leftNodeId); inode <= M_rightNodeId ; ++inode )
            {
                (*solution["W1"])[inode] = value1 *
                  ( 1 + M_physics -> Data() -> multiplier() * std::exp( exponent[inode - 1] ) );
                (*solution["W2"])[inode] = value2 *
                  ( 1 + M_physics -> Data() -> multiplier() * std::exp( exponent[inode - 1] ) );

                M_physics -> U_from_W( (*solution["A"])[inode],  (*solution["Q"])[inode],
                                     (*solution["W1"])[inode], (*solution["W2"])[inode], inode - 1);
            }

            break;

        case OneD_InitializeRiemann2:
            Debug( 6310 ) << "[initialize] 0- OneDInitRiemann2\n";

            value1 = M_physics -> Data() -> initialValue();
            value2 = - value1;

            std::cout << "[initialize] WARNING! Initializing W1 = - W2" << " (assuming Q = 0)" << std::endl;

            for (UInt inode = M_leftNodeId; inode <= M_rightNodeId; ++inode )
            {
                (*solution["W2"])[inode] = value1 *
                  ( 1 + M_physics -> Data() -> multiplier() * std::exp( exponent[inode - 1] ) );
                (*solution["W1"])[inode] = value2 *
                  ( 1 + M_physics -> Data() -> multiplier() * std::exp( exponent[inode - 1] ) );

                M_physics -> U_from_W( (*solution["A"])[inode],  (*solution["Q"])[inode],
                                     (*solution["W1"])[inode], (*solution["W2"])[inode], inode -1 );
            }

            break;

        case OneD_InitializePressure:
            Debug( 6310 ) << "[initialize] 0- OneDInitPressure\n";

            value1 = M_physics -> Data() -> restValue();
            value2 = 0; // HYPOTHESIS: when initializing pressure, flux is imposed constant = 0

            *solution["Q"] = value2;

            for ( UInt inode = M_leftNodeId; inode <= M_rightNodeId ; ++inode )
            {
                value2 = value1*( 1 + M_physics -> Data() -> multiplier()*std::exp( exponent[inode - 1] ) );

                (*solution["A"])[inode] = M_physics -> A_from_P( value2 );

                M_physics -> W_from_U( (*solution["W1"])[inode], (*solution["W2"])[inode],
                                     (*solution["A"])[inode],  (*solution["Q"])[inode], inode - 1);

            }
            break;

        default:
            ERROR_MSG("No such initializing option.");
    }

    Debug( 6310 ) << "[initialize]\t\tvalue1         = " << value1 << "\n";
    Debug( 6310 ) << "[initialize]\t\tvalue1_step    = " << value1 * M_physics -> Data() -> multiplier() << "\n";
    Debug( 6310 ) << "[initialize]\t\tvalue2         = " << value2 << "\n";

    // Compute A/A0 from A
    computeAreaRatio( solution );

    // Compute initial pressure (taking into account the viscoelastic wall)
    computePressure( solution, M_physics -> Data() -> dataTime() -> getTimeStep() );

    // Prepare the buffers
    //openFileBuffers( solution );
}

void
OneDimensionalModel_Solver::initialize( solution_Type& solution, const Real& u10, const Real& u20, const std::string& var )
{
    Debug( 6310 ) << "[OneDimensionalModel_Solver::initialize] O- Initialize: var " << var << "\n";

    if( var == "physical")
    {
        Debug( 6310 ) << "[OneDimensionalModel_Solver::initialize] O- Imposing real values ... \n";
        Debug( 6310 ) << "[OneDimensionalModel_Solver::initialize] O- A0 = " << u10 << "\n";
        Debug( 6310 ) << "[OneDimensionalModel_Solver::initialize] O- Q0 = " << u20 << "\n";

        for (UInt ielem = 0; ielem < M_FESpace -> dim() ; ++ielem )
        {
            (*solution["A"])[ielem + 1] = u10;
            (*solution["Q"])[ielem + 1] = u20;

            M_physics -> W_from_U( (*solution["W1"])[ielem + 1], (*solution["W2"])[ielem + 1],
                                 (*solution["A"])[ielem + 1], (*solution["Q"])[ielem + 1],
                                 ielem );
        }
        Debug( 6310 ) << "[OneDimensionalModel_Solver::initialize] O- ok\n";
    }
    else if( var == "Riemann" )
    {
        (*solution["W1"]) = vector_Type( M_FESpace -> map() );
        (*solution["W1"]) = u10;
        (*solution["W2"]) = vector_Type( M_FESpace -> map() );
        (*solution["W2"]) = u20;

        for (UInt ielem = 0; ielem < M_FESpace -> dim() ; ++ielem )
            M_physics -> U_from_W( (*solution["A"])[ielem + 1],  (*solution["Q"])[ielem + 1],
                                 (*solution["W1"])[ielem + 1], (*solution["W2"])[ielem + 1],
                                   ielem + 1 ); // WARNING the +1 is not debugged yet (GF 12/2009)
    }
    else
    {
        std::cout << "[initialize] trying to initialize " << var << " variables!" << std::endl;
        abort();
    }

    // Compute A/A0
    computeAreaRatio( solution );

    // Compute initial pressure (taking into account the viscoelastic wall)
    computePressure( solution, M_physics -> Data() -> dataTime() -> getTimeStep() );

    // Prepare the buffers
    //openFileBuffers( solution );
}

void
OneDimensionalModel_Solver::initialize( solution_Type& solution, const vector_Type& u10, const vector_Type& u20 )
{
    (*solution["A"]) = u10;
    (*solution["Q"]) = u20;

    for (UInt ielem=0; ielem <= M_FESpace -> dim() ; ++ielem )
    {
        M_physics -> W_from_U( (*solution["W1"])[ielem], (*solution["W2"])[ielem],
                             (*solution["A"])[ielem], (*solution["Q"])[ielem], ielem );
    }

    // Compute A/A0
    computeAreaRatio( solution );

    // Compute initial pressure (taking into account the viscoelastic wall)
    computePressure( solution, M_physics -> Data() -> dataTime() -> getTimeStep() );

    // Prepare the buffers
    //openFileBuffers( solution );
}

void
OneDimensionalModel_Solver::initialize( solution_Type& solution, const Real& u20 )
{
    (*solution["Q"]) = vector_Type( M_FESpace -> map() );
    (*solution["Q"]) = u20;

    for (UInt ielem = 0; ielem <= M_FESpace -> dim() ; ++ielem )
    {
        (*solution["A"])[ielem]=M_physics -> Data() -> Area0(ielem);
        M_physics -> W_from_U( (*solution["W1"])[ielem], (*solution["W2"])[ielem],
                             (*solution["A"])[ielem], (*solution["Q"])[ielem], ielem );
    }

    // Compute A/A0
    computeAreaRatio( solution );

    // Compute initial pressure (taking into account the viscoelastic wall)
    computePressure( solution, M_physics -> Data() -> dataTime() -> getTimeStep() );

    // Prepare the buffers
    //openFileBuffers( solution );
}
*/

void
OneDimensionalModel_Solver::updateRHS( const solution_Type& solution, const Real& TimeStep )
{
    updateFluxDer( solution );   // Update the vector containing the values of the flux at the nodes and its jacobian
    updateSourceDer( solution ); // Update the vector containing the values of the source term at the nodes and its jacobian
    updateMatrices();            // Update the matrices for the non-linear terms

    // Taylor-Galerkin scheme: (explicit, U = [U1,U2]^T, with U1=A, U2=Q )
    // (Un+1, phi) =          (               Un,     phi     ) ->  massFactor^{-1} * Un+1 = mass * U
    //             + dt     * (           Fh(Un),     dphi/dz ) ->             grad * F(U)
    //             - dt^2/2 * (diffFh(Un) Sh(Un),     dphi/dz ) ->  gradDiffFlux(U) * S(U)
    //             + dt^2/2 * (diffSh(Un) dFh/dz(Un), phi     ) ->    divDiffSrc(U) * F(U)
    //             - dt^2/2 * (diffFh(Un) dFh/dz(Un), dphi/dz ) -> stiffDiffFlux(U) * F(U)
    //             - dt     * (           Sh(Un),     phi     ) ->             mass * S(U)
    //             + dt^2/2 * (diffSh(Un) Sh(Un),     phi     ) ->   massDiffSrc(U) * S(U)

    Real dt2over2 = TimeStep * TimeStep * 0.5;

    // rhs = mass * Un
    M_rhs[0]  = ( *M_massMatrix ) * *solution.find("A") -> second;
    M_rhs[1]  = ( *M_massMatrix ) * *solution.find("Q") -> second;
    for ( UInt i(0); i < 2; ++i )
    {
        // rhs = rhs + dt * grad * F(Un)
        M_rhs[i] += ( *M_gradMatrix ) * ( TimeStep * M_fluxVector[i] );

        // rhs = rhs - dt * mass * S(Un)
        M_rhs[i] += ( *M_massMatrix ) * ( - TimeStep * M_sourceVector[i] );

        for ( UInt j(0); j < 2; ++j )
        {
            // rhs = rhs - dt^2/2 * gradDiffFlux * S(Un)
            M_gradMatrixDiffFlux[2*i + j] -> GlobalAssemble();
            M_rhs[i] += *M_gradMatrixDiffFlux[2*i + j] * ( -dt2over2 * M_sourceVector[j] );

            // rhs = rhs + dt^2/2 * divDiffSrc * F(Un)
            M_divMatrixDiffSrc[2*i + j] -> GlobalAssemble();
            M_rhs[i] += *M_divMatrixDiffSrc[2*i + j] * ( dt2over2 * M_fluxVector[j] );

            // rhs = rhs - dt^2/2 * stiffDiffFlux * F(Un)
            M_stiffMatrixDiffFlux[2*i + j] -> GlobalAssemble();
            M_rhs[i] += *M_stiffMatrixDiffFlux[2*i + j]*(-dt2over2*M_fluxVector[j]);

            // rhs = rhs + dt^2/2 * massDiffSrc * S(Un)
            M_massMatrixDiffSrc[2*i + j] -> GlobalAssemble();
            M_rhs[i] += *M_massMatrixDiffSrc[2*i + j] * ( dt2over2 * M_sourceVector[j] );
        }
    }
}

void
OneDimensionalModel_Solver::iterate( OneDimensionalModel_BCHandler& bcH, solution_Type& solution, const Real& Time, const Real& TimeStep )
{
    // Compute BC
    Debug( 6310 ) << "[timeAdvance] \tcompute BC\n";
    bcH.applyBC( Time, TimeStep, solution, M_flux, M_bcDirLeft, M_bcDirRight );

    // Apply BC to RHS
    Debug( 6310 ) << "[timeAdvance] \tcompute BC dirichlet vector\n";
    updateBCDirichletVector();

    if ( M_physics -> Data() -> UW() )
    {
        Real Ainode, Qinode;

        Real lambda1_plus  = 0.;
        Real lambda2_plus  = 0.;

        Real lambda1_minus = 0.;
        Real lambda2_minus = 0.;

        container2D_Type eigenvalues;
        container2D_Type leftEigenvector1;
        container2D_Type leftEigenvector2;

        //Real deltaX = M_FESpace -> mesh() -> edgeList( 1 ).point(1).x() - M_FESpace -> mesh() -> edgeList( 1 ).point(2).x();
        Real delta  = - std::sqrt( ( ( M_FESpace -> mesh() -> edgeList[ 1 ].point( 2 ) ).x() - ( M_FESpace -> mesh() -> edgeList[ 1 ].point( 1 ) ).x() ) *
                                   ( ( M_FESpace -> mesh() -> edgeList[ 1 ].point( 2 ) ).x() - ( M_FESpace -> mesh() -> edgeList[ 1 ].point( 1 ) ).x() ) +
                                   ( ( M_FESpace -> mesh() -> edgeList[ 1 ].point( 2 ) ).y() - ( M_FESpace -> mesh() -> edgeList[ 1 ].point( 1 ) ).y() ) *
                                   ( ( M_FESpace -> mesh() -> edgeList[ 1 ].point( 2 ) ).y() - ( M_FESpace -> mesh() -> edgeList[ 1 ].point( 1 ) ).y() ) +
                                   ( ( M_FESpace -> mesh() -> edgeList[ 1 ].point( 2 ) ).z() - ( M_FESpace -> mesh() -> edgeList[ 1 ].point( 1 ) ).z() ) *
                                   ( ( M_FESpace -> mesh() -> edgeList[ 1 ].point( 2 ) ).z() - ( M_FESpace -> mesh() -> edgeList[ 1 ].point( 1 ) ).z() )
                                 );

        // working on riemann invariants
        vector_Type W1_UW(M_FESpace -> map());
        vector_Type W2_UW(M_FESpace -> map());

        // converting boundary conditions on physical variables
        // in boundary conditions on characteristic variables
        M_physics -> W_from_U( W1_UW(0), W2_UW(0), M_rhs[0][0], M_rhs[1][0], 0 );
        M_physics -> W_from_U( W1_UW(M_FESpace -> dim()-1), W2_UW(M_FESpace -> dim()-1),
                             M_rhs[0][M_FESpace -> dim()-1], M_rhs[1][M_FESpace -> dim()-1],
                             M_FESpace -> dim() - 1 );

        for ( UInt ii=1; ii < (M_FESpace -> dim()-1) ; ++ii )
        {
            // compute the eigenvalues at node
            Ainode = (*solution["A"])( ii );
            Qinode = (*solution["Q"])( ii );
            M_flux -> eigenValuesEigenVectors( Ainode, Qinode,
                                             eigenvalues,
                                             leftEigenvector1,
                                             leftEigenvector2,
                                             ii - 1 );

            lambda1_plus  = std::max<Real>( eigenvalues[0], 0. );
            lambda1_minus = std::min<Real>( eigenvalues[0], 0. );
            lambda2_plus  = std::max<Real>( eigenvalues[1], 0. );
            lambda2_minus = std::min<Real>( eigenvalues[1], 0. );

            // update the solution for the next time step
            W1_UW[ii] = ( *solution["W1"] ) [ii]
                        - ( TimeStep / delta )  * lambda1_plus  * ( ( *solution["W1"] ) [ii]   - ( *solution["W1"] ) [ii-1] )
                        - ( TimeStep / delta )  * lambda1_minus * ( ( *solution["W1"] ) [ii+1] - ( *solution["W1"] ) [ii] )
                        -  TimeStep * ( leftEigenvector1[0] * M_sourceVector[0][ii] + leftEigenvector1[1] * M_sourceVector[1][ii] );
            W2_UW[ii] = ( *solution["W2"] ) [ii]
                        - ( TimeStep / delta )  * lambda2_plus  * ( ( *solution["W2"] ) [ii]   - ( *solution["W2"] ) [ii-1] )
                        - ( TimeStep / delta )  * lambda2_minus * ( ( *solution["W2"] ) [ii+1] - ( *solution["W2"] ) [ii] )
                        -   TimeStep * ( leftEigenvector2[0] * M_sourceVector[0][ii] + leftEigenvector2[1] * M_sourceVector[1][ii] );
        }

        (*solution["W1"]) = W1_UW;
        (*solution["W2"]) = W2_UW;

        for (UInt ielem=0; ielem <= M_FESpace -> dim() ; ++ielem )
            M_physics -> U_from_W( ( *solution["A"] ) [ielem], ( *solution["Q"] ) [ielem],
                                   ( *solution["W1"] ) [ielem], ( *solution["W2"] ) [ielem],
                                     ielem );
    }
    else
    {
        // Compute A^n+1
        vector_Type area( M_rhs[0] );
        M_linearSolver -> solveSystem( M_rhs[0], area, M_massMatrix );

        // Compute Q^n+1
        vector_Type flowrate( M_rhs[1] );
        M_linearSolver -> solveSystem( M_rhs[1], flowrate, M_massMatrix );

        // Correct flux with inertial, viscoelastic and longitudinal terms
        if ( M_physics -> Data() -> inertialWall() )
        {
            *solution["Q_inert"] = inertialFluxCorrection( flowrate );
            flowrate += *solution["Q_inert"];
        }

        if ( M_physics -> Data() -> viscoelasticWall() )
        {
            *solution["Q_visc"] = viscoelasticFluxCorrection( flowrate, TimeStep );
            flowrate += *solution["Q_visc"];
        }

        if ( M_physics -> Data() -> longitudinalWall() )
        {
            *solution["Q_long"] = longitudinalFluxCorrection();
            flowrate += *solution["Q_long"];
        }

        // compute L2 projection of d2Q_dx2
        //    if( M_physics -> Data() -> fluxSecondDer() )
        //      M_d2_U2_dx2 = _compute_d2Q_dx2( flowrate );

        // Update the solution container
        *solution["A"] = area;
        *solution["Q"] = flowrate;

        // Compute W1 and W2 from A and Q
        computeW1W2( solution );

        // Compute A/A0 from A
        computeAreaRatio( solution );

        // Update the pressure (taking into account the viscoelastic wall)
        computePressure( solution, TimeStep );
    }
}

Real
OneDimensionalModel_Solver::ComputeCFL( const solution_Type& solution, const Real& timeStep ) const
{
    Real lambdaMax = 0.;

    container2D_Type eigenvalues;
    container2D_Type leftEigenvector1;
    container2D_Type leftEigenvector2;

    for ( UInt inode(0); inode < M_FESpace -> dim() ; ++inode )
    {
        // compute the eigenvalues at node
        M_flux -> eigenValuesEigenVectors( ( *solution.find("A") -> second ) ( inode + 1 ),
                                           ( *solution.find("Q") -> second ) ( inode + 1 ),
                                             eigenvalues, leftEigenvector1, leftEigenvector2, inode );

        lambdaMax = std::max<Real>( std::max<Real>( std::fabs(eigenvalues[0]), std::fabs(eigenvalues[1]) ), lambdaMax );
    }

    return lambdaMax * timeStep / M_FESpace -> mesh() -> minH();
}

/*
void
OneDimensionalModel_Solver::savesol()
{
    M_UPreviousTime[0] = (*M_U_thistime)[0];
    M_UPreviousTime[1] = (*M_U_thistime)[1];
}

void
OneDimensionalModel_Solver::loadsol()
{
    container2D_Type U_leftbd, U_rightbd;
    for( UInt i=0; i<2; ++i )
    {
        U_leftbd[i] = (*M_U_thistime)[i]( M_leftNodeId );
        U_rightbd[i] = (*M_U_thistime)[i]( M_rightNodeId );
    }

    (*M_U_thistime)[0] = M_UPreviousTime[0];
    (*M_U_thistime)[1] = M_UPreviousTime[1];

    for( UInt i=0; i<2; ++i )
    {
        (*M_U_thistime)[i]( M_leftNodeId ) = U_leftbd[i];
        (*M_U_thistime)[i]( M_rightNodeId ) = U_rightbd[i];
    }

    for ( UInt inode=M_leftNodeId; inode <= M_rightNodeId ; ++inode )
        M_physics -> W_from_U( (*M_U_thistime)[2][inode], (*M_U_thistime)[3][inode],
                             (*M_U_thistime)[0][inode], (*M_U_thistime)[1][inode],
                             inode - 1 );
}
*/

void
OneDimensionalModel_Solver::resetOutput( const solution_Type& solution )
{
    std::ofstream outfile;
    for ( solutionConstIterator_Type i = solution.begin(); i != solution.end(); ++i )
    {
        std::string file = M_physics -> Data() -> postprocessingDirectory() + "/" + M_physics -> Data() -> postprocessingFile() + "_" + i -> first + ".m";
        outfile.open( file.c_str(), std::ios::trunc );
        outfile.close();
    }
}

void
OneDimensionalModel_Solver::postProcess( const solution_Type& solution )
{
    std::ofstream outfile;
    for ( solutionConstIterator_Type i = solution.begin(); i != solution.end(); ++i )
    {
        std::string file = M_physics -> Data() -> postprocessingDirectory() + "/" + M_physics -> Data() -> postprocessingFile() + "_" + i -> first + ".m";
        outfile.open( file.c_str(), std::ios::app );

        for ( UInt ii(0); ii < M_FESpace -> dim(); ++ii )
            outfile << (*i -> second)(ii + 1) << " ";

        outfile << std::endl;
        outfile.close();
    }
}

/*
void
OneDimensionalModel_Solver::create_movie_file()
{
    std::ofstream outfile;
    std::string file_output;
    file_output = M_physics -> Data() -> PostDirectory() + "/" + "areamovie"+M_physics -> Data() -> PostFile()+".m";
    outfile.open(file_output.c_str(), std::ios::out);
    outfile <<"Area"<<M_physics -> Data() -> PostFile()<<";\n"
            <<"[n,m]=size(A" << M_physics -> Data() -> PostFile() << ");\n"
            <<"Max=max(max(A" << M_physics -> Data() -> PostFile() << "));\n"
            <<"Min=min(min(A" << M_physics -> Data() -> PostFile() << "));\n"
            <<"for i=1:n\n"
            <<"  plot(A" << M_physics -> Data() -> PostFile() << "(i,:));\n"
            <<"  title((i-1)*"<<M_physics -> Data() -> dataTime() -> getTimeStep()<<"*"<<M_physics -> Data() -> verbose()<<");\n"
            <<"  axis([0 m Min-(Min/100) Max+(Min/100)]);\n"
            <<"  %pause;\n"
            <<"  F(i) = getframe;\n"
            <<"end\n"
            <<"movie(F)";

    outfile.close();
    file_output = M_physics -> Data() -> PostDirectory() + "/" + "portatamovie"+M_physics -> Data() -> PostFile()+".m";
    outfile.open(file_output.c_str(), std::ios::out);
    outfile <<"Portata"<<M_physics -> Data() -> PostFile()<<";\n"
            <<"[n,m]=size(Q" << M_physics -> Data() -> PostFile() << ");\n"
            <<"Max=max(max(Q" << M_physics -> Data() -> PostFile() << "));\n"
            <<"Min=min(min(Q" << M_physics -> Data() -> PostFile() << "));\n"
            <<"for i=1:n\n"
            <<"  plot(Q" << M_physics -> Data() -> PostFile() << "(i,:));\n"
            <<"  title((i-1)*"<< M_physics -> Data() -> dataTime() -> getTimeStep() <<"*"<<M_physics -> Data() -> verbose()<<");\n"
            <<"  axis([0 m Min-(Min/100) Max+(Min/100)]);\n"
            <<"  %pause;\n"
            <<"  F(i) = getframe;\n"
            <<"end\n"
            <<"movie(F)";
    outfile.close();
}

void
OneDimensionalModel_Solver::openFileBuffers( const solution_Type& solution )
{
    std::string file_output;
    boost::shared_ptr<std::ostringstream> buf;

    for ( solutionConstIterator_Type i = solution.begin(); i != solution.end(); ++i )
    {
        file_output = M_physics -> Data() -> PostDirectory() + "/" + M_physics -> Data() -> PostFile() + i -> first + ".m";
        Debug( 6310 ) << "[openFileBuffers] setting output for file " << file_output << "\n";

        buf.reset( new std::ostringstream() );
        buf -> setf(std::ios_base::scientific);
        buf -> precision(5);
        buf -> width(13);

        M_post_process_buffer.insert( std::map<std::string, boost::shared_ptr<std::ostringstream> >::value_type( file_output, buf ) );
        M_post_process_buffer_offset.insert( std::map<std::string, long>::value_type( file_output, buf -> tellp() ) );
    }
}

void
OneDimensionalModel_Solver::resetFileBuffers( const solution_Type& solution )
{
    closeFileBuffers();
    openFileBuffers( solution );
}

void
OneDimensionalModel_Solver::seekpFileBuffers()
{
    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator iter;

    for( iter = M_post_process_buffer.begin(); iter != M_post_process_buffer.end(); ++iter )
        (*iter).second -> seekp( M_post_process_buffer_offset[(*iter).first] );
}

void
OneDimensionalModel_Solver::tellpFileBuffers()
{
    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator iter;

    for( iter = M_post_process_buffer.begin(); iter != M_post_process_buffer.end(); ++iter )
        M_post_process_buffer_offset[(*iter).first] = (*iter).second -> tellp();
}

void
OneDimensionalModel_Solver::closeFileBuffers()
{
    // as I have a boost::shared_ptr, I expect the objects to be deallocated
    // now that the pointers are destroyed
    M_post_process_buffer.erase( M_post_process_buffer.begin(),
                                 M_post_process_buffer.end() );
    M_post_process_buffer_offset.erase( M_post_process_buffer_offset.begin(),
                                        M_post_process_buffer_offset.end() );
}
*/

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_Solver::setProblem( const physicsPtr_Type Physics,
                                        const fluxPtr_Type    Flux,
                                        const sourcePtr_Type  Source )
{
    M_physics = Physics;
    M_flux    = Flux;
    M_source  = Source;
}

void
OneDimensionalModel_Solver::setCommunicator( const commPtr_Type Comm )
{
    M_comm = Comm;
    M_displayer.setCommunicator( Comm );
}

void
OneDimensionalModel_Solver::setFESpace( const FESpacePtr_Type FESpace )
{
    M_FESpace = FESpace;

    //Id of left and right bc nodes
    M_leftNodeId          = 1;
    M_leftInternalNodeId  = M_leftNodeId + 1;
    M_rightNodeId         = M_FESpace -> dim();
    M_rightInternalNodeId = M_rightNodeId - 1;

    //Elementary Matrices
    M_elmatMass.reset ( new ElemMat ( M_FESpace -> fe().nbNode, 1, 1 ) );
    M_elmatStiff.reset( new ElemMat ( M_FESpace -> fe().nbNode, 1, 1 ) );
    M_elmatGrad.reset ( new ElemMat ( M_FESpace -> fe().nbNode, 1, 1 ) );
    M_elmatDiv.reset  ( new ElemMat ( M_FESpace -> fe().nbNode, 1, 1 ) );

    //Vectors
    M_rhs.resize(          2, vector_Type( M_FESpace -> map() ) );
    M_fluxVector.resize(   2, vector_Type( M_FESpace -> map() ) );
    M_sourceVector.resize( 2, vector_Type( M_FESpace -> map() ) );

    //Matrix
    M_massMatrix.reset( new matrix_Type( M_FESpace -> map() ) );
    M_gradMatrix.reset( new matrix_Type( M_FESpace -> map() ) );
}

void
OneDimensionalModel_Solver::setLinearSolver( const linearSolverPtr_Type linearSolver )
{
    M_linearSolver = linearSolver;
}

void
OneDimensionalModel_Solver::setBCValuesLeft( const Real& bcL1, const Real& bcL2 )
{
    M_bcDirLeft[0] = bcL1;
    M_bcDirLeft[1] = bcL2;
}

void
OneDimensionalModel_Solver::setBCValuesRight( const Real& bcR1, const Real& bcR2 )
{
    M_bcDirRight[0] = bcR1;
    M_bcDirRight[1] = bcR2;
}

// ===================================================
// Get Methods
// ===================================================
const OneDimensionalModel_Solver::physicsPtr_Type&
OneDimensionalModel_Solver::Physics() const
{
    return M_physics;
}

const OneDimensionalModel_Solver::fluxPtr_Type&
OneDimensionalModel_Solver::Flux() const
{
    return M_flux;
}

const OneDimensionalModel_Solver::sourcePtr_Type&
OneDimensionalModel_Solver::Source() const
{
    return M_source;
}

const UInt&
OneDimensionalModel_Solver::LeftNodeId() const
{
    return M_leftNodeId;
}

const UInt&
OneDimensionalModel_Solver::LeftInternalNodeId() const
{
    return M_leftInternalNodeId;
}

const UInt&
OneDimensionalModel_Solver::RightNodeId() const
{
    return M_rightNodeId;
}

const UInt&
OneDimensionalModel_Solver::RightInternalNodeId() const
{
    return M_rightInternalNodeId;
}

container2D_Type
OneDimensionalModel_Solver::BCValuesLeft( const solution_Type& solution ) const
{
    container2D_Type temp;

    temp[0] = (*solution.find("A") -> second)( LeftNodeId() );
    temp[1] = (*solution.find("Q") -> second)( LeftNodeId() );

    return temp;
}

container2D_Type
OneDimensionalModel_Solver::BCValuesInternalLeft( const solution_Type& solution ) const
{
    container2D_Type temp;

    temp[0] = (*solution.find("A") -> second)( LeftInternalNodeId() );
    temp[1] = (*solution.find("Q") -> second)( LeftInternalNodeId() );

    return temp;
}

container2D_Type
OneDimensionalModel_Solver::BCValuesRight( const solution_Type& solution ) const
{
    container2D_Type temp;

    temp[0] = (*solution.find("A") -> second)( RightNodeId() );
    temp[1] = (*solution.find("Q") -> second)( RightNodeId() );

    return temp;
}

container2D_Type
OneDimensionalModel_Solver::BCValuesInternalRight( const solution_Type& solution ) const
{
    container2D_Type temp;

    temp[0] = (*solution.find("A") -> second)( RightInternalNodeId() );
    temp[1] = (*solution.find("Q") -> second)( RightInternalNodeId() );

    return temp;
}

Real
OneDimensionalModel_Solver::BoundaryValue( const solution_Type& solution, const OneD_BC& bcType, const OneD_BCSide& bcSide ) const
{
    UInt boundaryDof;

    switch ( bcSide )
    {
    case OneD_left:
        boundaryDof = 1;
        break;
    case OneD_right:
        boundaryDof = M_flux -> physics() -> Data() -> NumberOfElements() + 1;
        break;
    default:
        std::cout << "Warning: bcSide \"" << bcSide << "\" not available!" << std::endl;
        return 0;
        break;
    }

    switch ( bcType )
    {
    case OneD_A:
        return (*solution.find("A") -> second)( boundaryDof );
    case OneD_Q:
        // Flow rate is positive with respect to the outgoing normal
        return (*solution.find("Q") -> second)( boundaryDof ) * ( ( bcSide == OneD_left ) ? -1. : 1. );
    case OneD_W1:
        return (*solution.find("W1") -> second)( boundaryDof );
    case OneD_W2:
        return (*solution.find("W2") -> second)( boundaryDof );
    case OneD_P:
        return (*solution.find("P") -> second)( boundaryDof );
    default:
        std::cout << "Warning: bcType \"" << bcType << "\"not available!" << std::endl;
        return 0.;
    }
}

void
OneDimensionalModel_Solver::BoundaryEigenValuesEigenVectors( const OneD_BCSide& bcSide,
                                                             const solution_Type& solution,
                                                             container2D_Type& eigenvalues,
                                                             container2D_Type& leftEigenvector1,
                                                             container2D_Type& leftEigenvector2 )
{
    UInt boundaryDof;

    switch ( bcSide )
    {
    case OneD_left:
        boundaryDof = 0;
        break;
    case OneD_right:
        boundaryDof = M_flux -> physics() -> Data() -> NumberOfElements();
        break;
    default:
        std::cout << "Warning: bcSide \"" << bcSide << "\" not available!" << std::endl;
        return;
    }

    M_flux -> eigenValuesEigenVectors( (*solution.find("A") -> second)( boundaryDof + 1 ),
                                     (*solution.find("Q") -> second)( boundaryDof + 1 ),
                                     eigenvalues, leftEigenvector1, leftEigenvector2,
                                     boundaryDof );
}

// ===================================================
// Private Methods
// ===================================================
void
OneDimensionalModel_Solver::updateFlux( const solution_Type& solution )
{
    Real Aii, Qii;

    for ( UInt ii = 1; ii <= M_FESpace -> dim() ; ++ii )
    {
        Aii = (*solution.find("A") -> second)( ii );
        Qii = (*solution.find("Q") -> second)( ii );
        M_fluxVector[0]( ii ) = ( *M_flux )( Aii, Qii, 1, ii - 1 );
        M_fluxVector[1]( ii ) = ( *M_flux )( Aii, Qii, 2, ii - 1 );
    }
}

void
OneDimensionalModel_Solver::updateFluxDer( const solution_Type& solution )
{
    // first update the Flux vector
    updateFlux( solution );

    // then update the derivative of the Flux vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;
    UInt ii, iip1;

    for ( UInt ielem = 0; ielem < M_FESpace -> dim() - 1; ++ielem )
    {
        // for P1Seg and appropriate mesh only!
        ii    = ielem;      // left node of current element
        iip1  = ielem + 1;  // right node of current element

        Aii   = (*solution.find("A") -> second)( ielem + 1 );
        Qii   = (*solution.find("Q") -> second)( ielem + 1 );

        Aiip1 = (*solution.find("A") -> second)( ielem + 2 );
        Qiip1 = (*solution.find("Q") -> second)( ielem + 2 );

        for ( UInt ii=1; ii<3; ++ii )
        {
            for ( UInt jj=1; jj<3; ++jj )
            {
                tmp  = M_flux -> diff(   Aii,   Qii, ii, jj, ii );
                tmp += M_flux -> diff( Aiip1, Qiip1, ii, jj, iip1 );

                M_diffFlux[ 2*(ii - 1) + jj - 1 ]( ielem ) = 0.5*tmp;
            }
        }
    }
}

void
OneDimensionalModel_Solver::updateSource( const solution_Type& solution )
{
    Real Aii, Qii;

    for ( UInt ii = 1; ii <= M_FESpace -> dim() ; ++ii )
    {
        Aii = (*solution.find("A") -> second)( ii );
        Qii = (*solution.find("Q") -> second)( ii );
        for (UInt k=0; k<2; ++k)
            M_sourceVector[k]( ii ) = ( *M_source )( Aii, Qii, k+1, ii - 1);
    }
}

void
OneDimensionalModel_Solver::updateSourceDer( const solution_Type& solution )
{
    // first update the Source vector
    updateSource( solution );

    // then update the derivative of the Source vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;
    UInt ii, iip1;

    for ( UInt ielem=0; ielem < M_FESpace -> dim() - 1 ; ++ielem )
    {
        // for P1Seg and appropriate mesh only!
        ii = ielem;        // left node of current element
        iip1 = ielem + 1;  // right node of current element

        Aii   = (*solution.find("A") -> second)( ielem + 1);
        Qii   = (*solution.find("Q") -> second)( ielem + 1);
        Aiip1 = (*solution.find("A") -> second)( ielem + 2 );
        Qiip1 = (*solution.find("Q") -> second)( ielem + 2 );

        for ( UInt ii=1; ii<3; ++ii )
        {
            for ( UInt jj=1; jj<3; ++jj )
            {
                tmp =  M_source -> diff(   Aii,   Qii, ii, jj, ii );
                tmp += M_source -> diff( Aiip1, Qiip1, ii, jj, iip1 );
                M_diffSrc[ 2*(ii - 1) + jj - 1 ]( ielem ) = 0.5 * tmp;
            }
        }
    }
}

void
OneDimensionalModel_Solver::updateMatrices()
{
    // Matrices initialization
    for ( UInt i = 0; i < 4; ++i )
    {
        M_massMatrixDiffSrc[i].reset  ( new matrix_Type( M_FESpace -> map() ));
        M_stiffMatrixDiffFlux[i].reset( new matrix_Type( M_FESpace -> map() ));
        M_gradMatrixDiffFlux[i].reset ( new matrix_Type( M_FESpace -> map() ));
        M_divMatrixDiffSrc[i].reset   ( new matrix_Type( M_FESpace -> map() ));
    }

    // Elementary computation and matrix assembling
    for (UInt iedge = 1; iedge <= M_FESpace -> mesh() -> numEdges(); ++iedge)
    {
        // update the current element
        M_FESpace -> fe().updateFirstDerivQuadPt(M_FESpace -> mesh() -> edgeList(iedge));

        for (UInt ii = 1; ii <= 2; ++ii)
        {
            for (UInt jj = 1; jj <= 2; ++jj)
            {
                // update the M_coeff*
                updateMatrixCoefficients( ii , jj, iedge);

                // update the M_elmat*
                updateElemMatrices();

                // assemble the global matrices
                matrixAssemble( ii, jj );
            }
        }
    }
}

void
OneDimensionalModel_Solver::updateMatrixCoefficients( const UInt& ii, const UInt& jj , const UInt& iedge )
{
    Real dFluxdUelem(0), dSrcdUelem(0);

    ASSERT_BD( 0 < ii && ii < 3 && 0 < jj && jj < 3 );

    dFluxdUelem = M_diffFlux[ 2*(ii-1) + jj-1 ]( iedge - 1 ); // iedge starts from 1...
    dSrcdUelem  = M_diffSrc [ 2*(ii-1) + jj-1 ]( iedge - 1 );

    M_coeffGrad  = dFluxdUelem; // term gradDiffFlux(U) [* S(U)]
    M_coeffDiv   = dSrcdUelem;  // term  divDiffSrc(U) [* F(U)]
    M_coeffStiff = dFluxdUelem; // term stiffDiffFlux(U) [* F(U)]
    M_coeffMass  = dSrcdUelem;  // term  massDiffSrc(U) [* S(U)]
}

void
OneDimensionalModel_Solver::updateElemMatrices()
{
    // set the elementary matrices to 0.
    M_elmatMass -> zero();
    M_elmatStiff -> zero();
    M_elmatGrad -> zero();
    M_elmatDiv -> zero();

    // update the mass matrix
    mass( M_coeffMass, *M_elmatMass, M_FESpace -> fe(),0, 0 );
    //  M_elmatMass -> showMe( std::cout );

    // update the stiffness matrix
    stiff( M_coeffStiff, *M_elmatStiff, M_FESpace -> fe(),0 ,0 );
    // std::cout << "Elem Stiff matrix :" << std::endl;
    // M_elmatStiff -> showMe( std::cout );

    /*! update the gradient matrix
      gradient operator:
      grad_{ij} = \int_{fe} coeff \phi_j \frac{d \phi_i}{d x}

      BEWARE :
      \param 0: the first argument "0" corresponds to the first
      and only coordinate (1D!), and HERE it starts from 0... (Damm'!)

      \param - M_coeffGrad: the sign "-" in the second argument
      is added to correspond to the described operator.
      (There is a minus in the elemOper implementation).
    */
    grad( 0 , - M_coeffGrad, *M_elmatGrad, M_FESpace -> fe(), M_FESpace -> fe(), 0, 0 );
    //  std::cout << "Elem Grad matrix :" << std::endl;
    //  M_elmatGrad -> showMe( std::cout );

    /*! update the divergence matrix
      divergence operator: (transpose of the gradient)
      div_{ij} = \int_{fe} coeff \frac{d \phi_j}{d x} \phi_i

      \note formally this M_elmatDiv is not necessary
      as it is the transpose of the M_elmatGrad.
      But for the sake of clarity, I prefer to keep it. (low cost!)

      BEWARE : same remarks as grad (see above).
    */
    div( 0 , - M_coeffDiv, *M_elmatDiv, M_FESpace -> fe(), M_FESpace -> fe(), 0, 0 );
    //  std::cout << "Elem Div matrix :" << std::endl;
    //  M_elmatDiv -> showMe( std::cout );
}

void
OneDimensionalModel_Solver::matrixAssemble( const UInt& ii, const UInt& jj )
{
    ASSERT_BD( 0 < ii && ii < 3 && 0 < jj && jj < 3 );

    // assemble the mass matrix
    assemb_mat( *M_massMatrixDiffSrc[ 2*(ii-1) + jj-1 ]  , *M_elmatMass, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );

    // assemble the stiffness matrix
    assemb_mat( *M_stiffMatrixDiffFlux[ 2*(ii-1) + jj-1 ], *M_elmatStiff, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );

    // assemble the gradient matrix
    assemb_mat( *M_gradMatrixDiffFlux[ 2*(ii-1) + jj-1 ] , *M_elmatGrad, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );

    // assemble the divergence matrix
    assemb_mat( *M_divMatrixDiffSrc[ 2*(ii-1) + jj-1 ]   , *M_elmatDiv, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );
}

void
OneDimensionalModel_Solver::updateBCDirichletVector()
{
    Debug( 6310 ) << "[updateBCDirichletVector] \t firstDof = " << M_leftNodeId
    << " lastDof = " << M_rightNodeId << ";\n";

    Debug( 6310 ) << "[updateBCDirichletVector] \t bcDirLeft[0] = " << M_bcDirLeft[0]
    << " bcDirLeft[1] = " << M_bcDirLeft[1] << ";\n";

    Debug( 6310 ) << "[updateBCDirichletVector] \t bcDirRight[0] = " << M_bcDirRight[0]
    << " bcDirRight[1] = " << M_bcDirRight[1] << ";\n";

    // unsymmetric treatment (LU must be used!)
    // first row modified
    M_rhs[0]( M_leftNodeId ) = M_bcDirLeft[0];
    M_rhs[1]( M_leftNodeId ) = M_bcDirLeft[1];

    // last row modified
    M_rhs[0]( M_rightNodeId ) = M_bcDirRight[0];
    M_rhs[1]( M_rightNodeId ) = M_bcDirRight[1];

    // symmetric treatment (cholesky can be used)
//    for(UInt i=0; i<2; ++i)
//    {
//        // first row modified (Dirichlet)
//        M_rhs[i]( firstDof ) = M_bcDirLeft[i];
//        // second row modified (for symmetry)
//        M_rhs[i]( firstDof + 1 ) += - M_massMatrix.LowDiag()( firstDof ) * M_bcDirLeft[i];
//        // last row modified (Dirichlet)
//        M_rhs[i]( lastDof ) = M_bcDirRight[i];
//        // penultimate row modified (for symmetry)
//        M_rhs[i]( lastDof - 1 ) += - M_massMatrix.UpDiag()( lastDof - 1 ) * M_bcDirRight[i];
//    }
}

void
OneDimensionalModel_Solver::updateBCDirichletMatrix( matrix_Type& mat )
{

    // unsymmetric treatment (LU must be used!)
    mat.diagonalize(M_leftNodeId - 1, 1, 0);
    mat.diagonalize(M_rightNodeId - 1, 1, 0);

    // symmetric treatment (cholesky can be used)
//     mat.Diag()( firstDof )    = 1.;
//     mat.UpDiag()( firstDof )  = 0.;
//     mat.LowDiag()( firstDof ) = 0.; //and second row

//     mat.Diag()( lastDof )      = 1.;
//     mat.UpDiag()( lastDof-1 )  = 0.; //and penultimate row
//     mat.LowDiag()( lastDof-1 ) = 0.;
}

OneDimensionalModel_Solver::vector_Type
OneDimensionalModel_Solver::inertialFluxCorrection( const vector_Type& flux )
{
    matrix_Type _matrixLHS(M_FESpace -> map());
    matrix_Type _stiffRHS (M_FESpace -> map());

    ElemMat _elmatMassLHS  (M_FESpace -> fe().nbNode, 1, 1);
    ElemMat _elmatStiffLHS (M_FESpace -> fe().nbNode, 1, 1);
    ElemMat _elmatStiffRHS (M_FESpace -> fe().nbNode, 1, 1);

    vector_Type _rhs(M_FESpace -> map());

    Real _coeffMass;
    Real _coeffStiff;

    Real m, meanA0;

    //    std::ostringstream output;

    _matrixLHS *= 0.;
    _stiffRHS  *= 0.;

    // Elementary computation and matrix assembling
    // Loop on elements
    for (UInt iedge = 1; iedge <= M_FESpace -> mesh() -> numEdges(); ++iedge)
    {
        // set the elementary matrices to 0.
        _elmatMassLHS. zero();
        _elmatStiffLHS.zero();
        _elmatStiffRHS.zero();

        _coeffMass  = M_rhs[0]( iedge - 1 ) + M_rhs[0]( iedge );
        _coeffMass /= 2;
        _coeffMass  = 1./_coeffMass;

        meanA0  = M_physics -> Data() -> Area0(iedge - 1) + M_physics -> Data() -> Area0(iedge);
        meanA0 /= 2;

        m = M_physics -> Data() -> DensityWall()*M_physics -> Data() -> Thickness(iedge)/
            ( 2*std::sqrt(4*std::atan(1))*std::sqrt(meanA0) );

        _coeffStiff = m/M_physics -> Data() -> DensityRho();

        // update the current element
        M_FESpace -> fe().updateFirstDerivQuadPt(M_FESpace -> mesh() -> edgeList(iedge));

        mass (   _coeffMass,  _elmatMassLHS,  M_FESpace -> fe(), 0, 0 );
        stiff(   _coeffStiff, _elmatStiffLHS, M_FESpace -> fe(), 0, 0 );
        stiff( - _coeffStiff, _elmatStiffRHS, M_FESpace -> fe(), 0, 0 );

        // assemble the mass and grad matrices
        assemb_mat( _matrixLHS, _elmatMassLHS, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );
        assemb_mat( _matrixLHS, _elmatStiffLHS, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );

        assemb_mat( _stiffRHS, _elmatStiffRHS, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );

        Debug( 6310 ) << "\n\tm = "           << m
        << "\n\t_coeffMass = "  << _coeffMass
        << "\n\t_coeffStiff = " << _coeffStiff << "\n";
    } // end loop on elements

    // update rhs
    //_stiffRHS.Axpy( 1., flux , 0., _rhs );

    _rhs = _stiffRHS*flux;

    UInt firstDof = 1;
    UInt lastDof  = _rhs.size();

    // symmetric treatment (cholesky can be used)
    // first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    // second row modified (for symmetry)
    //  _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    // last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    // penultimate row modified (for symmetry)
    //  _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    updateBCDirichletMatrix(_matrixLHS);

    //_tridiagsolver.Factor( _matrixLHS );

    // cholesky or lapack lu solve
    // solve the system: rhs1 = massFactor^{-1} * rhs1
    //_tridiagsolver.Solve( _matrixLHS, _rhs );

    vector_Type _sol(_rhs);

    M_linearSolver -> setMatrix(_matrixLHS);
    //@int numIter = M_linearSolver -> solveSystem( _rhs, _sol, _matrixLHS, true);

    //std::cout <<" iterations number :  " << numIter << std::endl;

    return _sol;
}

OneDimensionalModel_Solver::vector_Type
OneDimensionalModel_Solver::viscoelasticFluxCorrection( const vector_Type& flux, const Real& TimeStep )
{
    matrix_Type _matrixLHS(M_FESpace -> map());
    matrix_Type _stiffRHS (M_FESpace -> map());

    //TriDiagCholesky< Real, matrix_Type, Vector > _tridiagsolver(M_FESpace -> dim());

    ElemMat _elmatMassLHS  (M_FESpace -> fe().nbNode,1,1);
    ElemMat _elmatStiffLHS (M_FESpace -> fe().nbNode,1,1);
    ElemMat _elmatStiffRHS (M_FESpace -> fe().nbNode,1,1);

    vector_Type _rhs(M_FESpace -> map());

    Real _coeffMass;
    Real _coeffStiff;

    Real gamma, meanA0(1.);

    _matrixLHS *= 0.;
    _stiffRHS  *= 0.;

    // Elementary computation and matrix assembling
    // Loop on elements
    for (UInt iedge = 1; iedge <= M_FESpace -> mesh() -> numEdges(); ++iedge)
    {
        // set the elementary matrices to 0.
        _elmatMassLHS.zero();
        _elmatStiffLHS.zero();
        _elmatStiffRHS.zero();

        // this comes from the exact derivation of generalized string model
        // + Voigt viscoelasticity
        _coeffMass = M_rhs[0]( iedge - 1 ) + M_rhs[0]( iedge );
        _coeffMass *= 0.5;
        _coeffMass = 1./std::sqrt(_coeffMass);

        if (M_physics -> Data() -> linearizeStringModel())
        {
            // this is to recover the linearized version (_coeffMass = 1/A)
            _coeffMass *= _coeffMass;

            meanA0  = M_physics -> Data() -> Area0( iedge - 1 ) + M_physics -> Data() -> Area0( iedge );
            meanA0 *= 0.5;

            if (M_physics -> Data() -> linearizeEquations())
            {
                // when using linearized equations, A \simeq A0
                _coeffMass = 1./ meanA0;
            }
        }
        gamma = M_physics -> Data() -> ViscoelasticModulus() / ( 2 * std::sqrt(4*std::atan(1)) );
        gamma *= 1 / std::sqrt( meanA0 );
        _coeffStiff = TimeStep * gamma / M_physics -> Data() -> DensityRho();

        // update the current element
        M_FESpace -> fe().updateFirstDerivQuadPt(M_FESpace -> mesh() -> edgeList(iedge));
        //std::cout << M_FESpace -> fe().currentId() << std::endl;

        mass (   _coeffMass,  _elmatMassLHS,  M_FESpace -> fe(),0, 0 );
        stiff(   _coeffStiff, _elmatStiffLHS, M_FESpace -> fe(),0, 0 );
        stiff( - _coeffStiff, _elmatStiffRHS, M_FESpace -> fe(),0, 0 );

        // assemble the mass and grad matrices
        assemb_mat( _matrixLHS, _elmatMassLHS, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );
        assemb_mat( _matrixLHS, _elmatStiffLHS, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );

        assemb_mat( _stiffRHS, _elmatStiffRHS, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );

        Debug( 6310 ) << "\n\tgamma = " << gamma
        << "\n\t_coeffMass = " << _coeffMass
        << "\n\t_coeffStiff = " << _coeffStiff << "\n";

    } // end loop on elements

    // update rhs
    // rhs = _stiffRHS * rhs
    _rhs = _stiffRHS*flux;

    // NOTE I should add to rhs the value of boundary integral ("neumann-like")
    // BUT it's useless since I'm going to impose dirichlet conditions on all boundaries

    UInt firstDof = 1;
    UInt lastDof  = _rhs.size();

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

    updateBCDirichletMatrix(_matrixLHS);

    vector_Type _sol(_rhs);

    M_linearSolver -> setMatrix(_matrixLHS);
    //int numIter = M_linearSolver -> solveSystem( _rhs, _sol, _matrixLHS, true);

    return _sol;
}

OneDimensionalModel_Solver::vector_Type
OneDimensionalModel_Solver::longitudinalFluxCorrection()
{
    matrix_Type _massLHS(M_FESpace -> map());
    matrix_Type _massRHS(M_FESpace -> map());

    //TriDiagCholesky< Real, matrix_Type, Vector > _tridiagsolver(M_FESpace -> dim());

    ElemMat _elmatMassLHS (M_FESpace -> fe().nbNode,1,1);
    ElemMat _elmatMassRHS (M_FESpace -> fe().nbNode,1,1);

    vector_Type _rhs(M_FESpace -> map());
    // let g = sqrt(A) - sqrt(A0)
    // now f = _d3g_dz3
    vector_Type _g(M_FESpace -> map());
    vector_Type _f(M_FESpace -> map());

    //          _g = M_rhs[0];
    for ( UInt i=0; i<M_FESpace -> dim(); ++i )
        _g(i) = std::sqrt(M_rhs[0](i)) - std::sqrt(M_physics -> Data() -> Area0(i));

    UInt inode;

    Real _coeffMassLHS;
    Real _coeffMassRHS;

    Real _a;
    //    std::ostringstream output;
    _massLHS *= 0.;
    _massRHS *= 0.;

    Real _h( M_physics -> Data() -> Length() / static_cast<Real>(M_physics -> Data() -> NumberOfElements() - 1) );

    // Elementary computation and matrix assembling
    // Loop on elements
    for (UInt iedge = 1; iedge <= M_FESpace -> mesh() -> numEdges(); ++iedge)
    {
        inode = iedge - 1;

        // set the elementary matrices to 0.
        _elmatMassLHS.zero();
        _elmatMassRHS.zero();

        // coeff (1/A) (average values over the element)
        _coeffMassLHS = M_rhs[0]( inode ) + M_rhs[0]( inode+1 );
        _coeffMassLHS /= 2;
        _coeffMassLHS = 1./_coeffMassLHS;

        _a = M_physics -> Data() -> InertialModulus() / std::sqrt(4*std::atan(1));
        _coeffMassRHS = M_physics -> Data() -> dataTime() -> getTimeStep() * _a / M_physics -> Data() -> DensityRho();

        // backward differentiation when near to the left boundary
        // d3 A (xi)/ dz3 = (1/h^3) * ( -A(xi) + A(xi+3) - 3A(xi+2) + 3A(xi+1) )

        // forward differentiation when near to the right boundary
        // d3 A (xi)/ dz3 = (1/h^3) * ( A(xi) - A(xi-3) + 3A(xi-2) - 3A(xi-1) )

        // central differentiation otherwise
        // d3 A (xi)/ dz3 = (1/h^3) * ( -A(xi-2) + 2A(xi-1) - 2A(xi+1) + A(xi+2) )

        Debug( 6310 ) << "\ninode = " << inode << "\n";
        if (inode<2)
        { // backward differentiation
            _f( inode ) = -_g( inode ) + _g( inode+3 )
                          - 3*_g( inode+2 ) + 3*_g( inode+1 );
            Debug( 6310 ) << "\n\tbackward differentiation = " << _coeffMassLHS << "\n";
        }
        else if (inode>M_FESpace -> mesh() -> numEdges()-2)
        { // forward differentiation
            _f( inode ) = _g( inode ) - _g( inode-3 )
                          + 3*_g( inode-2 ) - 3*_g( inode-1 );
            Debug( 6310 ) << "\n\forward differentiation = " << _coeffMassLHS << "\n";
        }
        else
        { // central differentiation
            _f( inode ) = -_g( inode - 2 ) + 2*_g( inode-1 )
                          - 2*_g( inode+1 ) + _g( inode+2 );
            Debug( 6310 ) << "\n\tcentral differentiation = " << _coeffMassLHS << "\n";
        }

        _f(inode) *= 1/(2*std::pow(_h,3));

        // update the current element
        M_FESpace -> fe().updateFirstDerivQuadPt(M_FESpace -> mesh() -> edgeList(iedge));

        mass( _coeffMassLHS, _elmatMassLHS, M_FESpace -> fe(),0, 0 );
        mass( _coeffMassRHS, _elmatMassRHS, M_FESpace -> fe(),0, 0 );

        // assemble the mass and grad matrices
        assemb_mat( _massLHS, _elmatMassLHS, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );
        assemb_mat( _massRHS, _elmatMassRHS, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );

        Debug( 6310 ) << "\n\t_coeffMassLHS = " << _coeffMassLHS << "\n";
        Debug( 6310 ) << "\n\t_coeffMassRHS = " << _coeffMassRHS << "\n";

    } // end loop on elements

    // update rhs
    //    _massLHS.Axpy( 1., (*solution["P"]) , 0., _rhs );
    _rhs = _massRHS*_f;

    UInt firstDof = 1;
    UInt lastDof  = _rhs.size();

    // symmetric treatment (cholesky can be used)
    // first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    // second row modified (for symmetry)
    //    _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    // last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    // penultimate row modified (for symmetry)
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    updateBCDirichletMatrix(_massLHS);

    //@_tridiagsolver.Factor( _massLHS );

    // cholesky or lapack lu solve
    // solve the system: rhs1 = massFactor^{-1} * rhs1
    //@_tridiagsolver.Solve( _massLHS, _rhs );

    return _rhs;
}

/*
template< class Params, class Flux, class Source >
ScalVec
OneDimensionalModel_Solver::_compute_d2Q_dx2( const ScalVec& flux )
{
    matrix_Type _massLHS (M_FESpace -> dim());
    matrix_Type _stiffRHS(M_FESpace -> dim());

    TriDiagCholesky< Real, matrix_Type, Vector > _tridiagsolver(M_FESpace -> dim());

    ElemMat _elmatMassLHS (M_FESpace -> fe().nbNode,1,1);
    ElemMat _elmatStiffRHS (M_FESpace -> fe().nbNode,1,1);

    ScalVec _rhs(M_FESpace -> dim());

    _massLHS  *= 0.;
    _stiffRHS *= 0.;

    // Elementary computation and matrix assembling
    // Loop on elements
    for(UInt iedge = 1; iedge <= M_FESpace -> mesh() -> numEdges(); ++iedge){

        // set the elementary matrices to 0.
        _elmatMassLHS *= 0.;
        _elmatStiffRHS *= 0.;

        // update the current element
        M_FESpace -> fe().updateFirstDerivQuadPt(M_FESpace -> mesh() -> edgeList(iedge));
        //std::cout << M_FESpace -> fe().currentId() << std::endl;

        mass( 1., _elmatMassLHS, M_FESpace -> fe(),0, 0 );
        stiff( -1., _elmatStiffRHS, M_FESpace -> fe(),0, 0 );

        // assemble the mass and grad matrices
        assemb_mat( _massLHS, _elmatMassLHS, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );
        assemb_mat( _stiffRHS, _elmatStiffRHS, M_FESpace -> fe(), M_FESpace -> dof() , 0, 0 );

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
