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
 *  @brief File containing some functions for the boundary conditions of 1D models.
 *
 *  @version 1.0
 *  @author Lucia Mirabella  <lucia.mirabella@gmail.com>
 *  @date 01-08-2006
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 20-04-2010
 *  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 */

#include <lifemc/lifefem/OneDimensionalModel_BCFunction_Default.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BCFunction_Default::OneDimensionalModel_BCFunction_Default( const OneD_BCSide&     bcSide,
                                                                                const OneD_BC&         bcType ):
        M_Flux                          (),
        M_Source                        (),
        M_bcNode                        (),
        M_bcSide                        ( bcSide ),
        M_bcType                        ( bcType )
{
}

OneDimensionalModel_BCFunction_Default::OneDimensionalModel_BCFunction_Default( const OneDimensionalModel_BCFunction_Default& BCF_Default ) :
        M_Flux                          ( BCF_Default.M_Flux ),        // Ptr copy
        M_Source                        ( BCF_Default.M_Source ),      // Ptr copy
        M_Solution                      ( BCF_Default.M_Solution ),    // Ptr copy
        M_bcNode                        ( BCF_Default.M_bcNode ),
        M_bcType                        ( BCF_Default.M_bcType )
{}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalModel_BCFunction_Default::operator() ( const Real& /*time*/, const Real& /*timeStep*/ )
{
    assert( false );
    return 0.;
}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_BCFunction_Default::setFluxSource( const Flux_PtrType& flux, const Source_PtrType& source )
{
    M_Flux   = flux;
    M_Source = source;

    this->setupNode();
}

// ===================================================
// Protected Methods
// ===================================================
void
OneDimensionalModel_BCFunction_Default::setupNode()
{
    switch ( M_bcSide )
    {
    case OneD_left:
        M_bcNode = 1;
        break;

    case OneD_right:
        M_bcNode = M_Flux->physics()->data()->numberOfNodes();
        break;

    default:
        std::cout << "Warning: bcSide \"" << M_bcSide << "\" not available!" << std::endl;
    }
}



// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BCFunction_Riemann::OneDimensionalModel_BCFunction_Riemann( const OneD_BCSide&     bcSide,
                                                                                const OneD_BC&         bcType ) :
        super                           ( bcSide, bcType ),
        M_bcU                           (),
        M_bcW                           ()
{}

OneDimensionalModel_BCFunction_Riemann::OneDimensionalModel_BCFunction_Riemann( const OneDimensionalModel_BCFunction_Riemann& BCF_Riemann ) :
        super                           ( BCF_Riemann ),
        M_bcU                           ( BCF_Riemann.M_bcU ),
        M_bcW                           ( BCF_Riemann.M_bcW )
{}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalModel_BCFunction_Riemann::operator()( const Real& /*time*/, const Real& /*timeStep*/ )
{
    updateBCVariables();

    return ( ( M_bcType == OneD_W1 ) ? M_bcW[0] : M_bcW[1] );
}

// ===================================================
// Protected Methods
// ===================================================
void
OneDimensionalModel_BCFunction_Riemann::updateBCVariables()
{
    M_bcU[0] = (*(*M_Solution)["A"])(M_bcNode);
    M_bcU[1] = (*(*M_Solution)["Q"])(M_bcNode);
    M_bcW[0] = (*(*M_Solution)["W1"])(M_bcNode);
    M_bcW[1] = (*(*M_Solution)["W2"])(M_bcNode);
}



// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BCFunction_Compatibility::OneDimensionalModel_BCFunction_Compatibility( const OneD_BCSide&     bcSide,
        const OneD_BC&         bcType ):
        super                           ( bcSide, bcType ),
        M_bcInternalNode                (),
        M_boundaryPoint                 (),
        M_internalBdPoint               (),
        M_eigenvalues                   (),
        M_deltaEigenvalues              (),
        M_leftEigenvector1              (),
        M_leftEigenvector2              (),
        M_deltaLeftEigenvector1         (),
        M_deltaLeftEigenvector2         ()
{
}

OneDimensionalModel_BCFunction_Compatibility::OneDimensionalModel_BCFunction_Compatibility( const OneDimensionalModel_BCFunction_Compatibility& BCF_Compatibility ) :
        super                           ( BCF_Compatibility ),
        M_bcInternalNode                ( BCF_Compatibility.M_bcInternalNode ),
        M_boundaryPoint                 ( BCF_Compatibility.M_boundaryPoint ),
        M_internalBdPoint               ( BCF_Compatibility.M_internalBdPoint ),
        M_eigenvalues                   ( BCF_Compatibility.M_eigenvalues ),
        M_deltaEigenvalues              ( BCF_Compatibility.M_deltaEigenvalues ),
        M_leftEigenvector1              ( BCF_Compatibility.M_leftEigenvector1 ),
        M_leftEigenvector2              ( BCF_Compatibility.M_leftEigenvector2 ),
        M_deltaLeftEigenvector1         ( BCF_Compatibility.M_deltaLeftEigenvector1 ),
        M_deltaLeftEigenvector2         ( BCF_Compatibility.M_deltaLeftEigenvector2 )
{
}

// ===================================================
// Protected Methods
// ===================================================
void
OneDimensionalModel_BCFunction_Compatibility::setupNode()
{
    super::setupNode();

    Mesh_Type::EdgeType boundaryEdge;
    switch ( M_bcSide )
    {
    case OneD_left:
        M_bcInternalNode      = M_bcNode + 1;
        boundaryEdge          = M_Flux->physics()->data()->mesh()->edgeList(1);
        M_boundaryPoint[0]    = boundaryEdge.point(1).x();
        M_boundaryPoint[1]    = boundaryEdge.point(1).y();
        M_boundaryPoint[2]    = boundaryEdge.point(1).z();
        M_internalBdPoint[0]  = boundaryEdge.point(2).x();
        M_internalBdPoint[1]  = boundaryEdge.point(2).y();
        M_internalBdPoint[2]  = boundaryEdge.point(2).z();
        break;

    case OneD_right:
        M_bcInternalNode      = M_bcNode - 1;
        boundaryEdge          = M_Flux->physics()->data()->mesh()->edgeList(M_bcNode - 1);
        M_boundaryPoint[0]    = boundaryEdge.point(2).x();
        M_boundaryPoint[1]    = boundaryEdge.point(2).y();
        M_boundaryPoint[2]    = boundaryEdge.point(2).z();
        M_internalBdPoint[0]  = boundaryEdge.point(1).x();
        M_internalBdPoint[1]  = boundaryEdge.point(1).y();
        M_internalBdPoint[2]  = boundaryEdge.point(1).z();
        break;

    default:
        std::cout << "Warning: bcSide \"" << M_bcSide << "\" not available!" << std::endl;
    }
}

Real
OneDimensionalModel_BCFunction_Compatibility::computeRHS( const Real& timeStep )
{
    updateBCVariables();
    computeEigenValuesVectors();

    switch ( M_bcType )
    {
    case OneD_W1:
        return evaluateRHS( M_eigenvalues[0], M_leftEigenvector1, M_deltaLeftEigenvector1, timeStep );
        break;

    case OneD_W2:
        return evaluateRHS( M_eigenvalues[1], M_leftEigenvector2, M_deltaLeftEigenvector2, timeStep );
        break;

    default:
        std::cout << "Warning: bcType \"" << M_bcType << "\"not available!" << std::endl;
        return 0.;
    }
}

void
OneDimensionalModel_BCFunction_Compatibility::computeEigenValuesVectors()
{
    M_Flux->eigenValuesEigenVectors( M_bcU[0], M_bcU[1],
                                     M_eigenvalues, M_leftEigenvector1, M_leftEigenvector2,
                                     M_bcNode - 1 );

    M_Flux->deltaEigenValuesEigenVectors( M_bcU[0], M_bcU[1],
                                          M_deltaEigenvalues, M_deltaLeftEigenvector1, M_deltaLeftEigenvector2,
                                          M_bcNode - 1 );
}

Real
OneDimensionalModel_BCFunction_Compatibility::evaluateRHS( const Real& eigenvalue, const container2D_Type& eigenvector, const container2D_Type& deltaEigenvector, const Real& timeStep )
{
    Real cfl = computeCFL( eigenvalue, timeStep );

    container2D_Type U_interpolated;
    U_interpolated[0] = ( 1 - cfl ) * M_bcU[0]  + cfl * (*(*M_Solution)["A"])( M_bcInternalNode );
    U_interpolated[1] = ( 1 - cfl ) * M_bcU[1]  + cfl * (*(*M_Solution)["Q"])( M_bcInternalNode );

    container2D_Type U;

    container2D_Type bcNodes;
    switch ( M_bcSide )
    {
    case OneD_left:
        bcNodes[0] = M_bcNode - 1; // Boundary node
        bcNodes[1] = M_bcNode;     // Inner node
        break;

    case OneD_right:
        bcNodes[0] = M_bcNode - 1; // Boundary node
        bcNodes[1] = M_bcNode - 2; // Inner node
        break;

    default:
        bcNodes[0] = M_bcNode - 1; // Boundary node
        bcNodes[1] = M_bcNode - 1; // Inner node = boundary node -> no interpolation
        std::cout << "Warning: bcSide \"" << M_bcSide << "\" not available!" << std::endl;
    }

    U[0] = U_interpolated[0] - timeStep * M_Source->interpolatedQuasiLinearSource( U_interpolated[0], U_interpolated[1], 1, bcNodes, cfl );
    U[1] = U_interpolated[1] - timeStep * M_Source->interpolatedQuasiLinearSource( U_interpolated[0], U_interpolated[1], 2, bcNodes, cfl );

    return dot( eigenvector, U ) + timeStep * eigenvalue * dot( deltaEigenvector, U_interpolated );
}

Real
OneDimensionalModel_BCFunction_Compatibility::computeCFL( const Real& eigenvalue, const Real& timeStep ) const
{
    Real deltaX(0);
    switch ( M_bcSide )
    {
    case OneD_left:
        deltaX = M_Flux->physics()->data()->mesh()->edgeLength( 0 );
        break;

    case OneD_right:
        deltaX = M_Flux->physics()->data()->mesh()->edgeLength( M_bcNode - 2 );
        break;

    default:
        std::cout << "Warning: bcSide \"" << M_bcSide << "\" not available!" << std::endl;
    }

    Real cfl = eigenvalue * timeStep / deltaX;

    if ( M_bcInternalNode == 2 ) // the edge is on the left of the domain
    {
        ASSERT( -1. < cfl && cfl < 0. , "This characteristics is wrong!\nEither it is not outcoming (eigenvalue>0 at the left of the domain),\n or CFL is too high.");

    }
    else                         // the edge is on the right of the domain
    {
        ASSERT( 0. < cfl && cfl < 1. , "This characteristics is wrong!\nEither it is not outcoming (eigenvalue<0 at the right of the domain),\n or CFL is too high.");
    }

    return std::abs(cfl);
}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalModel_BCFunction_Absorbing::operator()( const Real& /*time*/, const Real& timeStep )
{
    updateBCVariables();
    computeEigenValuesVectors();

    Real W_out(0.), W_out_old(0.);
    switch ( M_bcType )
    {
    case OneD_W1:
        W_out = M_bcW[1] + evaluateRHS( M_eigenvalues[1], M_leftEigenvector2, M_deltaLeftEigenvector2, timeStep ) - dot( M_leftEigenvector2, M_bcU );
        W_out_old = -M_bcW[0] + dot( M_leftEigenvector1, M_bcU );
        break;

    case OneD_W2:
        W_out = M_bcW[0] + evaluateRHS( M_eigenvalues[0], M_leftEigenvector1, M_deltaLeftEigenvector1, timeStep ) - dot( M_leftEigenvector1, M_bcU );
        W_out_old = -M_bcW[1] + dot( M_leftEigenvector2, M_bcU );
        break;
    default:
        std::cout << "Warning: bcType \"" << M_bcType  << "\"not available!" << std::endl;
    }

    Real a1, a2, a11, a22, b1, b2, c1, c2;
    a1 = M_Flux->physics()->elasticPressure(M_bcU[0], M_bcNode - 1) - M_Flux->physics()->data()->externalPressure(); // pressure at previous time step
    a2 = M_bcU[1]; // flux at previous time step

    b1 = M_Flux->physics()->dPdW( M_bcW[0], M_bcW[1], 1, M_bcNode - 1);  // dP / dW1
    b2 = M_bcU[0] / 2; // dQ / dW1

    c1 = M_Flux->physics()->dPdW( M_bcW[0], M_bcW[1], 2, M_bcNode - 1);  // dP / dW2
    c2 = b2; // dQ / dW2

    a11 = a1 - b1*M_bcW[0] - c1*M_bcW[1];
    a22 = a2 - b2*M_bcW[0] - c2*M_bcW[1];

    Real resistance(b1 / b2);

    this->resistance( resistance );

    return W_out_old + W_out * (b2*resistance-b1)/(c1-c2*resistance) + (a22*resistance-a11)/(c1-c2*resistance);
}

// ===================================================
// Protected Methods
// ===================================================
/*void
OneDimensionalModel_BCFunction_Absorbing::resistance( Real& resistance )
{
    //Do nothing => absorbing!
}*/



// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BCFunction_Resistance::OneDimensionalModel_BCFunction_Resistance( const OneD_BCSide&     bcSide,
        const OneD_BC&         bcType,
        const Real&            resistance ):
        super                           ( bcSide, bcType ),
        M_resistance                    ( resistance )
{}

OneDimensionalModel_BCFunction_Resistance::OneDimensionalModel_BCFunction_Resistance( const OneDimensionalModel_BCFunction_Resistance& BCF_Resistance ) :
        super                           ( BCF_Resistance ),
        M_resistance                    ( BCF_Resistance.M_resistance )
{}
// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BCFunction_Windkessel3::OneDimensionalModel_BCFunction_Windkessel3( const OneD_BCSide&     bcSide,
        const OneD_BC&         bcType,
        const Real&            resistance1,
        const Real&            resistance2,
        const Real&            compliance,
        const bool&            absorbing1,
        const Real&            venousPressure ):
        super                           ( bcSide, bcType ),
        M_resistance1                   ( resistance1 ),
        M_resistance2                   ( resistance2 ),
        M_compliance                    ( compliance ),
        M_absorbing1                    ( absorbing1 ),
        M_venousPressure                ( venousPressure ),
        M_P0                            ( 0. ),
        M_Q_tn                          ( 0. ),
        M_dQdt_tn                       ( 0. ),
        M_integral_tn                   ( 0. )
{}

OneDimensionalModel_BCFunction_Windkessel3::OneDimensionalModel_BCFunction_Windkessel3( const OneDimensionalModel_BCFunction_Windkessel3& BCF_Windkessel3 ) :
        super                           ( BCF_Windkessel3 ),
        M_resistance1                   ( BCF_Windkessel3.M_resistance1 ),
        M_resistance2                   ( BCF_Windkessel3.M_resistance2 ),
        M_compliance                    ( BCF_Windkessel3.M_compliance ),
        M_absorbing1                    ( BCF_Windkessel3.M_absorbing1 ),
        M_venousPressure                ( BCF_Windkessel3.M_venousPressure ),
        M_P0                            ( BCF_Windkessel3.M_P0 ),
        M_Q_tn                          ( BCF_Windkessel3.M_Q_tn ),
        M_dQdt_tn                       ( BCF_Windkessel3.M_dQdt_tn ),
        M_integral_tn                   ( BCF_Windkessel3.M_integral_tn )
{}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalModel_BCFunction_Windkessel3::operator()( const Real& time, const Real& timeStep )
{
    UInt W_outID;
    Real W_out;

    switch ( M_bcType )
    {
    case OneD_W1:
        W_outID = 2;
        W_out = computeRHS( timeStep );
        break;
    case OneD_W2:
        W_outID = 1;
        W_out = computeRHS( timeStep );
        break;
    default:
        std::cout << "Warning: bcType \"" << M_bcType  << "\"not available!" << std::endl;
    }

    Real A( M_bcU[0] );
    Real Q( M_bcU[1] );

    if ( M_absorbing1 )
    {
        Real b1( M_Flux->physics()->dPdW( A, Q, 1, M_bcNode) );  // dP / dW1
        Real b2( A / 2 ); // dQ / dW1
        M_resistance1 = b1 / b2;
    }

    //  Trapezoidal rule for given integral
    //
    //  \int_0^t ( pv/(R2*C) + (R1+R2)/(R2*C) Q(s) + R1 dQ(s)/ds ) exp( s / (R2*C) ) ds
    //
    //  this gives: \int_0^t(n+1) = \int_0^t(n) + \int_t(n)^t(n+1)
    //
    //  \int_0^t(n) is stored from previous time step
    //  \int_t(n)^t(n+1) f(s) exp( a * s ) ds =
    //       =  dt/2 * ( f(t(n+1)) exp( a * t(n+1) ) + f(t(n)) exp( a * t(n) ) ) + O( dt^2 )
    //       =  dt/2 * [ exp(a*t(n)) * ( f(t(n+1))*exp(a*dt) + f(t(n)) ) ] + O( dt^2 )
    Real a        ( 1 / ( M_compliance * M_resistance2 ) );
    Real dQdt     ( ( Q - M_Q_tn ) / timeStep );
    Real F        ( a * M_venousPressure + a * (M_resistance1 + M_resistance2) * Q + M_resistance1 * dQdt );
    Real Fn       ( a * M_venousPressure + a * (M_resistance1 + M_resistance2) * M_Q_tn + M_resistance1 * M_dQdt_tn );
    Real EXP      ( std::exp( a * time     ) );
    Real EXPdt    ( std::exp( a * timeStep ) );
    Real integral ( M_integral_tn + ( timeStep / 2 ) * ( EXP * ( F * EXPdt + Fn ) ) );

    // Compute the solution of circuital ODE
    // P(t) = P(0) + [ \int_0^t ( pv/(R2*C) + (R1+R2)/(R2*C) Q(s)
    //                 + R1 dQ(s)/ds ) exp( s / (R2*C) ) ds ] * exp( - t / (R2*C) )
    Real P        ( M_P0 + integral / EXP );

    // Update variable for the next time step
    M_integral_tn = integral;
    M_dQdt_tn     = dQdt;
    M_Q_tn        = Q;

    return  M_Flux->physics()->fromPToW( P, W_out, W_outID, M_bcNode ); // W_in
}

}
