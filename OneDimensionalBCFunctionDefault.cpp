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
 *  @date 01-08-2006
 *  @author Lucia Mirabella  <lucia.mirabella@gmail.com>
 *
 *  @version 2.0
 *  @date 20-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifefem/OneDimensionalModel_BCFunction_Default.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalBCFunctionDefault::OneDimensionalBCFunctionDefault( const bcSide_Type& bcSide, const bcType_Type& bcType ):
        M_flux                          (),
        M_source                        (),
        M_bcNode                        (),
        M_bcSide                        ( bcSide ),
        M_bcType                        ( bcType )
{
}

OneDimensionalBCFunctionDefault::OneDimensionalBCFunctionDefault( const OneDimensionalBCFunctionDefault& bcFunctionDefault ) :
        M_flux                          ( bcFunctionDefault.M_flux ),        // Ptr copy
        M_source                        ( bcFunctionDefault.M_source ),      // Ptr copy
        M_solution                      ( bcFunctionDefault.M_solution ),    // Ptr copy
        M_bcNode                        ( bcFunctionDefault.M_bcNode ),
        M_bcType                        ( bcFunctionDefault.M_bcType )
{}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalBCFunctionDefault::operator() ( const Real& /*time*/, const Real& /*timeStep*/ )
{
#ifdef HAVE_LIFEV_DEBUG
    assert( false );
#endif
    return 0.;
}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalBCFunctionDefault::setFluxSource( const fluxPtr_Type& flux, const sourcePtr_Type& source )
{
    M_flux   = flux;
    M_source = source;

    this->setupNode();
}

// ===================================================
// Protected Methods
// ===================================================
void
OneDimensionalBCFunctionDefault::setupNode()
{
    switch ( M_bcSide )
    {
    case OneDimensional::left:
        M_bcNode = 1;
        break;

    case OneDimensional::right:
        M_bcNode = M_flux->physics()->data()->numberOfNodes();
        break;

    default:
        std::cout << "Warning: bcSide \"" << M_bcSide << "\" not available!" << std::endl;
    }
}



// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalBCFunctionRiemann::OneDimensionalBCFunctionRiemann( const bcSide_Type& bcSide, const bcType_Type& bcType ) :
        super                           ( bcSide, bcType ),
        M_bcU                           (),
        M_bcW                           ()
{}

OneDimensionalBCFunctionRiemann::OneDimensionalBCFunctionRiemann( const OneDimensionalBCFunctionRiemann& bcFunctionRiemann ) :
        super                           ( bcFunctionRiemann ),
        M_bcU                           ( bcFunctionRiemann.M_bcU ),
        M_bcW                           ( bcFunctionRiemann.M_bcW )
{}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalBCFunctionRiemann::operator()( const Real& /*time*/, const Real& /*timeStep*/ )
{
    updateBCVariables();

    return ( ( M_bcType == OneDimensional::W1 ) ? M_bcW[0] : M_bcW[1] );
}

// ===================================================
// Protected Methods
// ===================================================
void
OneDimensionalBCFunctionRiemann::updateBCVariables()
{
    M_bcU[0] = (*(*M_solution)["A"])(M_bcNode);
    M_bcU[1] = (*(*M_solution)["Q"])(M_bcNode);
    M_bcW[0] = (*(*M_solution)["W1"])(M_bcNode);
    M_bcW[1] = (*(*M_solution)["W2"])(M_bcNode);
}



// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalBCFunctionCompatibility::OneDimensionalBCFunctionCompatibility( const bcSide_Type& bcSide, const bcType_Type& bcType ):
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

OneDimensionalBCFunctionCompatibility::OneDimensionalBCFunctionCompatibility( const OneDimensionalBCFunctionCompatibility& bcFunctionCompatibility ) :
        super                           ( bcFunctionCompatibility ),
        M_bcInternalNode                ( bcFunctionCompatibility.M_bcInternalNode ),
        M_boundaryPoint                 ( bcFunctionCompatibility.M_boundaryPoint ),
        M_internalBdPoint               ( bcFunctionCompatibility.M_internalBdPoint ),
        M_eigenvalues                   ( bcFunctionCompatibility.M_eigenvalues ),
        M_deltaEigenvalues              ( bcFunctionCompatibility.M_deltaEigenvalues ),
        M_leftEigenvector1              ( bcFunctionCompatibility.M_leftEigenvector1 ),
        M_leftEigenvector2              ( bcFunctionCompatibility.M_leftEigenvector2 ),
        M_deltaLeftEigenvector1         ( bcFunctionCompatibility.M_deltaLeftEigenvector1 ),
        M_deltaLeftEigenvector2         ( bcFunctionCompatibility.M_deltaLeftEigenvector2 )
{
}

// ===================================================
// Protected Methods
// ===================================================
void
OneDimensionalBCFunctionCompatibility::setupNode()
{
    super::setupNode();

    mesh_Type::EdgeType boundaryEdge;
    switch ( M_bcSide )
    {
    case OneDimensional::left:
        M_bcInternalNode      = M_bcNode + 1;
        boundaryEdge          = M_flux->physics()->data()->mesh()->edgeList(1);
        M_boundaryPoint[0]    = boundaryEdge.point(1).x();
        M_boundaryPoint[1]    = boundaryEdge.point(1).y();
        M_boundaryPoint[2]    = boundaryEdge.point(1).z();
        M_internalBdPoint[0]  = boundaryEdge.point(2).x();
        M_internalBdPoint[1]  = boundaryEdge.point(2).y();
        M_internalBdPoint[2]  = boundaryEdge.point(2).z();
        break;

    case OneDimensional::right:
        M_bcInternalNode      = M_bcNode - 1;
        boundaryEdge          = M_flux->physics()->data()->mesh()->edgeList(M_bcNode - 1);
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
OneDimensionalBCFunctionCompatibility::computeRHS( const Real& timeStep )
{
    updateBCVariables();
    computeEigenValuesVectors();

    switch ( M_bcType )
    {
    case OneDimensional::W1:
        return evaluateRHS( M_eigenvalues[0], M_leftEigenvector1, M_deltaLeftEigenvector1, timeStep );
        break;

    case OneDimensional::W2:
        return evaluateRHS( M_eigenvalues[1], M_leftEigenvector2, M_deltaLeftEigenvector2, timeStep );
        break;

    default:
        std::cout << "Warning: bcType \"" << M_bcType << "\"not available!" << std::endl;
        return 0.;
    }
}

void
OneDimensionalBCFunctionCompatibility::computeEigenValuesVectors()
{
    M_flux->eigenValuesEigenVectors( M_bcU[0], M_bcU[1],
                                     M_eigenvalues, M_leftEigenvector1, M_leftEigenvector2,
                                     M_bcNode - 1 );

    M_flux->deltaEigenValuesEigenVectors( M_bcU[0], M_bcU[1],
                                          M_deltaEigenvalues, M_deltaLeftEigenvector1, M_deltaLeftEigenvector2,
                                          M_bcNode - 1 );
}

Real
OneDimensionalBCFunctionCompatibility::evaluateRHS( const Real& eigenvalue, const container2D_Type& eigenvector,
                                                           const container2D_Type& deltaEigenvector, const Real& timeStep )
{
    Real cfl = computeCFL( eigenvalue, timeStep );

    container2D_Type U_interpolated;
    U_interpolated[0] = ( 1 - cfl ) * M_bcU[0]  + cfl * (*(*M_solution)["A"])( M_bcInternalNode );
    U_interpolated[1] = ( 1 - cfl ) * M_bcU[1]  + cfl * (*(*M_solution)["Q"])( M_bcInternalNode );

    container2D_Type U;

    container2D_Type bcNodes;
    switch ( M_bcSide )
    {
    case OneDimensional::left:
        bcNodes[0] = M_bcNode - 1; // Boundary node
        bcNodes[1] = M_bcNode;     // Inner node
        break;

    case OneDimensional::right:
        bcNodes[0] = M_bcNode - 1; // Boundary node
        bcNodes[1] = M_bcNode - 2; // Inner node
        break;

    default:
        bcNodes[0] = M_bcNode - 1; // Boundary node
        bcNodes[1] = M_bcNode - 1; // Inner node = boundary node -> no interpolation
        std::cout << "Warning: bcSide \"" << M_bcSide << "\" not available!" << std::endl;
    }

    U[0] = U_interpolated[0] - timeStep * M_source->interpolatedQuasiLinearSource( U_interpolated[0], U_interpolated[1], 1, bcNodes, cfl );
    U[1] = U_interpolated[1] - timeStep * M_source->interpolatedQuasiLinearSource( U_interpolated[0], U_interpolated[1], 2, bcNodes, cfl );

    return scalarProduct( eigenvector, U ) + timeStep * eigenvalue * scalarProduct( deltaEigenvector, U_interpolated );
}

Real
OneDimensionalBCFunctionCompatibility::computeCFL( const Real& eigenvalue, const Real& timeStep ) const
{
    Real deltaX(0);
    switch ( M_bcSide )
    {
    case OneDimensional::left:
        deltaX = M_flux->physics()->data()->mesh()->edgeLength( 0 );
        break;

    case OneDimensional::right:
        deltaX = M_flux->physics()->data()->mesh()->edgeLength( M_bcNode - 2 );
        break;

    default:
        std::cout << "Warning: bcSide \"" << M_bcSide << "\" not available!" << std::endl;
    }

    Real cfl = eigenvalue * timeStep / deltaX;

#ifdef HAVE_LIFEV_DEBUG
    if ( M_bcInternalNode == 2 ) // the edge is on the left of the domain
    {
        ASSERT( -1. < cfl && cfl < 0. , "This characteristics is wrong!\nEither it is not outcoming (eigenvalue>0 at the left of the domain),\n or CFL is too high.");

    }
    else                         // the edge is on the right of the domain
    {
        ASSERT( 0. < cfl && cfl < 1. , "This characteristics is wrong!\nEither it is not outcoming (eigenvalue<0 at the right of the domain),\n or CFL is too high.");
    }
#endif

    return std::abs(cfl);
}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalBCFunctionAbsorbing::operator()( const Real& /*time*/, const Real& timeStep )
{
    updateBCVariables();
    computeEigenValuesVectors();

    Real W_out(0.), W_out_old(0.);
    switch ( M_bcType )
    {
    case OneDimensional::W1:
        W_out = M_bcW[1] + evaluateRHS( M_eigenvalues[1], M_leftEigenvector2, M_deltaLeftEigenvector2, timeStep ) - scalarProduct( M_leftEigenvector2, M_bcU );
        W_out_old = -M_bcW[0] + scalarProduct( M_leftEigenvector1, M_bcU );
        break;

    case OneDimensional::W2:
        W_out = M_bcW[0] + evaluateRHS( M_eigenvalues[0], M_leftEigenvector1, M_deltaLeftEigenvector1, timeStep ) - scalarProduct( M_leftEigenvector1, M_bcU );
        W_out_old = -M_bcW[1] + scalarProduct( M_leftEigenvector2, M_bcU );
        break;
    default:
        std::cout << "Warning: bcType \"" << M_bcType  << "\"not available!" << std::endl;
    }

    Real a1, a2, a11, a22, b1, b2, c1, c2;
    a1 = M_flux->physics()->elasticPressure(M_bcU[0], M_bcNode - 1) - M_flux->physics()->data()->externalPressure(); // pressure at previous time step
    a2 = M_bcU[1]; // flux at previous time step

    b1 = M_flux->physics()->dPdW( M_bcW[0], M_bcW[1], 1, M_bcNode - 1);  // dP / dW1
    b2 = M_bcU[0] / 2; // dQ / dW1

    c1 = M_flux->physics()->dPdW( M_bcW[0], M_bcW[1], 2, M_bcNode - 1);  // dP / dW2
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
OneDimensionalBCFunctionAbsorbing::resistance( Real& resistance )
{
    //Do nothing => absorbing!
}*/



// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalBCFunctionResistance::OneDimensionalBCFunctionResistance( const bcSide_Type& bcSide, const bcType_Type& bcType, const Real& resistance ):
        super                           ( bcSide, bcType ),
        M_resistance                    ( resistance )
{}

OneDimensionalBCFunctionResistance::OneDimensionalBCFunctionResistance( const OneDimensionalBCFunctionResistance& bcFunctionResistance ) :
        super                           ( bcFunctionResistance ),
        M_resistance                    ( bcFunctionResistance.M_resistance )
{}
// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalBCFunctionWindkessel3::OneDimensionalBCFunctionWindkessel3( const bcSide_Type& bcSide,
                                                                          const bcType_Type& bcType,
                                                                          const Real&        resistance1,
                                                                          const Real&        resistance2,
                                                                          const Real&        compliance,
                                                                          const bool&        absorbing1,
                                                                          const Real&        venousPressure ):
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

OneDimensionalBCFunctionWindkessel3::OneDimensionalBCFunctionWindkessel3( const OneDimensionalBCFunctionWindkessel3& bcFunctionWindkessel3 ) :
        super                           ( bcFunctionWindkessel3 ),
        M_resistance1                   ( bcFunctionWindkessel3.M_resistance1 ),
        M_resistance2                   ( bcFunctionWindkessel3.M_resistance2 ),
        M_compliance                    ( bcFunctionWindkessel3.M_compliance ),
        M_absorbing1                    ( bcFunctionWindkessel3.M_absorbing1 ),
        M_venousPressure                ( bcFunctionWindkessel3.M_venousPressure ),
        M_P0                            ( bcFunctionWindkessel3.M_P0 ),
        M_Q_tn                          ( bcFunctionWindkessel3.M_Q_tn ),
        M_dQdt_tn                       ( bcFunctionWindkessel3.M_dQdt_tn ),
        M_integral_tn                   ( bcFunctionWindkessel3.M_integral_tn )
{}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalBCFunctionWindkessel3::operator()( const Real& time, const Real& timeStep )
{
    UInt W_outID;
    Real W_out;

    switch ( M_bcType )
    {
    case OneDimensional::W1:
        W_outID = 2;
        W_out = computeRHS( timeStep );
        break;
    case OneDimensional::W2:
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
        Real b1( M_flux->physics()->dPdW( A, Q, 1, M_bcNode) );  // dP / dW1
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

    return  M_flux->physics()->fromPToW( P, W_out, W_outID, M_bcNode ); // W_in
}

}
