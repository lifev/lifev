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

#include <lifev/one_d_fsi/function/OneDFSIFunctionSolverDefined.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
OneDFSIFunctionSolverDefined::OneDFSIFunctionSolverDefined ( const bcSide_Type& bcSide, const bcType_Type& bcType ) :
    M_fluxPtr                       (),
    M_sourcePtr                     (),
    M_solutionPtr                   (),
    M_bcNode                        (),
    M_bcSide                        ( bcSide ),
    M_bcType                        ( bcType )
{
}

OneDFSIFunctionSolverDefined::OneDFSIFunctionSolverDefined ( const OneDFSIFunctionSolverDefined& bcFunctionDefault ) :
    M_fluxPtr                       ( bcFunctionDefault.M_fluxPtr ),        // Ptr copy
    M_sourcePtr                     ( bcFunctionDefault.M_sourcePtr ),      // Ptr copy
    M_solutionPtr                   ( bcFunctionDefault.M_solutionPtr ),    // Ptr copy
    M_bcNode                        ( bcFunctionDefault.M_bcNode ),
    M_bcType                        ( bcFunctionDefault.M_bcType )
{}

// ===================================================
// Methods
// ===================================================
Real
OneDFSIFunctionSolverDefined::operator() ( const Real& /*time*/, const Real& /*timeStep*/ )
{
#ifdef HAVE_LIFEV_DEBUG
    assert ( false );
#endif
    return 0.;
}

// ===================================================
// Set Methods
// ===================================================
void
OneDFSIFunctionSolverDefined::setFluxSource ( const fluxPtr_Type& fluxPtr, const sourcePtr_Type& sourcePtr )
{
    M_fluxPtr   = fluxPtr;
    M_sourcePtr = sourcePtr;

    this->setupNode();
}

// ===================================================
// Protected Methods
// ===================================================
void
OneDFSIFunctionSolverDefined::setupNode()
{
    ( M_bcSide == OneDFSI::left ) ? M_bcNode = 0 : M_bcNode = M_fluxPtr->physics()->data()->numberOfNodes() - 1;
}



// ===================================================
// Constructors & Destructor
// ===================================================
OneDFSIFunctionSolverDefinedRiemann::OneDFSIFunctionSolverDefinedRiemann ( const bcSide_Type& bcSide, const bcType_Type& bcType ) :
    super                           ( bcSide, bcType ),
    M_bcU                           (),
    M_bcW                           ()
{}

OneDFSIFunctionSolverDefinedRiemann::OneDFSIFunctionSolverDefinedRiemann ( const OneDFSIFunctionSolverDefinedRiemann& bcFunctionRiemann ) :
    super                           ( bcFunctionRiemann ),
    M_bcU                           ( bcFunctionRiemann.M_bcU ),
    M_bcW                           ( bcFunctionRiemann.M_bcW )
{}

// ===================================================
// Methods
// ===================================================
Real
OneDFSIFunctionSolverDefinedRiemann::operator() ( const Real& /*time*/, const Real& /*timeStep*/ )
{
    updateBCVariables();

    return ( ( M_bcType == OneDFSI::W1 ) ? M_bcW[0] : M_bcW[1] );
}

// ===================================================
// Protected Methods
// ===================================================
void
OneDFSIFunctionSolverDefinedRiemann::updateBCVariables()
{
    M_bcU[0] = (* (*M_solutionPtr) ["A"]) (M_bcNode);
    M_bcU[1] = (* (*M_solutionPtr) ["Q"]) (M_bcNode);
    M_bcW[0] = (* (*M_solutionPtr) ["W1"]) (M_bcNode);
    M_bcW[1] = (* (*M_solutionPtr) ["W2"]) (M_bcNode);
}



// ===================================================
// Constructors & Destructor
// ===================================================
OneDFSIFunctionSolverDefinedCompatibility::OneDFSIFunctionSolverDefinedCompatibility ( const bcSide_Type& bcSide, const bcType_Type& bcType ) :
    super                           ( bcSide, bcType ),
    M_bcElement                     (),
    M_bcInternalNode                (),
    M_eigenvalues                   (),
    M_deltaEigenvalues              (),
    M_leftEigenvector1              (),
    M_leftEigenvector2              (),
    M_deltaLeftEigenvector1         (),
    M_deltaLeftEigenvector2         ()
{
}

OneDFSIFunctionSolverDefinedCompatibility::OneDFSIFunctionSolverDefinedCompatibility ( const OneDFSIFunctionSolverDefinedCompatibility& bcFunctionCompatibility ) :
    super                           ( bcFunctionCompatibility ),
    M_bcElement                     ( bcFunctionCompatibility.M_bcElement ),
    M_bcInternalNode                ( bcFunctionCompatibility.M_bcInternalNode ),
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
OneDFSIFunctionSolverDefinedCompatibility::setupNode()
{
    super::setupNode();

    mesh_Type::edge_Type boundaryEdge;
    switch ( M_bcSide )
    {
        case OneDFSI::left:
            M_bcElement       = 0;
            M_bcInternalNode  = M_bcNode + 1;
            break;

        case OneDFSI::right:
            M_bcElement       = M_bcNode - 1;
            M_bcInternalNode  = M_bcNode - 1;
            break;

        default:
            std::cout << "Warning: bcSide \"" << M_bcSide << "\" not available!" << std::endl;
            break;
    }
}

Real
OneDFSIFunctionSolverDefinedCompatibility::computeRHS ( const Real& timeStep )
{
    updateBCVariables();
    computeEigenValuesVectors();

    switch ( M_bcType )
    {
        case OneDFSI::W1:
            return evaluateRHS ( M_eigenvalues[0], M_leftEigenvector1, M_deltaLeftEigenvector1, timeStep );
            break;

        case OneDFSI::W2:
            return evaluateRHS ( M_eigenvalues[1], M_leftEigenvector2, M_deltaLeftEigenvector2, timeStep );
            break;

        default:
            std::cout << "Warning: bcType \"" << M_bcType << "\"not available!" << std::endl;
            return 0.;
    }
}

void
OneDFSIFunctionSolverDefinedCompatibility::computeEigenValuesVectors()
{
    M_fluxPtr->eigenValuesEigenVectors ( M_bcU[0], M_bcU[1],
                                         M_eigenvalues, M_leftEigenvector1, M_leftEigenvector2,
                                         M_bcNode );

    M_fluxPtr->deltaEigenValuesEigenVectors ( M_bcU[0], M_bcU[1],
                                              M_deltaEigenvalues, M_deltaLeftEigenvector1, M_deltaLeftEigenvector2,
                                              M_bcNode );
}

Real
OneDFSIFunctionSolverDefinedCompatibility::evaluateRHS ( const Real& eigenvalue, const container2D_Type& eigenvector,
                                                         const container2D_Type& deltaEigenvector, const Real& timeStep )
{
    Real cfl = computeCFL ( eigenvalue, timeStep );

    container2D_Type U_interpolated;
    Real Qvisco_interpolated;
    U_interpolated[0]   = ( 1 - cfl ) * M_bcU[0]  + cfl * (* (*M_solutionPtr) ["A"]) ( M_bcInternalNode );
    U_interpolated[1]   = ( 1 - cfl ) * M_bcU[1]  + cfl * (* (*M_solutionPtr) ["Q"]) ( M_bcInternalNode );

    // The second condition detects if there is a viscoelastic flow on the bondary.
    if ( M_fluxPtr->physics()->data()->viscoelasticWall() && (* (*M_solutionPtr) ["Q_visc"]) (M_bcNode) > 1e-10 )
    {
        Qvisco_interpolated = ( 1 - cfl ) * (* (*M_solutionPtr) ["Q_visc"]) (M_bcNode)  + cfl * (* (*M_solutionPtr) ["Q_visc"]) ( M_bcInternalNode );
    }
    else
    {
        Qvisco_interpolated = 0;
    }

    container2D_Type U;

    container2D_Type bcNodes;
    bcNodes[0] = M_bcNode;
    bcNodes[1] = M_bcInternalNode;

#ifdef OLD_COMPATIBILITY
    U[0] = U_interpolated[0] - timeStep * M_sourcePtr->interpolatedNonConservativeSource ( U_interpolated[0], U_interpolated[1], 0, bcNodes, cfl );
    U[1] = U_interpolated[1] - timeStep * M_sourcePtr->interpolatedNonConservativeSource ( U_interpolated[0], U_interpolated[1], 1, bcNodes, cfl );
#else
    container2D_Type U0_interpolated;
    U0_interpolated[0] = ( 1 - cfl ) * M_fluxPtr->physics()->data()->area0 (bcNodes[0]) + cfl * M_fluxPtr->physics()->data()->area0 (bcNodes[1]);
    U0_interpolated[1] = 0;

    U[0] = U_interpolated[0]
           - U0_interpolated[0] - timeStep * ( M_sourcePtr->interpolatedNonConservativeSource ( U_interpolated[0], U_interpolated[1], 0, bcNodes, cfl ) -
                                               M_sourcePtr->interpolatedNonConservativeSource ( U0_interpolated[0], U0_interpolated[1], 0, bcNodes, cfl ) );
    U[1] = (U_interpolated[1] - Qvisco_interpolated) // We consider just the elastic component
           - U0_interpolated[1] - timeStep * ( M_sourcePtr->interpolatedNonConservativeSource ( U_interpolated[0], U_interpolated[1], 1, bcNodes, cfl ) -
                                               M_sourcePtr->interpolatedNonConservativeSource ( U0_interpolated[0], U0_interpolated[1], 1, bcNodes, cfl ) );

    U_interpolated[0] -= U0_interpolated[0];
    U_interpolated[1] -= U0_interpolated[1];

    // Adding Z_0
    U[0] += M_fluxPtr->physics()->data()->area0 (bcNodes[0]);
    U[1] += 0;
#endif
    return scalarProduct ( eigenvector, U ) + timeStep * eigenvalue * scalarProduct ( deltaEigenvector, U_interpolated );
}

Real
OneDFSIFunctionSolverDefinedCompatibility::computeCFL ( const Real& eigenvalue, const Real& timeStep ) const
{
    Real cfl = eigenvalue * timeStep / edgeLength (M_fluxPtr->physics()->data()->mesh()->edge ( M_bcElement ) );

#ifdef HAVE_LIFEV_DEBUG
    if ( M_bcInternalNode == 1 ) // the edge is on the left of the domain
    {
        ASSERT ( -1. < cfl && cfl < 0. , "This characteristics is wrong!\nEither it is not outcoming (eigenvalue>0 at the left of the domain),\n or CFL is too high.");

    }
    else                         // the edge is on the right of the domain
    {
        ASSERT ( 0. < cfl && cfl < 1. , "This characteristics is wrong!\nEither it is not outcoming (eigenvalue<0 at the right of the domain),\n or CFL is too high.");
    }
#endif

    return std::abs (cfl);
}



// ===================================================
// Methods
// ===================================================
Real
OneDFSIFunctionSolverDefinedAbsorbing::operator() ( const Real& /*time*/, const Real& timeStep )
{
    updateBCVariables();
    computeEigenValuesVectors();

    Real W_out (0.), W_out_old (0.);
    switch ( M_bcType )
    {
        case OneDFSI::W1:
            W_out = M_bcW[1] + evaluateRHS ( M_eigenvalues[1], M_leftEigenvector2, M_deltaLeftEigenvector2, timeStep ) - scalarProduct ( M_leftEigenvector2, M_bcU );
            W_out_old = -M_bcW[0] + scalarProduct ( M_leftEigenvector1, M_bcU );
            break;

        case OneDFSI::W2:
            W_out = M_bcW[0] + evaluateRHS ( M_eigenvalues[0], M_leftEigenvector1, M_deltaLeftEigenvector1, timeStep ) - scalarProduct ( M_leftEigenvector1, M_bcU );
            W_out_old = -M_bcW[1] + scalarProduct ( M_leftEigenvector2, M_bcU );
            break;
        default:
            std::cout << "Warning: bcType \"" << M_bcType  << "\"not available!" << std::endl;
            break;
    }

    Real a1, a2, a11, a22, b1, b2, c1, c2;
    a1 = M_fluxPtr->physics()->pressure ( M_bcU[0], timeStep, M_bcNode ) - this->venousPressure(); // pressure at previous time step
    a2 = M_bcU[1]; // flux at previous time step

    b1 = M_fluxPtr->physics()->dPdW ( M_bcW[0], M_bcW[1], 0, M_bcNode ); // dP / dW1
    b2 = M_bcU[0] / 2; // dQ / dW1

    c1 = M_fluxPtr->physics()->dPdW ( M_bcW[0], M_bcW[1], 1, M_bcNode ); // dP / dW2
    c2 = b2; // dQ / dW2

    a11 = a1 - b1 * M_bcW[0] - c1 * M_bcW[1];
    a22 = a2 - b2 * M_bcW[0] - c2 * M_bcW[1];

    Real resistance (b1 / b2);
    this->resistance ( resistance );

    return W_out_old + W_out * (b2 * resistance - b1) / (c1 - c2 * resistance) + (a22 * resistance - a11) / (c1 - c2 * resistance);
}



// ===================================================
// Constructors & Destructor
// ===================================================
OneDFSIFunctionSolverDefinedResistance::OneDFSIFunctionSolverDefinedResistance ( const bcSide_Type& bcSide, const bcType_Type& bcType, const Real& resistance ) :
    super                           ( bcSide, bcType ),
    M_resistance                    ( resistance )
{}

OneDFSIFunctionSolverDefinedResistance::OneDFSIFunctionSolverDefinedResistance ( const OneDFSIFunctionSolverDefinedResistance& bcFunctionResistance ) :
    super                           ( bcFunctionResistance ),
    M_resistance                    ( bcFunctionResistance.M_resistance )
{}
// ===================================================
// Constructors & Destructor
// ===================================================
OneDFSIFunctionSolverDefinedWindkessel3::OneDFSIFunctionSolverDefinedWindkessel3 ( const bcSide_Type& bcSide,
        const bcType_Type& bcType,
        const Real&        resistance1,
        const Real&        resistance2,
        const Real&        compliance,
        const bool&        absorbing,
        const Real&        venousPressure ) :
    super                           ( bcSide, bcType ),
    M_resistance1                   ( resistance1 ),
    M_resistance2                   ( resistance2 ),
    M_compliance                    ( compliance ),
    M_absorbing                     ( absorbing ),
    M_venousPressure                ( venousPressure ),
    M_P0                            ( 0. ),
    M_Q_tn                          ( 0. ),
    M_dQdt_tn                       ( 0. ),
    M_integral_tn                   ( 0. )
{}

OneDFSIFunctionSolverDefinedWindkessel3::OneDFSIFunctionSolverDefinedWindkessel3 ( const OneDFSIFunctionSolverDefinedWindkessel3& bcFunctionWindkessel3 ) :
    super                           ( bcFunctionWindkessel3 ),
    M_resistance1                   ( bcFunctionWindkessel3.M_resistance1 ),
    M_resistance2                   ( bcFunctionWindkessel3.M_resistance2 ),
    M_compliance                    ( bcFunctionWindkessel3.M_compliance ),
    M_absorbing                     ( bcFunctionWindkessel3.M_absorbing ),
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
OneDFSIFunctionSolverDefinedWindkessel3::operator() ( const Real& time, const Real& timeStep )
{
    UInt W_outID;
    Real W_out;

    switch ( M_bcType )
    {
        case OneDFSI::W1:
            W_outID = 1;
            W_out = computeRHS ( timeStep );
            break;
        case OneDFSI::W2:
            W_outID = 0;
            W_out = computeRHS ( timeStep );
            break;
        default:
            std::cout << "Warning: bcType \"" << M_bcType  << "\"not available!" << std::endl;
            break;
    }

    Real A ( M_bcU[0] );
    Real Q ( M_bcU[1] );

    if ( M_absorbing )
    {
        Real b1 ( M_fluxPtr->physics()->dPdW ( A, Q, 0, M_bcNode ) ); // dP / dW1 - Missing W_outID ???
        Real b2 ( A / 2 ); // dQ / dW1
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
    Real EXP      ( std::exp ( a * time     ) );
    Real EXPdt    ( std::exp ( a * timeStep ) );
    Real integral ( M_integral_tn + ( timeStep / 2 ) * ( EXP * ( F * EXPdt + Fn ) ) );

    // Compute the solution of circuital ODE
    // P(t) = P(0) + [ \int_0^t ( pv/(R2*C) + (R1+R2)/(R2*C) Q(s)
    //                 + R1 dQ(s)/ds ) exp( s / (R2*C) ) ds ] * exp( - t / (R2*C) )
    Real P        ( M_P0 + integral / EXP );

    // Update variable for the next time step
    M_integral_tn = integral;
    M_dQdt_tn     = dQdt;
    M_Q_tn        = Q;

    // Remove this call to fromPToW it is not compatible with viscoelasticity!
    return  M_fluxPtr->physics()->fromPToW ( P, W_out, W_outID, M_bcNode ); // W_in
}

}
