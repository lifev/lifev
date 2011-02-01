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
 *  @brief File containing a base class providing physical operations for the 1D model data.
 *
 *  @version 1.0
 *  @date 01-07-2004
 *  @author Vincent Martin
 *
 *  @version 2.0
 *  @date 13-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributor Simone Rossi <simone.rossi@epfl.ch>
 *  @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDimensionalPhysics_H
#define OneDimensionalPhysics_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalData.hpp>

namespace LifeV
{

//! OneDimensionalPhysics - Base class providing physical operations for the 1D model data.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 */
class OneDimensionalPhysics
{
public :
    //! @name Type definitions and Enumerators
    //@{

    typedef FactorySingleton< Factory< OneDimensionalPhysics, OneDimensional::physicsType_Type > > factoryPhysics_Type;

    typedef OneDimensionalData                    data_Type;
    typedef boost::shared_ptr< data_Type >        dataPtr_Type;

    typedef VectorEpetra                          vector_Type;
    typedef boost::shared_ptr< vector_Type >      vectorPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit OneDimensionalPhysics() : M_data(), M_area_tn() {}

    explicit OneDimensionalPhysics( const dataPtr_Type data ) : M_data ( data ), M_area_tn() {}

    //! Destructor
    virtual ~OneDimensionalPhysics() {}

    //@}


    //! @name Conversion methods
    //@{

    //! Compute U from W
    virtual void fromWToU( Real& U1, Real& U2, const Real& W1, const Real& W2, const UInt& iNode ) const = 0;

    //! Compute W from U
    virtual void fromUToW( Real& W1, Real& W2, const Real& U1, const Real& U2, const UInt& iNode ) const = 0;

    //! Compute the pressure as a function of W1, W2:
    virtual Real fromWToP( const Real& W1, const Real& W2, const UInt& iNode ) const = 0;

    //! Compute W1 or W2 given the pressure:
    virtual Real fromPToW( const Real& P, const Real& W, const ID& i, const UInt& iNode ) const = 0;

    //! Compute W1 or W2 given the flux
    virtual Real fromQToW( const Real& Q, const Real& W_n, const Real& W, const ID& i, const UInt& iNode ) const = 0;

    //! Compute area given the elastic pressure.
    /*!
     *  To be used in initialization, when time derivative of A is supposed null
     *  @return A = A0 * ( (P - Pext) / beta0 + 1 )^(1/beta1)
     */
    Real fromPToA( const Real& P, const Real& timeStep, const UInt& iNode ) const;

    //@}


    //! @name Derivatives methods
    //@{

    //! Compute the derivative of the area with respect to the time.
    /*!
     * @return dA(t)/dt
     */
    Real dAdt( const Real& Anp1, const Real& timeStep, const UInt& iNode ) const;

    //! Compute the derivative of pressure with respect to W1 and W2
    virtual Real dPdW( const Real& W1, const Real& W2, const ID& i, const UInt& iNode ) const = 0;

    //! Compute the derivative of the pressure with respect to A
    /*!
     * @return dP(A)/dA = dPelastic(A)/dA + dPviscoelastic(A)/dA
     */
    Real dPdA( const Real& A, const Real& timeStep, const UInt& iNode ) const;

    //! Compute the derivative of the elastic pressure with respect to A
    /*!
     * @return dP(A)/dA = beta1 * beta0 * ( A / Area0 )^beta1 / A
     */
    Real dPdAelastic( const Real& A, const UInt& iNode ) const;

    //! Compute the derivative of the viscoelastic pressure with respect to A
    /*!
     * @return dP(A)/dA = gamma / ( A^(3/2) ) * ( 1 / deltaT - 3 * dA/dT / ( 2 * A ) )
     */
    Real dPdAviscoelastic( const Real& A, const Real& timeStep, const UInt& iNode ) const;

    //! Compute the derivative of the elastic pressure with respect to A
    /*!
     * @return dA(A)/dP = A0 / ( beta0 * beta1 ) * ( 1 + ( P - Pext )/ beta0 )^(1/beta1 - 1)
     */
    Real dAdP( const Real& P, const Real& timeStep, const UInt& iNode ) const;

    //! Compute the derivative of total pressure (P is the elastic pressure) with respect to A and Q.
    /*!
     * @return dPt/dU_ii = dP/dU_ii + rho/2 * d(Q/A)^2/dU_ii
     */
    Real dPTdU( const Real& A, const Real& Q, const Real& timeStep, const ID& id, const UInt& iNode ) const;

    //@}


    //! @name Methods
    //@{

    Real celerity0( const UInt& iNode ) const;

    //! Compute the pressure.
    /*!
     * Includes the contribution of the external, elastic and viscoelastic pressure.
     * @return P = beta0 * ( ( A / Area0 )^beta1 - 1 ) + Pext
     */
    Real pressure( const Real& A, const Real& timeStep, const UInt& iNode ) const;

    //! Compute the elastic pressure.
    /*!
     * It includes the contribution of the external pressure.
     * @return P = beta0 * ( ( A / Area0 )^beta1 - 1 ) + Pext
     */
    Real elasticPressure( const Real& A, const UInt& iNode ) const;

    //! Compute the viscoelastic pressure.
    /*!
     * @return P = gamma * 1/(2*sqrt(pi*A)) * dA / dt
     */
    Real viscoelasticPressure( const Real& A, const Real& timeStep, const UInt& iNode ) const;

    //! Compute the total pressure (P is the elastic pressure)
    /*!
     * @return Pt = P + rho/2 * (Q/A)^2
     */
    Real totalPressure( const Real& A, const Real& Q, const UInt& iNode ) const;

    //! Make the vessel stiffer on the left side of interval [xl, xr]
    /*!
     *  These routines change the elastic modulus along the vessel
     *
     *  When x < alpha - delta/2, the Young modulus is E * factor
     *  When x > alpha + delta/2, the Young modulus is E
     *  When alpha - delta/2 < x < alpha + delta/2, the Young modulus changes
     *  smoothly from the larger to the smaller value, according to a
     *  polynomial law of order n.
     *
     *  The grid size can be adapted (yesadaptive=1) in the nieghborhood of alpha,
     *  where the spatial derivative of the parameter will be maximum.
     *  However, the grid size is not allowed to be smaller than min_deltax
     */
    void stiffenVesselLeft( const Real& xl,          const Real& xr,
                            const Real& factor,      const Real& alpha,
                            const Real& delta,       const Real& n,
                            const Real& minDeltaX=1, const UInt& yesAdaptive=0 );

    //! Make the vessel stiffer on the right side of interval [xl, xr]
    /*!
     * \sa stiffenVesselLeft
     */
    void stiffenVesselRight( const Real& xl,          const Real& xr,
                             const Real& factor,      const Real& alpha,
                             const Real& delta,       const Real& n,
                             const Real& minDeltaX=1, const UInt& yesAdaptive=0  );

    //@}


    //! @name Set Methods
    //@{

    void setData( const dataPtr_Type& data ) { M_data = data; }

    void setArea_tn( const vector_Type& area_tn ) { M_area_tn.reset( new vector_Type ( area_tn ) ); }

    //@}

    //! @name Get Methods
    //@{

    dataPtr_Type data() const { return M_data; }

    //@}

protected:

    dataPtr_Type                      M_data;

private:

    //! @name Unimplemented Methods
    //@{

    OneDimensionalPhysics& operator=( const dataPtr_Type data );

    //@}

    vectorPtr_Type                    M_area_tn;
};

// ===================================================
// Inline conversion methods
// ===================================================
inline Real
OneDimensionalPhysics::fromPToA( const Real& P, const Real& timeStep, const UInt& iNode ) const
{
    if ( M_data->viscoelasticWall() )
    {
        // Newton method to solve the non linear equation
        Real tolerance(1e-8);
        Real maxIT(100);
        UInt i(0);

        Real A( M_data->area0( iNode ) );
        for ( ; i < maxIT ; ++i )
        {
            if ( std::abs( pressure( A, timeStep, iNode ) - P ) < tolerance )
                break;
            A -= ( pressure( A, timeStep, iNode ) - P ) / dPdA( A, timeStep, iNode );
        }
        if ( i == maxIT )
            std::cout << "!!! Warning: conversion fromPToA below tolerance !!! " << std::endl;

        return A;
    }
    else
        return ( M_data->area0( iNode ) * OneDimensional::pow20( ( P - M_data->externalPressure() ) / M_data->beta0( iNode ) + 1, 1 / M_data->beta1( iNode ) )  );
}

// ===================================================
// Inline derivatives methods
// ===================================================
inline Real
OneDimensionalPhysics::dAdt( const Real& Anp1, const Real& timeStep, const UInt& iNode ) const
{
    return ( Anp1 - (*M_area_tn)[iNode+1] ) / timeStep;
}

inline Real
OneDimensionalPhysics::dPdA( const Real& A, const Real& timeStep, const UInt& iNode ) const
{
    return dPdAelastic( A, iNode ) + dPdAviscoelastic( A, timeStep, iNode );
}

inline Real
OneDimensionalPhysics::dPdAelastic( const Real& A, const UInt& iNode ) const
{
    return M_data->beta0( iNode ) * M_data->beta1( iNode ) * OneDimensional::pow05( A / M_data->area0( iNode ), M_data->beta1( iNode ) ) / A;
}

inline Real
OneDimensionalPhysics::dPdAviscoelastic( const Real& A, const Real& timeStep, const UInt& iNode ) const
{
    return M_data->viscoelasticCoefficient( iNode ) / ( A * std::sqrt( A ) ) * ( 1 / timeStep - 3 * dAdt( A, timeStep, iNode ) / ( 2 * A ) );
}

inline Real
OneDimensionalPhysics::dAdP( const Real& P, const Real& timeStep, const UInt& iNode ) const
{
    if ( M_data->viscoelasticWall() )
    {
        // Finite difference approach
        return ( fromPToA( P + M_data->jacobianPerturbationPressure(), timeStep, iNode ) - fromPToA( P, timeStep, iNode ) ) / M_data->jacobianPerturbationPressure();
    }
    else
        return M_data->area0( iNode ) / ( M_data->beta0( iNode ) * M_data->beta1( iNode ) )
                                      * OneDimensional::pow10( 1 + ( P - M_data->externalPressure() )
                                      / M_data->beta0( iNode ), 1 / M_data->beta1( iNode ) - 1 );
}

inline Real
OneDimensionalPhysics::dPTdU( const Real& A, const Real& Q, const Real& timeStep, const ID& id, const UInt& iNode ) const
{
    if ( id == 0 ) // dPt/dA
        return dPdA( A, timeStep, iNode ) - M_data->densityRho() * Q * Q / ( A * A * A );

    if ( id == 1 ) // dPt/dQ
        return M_data->densityRho() * Q / ( A * A );

    ERROR_MSG("Total pressure's differential function has only 2 components.");
    return -1.;
}

// ===================================================
// Inline methods
// ===================================================
inline Real
OneDimensionalPhysics::celerity0( const UInt& iNode ) const
{
    return std::sqrt( M_data->beta0( iNode ) * M_data->beta1( iNode ) / M_data->densityRho() );
}

inline Real
OneDimensionalPhysics::pressure( const Real& A, const Real& timeStep, const UInt& iNode ) const
{
    return elasticPressure( A, iNode ) + viscoelasticPressure( A, timeStep, iNode );
}

inline Real
OneDimensionalPhysics::elasticPressure( const Real& A, const UInt& iNode ) const
{
    return ( M_data->beta0( iNode ) * ( OneDimensional::pow05( A/M_data->area0( iNode ), M_data->beta1( iNode ) ) - 1 ) ) + M_data->externalPressure();
}

inline Real
OneDimensionalPhysics::viscoelasticPressure( const Real& A, const Real& timeStep, const UInt& iNode ) const
{
    if ( M_data->viscoelasticWall() )
        return M_data->viscoelasticCoefficient( iNode ) / ( A * std::sqrt( A ) ) * dAdt( A, timeStep, iNode );
    else
        return 0.;
}

inline Real
OneDimensionalPhysics::totalPressure( const Real& A, const Real& Q, const UInt& iNode ) const
{
    return elasticPressure( A, iNode ) + M_data->densityRho() / 2 * Q * Q / ( A * A );
}

}

#endif // OneDimensionalPhysics_H
