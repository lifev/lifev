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

#ifndef ONEDIMENSIONALMODEL_PHYSICS_H
#define ONEDIMENSIONALMODEL_PHYSICS_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Data.hpp>

namespace LifeV
{

//! OneDimensionalModel_Physics - Base class providing physical operations for the 1D model data.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 */
class OneDimensionalModel_Physics
{
public :

    //! @name Type definitions and Enumerators
    //@{

    typedef singleton< factory< OneDimensionalModel_Physics, OneDimensional::physicsType_Type > > factoryPhysics_Type;

    typedef OneDimensionalModel_Data              data_Type;
    typedef boost::shared_ptr< data_Type >        dataPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit OneDimensionalModel_Physics() : M_data () {}

    explicit OneDimensionalModel_Physics( const dataPtr_Type data ) : M_data ( data ) {}

    //! Destructor
    virtual ~OneDimensionalModel_Physics() {}

    //@}


    //! @name Conversion methods
    //@{

    //! Compute U from W
    virtual void fromWToU( Real& U1, Real& U2, const Real& W1, const Real& W2, const UInt& indz ) const = 0;

    //! Compute W from U
    virtual void fromUToW( Real& W1, Real& W2, const Real& U1, const Real& U2, const UInt& indz ) const = 0;

    //! Compute the pressure as a function of W1, W2:
    virtual Real fromWToP( const Real& W1, const Real& W2, const UInt& indz = 0 ) const = 0;

    //! Compute W1 or W2 given the pressure:
    virtual Real fromPToW( const Real& P, const Real& W, const ID& i, const UInt& indz ) const = 0;

    //! Compute W1 or W2 given the flux
    virtual Real fromQToW( const Real& Q, const Real& W_n, const Real& W, const ID& i, const UInt& indz ) const = 0;

    //! Compute area given the elastic pressure.
    /*!
     *  To be used in initialization, when time derivative of A is supposed null
     *  @return A = A0 * ( (P - Pext) / beta0 + 1 )^(1/beta1)
     */
    Real fromPToA( const Real& P, const UInt& i=0 ) const;

    //@}


    //! @name Derivatives methods
    //@{

    //! Compute the derivative of the area with respect to the time.
    /*!
     * @return dA(t)/dt
     */
    Real dAdt( const Real& Anp1, const Real& An, const Real& Anm1, const Real& timeStep ) const;

    //! Compute the derivative of pressure with respect to W1 and W2
    virtual Real dPdW( const Real& W1, const Real& W2, const ID& i, const UInt& indz = 0 ) const = 0;

    //! Compute the derivative of the elastic pressure with respect to A
    /*!
     * @return dP(A)/dA = beta1 * beta0 * ( A / Area0 )^beta1 / A
     */
    Real dPdA( const Real& A, const UInt& i = 0 ) const;

    //! Compute the derivative of the elastic pressure with respect to A
    /*!
     * @return dA(A)/dP = A0 / ( beta0 * beta1 ) * ( 1 + ( P - Pext )/ beta0 )^(1/beta1 - 1)
     */
    Real dAdP( const Real& P, const UInt& i = 0 ) const;

    //! Compute the derivative of total pressure (P is the elastic pressure) with respect to A and Q.
    /*!
     * @return dPt/dU_ii = dP/dU_ii + rho/2 * d(Q/A)^2/dU_ii
     */
    Real dPTdU( const Real& A, const Real& Q, const ID& id,  const UInt& i = 0 ) const;

    //@}


    //! @name Methods
    //@{

    Real celerity0( const UInt& i ) const;

    //! Compute the elastic pressure.
    /*!
     * It includes the contribution of the external pressure.
     * @return P = beta0 * ( ( A / Area0 )^beta1 - 1 ) + Pext
     */
    Real elasticPressure( const Real& A, const UInt& indz = 0 ) const;

    //! Compute the viscoelastic pressure.
    /*!
     * @return P = gamma * 1/(2*sqrt(pi*A)) * dA / dt
     */
    Real viscoelasticPressure( const Real& Anp1, const Real& An, const Real& Anm1, const Real& timeStep, const UInt& i ) const;

    //! Compute the total pressure (P is the elastic pressure)
    /*!
     * @return Pt = P + rho/2 * (Q/A)^2
     */
    Real totalPressure( const Real& A, const Real& Q, const UInt& i = 0 ) const;

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
    void stiffenVesselLeft( const Real& xl,           const Real& xr,
                            const Real& factor,       const Real& alpha,
                            const Real& delta,        const Real& n,
                            const Real& minDeltaX=1, const UInt& yesAdaptive=0 );

    //! Make the vessel stiffer on the right side of interval [xl, xr]
    /*!
     * \sa stiffenVesselLeft
     */
    void stiffenVesselRight( const Real& xl,           const Real& xr,
                             const Real& factor,       const Real& alpha,
                             const Real& delta,        const Real& n,
                             const Real& minDeltaX=1, const UInt& yesAdaptive=0  );

    //@}


    //! @name Set Methods
    //@{

    void setData( const dataPtr_Type& data ) { M_data = data; }

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

    OneDimensionalModel_Physics& operator=( const dataPtr_Type data );

    //@}
};

// ===================================================
// Inline conversion methods
// ===================================================
inline Real
OneDimensionalModel_Physics::fromPToA( const Real& P, const UInt& i ) const
{
    return ( M_data->area0(i) * std::pow( ( P - M_data->externalPressure() )
                                          / M_data->beta0(i) + 1, 1/M_data->beta1(i) )  );
}

// ===================================================
// Inline derivatives methods
// ===================================================
inline Real
OneDimensionalModel_Physics::dAdt( const Real& Anp1, const Real& An, const Real& Anm1, const Real& timeStep ) const
{
    if ( M_data->dPdtSteps() == 0 )
        return ( Anp1 - An ) / timeStep;
    else
        return ( 3 / 2 * Anp1 - 2 * An + 1 / 2 * Anm1 ) / timeStep;
}

inline Real
OneDimensionalModel_Physics::dPdA( const Real& A, const UInt& i ) const
{
    return M_data->beta0(i) * M_data->beta1(i)
                              * std::pow( A / M_data->area0(i), M_data->beta1(i) ) / A;
}

inline Real
OneDimensionalModel_Physics::dAdP( const Real& P, const UInt& i ) const
{
    return M_data->area0(i) / ( M_data->beta0(i) * M_data->beta1(i) )
                            * std::pow( 1 + ( P - M_data->externalPressure() )
                            / M_data->beta0(i), 1 / M_data->beta1(i) - 1 );
}

inline Real
OneDimensionalModel_Physics::dPTdU( const Real& A, const Real& Q, const ID& id, const UInt& i) const
{
    if ( id == 1 ) // dPt/dA
        return dPdA( A, i ) - M_data->densityRho() * Q * Q / ( A * A * A );

    if ( id == 2 ) // dPt/dQ
        return M_data->densityRho() * Q / ( A * A);

    ERROR_MSG("Total pressure's differential function has only 2 components.");
    return -1.;
}

// ===================================================
// Inline methods
// ===================================================
inline Real
OneDimensionalModel_Physics::celerity0( const UInt& i ) const
{
    return std::sqrt( M_data->beta0(i) * M_data->beta1(i) / M_data->densityRho() );
}

inline Real
OneDimensionalModel_Physics::elasticPressure( const Real& A, const UInt& i ) const
{
    return ( M_data->beta0(i) * ( std::pow( A/M_data->area0(i), M_data->beta1(i) ) - 1 ) ) + M_data->externalPressure();
}

inline Real
OneDimensionalModel_Physics::viscoelasticPressure( const Real& Anp1, const Real& An, const Real& Anm1, const Real& timeStep, const UInt& i ) const
{
    Real area(Anp1);
    if ( M_data->linearizeStringModel() )
        area = M_data->area0(i);

    return M_data->viscoelasticModulus() / ( 2*sqrt( Pi*area ) ) * dAdt(Anp1, An, Anm1, timeStep);
}

inline Real
OneDimensionalModel_Physics::totalPressure( const Real& A, const Real& Q, const UInt& i ) const
{
    return elasticPressure( A, i ) + M_data->densityRho() / 2 * Q * Q / ( A * A );
}

}

#endif //ONEDIMENSIONALMODEL_PHYSICS_H
