//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing a base class providing physical operations for the 1D model data.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @date 01-07-2004
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 13-04-2010
 */

#ifndef ONEDIMENSIONALMODEL_PHYSICS_H
#define ONEDIMENSIONALMODEL_PHYSICS_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Data.hpp>

namespace LifeV {

//! OneDimensionalModel_Physics - Base class providing physical operations for the 1D model data.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 */
class OneDimensionalModel_Physics
{
public :

    typedef OneDimensionalModel_Data              Data_Type;
    typedef boost::shared_ptr< Data_Type >        Data_PtrType;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_Physics();

    OneDimensionalModel_Physics( const Data_PtrType Data );

    //! Destructor
    virtual ~OneDimensionalModel_Physics() {}

    //@}


    //! @name Pure Virtual Methods
    //@{

    //! Compute U from W
    virtual void U_from_W( Real& _U1, Real& _U2, const Real& _W1, const Real& _W2, const UInt& indz ) const = 0;

    //! Compute W from U
    virtual void W_from_U( Real& _W1, Real& _W2, const Real& _U1, const Real& _U2, const UInt& indz ) const = 0;

    //! Compute the pressure as a function of W1, W2:
    virtual Real pressure_W( const Real& _W1, const Real& _W2, const UInt& indz = 0 ) const = 0;

    //! Compute the derivative of pressure with respect to W1 and W2
    virtual Real pressure_WDiff( const Real& _W1, const Real& _W2, const ID& i, const UInt& indz = 0 ) const = 0;

    //! Compute W1 or W2 given the pressure:
    virtual Real W_from_P( const Real& _P, const Real& _W, const ID& i, const UInt& indz ) const = 0;

    //! Compute W1 or W2 given the flux
    virtual Real W_from_Q( const Real& _Q, const Real& _W_n, const Real& _W, const ID& i, const UInt& indz ) const = 0;

    //@}


    //! @name Methods
    //@{

    Real Celerity0( const UInt& i ) const;

    Real Length() const;

    //! Compute the pressure (with viscoelastic term)
    /*!
     * @return P = beta0 * ( ( _A / Area0 )^beta1 - 1 ) + 1/(2*sqrt(pi*A)) * gamma * dA / dt
     */
    ScalVec pressure( const Real& _A,   const Real& _A_n, const Real& _A_nm1,
                      const Real& dt,   const UInt& indz, const UInt& steps = 1,
                      const bool& visco = 1, const bool& linearized = 1 ) const;

    //! compute the pressure : beta0 * ( ( _A / Area0 )^beta1 - 1 )
    Real pressure( const Real& _A, const UInt& indz = 0 ) const;

    //! compute the derivative of the pressure with respect to A
    /*!
     * @return dP(A)/dA = beta1 * beta0 * ( _A / Area0 )^beta1 / A
     */
    Real pressureDiff( const Real& _A, const UInt& indz = 0 ) const;

    //! Compute the total pressure
    /*!
     * @return Pt = P + rho/2 * (Q/A)^2
     */
    Real totalPressure( const Real& _A, const Real& _Q, const UInt& indz = 0 ) const;

    //! Compute the derivative of total pressure with respect to A and Q (used for interface conditions)
    /*!
     * @return dPt/dU_ii = dP/dU_ii + rho/2 * d(Q/A)^2/dU_ii
     */
    Real totalPressureDiff( const Real& _A, const Real& _Q,
                            const ID& i,
                            const UInt& indz = 0 ) const;

    //! Compute area given the pressure: A = A0 * ( P / beta0 + 1 )^(1/beta1)
    /*!
     *  To be used in initialization! when time derivative of A is supposed null
     *  @return A = A0 * ( P / beta0 + 1 )^(1/beta1)
     */
    Real A_from_P( const Real& P, const UInt& indz=0 ) const;

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
    void stiffenVesselLeft( const Real& xl, const Real& xr,
                            const Real& factor, const Real& alpha,
                            const Real& delta, const Real& n,
                            const Real& min_deltax=1, const UInt& yesAdaptive=0 );
    //! Make the vessel stiffer on the right side of interval [xl, xr]
    /*!
     * \sa stiffenVesselLeft
     */
    void stiffenVesselRight( const Real& xl, const Real& xr,
                             const Real& factor, const Real& alpha,
                             const Real& delta, const Real& n,
                             const Real& min_deltax=1, const UInt& yesAdaptive=0  );

    //@}


    //! @name Set Methods
    //@{

    void SetData( const Data_PtrType Data );

    //@}

    //! @name Get Methods
    //@{

    Data_PtrType Data() const ;

    //@}

protected:

    Data_PtrType                      M_Data;

};

}

#endif //ONEDIMENSIONALMODEL_PHYSICS_H
