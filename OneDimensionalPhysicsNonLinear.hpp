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
 *  @brief File containing a class providing non linear physical operations for the 1D model data.
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

#ifndef OneDimensionalPhysicsNonLinear_H
#define OneDimensionalPhysicsNonLinear_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalPhysics.hpp>

namespace LifeV
{

//! OneDimensionalPhysicsNonLinear - Class providing non linear physical operations for the 1D model data.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 *
 *  Parameters:
 *  Area0, alpha, beta0, beta1, Kr, rho.
 *
 *  Euler equations
 *  dA/dt + dQ/dz = 0
 *  dQ/dt + d/dz(alpha * Q^2/A) + A/rho * dP/dz + Kr * Q/A = 0
 *
 *  with
 *  P - P_ext = beta0 [ ( A / Area0 )^{beta1} - 1 ]
 */
class OneDimensionalPhysicsNonLinear : public OneDimensionalPhysics
{
public :

    //! @name Type definitions and Enumerators
    //@{

    typedef OneDimensionalPhysics           super;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit OneDimensionalPhysicsNonLinear() : super() {}

    explicit OneDimensionalPhysicsNonLinear( const dataPtr_Type data ) : super( data ) {}

    //! Destructor
    virtual ~OneDimensionalPhysicsNonLinear() {}

    //@}


    //! @name Conversion methods
    //@{

    //! Compute W from U
    /*!
     *  Riemann Invariants corresponding to data (Q, A) at node indz
     *  W1,2 = (Q / A) +- (2 / beta1) * sqrt(chi) * (celerity - celerity0)
     *  being chi the correction coefficient proposed by A. Robertson and H. Zakaria
     */
    void fromUToW( Real& W1, Real& W2, const Real& U1, const Real& U2, const UInt& indz ) const;

    //! Compute U from W
    /*!
     *  Physical variables corresponding to (W1, W2) at node indz
     *  A = A0 * ( rho / (beta0 * beta1) )^(1/beta1)
     *    * ( beta1 / (4 * sqrt(chi) ) * (W1 - W2) + celerity0 )^(2/beta1)
     *
     *  Q = A (W1 + W2) / 2
     */
    void fromWToU( Real& U1, Real& U2, const Real& W1, const Real& W2, const UInt& indz ) const;

    //! Compute the pressure as a fuction of W1, W2
    /*!
     *  P = beta0 * ( rho / (beta0 * beta1) * ( beta1 / (4 * sqrt(chi)) * (W1 - W2) + celerity0 )^2 - 1 )
     */
    Real fromWToP( const Real& W1, const Real& W2, const UInt& indz = 0 ) const;

    //! Compute W1 or W2 given the pressure:
    /*!
     *  W1 - W2 = (4 * sqrt(chi) / beta1) * sqrt( beta0 * beta1 / rho ) ) * ( sqrt( P / beta0 + 1 ) - 1 )
     *  W1 - W2 = 4 * sqrt( beta0 / (beta1 * rho ) ) * ( sqrt( P / beta0 + 1 ) - 1 )
     */
    Real fromPToW( const Real& P, const Real& W, const ID& i, const UInt& indz ) const;

    //! Compute W1 or W2 given the flux: fixed point problem
    /*!
     *  ( W1 - W2 + celerity0/K0 )^(2/beta1) * ( W1 + W2 ) = Q / K1
     *
     *  where
     *
     *  K0 = beta1 / ( 4 * sqrt(chi) )
     *  K1 = A0 / 2 * ( rho / (beta0*beta1) )^(1/beta1) * K0^(2/beta1)
     */
    Real fromQToW( const Real& Q, const Real& W_n, const Real& W, const ID& i, const UInt& indz ) const;

    //@}


    //! @name Derivatives methods
    //@{

    //! Compute the derivative of pressure with respect to W1 and W2
    /*!
     *  Derivative of pressure as a function of (W1, W2)
     *  dP(W1,W2)/dW_1 = rho / (2 * sqrt(chi)) * ( beta1 / (4 * sqrt(chi)) * (W1 - W2) + celerity0 )
     *  dP(W1,W2)/dW_2 = - rho / (2 * sqrt(chi)) * ( beta1 / (4 * sqrt(chi)) * (W1 - W2) + celerity0 )
     */
    Real dPdW( const Real& W1, const Real& W2, const ID& i, const UInt& indz = 0) const;

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    OneDimensionalPhysicsNonLinear& operator=( const dataPtr_Type data );

    //@}

};

//! Factory create function
inline OneDimensionalPhysics* createOneDimensionalPhysicsNonLinear()
{
    return new OneDimensionalPhysicsNonLinear();
}

}

#endif // OneDimensionalPhysicsNonLinear_H
