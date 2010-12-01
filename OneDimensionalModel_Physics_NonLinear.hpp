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
 *  @brief File containing a class providing non linear physical operations for the 1D model data.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @date 01-07-2004
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 13-04-2010
 */

#ifndef ONEDIMENSIONALMODEL_PHYSICS_NONLINEAR_H
#define ONEDIMENSIONALMODEL_PHYSICS_NONLINEAR_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Physics.hpp>

namespace LifeV
{

//! OneDimensionalModel_Physics_NonLinear - Class providing non linear physical operations for the 1D model data.
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
class OneDimensionalModel_Physics_NonLinear : public OneDimensionalModel_Physics
{
public :

    typedef OneDimensionalModel_Physics           super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_Physics_NonLinear();

    OneDimensionalModel_Physics_NonLinear( const Data_PtrType Data );

    //! Destructor
    ~OneDimensionalModel_Physics_NonLinear() {}

    //@}


    //! @name Methods
    //@{

    //! Compute W from U
    /*!
     *  Riemann Invariants corresponding to data (Q, A) at node indz
     *  W1,2 = (Q / A) +- (2 / beta1) * sqrt(chi) * (celerity - celerity0)
     *  being chi the correction coefficient proposed by A. Robertson and H. Zakaria
     */
    void W_from_U( Real& _W1, Real& _W2,
                   const Real& _U1, const Real& _U2,
                   const UInt& indz ) const;

    //! Compute U from W
    /*!
     *  Physical variables corresponding to (W1, W2) at node indz
     *  A = A0 * ( rho / (beta0 * beta1) )^(1/beta1)
     *    * ( beta1 / (4 * sqrt(chi) ) * (W1 - W2) + celerity0 )^(2/beta1)
     *
     *  Q = A (W1 + W2) / 2
     */
    void U_from_W( Real& _U1, Real& _U2,
                   const Real& _W1, const Real& _W2,
                   const UInt& indz ) const;

    //! Compute the pressure as a fuction of W1, W2
    /*!
     *  P = beta0 * ( rho / (beta0 * beta1) * ( beta1 / (4 * sqrt(chi)) * (W1 - W2) + celerity0 )^2 - 1 )
     */
    Real pressure_W( const Real& _W1, const Real& _W2, const UInt& indz = 0 ) const;

    //! Compute the derivative of pressure with respect to W1 and W2
    /*!
     *  Derivative of pressure as a function of (W1, W2)
     *  dP(W1,W2)/dW_1 = rho / (2 * sqrt(chi)) * ( beta1 / (4 * sqrt(chi)) * (W1 - W2) + celerity0 )
     *  dP(W1,W2)/dW_2 = - rho / (2 * sqrt(chi)) * ( beta1 / (4 * sqrt(chi)) * (W1 - W2) + celerity0 )
     */
    Real pressure_WDiff( const Real& _W1, const Real& _W2,
                         const ID& i,
                         const UInt& indz = 0) const;

    //! Compute W1 or W2 given the pressure:
    /*!
     *  W1 - W2 = (4 * sqrt(chi) / beta1) * sqrt( beta0 * beta1 / rho ) ) * ( sqrt( P / beta0 + 1 ) - 1 )
     *  W1 - W2 = 4 * sqrt( beta0 / (beta1 * rho ) ) * ( sqrt( P / beta0 + 1 ) - 1 )
     */
    Real W_from_P( const Real& _P, const Real& _W, const ID& i, const UInt& indz ) const;

    //! Compute W1 or W2 given the flux: fixed point problem
    /*!
     *  ( W1 - W2 + celerity0/K0 )^(2/beta1) * ( W1 + W2 ) = Q / K1
     *
     *  where
     *
     *  K0 = beta1 / ( 4 * sqrt(chi) )
     *  K1 = A0 / 2 * ( rho / (beta0*beta1) )^(1/beta1) * K0^(2/beta1)
     */
    Real W_from_Q( const Real& _Q, const Real& _W_n, const Real& _W, const ID& i, const UInt& indz ) const;

    //@}

};

//! Factory create function
inline OneDimensionalModel_Physics* Create_OneDimensionalModel_Physics_NonLinear()
{
    return new OneDimensionalModel_Physics_NonLinear();
}

}

#endif // ONEDIMENSIONALMODEL_PHYSICS_NONLINEAR_H
