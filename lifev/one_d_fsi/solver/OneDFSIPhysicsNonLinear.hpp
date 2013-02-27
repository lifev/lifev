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
 *  @maintainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDFSIPhysicsNonLinear_H
#define OneDFSIPhysicsNonLinear_H

#include <lifev/one_d_fsi/solver/OneDFSIPhysics.hpp>

namespace LifeV
{

//! OneDFSIPhysicsNonLinear - Class providing non linear physical operations for the 1D model data.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 *  @see Equations and networks of 1-D models \cite FormaggiaLamponi2003
 *  @see Geometrical multiscale coupling of 1-D models \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite BonnemainMalossi2012LVAD
 *
 *  It contains the following methods:
 *  <ol>
 *      <li> utilities for converting Riemann variables to physical quantities (and viceversa);
 *      <li> utilities to compute the different pressure components (and derivatives).
 *  </ol>
 */
class OneDFSIPhysicsNonLinear : public OneDFSIPhysics
{
public :

    //! @name Type definitions and Enumerators
    //@{

    typedef OneDFSIPhysics           super;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    explicit OneDFSIPhysicsNonLinear() : super() {}

    //! Constructor
    /*!
     * @param dataPtr pointer to the data container of the problem
     */
    explicit OneDFSIPhysicsNonLinear ( const dataPtr_Type dataPtr ) : super ( dataPtr ) {}

    //! Destructor
    virtual ~OneDFSIPhysicsNonLinear() {}

    //@}


    //! @name Conversion methods
    //@{

    //! Compute \f$\mathbf U\f$ from \f$\mathbf W\f$
    /*!
     *  \cond \TODO improve doxygen description with latex equation and other features \endcond
     *
     *  Riemann Invariants corresponding to data (Q, A) at node iNode
     *  W1,2 = (Q / A) +- (2 / beta1) * sqrt(chi) * (celerity - celerity0)
     *  being chi the correction coefficient proposed by A. Robertson and H. Zakaria
     *
     *  @param U1 first physical variable
     *  @param U2 second physical variable
     *  @param W1 first Riemann variable
     *  @param W2 second Riemann variable
     *  @param iNode node of the mesh
     */
    void fromUToW ( Real& W1, Real& W2, const Real& U1, const Real& U2, const UInt& iNode ) const;

    //! Compute \f$\mathbf W\f$ from \f$\mathbf U\f$
    /*!
     *  \cond \TODO improve doxygen description with latex equation and other features \endcond
     *
     *  Physical variables corresponding to (W1, W2) at node iNode
     *  A = A0 * ( rho / (beta0 * beta1) )^(1/beta1)
     *    * ( beta1 / (4 * sqrt(chi) ) * (W1 - W2) + celerity0 )^(2/beta1)
     *
     *  Q = A (W1 + W2) / 2
     *
     *  @param W1 first Riemann variable
     *  @param W2 second Riemann variable
     *  @param U1 first physical variable
     *  @param U2 second physical variable
     *  @param iNode node of the mesh
     */
    void fromWToU ( Real& U1, Real& U2, const Real& W1, const Real& W2, const UInt& iNode ) const;

    //! Compute \f$P\f$ from \f$\mathbf W\f$
    /*!
     *  \cond \TODO improve doxygen description with latex equation and other features \endcond
     *
     *  @param W1 first Riemann variable
     *  @param W2 second Riemann variable
     *  @param iNode node of the mesh
     *  @return P = beta0 * ( rho / (beta0 * beta1) * ( beta1 / (4 * sqrt(chi)) * (W1 - W2) + celerity0 )^2 - 1 )
     */
    Real fromWToP ( const Real& W1, const Real& W2, const UInt& iNode ) const;

    //! Compute \f$W_1\f$ or \f$W_2\f$ from \f$P\f$
    /*!
     *  \cond \TODO improve doxygen description with latex equation and other features \endcond
     *
     *  W1 - W2 = (4 * sqrt(chi) / beta1) * sqrt( beta0 * beta1 / rho ) ) * ( sqrt( P / beta0 + 1 ) - 1 )
     *  W1 - W2 = 4 * sqrt( beta0 / (beta1 * rho ) ) * ( sqrt( P / beta0 + 1 ) - 1 )
     *
     *  @param P pressure
     *  @param W Riemann variable
     *  @param iW Riemann variable ID (0 for \f$W_1\f$, 1 or \f$W_2\f$)
     *  @param iNode node of the mesh
     *  @return the other Riemann variable
     */
    Real fromPToW ( const Real& P, const Real& W, const ID& iW, const UInt& iNode ) const;

    //! Compute \f$W_1\f$ or \f$W_2\f$ from \f$Q\f$
    /*!
     *  \cond \TODO improve doxygen description with latex equation and other features \endcond
     *
     *  ( W1 - W2 + celerity0/K0 )^(2/beta1) * ( W1 + W2 ) = Q / K1
     *
     *  where
     *
     *  K0 = beta1 / ( 4 * sqrt(chi) )
     *  K1 = A0 / 2 * ( rho / (beta0*beta1) )^(1/beta1) * K0^(2/beta1)
     *
     *  @param Q pressure
     *  @param W_tn Riemann variable at time \f$t^n\f$
     *  @param W Riemann variable
     *  @param iW Riemann variable ID (0 for \f$W_1\f$, 1 or \f$W_2\f$)
     *  @param iNode node of the mesh
     *  @return the other Riemann variable
     */
    Real fromQToW ( const Real& Q, const Real& W_tn, const Real& W, const ID& iW, const UInt& iNode ) const;

    //@}


    //! @name Derivatives methods
    //@{

    //! Compute the derivative of pressure with respect to \f$ \mathbf W\f$
    /*!
     *  \cond \TODO improve doxygen description with latex equation and other features \endcond
     *
     *  dP(W1,W2)/dW_1 = rho / (2 * sqrt(chi)) * ( beta1 / (4 * sqrt(chi)) * (W1 - W2) + celerity0 )
     *  dP(W1,W2)/dW_2 = - rho / (2 * sqrt(chi)) * ( beta1 / (4 * sqrt(chi)) * (W1 - W2) + celerity0 )
     *
     *  @param W1 first Riemann variable
     *  @param W2 second Riemann variable
     *  @param iW Riemann derivative ID (0 for \f$\displaystyle\frac{dP}{dW_1}\f$, 1 or \f$\displaystyle\frac{dP}{dW_2}\f$)
     *  @param iNode node of the mesh
     *  @return \f$\displaystyle\frac{dP}{dW_1}\f$ or \f$\displaystyle\frac{dP}{dW_2}\f$
     */
    Real dPdW ( const Real& W1, const Real& W2, const ID& iW, const UInt& iNode ) const;

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    explicit OneDFSIPhysicsNonLinear ( const OneDFSIPhysicsNonLinear& physics );

    OneDFSIPhysicsNonLinear& operator= ( const OneDFSIPhysicsNonLinear& physics );

    //@}

};

//! Factory create function
inline OneDFSIPhysics* createOneDFSIPhysicsNonLinear()
{
    return new OneDFSIPhysicsNonLinear();
}

}

#endif // OneDFSIPhysicsNonLinear_H
