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
 *  @brief File containing a class for the non linear source term \f$\mathbf S\f$ of the 1D hyperbolic problem
 *
 *  @version 1.0
 *  @author Vincent Martin
 *
 *  @version 2.0
 *  @date 15-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributor Simone Rossi <simone.rossi@epfl.ch>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDFSISourceNonLinear_H
#define OneDFSISourceNonLinear_H

#include <lifev/one_d_fsi/solver/OneDFSISource.hpp>

namespace LifeV
{

//! OneDFSISourceNonLinear - Class for the non-linear source function B of the 1D hyperbolic problem.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 *  @see Equations and networks of 1-D models \cite FormaggiaLamponi2003
 *  @see Geometrical multiscale coupling of 1-D models \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite BonnemainMalossi2012LVAD
 *
 *  The conservative form of the generic hyperbolic problem is
 *
 *  \f[
 *  \frac{\partial \mathbf U}{\partial t} + \frac{\partial \mathbf F(\mathbf U)}{\partial z} + \mathbf S(\mathbf U) = 0,
 *  \f]
 *
 *  where \f$\mathbf U\f$ are the conservative variables, \f$\mathbf F\f$ the corresponding fluxes,
 *  and \f$\mathbf S\f$ represents the source terms.
 *
 *  In the present implementation we have:
 *
 *  \f[
 *  \mathbf F(\mathbf U) =
 *  \left[\begin{array}{c}
 *  Q \\[2ex]
 *  \alpha \displaystyle \frac{Q^2}{A} + \displaystyle \displaystyle\int_{0}^A \frac{A}{\rho}\frac{\partial \psi}{\partial A} dA
 *  \end{array}\right], \quad
 *  \mathbf S(\mathbf U) =  \mathbf B(\mathbf U) -
 *  \left[\begin{array}{c}
 *  0 \\[2ex]
 *  \displaystyle\frac{\partial}{\partial A^0}\displaystyle\int_{0}^A \displaystyle\frac{A}{\rho}\displaystyle\frac{\partial \psi}{\partial A} dA
 *  \displaystyle\frac{\partial A^0}{\partial z} + \displaystyle\frac{\partial}{\partial \beta_0}\displaystyle\int_{0}^A
 *  \displaystyle\frac{A}{\rho}\displaystyle\frac{\partial \psi}{\partial A} dA \displaystyle\frac{\partial \beta_0}{\partial z} +
 *  \displaystyle\frac{\partial}{\partial \beta_1}\displaystyle\int_{0}^A \displaystyle\frac{A}{\rho}\displaystyle\frac{\partial \psi}{\partial A} dA
 *  \displaystyle\frac{\partial \beta_1}{\partial z}
 *  \end{array}\right]
 *  \f]
 *
 *  where
 *
 *  \f[
 *  \mathbf B(\mathbf U) =
 *  \left[\begin{array}{c}
 *  0 \\[2ex]
 *  K_r \displaystyle\frac{Q}{A} +
 *  \displaystyle\frac{A}{\rho}\left(\displaystyle\frac{\partial \psi}{\partial A^0}\displaystyle\frac{\partial A^0}{\partial z} +
 *  \displaystyle\frac{\partial \psi}{\partial \beta_0}\displaystyle\frac{\partial \beta_0}{\partial z} +
 *  \displaystyle\frac{\partial \psi}{\partial \beta_1}\displaystyle\frac{\partial \beta_1}{\partial z}\right) +
 *  \displaystyle\frac{Q^2}{A}\displaystyle\frac{\partial \alpha}{\partial z}
 *  \end{array}\right]
 *  \f]
 *
 *  The assumed wall-law is
 *
 *  \f[
 *  P-P_\mathrm{ext} = \psi(A,A^0,\beta_0, \beta_1, \gamma) =
 *  \underbrace{\sqrt{\frac{\pi}{A^0}}\frac{h E}{1-\nu^2}}_{\beta_0} \left(\left(\frac{A}{A^0}\right)^{\beta_1}-1\right) +
 *  \underbrace{\frac{T \tan\phi}{4 \sqrt{\pi}}\frac{h E}{1-\nu^2}}_{\displaystyle\gamma} \frac{1}{A\sqrt{A}} \frac{\partial A}{\partial t}.
 *  \f]
 *
 *  This class implements all the interfaces required for the computation of \f$\mathbf S\f$ and its derivatives.
 */
class OneDFSISourceNonLinear : public OneDFSISource
{
public:

    //! @name Type definitions and Enumerators
    //@{

    typedef OneDFSISource         super;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    explicit OneDFSISourceNonLinear() : super() {}

    //! Constructor
    /*!
     * @param physicsPtr pointer to the physics of the problem
     */
    explicit OneDFSISourceNonLinear( const physicsPtr_Type physicsPtr ) : super( physicsPtr ) {}

    //! Do nothing destructor
    virtual ~OneDFSISourceNonLinear() {}

    //@}


    //! @name Methods
    //@{

    //! Evaluate the source term
    /*!
     *  \f[
     *  \begin{array}{rcl}
     *  \mathbf S(\mathbf U)_1 & = & 0,\\
     *  \mathbf S(\mathbf U)_2 & = &
     *  K_r\displaystyle\frac{Q}{A} -\displaystyle\frac{\beta_0 \beta_1}{\rho(\beta_1+1)}
     *  \left(\displaystyle\frac{A}{A^0}\right)^{\beta_1+1}\displaystyle\frac{\partial A^0}{\partial z}
     *  + \displaystyle\frac{1}{\rho}\left(\displaystyle\frac{A^0}{(\beta_1+1)}
     *  \left(\displaystyle\frac{A}{A^0}\right)^{\beta_1+1}-A\right)\displaystyle\frac{\partial \beta_0}{\partial z}\\[4ex]
     *  &+& \displaystyle\frac{A^0 \beta_0}{\rho(\beta_1+1)}\left(\ln\left(\displaystyle\frac{A}{A^0}\right)-
     *  \displaystyle\frac{1}{(\beta_1+1)}\right)\left(\displaystyle\frac{A}{A^0}\right)^{\beta_1+1}
     *  \displaystyle\frac{\partial \beta_1}{\partial z}+\displaystyle\frac{Q^2}{A}\displaystyle\frac{\partial \alpha}{\partial z},
     *  \end{array}
     *  \f]
     *
     *  @param A area
     *  @param Q flow rate
     *  @param row row of the source term
     *  @param iNode node of the mesh
     */
    Real source( const Real& A, const Real& Q, const ID& row, const UInt& iNode ) const ;

    //! Evaluate the derivative of the source term
    /*!
     *  \f[
     *  \displaystyle\frac{\partial \mathbf S}{\partial A} =
     *  \left[\begin{array}{c}
     *  0 \\[2ex]
     *  -K_r\displaystyle\frac{Q}{A^2}+\displaystyle\frac{1}{\rho}\left(\displaystyle\frac{\partial \psi}{\partial A^0}
     *  \displaystyle\frac{\partial A^0}{\partial z}
     *                    + \displaystyle\frac{\partial \psi}{\partial \beta_0}\displaystyle\frac{\partial \beta_0}{\partial z}
     *                    + \displaystyle\frac{\partial \psi}{\partial \beta_1}\displaystyle\frac{\partial \beta_1}{\partial z} \right)
     *                    -\left(\displaystyle\frac{Q}{A}\right)^2\displaystyle\frac{\partial \alpha}{\partial z}
     *  \end{array}\right],
     *  \quad
     *  \displaystyle\frac{\partial \mathbf S}{\partial Q} =
     *  \left[\begin{array}{c}
     *  0 \\[2ex]
     *  \displaystyle\frac{K_r}{A} + 2 \displaystyle\frac{Q}{A} \displaystyle\frac{\partial \alpha}{\partial z}
     *  \end{array}\right].
     *  \f]
     *
     *  @param A area
     *  @param Q flow rate
     *  @param row row of the derivative of the source term
     *  @param column column of the derivative of the source term
     *  @param iNode node of the mesh
     */
    Real dSdU( const Real& A, const Real& Q, const ID& row, const ID& column, const UInt& iNode ) const;

    //! Evaluate the non-conservative form of the source term at the foot of the outgoing characteristic.
    /*!
     *  This method is used for imposing the compatibility equations at the boundaries. It interpolates the value between to nodes.
     *
     *  @param A area
     *  @param Q flow rate
     *  @param row row of the source term
     *  @param bcNodes list of boundary nodes
     *  @param cfl cfl used to identify the foot of the characteristic
     */
    Real interpolatedNonConservativeSource( const Real& A, const Real& Q,
                                            const ID& row, const container2D_Type& bcNodes, const Real& cfl ) const ;

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    explicit OneDFSISourceNonLinear( const OneDFSISourceNonLinear& source );

    OneDFSISourceNonLinear& operator=( const OneDFSISourceNonLinear& source );

    //@}
};

//! Factory create function
inline OneDFSISource* createOneDFSISourceNonLinear()
{
    return new OneDFSISourceNonLinear();
}

}

#endif // OneDFSISourceNonLinear_H
