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
 *  @brief File containing a base class for non linear 1D model flux function.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *
 *  @version 2.0
 *  @date 15-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributors Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDFSIFluxNonLinear_H
#define OneDFSIFluxNonLinear_H

#include <lifev/one_d_fsi/solver/OneDFSIFlux.hpp>

namespace LifeV
{

//! OneDFSIFluxNonLinear - Class containing the non-linear flux term \f$\mathbf F\f$ of the 1D hyperbolic problem.
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
 *  \displaystyle\frac{\partial}{\partial A^0}\displaystyle\int_{0}^A
 *  \displaystyle\frac{A}{\rho}\displaystyle\frac{\partial \psi}{\partial A} dA \displaystyle\frac{\partial A^0}{\partial z} +
 *  \displaystyle\frac{\partial}{\partial \beta_0}\displaystyle\int_{0}^A
 *  \displaystyle\frac{A}{\rho}\displaystyle\frac{\partial \psi}{\partial A} dA \displaystyle\frac{\partial \beta_0}{\partial z} +
 *  \displaystyle\frac{\partial}{\partial \beta_1}\displaystyle\int_{0}^A
 *  \displaystyle\frac{A}{\rho}\displaystyle\frac{\partial \psi}{\partial A} dA \displaystyle\frac{\partial \beta_1}{\partial z}
 *  \end{array}\right]
 *  \f]
 *
 *  where
 *
 *  \f[
 *  \mathbf B(\mathbf U) =
 *  \left[\begin{array}{c}
 *  0 \\[2ex]
 *  K_r \displaystyle\frac{Q}{A} + \displaystyle\frac{A}{\rho}\left(\displaystyle\frac{\partial \psi}{\partial A^0}\displaystyle\frac{\partial A^0}{\partial z} +
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
 *  This class implements all the interfaces required for the computation of \f$\mathbf F\f$ and its derivatives.
 */
class OneDFSIFluxNonLinear : public OneDFSIFlux
{

public:

    //! @name Type definitions and Enumerators
    //@{

    typedef OneDFSIFlux           super;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    explicit OneDFSIFluxNonLinear() : super() {};

    //! Constructor
    /*!
     * @param physicsPtr pointer to the physics of the problem
     */
    explicit OneDFSIFluxNonLinear ( const physicsPtr_Type physicsPtr ) : super ( physicsPtr ) {};

    //! Do nothing destructor
    virtual ~OneDFSIFluxNonLinear() {}

    //@}


    //! @name Methods
    //@{

    //! Evaluate the flux term
    /*!
     *  \f[
     *  \mathbf F(\mathbf U) =
     *  \left[\begin{array}{c}
     *  Q \\[2ex]
     *  \alpha \displaystyle \frac{Q^2}{A} + \displaystyle\frac{\beta_0 \beta_1 A^0}{\rho(\beta_1+1)}\left(\displaystyle\frac{A}{A^0}\right)^{\beta_1+1}
     *  \end{array}\right]
     *  \f]
     *
     *  @param A area
     *  @param Q flow rate
     *  @param row row of the flux term
     *  @param iNode node of the mesh
     */
    Real flux ( const Real& A, const Real& Q, const ID& row, const UInt& iNode ) const;


    //! Evaluate the derivative of the flux term
    /*!
     *  \f[
     * \displaystyle\frac{\partial \mathbf F}{\partial \mathbf U} =
     *  \left[\begin{array}{cc}
     * 0 & 1 \\[2ex]
     * \displaystyle\frac{A}{\rho}\displaystyle\frac{\partial \psi}{\partial A} - \alpha \displaystyle\frac{Q^2}{A^2} & 2 \alpha \displaystyle\frac{Q}{A}
     *  \end{array}\right]
     *  \f]
     *
     *  @param A area
     *  @param Q flow rate
     *  @param row row of the derivative of the flux term
     *  @param column column of the derivative of the flux term
     *  @param iNode node of the mesh
     */
    Real dFdU ( const Real& A, const Real& Q, const ID& row, const ID& column, const UInt& iNode ) const;

    //! Eigenvalues and eigenvectors of the Jacobian matrix
    /*!
     *  \f[
     * \lambda_{1,2} = \alpha \displaystyle\frac{Q}{A} \pm \sqrt{\alpha (\alpha - 1)\left(\displaystyle\frac{Q}{A}\right)^2+
     * \displaystyle\frac{A}{\rho}\displaystyle\frac{\partial \psi}{\partial A}},
     *  \f]
     *
     *  \f[
     *  \displaystyle L =
     *  \varsigma
     *  \left[\begin{array}{cc}
     *  -\lambda_2 & 1\\
     *  -\lambda_1 & 1
     *  \end{array}\right]
     *  \f]
     *
     *  @param A area
     *  @param Q flow rate
     *  @param eigenvalues eigenvalues of the Jacobian matrix
     *  @param leftEigenvector1 first row of the left eigenvector matrix
     *  @param leftEigenvector2 second row of the left eigenvector matrix
     *  @param iNode node of the mesh
     */
    void eigenValuesEigenVectors ( const Real& A, const Real& Q,
                                   container2D_Type& eigenvalues,
                                   container2D_Type& leftEigenvector1,
                                   container2D_Type& leftEigenvector2,
                                   const UInt& iNode ) const;

    //! Derivatives of the eigenvalues and eigenvectors of the derivative of the Jacobian matrix
    /*!
     *
     *  \f[
     *  \begin{array}{@{}r@{}c@{}l}
     *  \displaystyle\frac{\partial \lambda_{1,2}}{\partial z} & = & \displaystyle\frac{\partial \lambda_{1,2}}{\partial A^0}\displaystyle\frac{\partial A^0}{\partial z}
     *  +   \displaystyle\frac{\partial \lambda_{1,2}}{\partial \beta_0}\displaystyle\frac{\partial \beta_0}{\partial z}
     *  +   \displaystyle\frac{\partial \lambda_{1,2}}{\partial \beta_1}\displaystyle\frac{\partial \beta_1}{\partial z}
     *  +   \displaystyle\frac{\partial \lambda_{1,2}}{\partial \alpha}\displaystyle\frac{\partial \alpha}{\partial z}\\[4ex]
     *  & = & \displaystyle\frac{Q}{A}\displaystyle\frac{\partial \alpha}{\partial z}
     *  \pm \displaystyle\frac{1}{2}\left(\alpha (\alpha - 1)\left(\displaystyle\frac{Q}{A}\right)^2+
     *  \displaystyle\frac{A}{\rho}\displaystyle\frac{\partial \psi}{\partial A}\right)^{-1/2}\left(\displaystyle\frac{A}{\rho}
     *  \left(\displaystyle\frac{\partial^2 \psi}{\partial A \partial A^0}\displaystyle\frac{\partial A^0}{\partial z}
     *  + \displaystyle\frac{\partial^2 \psi}{\partial A \partial \beta_0}\displaystyle\frac{\partial \beta_0}{\partial z}
     *  + \displaystyle\frac{\partial^2 \psi}{\partial A \partial \beta_1}\displaystyle\frac{\partial \beta_1}{\partial z}\right)
     *  + (2\alpha - 1)\left(\displaystyle\frac{Q}{A}\right)^2 \displaystyle\frac{\partial \alpha}{\partial z}\right).
     *  \end{array}
     *  \f]
     *
     *  \f[
     *  \displaystyle\frac{\partial L}{\partial z} =
     *  \varsigma
     *  \left[\begin{array}{cc}
     *  -\displaystyle\frac{\partial \lambda_2}{\partial z} & 0\\[4ex]
     *  -\displaystyle\frac{\partial \lambda_1}{\partial z} & 0
     *  \end{array}\right]
     *  \f]
     *
     *  @param A area
     *  @param Q flow rate
     *  @param deltaEigenvalues derivative of the eigenvalues of the derivative of the Jacobian matrix
     *  @param deltaLeftEigenvector1 derivative of the first row of the left eigenvector matrix
     *  @param deltaLeftEigenvector2 derivative of the second row of the left eigenvector matrix
     *  @param iNode node of the mesh
     */
    void deltaEigenValuesEigenVectors ( const Real& A, const Real& Q,
                                        container2D_Type& deltaEigenvalues,
                                        container2D_Type& deltaLeftEigenvector1,
                                        container2D_Type& deltaLeftEigenvector2,
                                        const UInt& iNode ) const;

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    explicit OneDFSIFluxNonLinear ( const OneDFSIFluxNonLinear& flux );

    OneDFSIFluxNonLinear& operator= ( const OneDFSIFluxNonLinear& flux );

    //@}
};

//! Factory create function
inline OneDFSIFlux* createOneDFSIFluxNonLinear()
{
    return new OneDFSIFluxNonLinear();
}

}

#endif // OneDFSIFluxNonLinear_H
