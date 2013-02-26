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
 *  @brief File containing a base class for linear 1D model flux function.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *
 *  @version 2.0
 *  @date 15-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributor Simone Rossi <simone.rossi@epfl.ch>
 *  @maintainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDFSIFluxLinear_H
#define OneDFSIFluxLinear_H

#include <lifev/one_d_fsi/solver/OneDFSIFlux.hpp>

namespace LifeV
{

//! OneDFSIFluxLinear - Class containing the linear flux term \f$\mathbf F\f$ of the 1D hyperbolic problem.
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
 *  \dots \\[2ex]
 *  \dots
 *  \end{array}\right], \quad
 *  \mathbf S(\mathbf U) =  \mathbf B(\mathbf U) -
 *  \left[\begin{array}{c}
 *  \dots \\[2ex]
 *  \dots
 *  \end{array}\right]
 *  \f]
 *
 *
 *  The assumed wall-law is
 *
 *  \f[
 *  P-P_\mathrm{ext} = \psi(A,A^0,\beta_0, \beta_1, \gamma) = \dots
 *  \f]
 *
 *  This class implements all the interfaces required for the computation of \f$\mathbf F\f$ and its derivatives.
 */
class OneDFSIFluxLinear : public OneDFSIFlux
{

public:

    //! @name Type definitions and Enumerators
    //@{

    typedef OneDFSIFlux           super;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    explicit OneDFSIFluxLinear() : super() {};

    //! Constructor
    /*!
     * @param physicsPtr pointer to the physics of the problem
     */
    explicit OneDFSIFluxLinear ( const physicsPtr_Type physicsPtr ) : super ( physicsPtr ) {};

    //! Do nothing destructor
    virtual ~OneDFSIFluxLinear() {}

    //@}


    //! @name Methods
    //@{

    //! Evaluate the source term
    /*!
     *  \f[
     *  \begin{array}{rcl}
     *  \mathbf F(\mathbf U)_1 & = & F_{11} U_1 + F_{12} U_2,\\
     *  \mathbf F(\mathbf U)_2 & = & F_{21} U_1 + F_{22} U_2
     *  \end{array}
     *  \f]]
     *
     *  @param U1 first unknown
     *  @param U2 second unknown
     *  @param row row of the source term
     *  @param iNode node of the mesh
     */
    Real flux ( const Real& U1, const Real& U2, const ID& row, const UInt& iNode ) const ;

    //! Evaluate the derivative of the flux term
    /*!
     *  @param A area
     *  @param Q flow rate
     *  @param row row of the derivative of the flux term
     *  @param column column of the derivative of the flux term
     *  @param iNode node of the mesh
     */
    Real dFdU ( const Real& U1, const Real& U2, const ID& row, const ID& column, const UInt& iNode ) const;


    //! Eigenvalues and eigenvectors of the Jacobian matrix
    /*!
     *  @param A area
     *  @param Q flow rate
     *  @param eigenvalues eigenvalues of the Jacobian matrix
     *  @param leftEigenvector1 first row of the left eigenvector matrix
     *  @param leftEigenvector2 second row of the left eigenvector matrix
     *  @param iNode node of the mesh
     */
    void eigenValuesEigenVectors ( const Real& U1, const Real& U2,
                                   container2D_Type& eigenvalues,
                                   container2D_Type& leftEigenvector1,
                                   container2D_Type& leftEigenvector2,
                                   const UInt& iNode ) const;

    //! Derivatives of the eigenvalues and eigenvectors of the derivative of the Jacobian matrix
    /*!
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

    explicit OneDFSIFluxLinear ( const OneDFSIFluxLinear& flux );

    OneDFSIFluxLinear& operator= ( const OneDFSIFluxLinear& flux );

    //@}

};

//! Factory create function
inline OneDFSIFlux* createOneDFSIFluxLinear()
{
    return new OneDFSIFluxLinear();
}

}

#endif // OneDFSIFluxLinear_H
