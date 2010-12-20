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
 *  @mantainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDimensionalFluxNonLinear_H
#define OneDimensionalFluxNonLinear_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Flux.hpp>

namespace LifeV
{

//! OneDimensionalFluxNonLinear - Class containing the non-linear flux function F of the 1D hyperbolic problem.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 *
 *  dU/dt + dF(U)/dz + B(U) = 0
 *  with U=[A,Q]^T
 */
class OneDimensionalFluxNonLinear : public OneDimensionalFlux
{

public:

    //! @name Type definitions and Enumerators
    //@{

    typedef OneDimensionalFlux           super;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit OneDimensionalFluxNonLinear() : super() {};

    explicit OneDimensionalFluxNonLinear( const physicsPtr_Type physics ) : super( physics ) {};

    //! Do nothing destructor
    virtual ~OneDimensionalFluxNonLinear() {}

    //@}


    //! @name Methods
    //@{

    //! operator()
    /*!
     *  F = [Q, alpha*Q^2/A + beta0*beta1/(rho*(beta1+1)*A0^beta1) * A^(beta1+1) ]
     *  \param indz : is the index position for the parameters
     *  when they are space dependent.
     *  This is NOT pretty. I should try to remove this dependency. VM 09/04
     */
    Real flux( const Real& A, const Real& Q, const ID& ii, const UInt& indz = 0 ) const ;


    //! Jacobian matrix
    /*!
     *  Hij = dFi/dxj
     *
     *  diff(1,1) = dF1/dx1    diff(1,2) = dF1/dx2
     *  diff(2,1) = dF2/dx1    diff(2,2) = dF2/dx2
     */
    Real dFdU( const Real& A, const Real& Q, const ID& ii, const ID& jj, const UInt& indz = 0 ) const;

    //! Second derivative tensor d2Fi/(dxj dxk)
    /*!
     *  diff2(1,1,1) = d2F1/dx1dx1    diff2(1,1,2) = d2F1/dx1dx2
     *  diff2(1,2,1) = d2F1/dx2dx1    diff2(1,2,2) = d2F1/dx2dx2
     *  diff2(2,1,1) = d2F2/dx1dx1    diff2(2,1,2) = d2F2/dx1dx2
     *  diff2(2,2,1) = d2F2/dx2dx1    diff2(2,2,2) = d2F2/dx2dx2
     *
     *  with d2Fi/dx1dx2 = d2Fi/dx2dx1
     */
//    Real diff2( const Real& A, const Real& Q,
//                const ID& ii,   const ID& jj, const ID& kk,
//                const UInt& indz = 0 ) const;

    //! Eigenvalues and eigenvectors of the Jacobian matrix dFi/dxj
    /*!
     * \param eigi is the ith eigen value of the matrix dF/dx (i=1,2).
     * \param lefteigvecij is the jth component of the left eigen vector associated to eigi. (i,j=1,2)
     */
    void eigenValuesEigenVectors( const Real& A, const Real& Q,
                                  container2D_Type& eigenvalues,
                                  container2D_Type& leftEigenvector1,
                                  container2D_Type& leftEigenvector2,
                                  const UInt& indz = 0 ) const;

    //! Compute the derivative of the eigenvalues and of the eigenvectors of the Jacobian matrix
    void deltaEigenValuesEigenVectors( const Real& A, const Real& Q,
                                       container2D_Type& deltaEigenvalues,
                                       container2D_Type& deltaLeftEigenvector1,
                                       container2D_Type& deltaLeftEigenvector2,
                                       const UInt& indz = 0 ) const;

    //@}
};

//! Factory create function
inline OneDimensionalFlux* createOneDimensionalFluxNonLinear()
{
    return new OneDimensionalFluxNonLinear();
}

}

#endif // OneDimensionalFluxNonLinear_H
