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
    @file
    @brief File containing a base class for linear 1D model flux function.

    @version 1.0
    @author Vincent Martin

    @version 2.0
    @date 15-04-2010
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>

    @contributor Simone Rossi <simone.rossi@epfl.ch>

    @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ONEDIMENSIONALMODEL_FLUX_LINEAR_H
#define ONEDIMENSIONALMODEL_FLUX_LINEAR_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Flux.hpp>

namespace LifeV
{

//! OneDimensionalModel_Flux_Linear - Class containing the linear flux function F of the 1D hyperbolic problem.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 *
 *  dU/dt + dF(U)/dz + S(U) = 0
 *  with U=[U1, U2]^T and F(U) = [F11 U1 + F12 U2 ; F21 U1 + F22 U2]
 *
 *  Fij are constant.
 */
class OneDimensionalModel_Flux_Linear : public OneDimensionalModel_Flux
{

public:

    typedef OneDimensionalModel_Flux           super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_Flux_Linear();

    OneDimensionalModel_Flux_Linear( const physicsPtr_Type physics );

    //! Do nothing destructor
    virtual ~OneDimensionalModel_Flux_Linear() {}

    //@}


    //! @name Methods
    //@{

    //! operator()
    /*!
     *  F = F(U) = [F11 U1 + F12 U2 , F21 U1 + F22 U2]^T
     *  \param indz : is the index position for the parameters
     *  when they are space dependent.
     *  This is NOT pretty. I should try to remove this dependency. VM 09/04
     */
    Real operator()( const Real& U1, const Real& U2,
                     const ID& ii,    const UInt& indz = 0) const ;

    //! Jacobian matrix
    /*!
     *  Hij = dFi/dxj
     *
     *  diff(1,1) = dF1/dx1    diff(1,2) = dF1/dx2
     *  diff(2,1) = dF2/dx1    diff(2,2) = dF2/dx2
     */
    Real diff( const Real& U1, const Real& U2,
               const ID& ii,    const ID& jj, const UInt& indz = 0) const;

    //! Second derivative tensor d2Fi/(dxj dxk)
    /*!
     *  diff2(1,1,1) = d2F1/dx1dx1    diff2(1,1,2) = d2F1/dx1dx2
     *  diff2(1,2,1) = d2F1/dx2dx1    diff2(1,2,2) = d2F1/dx2dx2
     *  diff2(2,1,1) = d2F2/dx1dx1    diff2(2,1,2) = d2F2/dx1dx2
     *  diff2(2,2,1) = d2F2/dx2dx1    diff2(2,2,2) = d2F2/dx2dx2
     *
     *  with d2Fi/dx1dx2 = d2Fi/dx2dx1
     */
//    Real diff2( const Real& U1, const Real& U2,
//                const ID& ii,    const ID& jj, const ID& kk,
//                const UInt& indz = 0 ) const;

    //! Eigenvalues and eigenvectors of the Jacobian matrix dFi/dxj
    /*!
     * \param eigi is the ith eigen value of the matrix dF/dx (i=1,2).
     * \param lefteigvecij is the jth component of the left eigen vector associated to eigi. (i,j=1,2)
     */
    void eigenValuesEigenVectors( const Real& U1, const Real& U2,
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
inline OneDimensionalModel_Flux* createOneDimensionalFluxLinear()
{
    return new OneDimensionalModel_Flux_Linear();
}

}

#endif // ONEDIMENSIONALMODEL_FLUX_LINEAR_H
