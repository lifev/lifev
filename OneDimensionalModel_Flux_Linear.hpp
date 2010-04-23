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
 *  @brief File containing a class for linear 1D model flux function.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @date
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 15-04-2010
 */

#ifndef ONEDIMENSIONALMODEL_FLUX_LINEAR_H
#define ONEDIMENSIONALMODEL_FLUX_LINEAR_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Flux.hpp>

namespace LifeV {

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

    OneDimensionalModel_Flux_Linear( const Physics_PtrType Physics );

    //! Do nothing destructor
    ~OneDimensionalModel_Flux_Linear() {}

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
    Real operator()( const Real& _U1, const Real& _U2,
                     const ID& ii,    const UInt& indz = 0) const ;

    //! Jacobian matrix
    /*!
     *  Hij = dFi/dxj
     *
     *  diff(1,1) = dF1/dx1    diff(1,2) = dF1/dx2
     *  diff(2,1) = dF2/dx1    diff(2,2) = dF2/dx2
     */
    Real diff( const Real& _U1, const Real& _U2,
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
    Real diff2( const Real& _U1, const Real& _U2,
                const ID& ii,    const ID& jj, const ID& kk,
                const UInt& indz = 0 ) const;

    //! Eigenvalues and eigenvectors of the Jacobian matrix dFi/dxj
    /*!
     * \param eigi is the ith eigen value of the matrix dF/dx (i=1,2).
     * \param lefteigvecij is the jth component of the left eigen vector associated to eigi. (i,j=1,2)
     */
    void jacobian_EigenValues_Vectors( const Real& _U1,    const Real& _U2,
                                             Real& eig1,         Real& eig2,
                                             Real& lefteigvec11, Real& lefteigvec12,
                                             Real& lefteigvec21, Real& lefteigvec22,
                                       const UInt& indz = 0 ) const;

    //@}

};

//! Factory create function
inline OneDimensionalModel_Flux* Create_OneDimensionalModel_Flux_Linear()
{
    return new OneDimensionalModel_Flux_Linear();
}

}

#endif // ONEDIMENSIONALMODEL_FLUX_LINEAR_H
