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
 *  @brief File containing a base class for 1D model flux function.
 *
 *  @version 1.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 15-04-2010
 */

#ifndef ONEDIMENSIONALMODEL_FLUX_H
#define ONEDIMENSIONALMODEL_FLUX_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Physics.hpp>

namespace LifeV {

//! OneDimensionalModel_Flux - Base class for the flux function F of the 1D hyperbolic problem.
/*!
 *  @author Cristiano Malossi
 */
class OneDimensionalModel_Flux
{

public:

    typedef OneDimensionalModel_Physics              Physics_Type;
    typedef boost::shared_ptr< Physics_Type >        Physics_PtrType;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_Flux();

    OneDimensionalModel_Flux( const Physics_PtrType Physics );

    //! Do nothing destructor
    virtual ~OneDimensionalModel_Flux() {}

    //@}


    //! @name Virtual Methods
    //@{

    //! operator()
    /*!
     *  F = [Q, alpha*Q^2/A + beta0*beta1/(rho*(beta1+1)*A0^beta1) * A^(beta1+1) ]
     *  \param indz : is the index position for the parameters
     *  when they are space dependent.
     *  This is NOT pretty. I should try to remove this dependency. VM 09/04
     */
    virtual Real operator()( const Real& A, const Real& Q,
                             const ID& ii,   const UInt& indz = 0 ) const = 0;


    //! Jacobian matrix
    /*!
     *  Hij = dFi/dxj
     *
     *  diff(1,1) = dF1/dx1    diff(1,2) = dF1/dx2
     *  diff(2,1) = dF2/dx1    diff(2,2) = dF2/dx2
     */
    virtual Real diff( const Real& A, const Real& Q,
                       const ID& ii,   const ID& jj, const UInt& indz = 0 ) const  = 0;

    //! Second derivative tensor d2Fi/(dxj dxk)
    /*!
     *  diff2(1,1,1) = d2F1/dx1dx1    diff2(1,1,2) = d2F1/dx1dx2
     *  diff2(1,2,1) = d2F1/dx2dx1    diff2(1,2,2) = d2F1/dx2dx2
     *  diff2(2,1,1) = d2F2/dx1dx1    diff2(2,1,2) = d2F2/dx1dx2
     *  diff2(2,2,1) = d2F2/dx2dx1    diff2(2,2,2) = d2F2/dx2dx2
     *
     *  with d2Fi/dx1dx2 = d2Fi/dx2dx1
     */
//    virtual Real diff2( const Real& A, const Real& Q,
//                        const ID& ii,   const ID& jj, const ID& kk,
//                        const UInt& indz = 0 ) const  = 0;

    //! Eigenvalues and eigenvectors of the Jacobian matrix dFi/dxj
    /*!
     * \param eigi is the ith eigen value of the matrix dF/dx (i=1,2).
     * \param lefteigvecij is the jth component of the left eigen vector associated to eigi. (i,j=1,2)
     */
    virtual void EigenValuesEigenVectors( const Real& A, const Real& Q,
                                                Container2D_Type& eigenvalues,
                                                Container2D_Type& leftEigenvector1,
                                                Container2D_Type& leftEigenvector2,
                                          const UInt& indz = 0 ) const = 0;

    //! Compute the derivative of the eigenvalues and of the eigenvectors of the Jacobian matrix
    virtual void deltaEigenValuesEigenVectors( const Real& A, const Real& Q,
                                                     Container2D_Type& deltaEigenvalues,
                                                     Container2D_Type& deltaLeftEigenvector1,
                                                     Container2D_Type& deltaLeftEigenvector2,
                                               const UInt& indz = 0 ) const = 0;

    //@}


    //! @name Set Methods
    //@{

    void SetPhysics( const Physics_PtrType& Physics );

    //@}


    //! @name Get Methods
    //@{

    Physics_PtrType Physics() const ;

    //@}

protected:

    Physics_PtrType                 M_Physics;
};

}

#endif // ONEDIMENSIONALMODEL_FLUX_H
