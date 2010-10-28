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
 *  @brief File containing a base class for the source function of the 1D hyperbolic problem.
 *
 *  @version 1.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 15-04-2010
 */

#ifndef ONEDIMENSIONALMODEL_SOURCE_H
#define ONEDIMENSIONALMODEL_SOURCE_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Physics.hpp>

namespace LifeV {

//! OneDimensionalModel_Source - Base class for the source function of the 1D hyperbolic problem.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 *
 *  dU/dt + dF(U)/dz + B(U) = 0
 *  with U=[A,Q]^T
 */
class OneDimensionalModel_Source
{
public:

    typedef OneDimensionalModel_Physics              Physics_Type;
    typedef boost::shared_ptr< Physics_Type >        Physics_PtrType;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_Source();

    OneDimensionalModel_Source( const Physics_PtrType Physics );

    //! Do nothing destructor
    ~OneDimensionalModel_Source() {}

    //@}


    //! @name Methods
    //@{

    //! B = [0, B2]^T
    /*!
     *  with B2 such that:
     *
     *  B2 = Kr*Q/A
     *     - beta1 * beta0/( rho*(beta1+1) ) * (A/A0)^(beta1+1)  * dA0/dz
     *     + A0/rho * [ 1/(beta1+1) * (A/A0)^(beta1+1) - A/A0 ]  * dbeta0/dz
     *     + A0    * beta0/( rho*(beta1+1) ) * (A/A0)^(beta1+1)
     *     * [ log(A/A0) - 1/(beta1+1) ]                         * dbeta1/dz
     *
     *  \param indz : is the index position for the parameter
     */
    virtual Real operator()( const Real& _A, const Real& _Q,
                             const ID& ii,   const UInt& indz = 0 ) const = 0;

    //! Jacobian matrix dBi/dxj
    virtual Real diff( const Real& _A, const Real& _Q,
                       const ID& ii,   const ID& jj,
                       const UInt& indz = 0 ) const = 0;

    //! Second derivative tensor d2Bi/(dxj dxk)
//    virtual Real diff2( const Real& _A, const Real& _Q,
//                        const ID& ii,   const ID& jj, const ID& kk,
//                        const UInt& indz = 0 ) const = 0;

    //! Sql = [Sql1, Sql2]^T
    /*!
     *  Sql source term of the equation under its quasi-linear formulation:
     *
     *  dU/dt + H(U) dU/dz + Sql(U) = 0
     *
     *  \param indz : is the index position for the parameter
     */
    virtual Real interpolatedQuasiLinearSource( const Real& _U1, const Real& _U2,
                                                const ID& ii,    const Container2D_Type& bcNodes, const Real& cfl ) const = 0;

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

#endif // ONEDIMENSIONALMODEL_SOURCE_H
