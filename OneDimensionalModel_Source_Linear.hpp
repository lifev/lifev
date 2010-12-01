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
 *  @brief File containing a class for the linear source function B of the 1D hyperbolic problem.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @date
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 15-04-2010
 */

#ifndef ONEDIMENSIONALMODEL_SOURCE_LINEAR_H
#define ONEDIMENSIONALMODEL_SOURCE_LINEAR_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Source.hpp>

namespace LifeV
{

//! OneDimensionalModel_Source_Linear - Class for the linear source function S of the 1D hyperbolic problem.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 *
 *  dU/dt + dF(U)/dz + S(U) = 0
 *  with U=[U1,U2]^T
 */
class OneDimensionalModel_Source_Linear : public OneDimensionalModel_Source
{

public:

    typedef OneDimensionalModel_Source         super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_Source_Linear();

    OneDimensionalModel_Source_Linear( const Physics_PtrType Physics );

    //! Do nothing destructor
    ~OneDimensionalModel_Source_Linear() {}

    //@}


    //! @name Methods
    //@{

    //! S = [S1, S2]^T
    /*!
     *
     * S1 = S10 + S11 U1 + S12 U2
     * S2 = S20 + S21 U1 + S22 U2
     *
     * \param indz : is the index position for the parameter
     */
    Real operator()( const Real& _U1, const Real& _U2,
                     const ID& ii,
                     const UInt& indz = 0 ) const ;

    //! Jacobian matrix dSi/dxj
    Real diff( const Real& _U1, const Real& _U2,
               const ID& ii,    const ID& jj,
               const UInt& indz = 0 ) const;

    //! Second derivative tensor d2Si/(dxj dxk)
//    Real diff2( const Real& _U1, const Real& _U2,
//                const ID& ii,    const ID& jj, const ID& kk,
//                const UInt& indz = 0 ) const;

    //! Sql = [Sql1, Sql2]^T
    /*!
     *  Sql source term of the equation under its quasi-linear formulation:
     *
     *  dU/dt + H(U) dU/dz + Sql(U) = 0
     *
     *  Here H is constant w.r. to U. And Sql = S(U), because there is no variation of
     *  the coefficients.
     *
     *  \param indz : is the index position for the parameter
     */
    Real interpolatedQuasiLinearSource( const Real& _U1, const Real& _U2,
                                        const ID& ii,    const Container2D_Type& bcNodes, const Real& cfl ) const ;

    //@}

};

//! Factory create function
inline OneDimensionalModel_Source* Create_OneDimensionalModel_Source_Linear()
{
    return new OneDimensionalModel_Source_Linear();
}

}

#endif // ONEDIMENSIONALMODEL_SOURCE_LINEAR_H
