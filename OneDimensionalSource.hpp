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
 *  @brief File containing a base class for the source function of the 1D hyperbolic problem.
 *
 *  @date 15-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributor Simone Rossi <simone.rossi@epfl.ch>
 *  @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDimensionalSource_H
#define OneDimensionalSource_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Physics.hpp>

namespace LifeV
{

//! OneDimensionalSource - Base class for the source function of the 1D hyperbolic problem.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 *
 *  dU/dt + dF(U)/dz + B(U) = 0
 *  with U=[A,Q]^T
 */
class OneDimensionalSource
{
public:

    //! @name Type definitions and Enumerators
    //@{

    typedef singleton< factory< OneDimensionalSource, OneDimensional::sourceTerm_Type > > factorySource_Type;

    typedef OneDimensionalPhysics                       physics_Type;
    typedef boost::shared_ptr< physics_Type >           physicsPtr_Type;

    typedef OneDimensionalData::container2D_Type        container2D_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit OneDimensionalSource() : M_physics() {}

    explicit OneDimensionalSource( const physicsPtr_Type physics ) : M_physics( physics ) {}

    //! Do nothing destructor
    virtual ~OneDimensionalSource() {}

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
    virtual Real source( const Real& A, const Real& Q, const ID& ii, const UInt& indz = 0 ) const = 0;

    //! Jacobian matrix dBi/dxj
    virtual Real dSdU( const Real& A, const Real& Q, const ID& ii, const ID& jj, const UInt& indz = 0 ) const = 0;

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
                                                const ID& ii, const container2D_Type& bcNodes, const Real& cfl ) const = 0;

    //@}


    //! @name Set Methods
    //@{

    void setPhysics( const physicsPtr_Type& physics ) { M_physics = physics; }

    //@}


    //! @name Get Methods
    //@{

    physicsPtr_Type physics() const { return M_physics; }

    //@}

protected:

    physicsPtr_Type                 M_physics;

private:

    //! @name Unimplemented Methods
    //@{

    OneDimensionalSource& operator=( const physicsPtr_Type physics );

    //@}
};

}

#endif // OneDimensionalSource_H
