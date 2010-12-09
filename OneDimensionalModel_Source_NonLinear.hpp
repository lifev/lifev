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
    @brief File containing a class for the non linear source function B of the 1D hyperbolic problem

    @version 1.0
    @author Vincent Martin

    @version 2.0
    @date 15-04-2010
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>

    @contributor Simone Rossi <simone.rossi@epfl.ch>

    @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ONEDIMENSIONALMODEL_SOURCE_NONLINEAR_H
#define ONEDIMENSIONALMODEL_SOURCE_NONLINEAR_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Source.hpp>

namespace LifeV
{

//! OneDimensionalModel_Source_NonLinear - Class for the non-linear source function B of the 1D hyperbolic problem.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 *
 *  dU/dt + dF(U)/dz + B(U) = 0
 *  with U=[A,Q]^T
 */
class OneDimensionalModel_Source_NonLinear : public OneDimensionalModel_Source
{
public:

    typedef OneDimensionalModel_Source         super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_Source_NonLinear();

    OneDimensionalModel_Source_NonLinear( const physicsPtr_Type physics );

    //! Do nothing destructor
    virtual ~OneDimensionalModel_Source_NonLinear() {}

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
    Real operator()(const Real& A, const Real& Q,
                    const ID& ii,
                    const UInt& indz = 0) const ;

    //! Jacobian matrix dBi/dxj
    Real diff(const Real& A, const Real& Q,
              const ID& ii, const ID& jj,
              const UInt& indz = 0) const;

    //! Second derivative tensor d2Bi/(dxj dxk)
//    Real diff2(const Real& _A, const Real& _Q,
//               const ID& ii, const ID& jj, const ID& kk,
//               const UInt& indz = 0) const;

    //! Sql = [Sql1, Sql2]^T
    /*!
     *  Sql source term of the equation under its quasi-linear formulation:
     *
     *  dU/dt + H(U) dU/dz + Sql(U) = 0
     *
     *  \param indz : is the index position for the parameter
     */
    Real interpolatedQuasiLinearSource( const Real& U1, const Real& U2,
                                        const ID& ii,    const container2D_Type& bcNodes, const Real& cfl ) const ;

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    OneDimensionalModel_Source_NonLinear& operator=( const physicsPtr_Type physics );

    //@}
};

//! Factory create function
inline OneDimensionalModel_Source* createOneDimensionalSourceNonLinear()
{
    return new OneDimensionalModel_Source_NonLinear();
}

}

#endif // ONEDIMENSIONALMODEL_SOURCE_NONLINEAR_H
