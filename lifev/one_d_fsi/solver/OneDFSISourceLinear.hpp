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
 *  @brief File containing a class for the linear source term \f$\mathbf S\f$ of the 1D hyperbolic problem
 *
 *  @version 1.0
 *  @author Vincent Martin
 *
 *  @version 2.0
 *  @date 15-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributor Simone Rossi <simone.rossi@epfl.ch>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDFSISourceLinear_H
#define OneDFSISourceLinear_H

#include <lifev/one_d_fsi/solver/OneDFSISource.hpp>

namespace LifeV
{

//! OneDFSISourceLinear - Class for the linear source function S of the 1D hyperbolic problem.
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
 *  This class implements all the interfaces required for the computation of \f$\mathbf S\f$ and its derivatives.
 */
class OneDFSISourceLinear : public OneDFSISource
{

public:

    //! @name Type definitions and Enumerators
    //@{

    typedef OneDFSISource         super;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    explicit OneDFSISourceLinear() : super() {}

    //! Constructor
    /*!
     * @param physicsPtr pointer to the physics of the problem
     */
    explicit OneDFSISourceLinear ( const physicsPtr_Type physicsPtr ) : super ( physicsPtr ) {}

    //! Do nothing destructor
    virtual ~OneDFSISourceLinear() {}

    //@}


    //! @name Methods
    //@{

    //! Evaluate the source term
    /*!
     *  \f[
     *  \begin{array}{rcl}
     *  \mathbf S(\mathbf U)_1 & = & S_{10} + S_{11} U_1 + S_{12} U_2,\\
     *  \mathbf S(\mathbf U)_2 & = & S_{20} + S_{21} U_1 + S_{22} U_2
     *  \end{array}
     *  \f]
     *
     *  @param U1 first unknown
     *  @param U2 second unknown
     *  @param row row of the source term
     *  @param iNode node of the mesh
     */
    Real source ( const Real& U1, const Real& U2, const ID& row, const UInt& iNode ) const ;

    //! Evaluate the derivative of the source term
    /*!
     *  @param U1 first unknown
     *  @param U2 second unknown
     *  @param row row of the derivative of the source term
     *  @param column column of the derivative of the source term
     *  @param iNode node of the mesh
     */
    Real dSdU ( const Real& U1, const Real& U2, const ID& row, const ID& colum, const UInt& iNode ) const;

    //! Evaluate the non-conservative form of the source term at the foot of the outgoing characteristic.
    /*!
     *  This method is used for imposing the compatibility equations at the boundaries.
     *  It interpolates the value between to nodes.
     *
     *  @param U1 first unknown
     *  @param U2 second unknown
     *  @param row row of the source term
     *  @param bcNodes list of boundary nodes
     *  @param cfl cfl used to identify the foot of the characteristic
     */
    Real interpolatedNonConservativeSource ( const Real& U1, const Real& U2,
                                             const ID& row, const container2D_Type& bcNodes, const Real& cfl ) const ;

    //@}
private:

    //! @name Unimplemented Methods
    //@{

    explicit OneDFSISourceLinear ( const OneDFSISourceLinear& source );

    OneDFSISourceLinear& operator= ( const OneDFSISourceLinear& source );

    //@}

};

//! Factory create function
inline OneDFSISource* createOneDFSISourceLinear()
{
    return new OneDFSISourceLinear();
}

}

#endif // OneDFSISourceLinear_H
