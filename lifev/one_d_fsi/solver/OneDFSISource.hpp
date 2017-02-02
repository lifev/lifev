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
 *  @maintainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDFSISource_H
#define OneDFSISource_H

#include <lifev/one_d_fsi/solver/OneDFSIPhysics.hpp>

namespace LifeV
{

//! OneDFSISource - Base class for the source term \f$\mathbf S\f$ of the 1D hyperbolic problem.
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
 *  This class implements all the interfaces required for the computation of \f$\mathbf S\f$ and its derivatives.
 */
class OneDFSISource
{
public:

    //! @name Type definitions and Enumerators
    //@{

    typedef FactorySingleton< Factory< OneDFSISource, OneDFSI::sourceTerm_Type > > factorySource_Type;

    typedef OneDFSIPhysics                              physics_Type;
    typedef std::shared_ptr< physics_Type >           physicsPtr_Type;

    typedef OneDFSIData::container2D_Type               container2D_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    explicit OneDFSISource() : M_physicsPtr() {}

    //! Constructor
    /*!
     * @param physicsPtr pointer to the physics of the problem
     */
    explicit OneDFSISource ( const physicsPtr_Type physicsPtr ) : M_physicsPtr ( physicsPtr ) {}

    //! Do nothing destructor
    virtual ~OneDFSISource() {}

    //@}


    //! @name Virtual methods
    //@{

    //! Evaluate the source term
    /*!
     *  @param A area
     *  @param Q flow rate
     *  @param row row of the source term
     *  @param iNode node of the mesh
     */
    virtual Real source ( const Real& A, const Real& Q, const ID& row, const UInt& iNode ) const = 0;

    //! Evaluate the derivative of the source term
    /*!
     *  @param A area
     *  @param Q flow rate
     *  @param row row of the derivative of the source term
     *  @param column column of the derivative of the source term
     *  @param iNode node of the mesh
     */
    virtual Real dSdU ( const Real& A, const Real& Q, const ID& row, const ID& column, const UInt& iNode ) const = 0;

    //! Evaluate the non-conservative form of the source term at the foot of the outgoing characteristic.
    /*!
     *  This method is used for imposing the compatibility equations at the boundaries.
     *  It interpolates the value between to nodes.
     *
     *  @param A area
     *  @param Q flow rate
     *  @param row row of the source term
     *  @param bcNodes list of boundary nodes
     *  @param cfl cfl used to identify the foot of the characteristic
     */
    virtual Real interpolatedNonConservativeSource ( const Real& A, const Real& Q,
                                                     const ID& row, const container2D_Type& bcNodes, const Real& cfl ) const = 0;

    //@}


    //! @name Set Methods
    //@{

    //! Set the physics of the problem.
    /*!
     * @param physicsPtr pointer to physics of the problem
     */
    void setPhysics ( const physicsPtr_Type& physicsPtr )
    {
        M_physicsPtr = physicsPtr;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Get the physics of the problem.
    /*!
     * @return physics of the problem
     */
    physicsPtr_Type physics() const
    {
        return M_physicsPtr;
    }

    //@}

protected:

    physicsPtr_Type                 M_physicsPtr;

private:

    //! @name Unimplemented Methods
    //@{

    explicit OneDFSISource ( const OneDFSISource& source );

    OneDFSISource& operator= ( const OneDFSISource& source );

    //@}
};

}

#endif // OneDFSISource_H
