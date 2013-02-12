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
 *  @brief File containing a base class for 1D model flux function.
 *
 *  @date 15-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributors Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 *  @maintainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDFSIFlux_H
#define OneDFSIFlux_H

#include <lifev/one_d_fsi/solver/OneDFSIPhysics.hpp>

namespace LifeV
{

//! OneDFSIFlux - Base class for the flux term \f$\mathbf F\f$ of the 1D hyperbolic problem.
/*!
 *  @author Cristiano Malossi
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
 *  This class implements all the interfaces required for the computation of \f$\mathbf F\f$ and its derivatives.
 */
class OneDFSIFlux
{

public:

    //! @name Type definitions and Enumerators
    //@{

    typedef FactorySingleton< Factory< OneDFSIFlux, OneDFSI::fluxTerm_Type > > factoryFlux_Type;

    typedef OneDFSIPhysics                              physics_Type;
    typedef boost::shared_ptr< physics_Type >           physicsPtr_Type;

    typedef OneDFSIData::container2D_Type               container2D_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    explicit OneDFSIFlux() : M_physicsPtr() {}

    //! Constructor
    /*!
     * @param physicsPtr pointer to the physics of the problem
     */
    explicit OneDFSIFlux( const physicsPtr_Type physicsPtr ) : M_physicsPtr( physicsPtr ) {}

    //! Do nothing destructor
    virtual ~OneDFSIFlux() {}

    //@}


    //! @name Virtual Methods
    //@{

    //! Evaluate the flux term
    /*!
     *  @param A area
     *  @param Q flow rate
     *  @param row row of the flux term
     *  @param iNode node of the mesh
     */
    virtual Real flux( const Real& A, const Real& Q, const ID& row, const UInt& iNode ) const = 0;


    //! Evaluate the derivative of the flux term
    /*!
     *  @param A area
     *  @param Q flow rate
     *  @param row row of the derivative of the flux term
     *  @param column column of the derivative of the flux term
     *  @param iNode node of the mesh
     */
    virtual Real dFdU( const Real& A, const Real& Q, const ID& row, const ID& column, const UInt& iNode ) const  = 0;

    //! Eigenvalues and eigenvectors of the Jacobian matrix
    /*!
     *  @param A area
     *  @param Q flow rate
     *  @param eigenvalues eigenvalues of the Jacobian matrix
     *  @param leftEigenvector1 first row of the left eigenvector matrix
     *  @param leftEigenvector2 second row of the left eigenvector matrix
     *  @param iNode node of the mesh
     */
    virtual void eigenValuesEigenVectors( const Real& A, const Real& Q,
                                          container2D_Type& eigenvalues,
                                          container2D_Type& leftEigenvector1,
                                          container2D_Type& leftEigenvector2,
                                          const UInt& iNode ) const = 0;

    //! Derivatives of the eigenvalues and eigenvectors of the derivative of the Jacobian matrix
    /*!
     *  @param A area
     *  @param Q flow rate
     *  @param deltaEigenvalues derivative of the eigenvalues of the derivative of the Jacobian matrix
     *  @param deltaLeftEigenvector1 derivative of the first row of the left eigenvector matrix
     *  @param deltaLeftEigenvector2 derivative of the second row of the left eigenvector matrix
     *  @param iNode node of the mesh
     */
    virtual void deltaEigenValuesEigenVectors( const Real& A, const Real& Q,
                                               container2D_Type& deltaEigenvalues,
                                               container2D_Type& deltaLeftEigenvector1,
                                               container2D_Type& deltaLeftEigenvector2,
                                               const UInt& iNode ) const = 0;

    //@}


    //! @name Set Methods
    //@{

    //! Set the physics of the problem.
    /*!
     * @param physicsPtr pointer to the physics of the problem
     */
    void setPhysics( const physicsPtr_Type& physicsPtr ) { M_physicsPtr = physicsPtr; }

    //@}


    //! @name Get Methods
    //@{

    //! Get the physics of the problem.
    /*!
     * @return physics of the problem
     */
    physicsPtr_Type physics() const { return M_physicsPtr; }

    //@}

protected:

    physicsPtr_Type                 M_physicsPtr;

private:

    //! @name Unimplemented Methods
    //@{

    explicit OneDFSIFlux( const OneDFSIFlux& flux );

    OneDFSIFlux& operator=( const OneDFSIFlux& flux );

    //@}

};

}

#endif // OneDFSIFlux_H
