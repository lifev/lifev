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
 *  @brief File containing a base class providing physical operations for the 1D model data.
 *
 *  @version 1.0
 *  @date 01-07-2004
 *  @author Vincent Martin
 *
 *  @version 2.0
 *  @date 13-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributor Simone Rossi <simone.rossi@epfl.ch>
 *  @maintainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDFSIPhysics_H
#define OneDFSIPhysics_H

#include <lifev/one_d_fsi/solver/OneDFSIData.hpp>

namespace LifeV
{

//! OneDFSIPhysics - Base class providing physical operations for the 1D model data.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 *  @see Equations and networks of 1-D models \cite FormaggiaLamponi2003
 *  @see Geometrical multiscale coupling of 1-D models \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite BonnemainMalossi2012LVAD
 *
 *  It contains the following methods:
 *  <ol>
 *      <li> utilities for converting Riemann variables to physical quantities (and viceversa);
 *      <li> utilities to compute the different pressure components (and derivatives).
 *  </ol>
 *
 *  \cond
 *  TODO: This class should be splitted in two separate classes
 *  <ol>
 *      <li> one class (no derivation) for the wall physics with a name like OneDFSIWallPhysics
 *      <li> a set of classes for the Riemann conversions like OneDFSIRiemannConverter and derived Linear/NonLinear versions.
 *  </ol>
 *  \endcond
 */
class OneDFSIPhysics
{
public :
    //! @name Type definitions and Enumerators
    //@{

    typedef FactorySingleton< Factory< OneDFSIPhysics, OneDFSI::physicsType_Type > > factoryPhysics_Type;

    typedef OneDFSIData                           data_Type;
    typedef std::shared_ptr< data_Type >        dataPtr_Type;

    typedef VectorEpetra                          vector_Type;
    typedef std::shared_ptr< vector_Type >      vectorPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    explicit OneDFSIPhysics() : M_dataPtr(), M_previousAreaPtr() {}

    //! Constructor
    /*!
     * @param dataPtr pointer to the data container of the problem
     */
    explicit OneDFSIPhysics ( const dataPtr_Type dataPtr ) : M_dataPtr ( dataPtr ), M_previousAreaPtr() {}

    //! Destructor
    virtual ~OneDFSIPhysics() {}

    //@}


    //! @name Conversion methods
    //@{

    //! Compute \f$\mathbf U\f$ from \f$\mathbf W\f$
    /*!
     *  @param U1 first physical variable
     *  @param U2 second physical variable
     *  @param W1 first Riemann variable
     *  @param W2 second Riemann variable
     *  @param iNode node of the mesh
     */
    virtual void fromWToU ( Real& U1, Real& U2, const Real& W1, const Real& W2, const UInt& iNode ) const = 0;

    //! Compute \f$\mathbf W\f$ from \f$\mathbf U\f$
    /*!
     *  @param W1 first Riemann variable
     *  @param W2 second Riemann variable
     *  @param U1 first physical variable
     *  @param U2 second physical variable
     *  @param iNode node of the mesh
     */
    virtual void fromUToW ( Real& W1, Real& W2, const Real& U1, const Real& U2, const UInt& iNode ) const = 0;

    //! Compute \f$P\f$ from \f$\mathbf W\f$
    /*!
     *  @param W1 first Riemann variable
     *  @param W2 second Riemann variable
     *  @param iNode node of the mesh
     *  @return pressure
     */
    virtual Real fromWToP ( const Real& W1, const Real& W2, const UInt& iNode ) const = 0;

    //! Compute \f$W_1\f$ or \f$W_2\f$ from \f$P\f$
    /*!
     *  @param P pressure
     *  @param W Riemann variable
     *  @param iW Riemann variable ID (0 for \f$W_1\f$, 1 or \f$W_2\f$)
     *  @param iNode node of the mesh
     *  @return the other Riemann variable
     */
    virtual Real fromPToW ( const Real& P, const Real& W, const ID& iW, const UInt& iNode ) const = 0;

    //! Compute \f$W_1\f$ or \f$W_2\f$ from \f$Q\f$
    /*!
     *  @param Q pressure
     *  @param W_tn Riemann variable at time \f$t^n\f$
     *  @param W Riemann variable
     *  @param iW Riemann variable ID (0 for \f$W_1\f$, 1 or \f$W_2\f$)
     *  @param iNode node of the mesh
     *  @return the other Riemann variable
     */
    virtual Real fromQToW ( const Real& Q, const Real& W_tn, const Real& W, const ID& iW, const UInt& iNode ) const = 0;

    //! Compute the area \f$A\f$ given the elastic pressure \f$P_\mathrm{elastic}\f$.
    /*!
     *  To be used in initialization, when time derivative of A is supposed equal to zero.
     *
     *  @param P elastic pressure
     *  @param timeStep the time step
     *  @param iNode node of the mesh
     *  @param elasticExternalNodes consider elastic the external nodes (neglect viscoelasticity)
     *  @return \f$ A = A^0 \left( \displaystyle\frac{P_\mathrm{elastic} - P_\mathrm{ext}}{\beta_0} + 1 \right)^{\left(\displaystyle\frac{1}{\beta_1}\right)} \f$
     */
#ifdef HAVE_NEUMANN_VISCOELASTIC_BC
    Real fromPToA ( const Real& P, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes = false ) const;
#else
    Real fromPToA ( const Real& P, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes = true ) const;
#endif

    //@}


    //! @name Derivatives methods
    //@{

    //! Compute the derivative of the area with respect to the time using a first order finite difference
    /*!
     *  @param Anp1 area at time \f$t^{n+1}\f$
     *  @param timeStep the time step
     *  @param iNode node of the mesh
     *  @return \f$\displaystyle\frac{dA(t)}{dt}\f$
     */
    Real dAdt ( const Real& Anp1, const Real& timeStep, const UInt& iNode ) const;

    //! Compute the derivative of pressure with respect to \f$ \mathbf W\f$
    /*!
     *  @param W1 first Riemann variable
     *  @param W2 second Riemann variable
     *  @param iW Riemann derivative ID (0 for \f$\displaystyle\frac{dP}{dW_1}\f$, 1 or \f$\displaystyle\frac{dP}{dW_2}\f$)
     *  @param iNode node of the mesh
     *  @return \f$\displaystyle\frac{dP}{dW_1}\f$ or \f$\displaystyle\frac{dP}{dW_2}\f$
     */
    virtual Real dPdW ( const Real& W1, const Real& W2, const ID& iW, const UInt& iNode ) const = 0;

    //! Compute the derivative of the pressure with respect to \f$A\f$
    /*!
     *  @param A area
     *  @param timeStep the time step
     *  @param iNode node of the mesh
     *  @param elasticExternalNodes consider elastic the external nodes (neglect viscoelasticity)
     *  @return \f$\displaystyle\frac{dP(A)}{dA} = \displaystyle\frac{dP_\mathrm{elastic}(A)}{dA} + \displaystyle\frac{dP_\mathrm{viscoelastic}(A)}{dA}\f$
     */
#ifdef HAVE_NEUMANN_VISCOELASTIC_BC
    Real dPdA ( const Real& A, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes = false ) const;
#else
    Real dPdA ( const Real& A, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes = true ) const;
#endif

    //! Compute the derivative of the elastic pressure with respect to \f$A\f$
    /*!
     *  @param A area
     *  @param iNode node of the mesh
     *  @return \f$\displaystyle\frac{dP_\mathrm{elastic}(A)}{dA} = \displaystyle\frac{\beta_1 \beta_0 ( \displaystyle\frac{A}{A^0} )^{\beta_1}}{A}\f$
     */
    Real dPdAelastic ( const Real& A, const UInt& iNode ) const;

    //! Compute the derivative of the viscoelastic pressure with respect to \f$A\f$
    /*!
     *  @param A area
     *  @param timeStep the time step
     *  @param iNode node of the mesh
     *  @param elasticExternalNodes consider elastic the external nodes (neglect viscoelasticity)
     *  @return \f$\displaystyle\frac{dP_\mathrm{viscoelastic}(A)}{dA} = \displaystyle\frac{\gamma}{A^{3/2}} \left( \displaystyle\frac{1}{\Delta t} -
     *  \displaystyle\frac{dA}{dt} \displaystyle\frac{3}{2A} \right)\f$
     */
#ifdef HAVE_NEUMANN_VISCOELASTIC_BC
    Real dPdAviscoelastic ( const Real& A, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes = false ) const;
#else
    Real dPdAviscoelastic ( const Real& A, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes = true ) const;
#endif

    //! Compute the derivative of the area with respect to \f$P\f$
    /*!
     *  @param P pressure
     *  @param timeStep the time step
     *  @param iNode node of the mesh
     *  @param elasticExternalNodes consider elastic the external nodes (neglect viscoelasticity)
     *  @return \f$\displaystyle\frac{dA(P)}{dP} = \displaystyle\frac{dA(P)}{dP_\mathrm{elastic}} +
     *  \displaystyle\frac{dA(P)}{dP_\mathrm{viscoelastic}}\f$, with \f$\displaystyle\frac{dA(P)}{dP_\mathrm{elastic}} = \displaystyle\frac{A^0}{\beta_0 \beta_1} \left( 1 +
     *  \displaystyle\frac{ P - P_\mathrm{ext} }{ \beta_0 }\right)^{\displaystyle\frac{1}{\beta_1} - 1}\f$
     */
#ifdef HAVE_NEUMANN_VISCOELASTIC_BC
    Real dAdP ( const Real& P, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes = false ) const;
#else
    Real dAdP ( const Real& P, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes = true ) const;
#endif

    //! Compute the derivative of total pressure with respect to \f$A\f$ or \f$Q\f$.
    /*!
     *  @param A area
     *  @param Q flow rate
     *  @param timeStep the time step
     *  @param variable ID (0 for \f$A\f$, 1 or \f$Q\f$)
     *  @param iNode node of the mesh
     *  @return \f$\displaystyle\frac{dP_t}{dA}\f$ or \f$\displaystyle\frac{dP_t}{dQ}\f$
     */
    Real dPTdU ( const Real& A, const Real& Q, const Real& timeStep, const ID& id, const UInt& iNode ) const;

    //@}


    //! @name Methods
    //@{

    //! Compute the reference celerity.
    /*!
     *  @param iNode node of the mesh
     *  @return reference celerity
     */
    Real celerity0 ( const UInt& iNode ) const;

    //! Compute the pressure.
    /*!
     *  Includes the contribution of the external, elastic and viscoelastic pressure.
     *  @param A area
     *  @param timeStep the time step
     *  @param iNode node of the mesh
     *  @param elasticExternalNodes consider elastic the external nodes (neglect viscoelasticity)
     *  @return \f$P = P_\mathrm{elastic} + P_\mathrm{viscoelastic} + P_\mathrm{external}\f$
     */
#ifdef HAVE_NEUMANN_VISCOELASTIC_BC
    Real pressure ( const Real& A, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes = false ) const;
#else
    Real pressure ( const Real& A, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes = true ) const;
#endif

    //! Return the external pressure.
    /*!
     * @return \f$P_\mathrm{external}\f$
     */
    const Real& externalPressure() const
    {
        return M_dataPtr->externalPressure();
    }

    //! Return the venous pressure.
    /*!
     * @return \f$P_\mathrm{venous}\f$
     */
    const Real& venousPressure() const
    {
        return M_dataPtr->venousPressure();
    }

    //! Compute the elastic pressure.
    /*!
     *  @param A area
     *  @param iNode node of the mesh
     *  @return \f$P_\mathrm{elastic} = \beta_0 \left( \left( \displaystyle\frac{A}{A^0} \right)^{\beta_1} - 1 \right)\f$
     */
    Real elasticPressure ( const Real& A, const UInt& iNode ) const;

    //! Compute the viscoelastic pressure.
    /*!
     *  @param A area
     *  @param timeStep the time step
     *  @param iNode node of the mesh
     *  @param elasticExternalNodes consider elastic the external nodes (neglect viscoelasticity)
     *  @return \f$P_\mathrm{viscoelastic} = \gamma \displaystyle\frac{1}{2\sqrt{\pi A}} \displaystyle\frac{dA}{dt}\f$
     */
#ifdef HAVE_NEUMANN_VISCOELASTIC_BC
    Real viscoelasticPressure ( const Real& A, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes = false ) const;
#else
    Real viscoelasticPressure ( const Real& A, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes = true ) const;
#endif

    //! Compute the total pressure
    /*!
     *  @param A area
     *  @param Q flow rate
     *  @param iNode node of the mesh
     *  @return \f$P_t = P + \displaystyle\frac{\rho}{2} \left(\displaystyle\frac{Q}{A}\right)^2\f$
     */
    Real totalPressure ( const Real& A, const Real& Q, const UInt& iNode ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the data container of the problem.
    /*!
     * @param dataPtr pointer to the data container of the problem
     */
    void setData ( const dataPtr_Type& dataPtr )
    {
        M_dataPtr = dataPtr;
    }

    //! Set the area at time \f$t^n\f$.
    /*!
     *  This parameter is required for computing the derivative of the area in time.
     * @param area_tn \f$A^{n}\f$
     */
    void setArea_tn ( const vector_Type& area_tn )
    {
        M_previousAreaPtr.reset ( new vector_Type ( area_tn ) );
    }

    //@}

    //! @name Get Methods
    //@{

    //! Get the data container of the problem.
    /*!
     * @return shared pointer to the data container of the problem
     */
    dataPtr_Type data() const
    {
        return M_dataPtr;
    }

    //@}

protected:

    dataPtr_Type                      M_dataPtr;

private:

    //! @name Unimplemented Methods
    //@{

    explicit OneDFSIPhysics ( const OneDFSIPhysics& physics );

    OneDFSIPhysics& operator= ( const OneDFSIPhysics& physics );

    //@}

    vectorPtr_Type                    M_previousAreaPtr;
};

// ===================================================
// Inline conversion methods
// ===================================================
inline Real
OneDFSIPhysics::fromPToA ( const Real& P, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes ) const
{
    if ( !M_dataPtr->viscoelasticWall() || ( ( iNode == 0 || iNode == M_dataPtr->numberOfNodes() - 1 ) && elasticExternalNodes ) )
    {
        return ( M_dataPtr->area0 ( iNode ) * OneDFSI::pow20 ( ( P - externalPressure() ) / M_dataPtr->beta0 ( iNode ) + 1, 1 / M_dataPtr->beta1 ( iNode ) )  );
    }
    else
    {
        // Newton method to solve the non linear equation
        Real tolerance (1e-6);
        Real maxIT (100);
        UInt i (0);

        Real A ( M_dataPtr->area0 ( iNode ) );
        Real newtonUpdate (0);
        for ( ; i < maxIT ; ++i )
        {
            if ( std::abs ( pressure ( A, timeStep, iNode, elasticExternalNodes ) - P ) < tolerance )
            {
                break;
            }

            newtonUpdate = ( pressure ( A, timeStep, iNode, elasticExternalNodes ) - P ) / dPdA ( A, timeStep, iNode, elasticExternalNodes );
            if ( A - newtonUpdate <= 0 )
            {
                A /= 2.0;    // Bisection
            }
            else
            {
                A -= newtonUpdate;    // Newton
            }
        }
        if ( i == maxIT )
        {
            std::cout << "!!! Warning: conversion fromPToA below tolerance !!! " << std::endl;
            std::cout << "Tolerance: " << tolerance << "; Residual: " << std::abs ( pressure ( A, timeStep, iNode, elasticExternalNodes ) - P ) << std::endl;
        }

        return A;
    }
}

// ===================================================
// Inline derivatives methods
// ===================================================
inline Real
OneDFSIPhysics::dAdt ( const Real& Anp1, const Real& timeStep, const UInt& iNode ) const
{
    return ( Anp1 - (*M_previousAreaPtr) [iNode] ) / timeStep;
}

inline Real
OneDFSIPhysics::dPdA ( const Real& A, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes ) const
{
    return dPdAelastic ( A, iNode ) + dPdAviscoelastic ( A, timeStep, iNode, elasticExternalNodes );
}

inline Real
OneDFSIPhysics::dPdAelastic ( const Real& A, const UInt& iNode ) const
{
    return M_dataPtr->beta0 ( iNode ) * M_dataPtr->beta1 ( iNode ) * OneDFSI::pow05 ( A / M_dataPtr->area0 ( iNode ), M_dataPtr->beta1 ( iNode ) ) / A;
}

inline Real
OneDFSIPhysics::dPdAviscoelastic ( const Real& A, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes ) const
{
    if ( !M_dataPtr->viscoelasticWall() || ( ( iNode == 0 || iNode == M_dataPtr->numberOfNodes() - 1 ) && elasticExternalNodes ) )
    {
        return 0;
    }
    else
    {
        return M_dataPtr->viscoelasticCoefficient ( iNode ) / ( A * std::sqrt ( A ) ) * ( 1 / timeStep - 3 * dAdt ( A, timeStep, iNode ) / ( 2 * A ) );
    }
}

inline Real
OneDFSIPhysics::dAdP ( const Real& P, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes ) const
{
    if ( !M_dataPtr->viscoelasticWall() || ( ( iNode == 0 || iNode == M_dataPtr->numberOfNodes() - 1 ) && elasticExternalNodes ) )
    {
        return M_dataPtr->area0 ( iNode ) / ( M_dataPtr->beta0 ( iNode ) * M_dataPtr->beta1 ( iNode ) )
               * OneDFSI::pow10 ( 1 + ( P - externalPressure() )
                                  / M_dataPtr->beta0 ( iNode ), 1 / M_dataPtr->beta1 ( iNode ) - 1 );
    }
    else
    {
        // Finite difference approach
        return ( fromPToA ( P + M_dataPtr->jacobianPerturbationStress(), timeStep, iNode, elasticExternalNodes ) - fromPToA ( P, timeStep, iNode, elasticExternalNodes ) )
               / M_dataPtr->jacobianPerturbationStress();
    }
}

inline Real
OneDFSIPhysics::dPTdU ( const Real& A, const Real& Q, const Real& timeStep, const ID& id, const UInt& iNode ) const
{
    if ( id == 0 ) // dPt/dA
    {
        return dPdA ( A, timeStep, iNode ) - M_dataPtr->densityRho() * Q * Q / ( A * A * A );
    }

    if ( id == 1 ) // dPt/dQ
    {
        return M_dataPtr->densityRho() * Q / ( A * A );
    }

    ERROR_MSG ("Total pressure's differential function has only 2 components.");
    return -1.;
}

// ===================================================
// Inline methods
// ===================================================
inline Real
OneDFSIPhysics::celerity0 ( const UInt& iNode ) const
{
    return std::sqrt ( M_dataPtr->beta0 ( iNode ) * M_dataPtr->beta1 ( iNode ) / M_dataPtr->densityRho() );
}

inline Real
OneDFSIPhysics::pressure ( const Real& A, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes ) const
{
    return elasticPressure ( A, iNode ) + viscoelasticPressure ( A, timeStep, iNode, elasticExternalNodes ) + externalPressure();
}

inline Real
OneDFSIPhysics::elasticPressure ( const Real& A, const UInt& iNode ) const
{
    return ( M_dataPtr->beta0 ( iNode ) * ( OneDFSI::pow05 ( A / M_dataPtr->area0 ( iNode ), M_dataPtr->beta1 ( iNode ) ) - 1 ) );
}

inline Real
OneDFSIPhysics::viscoelasticPressure ( const Real& A, const Real& timeStep, const UInt& iNode, const bool& elasticExternalNodes ) const
{
    if ( !M_dataPtr->viscoelasticWall() || ( ( iNode == 0 || iNode == M_dataPtr->numberOfNodes() - 1 ) && elasticExternalNodes ) )
    {
        return 0;
    }
    else
    {
        return M_dataPtr->viscoelasticCoefficient ( iNode ) / ( A * std::sqrt ( A ) ) * dAdt ( A, timeStep, iNode );
    }
}

inline Real
OneDFSIPhysics::totalPressure ( const Real& A, const Real& Q, const UInt& iNode ) const
{
    return elasticPressure ( A, iNode ) + M_dataPtr->densityRho() / 2 * Q * Q / ( A * A );
}

}

#endif // OneDFSIPhysics_H
