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
 *  @brief File containing some functions for the boundary conditions of 1D models.
 *
 *  @version 1.0
 *  @date 01-08-2006
 *  @author Lucia Mirabella  <lucia.mirabella@gmail.com>
 *
 *  @version 2.0
 *  @date 20-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDimensionalBCFunctionDefault_H
#define OneDimensionalBCFunctionDefault_H

// LIFEV - MATHCARD
#include <lifemc/lifefem/OneDimensionalModel_BCFunction.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Data.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Flux.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Source.hpp>

namespace LifeV
{

//! OneDimensionalModel_BCDefaultFunction - Base class for default BC functions
/*!
 *  @author Lucia Mirabella
 */
class OneDimensionalBCFunctionDefault
{
public:

    //! @name Type definitions and Enumerators
    //@{

    typedef OneDimensionalBCFunction                bcFunction_Type;
    typedef boost::shared_ptr<bcFunction_Type>      bcFunctionPtr_Type;

    typedef OneDimensionalFlux                      flux_Type;
    typedef boost::shared_ptr< flux_Type >          fluxPtr_Type;

    typedef OneDimensionalSource                    source_Type;
    typedef boost::shared_ptr< source_Type >        sourcePtr_Type;

    typedef OneDimensionalData                      data_Type;
    typedef data_Type::mesh_Type                    mesh_Type;

    typedef data_Type::container2D_Type             container2D_Type;

    typedef SolverAmesos                            linearSolver_Type;
    typedef linearSolver_Type::vector_type          vector_Type;
    typedef boost::shared_ptr< vector_Type >        vectorPtr_Type;

    typedef std::map< std::string, vectorPtr_Type > solution_Type;
    typedef boost::shared_ptr< solution_Type >      solutionPtr_Type;

    typedef OneDimensional::bcLine_Type             bcLine_Type;
    typedef OneDimensional::bcSide_Type             bcSide_Type;
    typedef OneDimensional::bcType_Type             bcType_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit OneDimensionalBCFunctionDefault( const bcSide_Type& bcSide, const bcType_Type& bcType );

    //! Copy constructor
    /*!
     * @param BCF_Default OneDimensionalBCFunctionDefault
     */
    explicit OneDimensionalBCFunctionDefault( const OneDimensionalBCFunctionDefault& bcFunctionDefault );

    //! Destructor
    virtual ~OneDimensionalBCFunctionDefault() {}

    //@}


    //! @name Methods
    //@{

    virtual Real operator() ( const Real& time, const Real& timeStep );

    //@}


    //! @name Set Methods
    //@{

    void setSolution( const solutionPtr_Type& solution ) { M_solution = solution; }

    void setFluxSource( const fluxPtr_Type& flux, const sourcePtr_Type& source );

    //@}

protected:

    //! @name Protected Methods
    //@{

    virtual void setupNode();

    //@}

    fluxPtr_Type                             M_flux;
    sourcePtr_Type                           M_source;
    solutionPtr_Type                         M_solution;

    UInt                                     M_bcNode;
    bcSide_Type                              M_bcSide;
    bcType_Type                              M_bcType;
};



//! Riemann - Class which implement Riemann BC for One Dimensional BC Functions
/*!
 *  @author Lucia Mirabella
 */
class OneDimensionalBCFunctionRiemann : public OneDimensionalBCFunctionDefault
{
public:

    //! @name Type definitions and Enumerators
    //@{

    typedef OneDimensionalBCFunctionDefault             super;
    typedef super::container2D_Type                     container2D_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit OneDimensionalBCFunctionRiemann( const bcSide_Type& bcSide, const bcType_Type& bcType );

    //! Copy constructor
    /*!
     * @param BCF_Riemann OneDimensionalBCFunctionRiemann
     */
    explicit OneDimensionalBCFunctionRiemann( const OneDimensionalBCFunctionRiemann& bcFunctionRiemann );

    //! Destructor
    virtual ~OneDimensionalBCFunctionRiemann() {}

    //@}


    //! @name Methods
    //@{

    virtual Real operator() ( const Real& time, const Real& timeStep );

    //@}

protected:

    //! @name Protected Methods
    //@{

    void updateBCVariables();

    //@}

    //! Value of U at the boundary
    container2D_Type                         M_bcU;

    //! Value of W at the boundary
    container2D_Type                         M_bcW;
};



//! Compatibility - Class which implement Compatibility BC for One Dimensional BC Functions
/*!
 *  @author Lucia Mirabella
 */
class OneDimensionalBCFunctionCompatibility : public OneDimensionalBCFunctionRiemann
{
public:

    //! @name Type definitions
    //@{

    typedef OneDimensionalBCFunctionRiemann              super;

    typedef super::fluxPtr_Type                          fluxPtr_Type;
    typedef super::sourcePtr_Type                        sourcePtr_Type;
    typedef super::solutionPtr_Type                      solutionPtr_Type;

    typedef super::mesh_Type                             mesh_Type;
    typedef super::container2D_Type                      container2D_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit OneDimensionalBCFunctionCompatibility( const bcSide_Type& bcSide,  const bcType_Type& bcType );

    //! Copy constructor
    /*!
     * @param BCF_Compatibility OneDimensionalBCFunctionCompatibility
     */
    explicit OneDimensionalBCFunctionCompatibility( const OneDimensionalBCFunctionCompatibility& bcFunctionCompatibility );

    //! Destructor
    virtual ~OneDimensionalBCFunctionCompatibility() {}

    //@}


    //! @name Methods
    //@{

    virtual Real operator() ( const Real& /*time*/, const Real& timeStep ) { return computeRHS( timeStep ); }

    //@}

protected:

    //! @name Protected Methods
    //@{

    void setupNode();

    Real computeRHS( const Real& timeStep );

    void computeEigenValuesVectors();

    Real evaluateRHS( const Real& eigenvalue, const container2D_Type& eigenvector, const container2D_Type& deltaEigenvector, const Real& timeStep );

    Real computeCFL( const Real& eigenvalue, const Real& timeStep ) const;

    //! Scalar product between 2D vectors
    Real scalarProduct( const container2D_Type& vector1, const container2D_Type& vector2 ) { return vector1[0]*vector2[0] + vector1[1]*vector2[1]; }

    //@}

    //! Dof of the internal node adjacent to the boundary
    UInt                                               M_bcInternalNode;

    //! Boundary point and internal boundary point
    boost::array< Real, NDIM >                         M_boundaryPoint;
    boost::array< Real, NDIM >                         M_internalBdPoint;

    //! Eigen values of the jacobian diffFlux (= dF/dU = H)
    container2D_Type                                   M_eigenvalues;
    container2D_Type                                   M_deltaEigenvalues;

    //! Left eigen vectors for the two eigen values
    container2D_Type                                   M_leftEigenvector1;
    container2D_Type                                   M_leftEigenvector2;
    container2D_Type                                   M_deltaLeftEigenvector1;
    container2D_Type                                   M_deltaLeftEigenvector2;
};



//! Absorbing - Class which implement Absorbing BC for One Dimensional BC Functions
/*!
 *  @author Lucia Mirabella
 */
class OneDimensionalBCFunctionAbsorbing : public OneDimensionalBCFunctionCompatibility
{
public:

    //! @name Type definitions
    //@{

    typedef OneDimensionalBCFunctionCompatibility        super;

    typedef super::fluxPtr_Type                          fluxPtr_Type;
    typedef super::sourcePtr_Type                        sourcePtr_Type;
    typedef super::solutionPtr_Type                      solutionPtr_Type;

    typedef super::mesh_Type                             mesh_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit OneDimensionalBCFunctionAbsorbing( const bcSide_Type& bcSide, const bcType_Type& bcType ) : super( bcSide, bcType ) {}

    //! Copy constructor
    /*!
     * @param BCF_Absorbing OneDimensionalBCFunctionAbsorbing
     */
    explicit OneDimensionalBCFunctionAbsorbing( const OneDimensionalBCFunctionAbsorbing& bcFunctionAbsorbing ) : super( bcFunctionAbsorbing ) {}

    //! Destructor
    virtual ~OneDimensionalBCFunctionAbsorbing() {}

    //@}


    //! @name Methods
    //@{

    Real operator() ( const Real& time, const Real& timeStep );

    //@}

protected:


   virtual void resistance( Real& /*resistance*/ ) { /*Do nothing => absorbing!*/ }

};



//! Resistance - Class which implement Resistance BC for One Dimensional BC Functions
/*!
 *  @author Lucia Mirabella
 */
class OneDimensionalBCFunctionResistance : public OneDimensionalBCFunctionAbsorbing
{
public:

    //! @name Type definitions
    //@{

    typedef OneDimensionalBCFunctionAbsorbing            super;

    typedef super::fluxPtr_Type                          fluxPtr_Type;
    typedef super::sourcePtr_Type                        sourcePtr_Type;
    typedef super::solutionPtr_Type                      solutionPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit OneDimensionalBCFunctionResistance( const bcSide_Type& bcSide,  const bcType_Type& bcType, const Real& resistance );

    //! Copy constructor
    /*!
     * @param BCF_Resistance OneDimensionalBCFunctionResistance
     */
    explicit OneDimensionalBCFunctionResistance( const OneDimensionalBCFunctionResistance& bcFunctionResistance );

    //! Destructor
    virtual ~OneDimensionalBCFunctionResistance() {}

    //@}

protected:

    void resistance( Real& resistance ) { resistance = M_resistance; }

    Real M_resistance;

};



//! Windkessel (3 elements)
/*!
 *  Q1D -> ---R1------R2---
 *         ^      |       ^
 *        P1D     C       Pv
 *         ^      |       ^
 *
 *  The holding equation is:
 *
 *  P1D + C R2 dP1D/dt = (R1 + R2 ) * Q1D + R1 R2 C dQ1D/dt + Pv
 *
 *  where the "venous" pressure Pv is taken constant and equal to 5mmHg (6666 dyn/cm^2).
 *
 *  You can solve this ODE in analytical fashion and obtain:
 *
 *  P1D(t) = P1D(0) + [ \int_0^t ( pv + (R1+R2) Q1D(s) + R1*R2*C dQ1D(s)/ds ) exp( s / (R2*C) ) ds ] * exp( - t / (R2*C) )
 *
 *  then you just need to exploit numerical integration rules.
 *  Here a simple trapezoidal rule is used, with a first order approximation of derivative
 *
 *  dQ1D(s)/ds = Q1D(t(n+1)) - Q1D(t(n)) * (1/dt)
 *
 *  @author Lucia Mirabella
 */
class OneDimensionalBCFunctionWindkessel3 : public OneDimensionalBCFunctionCompatibility
{
public:

    //! @name Type definitions
    //@{

    typedef OneDimensionalBCFunctionCompatibility        super;

    typedef super::fluxPtr_Type                          fluxPtr_Type;
    typedef super::sourcePtr_Type                        sourcePtr_Type;
    typedef super::solutionPtr_Type                      solutionPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit OneDimensionalBCFunctionWindkessel3( const bcSide_Type& bcSide, const bcType_Type& bcType,
                                                         const Real& resistance1, const Real& resistance2,
                                                         const Real& compliance,
                                                         const bool& absorbing1 = false,
                                                         const Real& venousPressure = 6666. );

    //! Copy constructor
    /*!
     * @param BCF_Resistance OneDimensionalBCFunctionResistance
     */
    explicit OneDimensionalBCFunctionWindkessel3( const OneDimensionalBCFunctionWindkessel3& bcFunctionWindkessel3 );

    //! Destructor
    virtual ~OneDimensionalBCFunctionWindkessel3() {}

    //@}


    //! @name Methods
    //@{

    Real operator() ( const Real& time, const Real& timeStep );

    //@}

protected:

    Real M_resistance1;
    Real M_resistance2;
    Real M_compliance;
    bool M_absorbing1;
    Real M_venousPressure;

    // Initial value of the pressure
    Real M_P0;

    // Q at previous time step
    Real M_Q_tn;

    // dQdt at the previous time step
    Real M_dQdt_tn;

    // Integral of the pressure at the previous time step
    Real M_integral_tn;
};

}

#endif // OneDimensionalBCFunctionDefault_H
