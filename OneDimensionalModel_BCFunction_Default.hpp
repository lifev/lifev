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
 *  @brief File containing some functions for the boundary conditions of 1D models.
 *
 *  @version 1.0
 *  @author Lucia Mirabella
 *  @date 01-08-2006
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 20-04-2010
 *
 */

#ifndef ONEDIMENSIONALMODEL_BCSOLVERFUNCTIONS_H
#define ONEDIMENSIONALMODEL_BCSOLVERFUNCTIONS_H

// LIFEV - MATHCARD
#include <lifemc/lifefem/OneDimensionalModel_BCFunction.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Data.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Flux.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Source.hpp>

namespace LifeV {

//! OneDimensionalModel_BCDefaultFunction - Base class for default BC functions
/*!
 *  @author Lucia Mirabella
 */
class OneDimensionalModel_BCFunction_Default
{
public:

    //! @name Type definitions and Enumerators
    //@{

    typedef OneDimensionalModel_BCFunction          BCFunction_Type;
    typedef boost::shared_ptr<BCFunction_Type>      BCFunction_PtrType;

    typedef OneDimensionalModel_Flux                Flux_Type;
    typedef boost::shared_ptr< Flux_Type >          Flux_PtrType;

    typedef OneDimensionalModel_Source              Source_Type;
    typedef boost::shared_ptr< Source_Type >        Source_PtrType;

    typedef OneDimensionalModel_Data                Data_Type;
    typedef Data_Type::Mesh_Type                    Mesh_Type;

    typedef SolverAmesos                            LinearSolver_Type;
    typedef LinearSolver_Type::vector_type          Vector_Type;
    typedef boost::shared_ptr< Vector_Type >        Vector_PtrType;

    typedef std::map< std::string, Vector_PtrType > Solution_Type;
    typedef boost::shared_ptr< Solution_Type >      Solution_PtrType;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_BCFunction_Default( const Flux_PtrType& flux, const Source_PtrType& source,
                                            const OneD_BCSide& side, const OneD_BC& bcType );

    //! Copy constructor
    /*!
     * @param BCF_Default OneDimensionalModel_BCFunction_Default
     */
    OneDimensionalModel_BCFunction_Default( const OneDimensionalModel_BCFunction_Default& BCF_Default );

    //! Destructor
    virtual ~OneDimensionalModel_BCFunction_Default() {}

    //@}


    //! @name Methods
    //@{

    virtual Real operator() ( const Real& /*time*/, const Real& /*timeStep*/ );

    //@}

    //! @name Set Methods
    //@{

    void setSolution( const Solution_PtrType& solution );

    //@}

protected:

    Flux_PtrType                             M_Flux;
    Source_PtrType                           M_Source;
    Solution_PtrType                         M_Solution;

    UInt                                     M_bcNode;
    OneD_BC                                  M_bcType;
};



//! Riemann - Class which implement Riemann BC for One Dimensional BC Functions
/*!
 *  @author Lucia Mirabella
 */
class OneDimensionalModel_BCFunction_Riemann : public OneDimensionalModel_BCFunction_Default
{
public:

    typedef OneDimensionalModel_BCFunction_Default      super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_BCFunction_Riemann( const Flux_PtrType& flux, const Source_PtrType& source,
                                            const OneD_BCSide& side, const OneD_BC& bcType );

    //! Copy constructor
    /*!
     * @param BCF_Riemann OneDimensionalModel_BCFunction_Riemann
     */
    OneDimensionalModel_BCFunction_Riemann( const OneDimensionalModel_BCFunction_Riemann& BCF_Riemann );

    //! Destructor
    virtual ~OneDimensionalModel_BCFunction_Riemann() {}

    //@}


    //! @name Methods
    //@{

    virtual Real operator() ( const Real& /*time*/, const Real& /*timeStep*/ );

    //@}

protected:

    //! @name Protected Methods
    //@{

    void updateBCVariables();

    //@}

    //! Value of U at the boundary
    Container2D_Type                         M_bcU;

    //! Value of W at the boundary
    Container2D_Type                         M_bcW;
};



//! Compatibility - Class which implement Compatibility BC for One Dimensional BC Functions
/*!
 *  @author Lucia Mirabella
 */
class OneDimensionalModel_BCFunction_Compatibility : public OneDimensionalModel_BCFunction_Riemann
{
public:

    //! @name Type definitions
    //@{

    typedef OneDimensionalModel_BCFunction_Riemann       super;

    typedef super::Flux_PtrType                          Flux_PtrType;
    typedef super::Source_PtrType                        Source_PtrType;
    typedef super::Solution_PtrType                      Solution_PtrType;

    typedef super::Mesh_Type                             Mesh_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_BCFunction_Compatibility( const Flux_PtrType& flux, const Source_PtrType& source,
                                                  const OneD_BCSide& side,  const OneD_BC& bcType );

    //! Copy constructor
    /*!
     * @param BCF_Compatibility OneDimensionalModel_BCFunction_Compatibility
     */
    OneDimensionalModel_BCFunction_Compatibility( const OneDimensionalModel_BCFunction_Compatibility& BCF_Compatibility );

    //! Destructor
    virtual ~OneDimensionalModel_BCFunction_Compatibility() {}

    //@}


    //! @name Methods
    //@{

    virtual Real operator() ( const Real& time, const Real& timeStep );

    //@}

protected:

    //! @name Protected Methods
    //@{

    Real extrapolateW( const Real& timeStep );

    void computeEigenValuesVectors();

    Real computeDeltaLU( const Real& eigenvalue, const Container2D_Type& eigenvector, const Real& timeStep );

    Container2D_Type interpolateU( const Real& eigenvalue, const Real& timeStep ) const;
    Container2D_Type interpolateW( const Real& eigenvalue, const Real& timeStep ) const;

    Real computeCFL( const Real& eigenvalue, const Real& timeStep ) const;

    //@}

    //! Dof of the internal node adjacent to the boundary
    UInt                                               M_bcInternalNode;

    //! Boundary point and internal boundary point
    boost::array< Real, NDIM >                         M_boundaryPoint;
    boost::array< Real, NDIM >                         M_internalBdPoint;

    //! Eigen values of the jacobian diffFlux (= dF/dU = H)
    Container2D_Type                                   M_eigenvalues;

    //! Left eigen vectors for the two eigen values
    Container2D_Type                                   M_leftEigenvector1;
    Container2D_Type                                   M_leftEigenvector2;
};



//! Absorbing - Class which implement Absorbing BC for One Dimensional BC Functions
/*!
 *  @author Lucia Mirabella
 */
class OneDimensionalModel_BCFunction_Absorbing : public OneDimensionalModel_BCFunction_Compatibility
{
public:

    //! @name Type definitions
    //@{

    typedef OneDimensionalModel_BCFunction_Compatibility super;

    typedef super::Flux_PtrType                          Flux_PtrType;
    typedef super::Source_PtrType                        Source_PtrType;
    typedef super::Solution_PtrType                      Solution_PtrType;

    typedef super::Mesh_Type                             Mesh_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_BCFunction_Absorbing( const Flux_PtrType& flux, const Source_PtrType& source,
                                              const OneD_BCSide& side,  const OneD_BC& bcType );

    //! Copy constructor
    /*!
     * @param BCF_Absorbing OneDimensionalModel_BCFunction_Absorbing
     */
    OneDimensionalModel_BCFunction_Absorbing( const OneDimensionalModel_BCFunction_Absorbing& BCF_Absorbing );

    //! Destructor
    virtual ~OneDimensionalModel_BCFunction_Absorbing() {}

    //@}


    //! @name Methods
    //@{

    Real operator() ( const Real& time, const Real& timeStep );

    //@}

protected:

   virtual void resistance( Real& /*resistance*/ );

};



//! Resistance - Class which implement Resistance BC for One Dimensional BC Functions
/*!
 *  @author Lucia Mirabella
 */
class OneDimensionalModel_BCFunction_Resistance : public OneDimensionalModel_BCFunction_Absorbing
{
public:

    //! @name Type definitions
    //@{

    typedef OneDimensionalModel_BCFunction_Absorbing     super;

    typedef super::Flux_PtrType                          Flux_PtrType;
    typedef super::Source_PtrType                        Source_PtrType;
    typedef super::Solution_PtrType                      Solution_PtrType;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_BCFunction_Resistance( const Flux_PtrType& flux, const Source_PtrType& source,
                                               const OneD_BCSide& side,  const OneD_BC& bcType,
                                               const Real& resistance );

    //! Copy constructor
    /*!
     * @param BCF_Resistance OneDimensionalModel_BCFunction_Resistance
     */
    OneDimensionalModel_BCFunction_Resistance( const OneDimensionalModel_BCFunction_Resistance& BCF_Resistance );

    //! Destructor
    ~OneDimensionalModel_BCFunction_Resistance() {}

    //@}

protected:

   void resistance( Real& resistance );

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
class OneDimensionalModel_BCFunction_Windkessel3 : public OneDimensionalModel_BCFunction_Compatibility
{
public:

    //! @name Type definitions
    //@{

    typedef OneDimensionalModel_BCFunction_Compatibility super;

    typedef super::Flux_PtrType                          Flux_PtrType;
    typedef super::Source_PtrType                        Source_PtrType;
    typedef super::Solution_PtrType                      Solution_PtrType;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_BCFunction_Windkessel3( const Flux_PtrType& flux,        const Source_PtrType& source,
                                                const OneD_BCSide&  side,        const OneD_BC&        bcType,
                                                const Real&         resistance1, const Real&           resistance2,
                                                const Real&         compliance,
                                                const bool&         absorbing1 = false,
                                                const Real&         venousPressure = 6666. );

    //! Copy constructor
    /*!
     * @param BCF_Resistance OneDimensionalModel_BCFunction_Resistance
     */
    OneDimensionalModel_BCFunction_Windkessel3( const OneDimensionalModel_BCFunction_Windkessel3& BCF_Windkessel3 );

    //! Destructor
    ~OneDimensionalModel_BCFunction_Windkessel3() {}

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

#endif // ONEDIMENSIONALMODEL_BCSOLVERFUNCTIONS_H
