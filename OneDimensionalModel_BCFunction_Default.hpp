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

    typedef OneDimensionalModel_BCFunction     BCFunction_Type;
    typedef boost::shared_ptr<BCFunction_Type> BCFunction_PtrType;

    typedef OneDimensionalModel_Flux           Flux_Type;
    typedef boost::shared_ptr< Flux_Type >     Flux_PtrType;

    typedef OneDimensionalModel_Source         Source_Type;
    typedef boost::shared_ptr< Source_Type >   Source_PtrType;

    typedef OneDimensionalModel_Data           Data_Type;
    typedef Data_Type::Mesh_Type               Mesh_Type;

    typedef SolverAmesos                       LinearSolver_Type;
    typedef LinearSolver_Type::vector_type     Vector_Type;

    typedef std::vector<Vector_Type>           Solution_Type;
    typedef boost::shared_ptr< Solution_Type > Solution_PtrType;

    //@}


    //! @name Constructors & Destructor
    //@{

    OneDimensionalModel_BCFunction_Default( const Flux_PtrType flux, const Source_PtrType source,
                                            const Solution_PtrType solution,
                                            const OneD_BCSide& side, const OneD_BC& bcType );

    virtual ~OneDimensionalModel_BCFunction_Default() {}

    //@}


    //! @name Methods
    //@{

    virtual Real operator() ( const Real& /*time*/ ) = 0;

    //@}

protected:

    Flux_PtrType                             M_Flux;
    Source_PtrType                           M_Source;
    Solution_PtrType                         M_Solution;

    UInt                                     M_boundaryDof;
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

    OneDimensionalModel_BCFunction_Riemann( const Flux_PtrType flux, const Source_PtrType source,
                                            const Solution_PtrType solution,
                                            const OneD_BCSide& side, const OneD_BC& bcType );

    virtual ~OneDimensionalModel_BCFunction_Riemann() {}

    //@}


    //! @name Methods
    //@{

    virtual Real operator() ( const Real& /*time*/ );

    //@}

protected:

    //! @name Protected Methods
    //@{

    void update_U_boundary();

    //@}

    //! Value of U at the boundary
    Container2D_Type                         M_U_boundary;

    //! Value of W at the boundary
    Container2D_Type                         M_W_boundary;
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

    OneDimensionalModel_BCFunction_Compatibility( const Flux_PtrType flux, const Source_PtrType source,
                                                  const Solution_PtrType solution,
                                                  const OneD_BCSide& side, const OneD_BC& bcType );

    virtual ~OneDimensionalModel_BCFunction_Compatibility() {}

    //@}


    //! @name Methods
    //@{

    virtual Real operator() ( const Real& time );

    //@}

protected:

    //! @name Protected Methods
    //@{

    Real extrapolate_W( OneD_BC const& W );

    void computeEigenValuesVectors();

    Real extrapolate_L_dot_U( Real const& eigval, Container2D_Type const& eigvec );

    Container2D_Type _interpolLinear( const Real&  deltaT,
                                      const Real&  eigenvalue,
                                      const Container2D_Type& U_bound ) const;

    //@}

    //! Dof of the internal node adjacent to the boundary
    UInt                                               M_internalBoundaryDof;

    //! Boundary point and internal boundary point
    boost::array< Real, NDIM >                         M_boundaryPoint;
    boost::array< Real, NDIM >                         M_internalBdPoint;

    //! Eigen values of the jacobian diffFlux (= dF/dU = H)
    Real                                               M_eigval1;
    Real                                               M_eigval2;

    //! Left eigen vectors for the eigen values eigval1 and eigval2
    Container2D_Type                                   M_left_eigvec1;
    Container2D_Type                                   M_left_eigvec2;
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

    OneDimensionalModel_BCFunction_Absorbing( const Flux_PtrType flux, const Source_PtrType source,
                                              const Solution_PtrType solution,
                                              const OneD_BCSide& side, const OneD_BC& bcType );

    virtual ~OneDimensionalModel_BCFunction_Absorbing() {}

    //@}


    //! @name Methods
    //@{

    Real operator() ( const Real& time );

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

    OneDimensionalModel_BCFunction_Resistance( const Flux_PtrType flux, const Source_PtrType source,
                                               const Solution_PtrType solution,
                                               const OneD_BCSide& side, const OneD_BC& bcType,
                                               const Real& resistance );

    ~OneDimensionalModel_BCFunction_Resistance() {}

    //@}

protected:

   void resistance( Real& resistance );

   Real M_resistance;

};

}

#endif // ONEDIMENSIONALMODEL_BCSOLVERFUNCTIONS_H
