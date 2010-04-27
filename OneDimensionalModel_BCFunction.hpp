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
 *  @brief File containing a class for the boundary conditions for 1D tubes.
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

#ifndef ONEDIMENSIONALMODEL_BCFUNCTION_H
#define ONEDIMENSIONALMODEL_BCFUNCTION_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Definitions.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Data.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Flux.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Source.hpp>

namespace LifeV {

//! OneDimensionalModel_BCFunction - Base class for One Dimensional BC Functions.
/*!
 *  @author Lucia Mirabella
 */
class OneDimensionalModel_BCFunction
{
public:

    //! @name Type definitions and Enumerators
    //@{

    //! Constructor

    typedef OneDimensionalModel_Flux          Flux_Type;
    typedef boost::shared_ptr< Flux_Type >    Flux_PtrType;

    typedef OneDimensionalModel_Source        Source_Type;
    typedef boost::shared_ptr< Source_Type >  Source_PtrType;

    typedef OneDimensionalModel_Data          Data_Type;
    typedef Data_Type::Mesh_Type              Mesh_Type;

    typedef FESpace< Mesh_Type, EpetraMap >   FESpace_Type;

    typedef SolverAmesos                      LinearSolver_Type;
    typedef LinearSolver_Type::vector_type    Vector_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_BCFunction() {}

    //! Destructor
    virtual ~OneDimensionalModel_BCFunction() {}

    //@}


    //! @name Methods
    //@{

    virtual Real evaluate( const Real& /*time*/ ) = 0;

    //@}
};





//! Riemann - Class which implement Riemann BC for One Dimensional BC Functions
/*!
 *  @author Lucia Mirabella
 */
class Riemann : public OneDimensionalModel_BCFunction
{
public:

    typedef OneDimensionalModel_BCFunction                     super;

    //! @name Constructors & Destructor
    //@{

    Riemann( const FESpace_Type&                   fespace,
             const Flux_PtrType                    fluxFun,
             const std::vector<Vector_Type>&       U_thistime,
             const std::string&                    border,
             const std::string&                    var );

    virtual ~Riemann() {}

    //@}

    //! @name Methods
    //@{

    virtual Real evaluate( const Real& time );

    //@}

protected:

    //! @name Protected Methods
    //@{

    void update_U_boundary();

    //@}

    //! Reference to the FESpace
    const FESpace_Type&                      M_FESpace;

    //! Reference to the solver non linear flux functions
    Flux_PtrType                             M_fluxFun;

    //! Reference to the solver current unknowns (U)
    const std::vector<Vector_Type>&          M_U_thistime;

    //! Number of Dof
    UInt                                     M_dimDof;

    //! Boundary Dof (right or left)
    UInt                                     M_boundaryDof;

    //! Value of U at the boundary
    Vec2D                                    M_U_boundary;

    //! Value of W at the boundary
    Vec2D                                    M_W_boundary;

    //! boundary
    std::string                              M_border;

    //! variable
    std::string                              M_var;

    std::map<std::string, OneDBCStringValue> M_oneDBCFunctionsMapStringValues;
};



//! Compatibility - Class which implement Compatibility BC for One Dimensional BC Functions
/*!
 *  @author Lucia Mirabella
 */
class Compatibility : public Riemann
{
public:

    //! @name Type definitions
    //@{

    typedef OneDimensionalModel_BCFunction                       super;

    typedef OneDimensionalModel_BCFunction::Flux_Type            Flux_Type;
    typedef OneDimensionalModel_BCFunction::Flux_PtrType         Flux_PtrType;

    typedef OneDimensionalModel_BCFunction::Source_Type          Source_Type;
    typedef OneDimensionalModel_BCFunction::Source_PtrType       Source_PtrType;

    typedef OneDimensionalModel_BCFunction::Data_Type            Data_Type;
    typedef OneDimensionalModel_BCFunction::Mesh_Type            Mesh_Type;

    typedef OneDimensionalModel_BCFunction::FESpace_Type         FESpace_Type;

    typedef OneDimensionalModel_BCFunction::LinearSolver_Type    LinearSolver_Type;
    typedef OneDimensionalModel_BCFunction::Vector_Type          Vector_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    Compatibility( const FESpace_Type&                   fespace,
                   const Flux_PtrType                    fluxFun,
                   const Source_PtrType                  sourceFun,
                   const std::vector<Vector_Type>&       U_thistime,
                   const Real&                           dt,
                   const std::string&                    border,
                   const std::string&                    var);

    virtual ~Compatibility() {}

    //@}


    //! @name Methods
    //@{

    virtual Real evaluate( const Real& time );

    //@}

protected:

    //! @name Protected Methods
    //@{

    void computeEigenValuesVectors();

    void update_U_internalBd();

    Real extrapolate_L_dot_U( Real const& eigval, Vec2D const& eigvec );

    Real extrapolate_W( OneDBCStringValue const& _W );

    Vec2D _interpolLinear( const Real&  deltaT,
                           const Real&  eigenvalue,
                           const Vec2D& U_bound,
                           const Vec2D& U_intern ) const;

    //@}

    //! Reference to the solver non linear source functions
    Source_PtrType                          M_sourceFun;

    //! Number of elements
    UInt                                    M_nb_elem;

    //! Time step
    const Real&                             M_time_step;

    //! Dof of the internal node adjacent to the boundary
    UInt                                    M_internalBoundaryDof;

    //! Boundary Edge
    Mesh_Type::EdgeType                     M_boundaryEdge;

    //! Boundary point and internal boundary point
    boost::array< Real, NDIM >              M_boundaryPoint;
    boost::array< Real, NDIM >              M_internalBdPoint;

    //! Eigen values of the jacobian diffFlux (= dF/dU = H)
    Real                                    M_eigval1;
    Real                                    M_eigval2;

    //! Left eigen vectors for the eigen values eigval1 and eigval2
    Vec2D                                   M_left_eigvec1;
    Vec2D                                   M_left_eigvec2;

    //! Value of U at the neighboring internal node
    Vec2D                                   M_U_internalBd;
};

}

#endif // ONEDIMENSIONALMODEL_BCFUNCTION_H
