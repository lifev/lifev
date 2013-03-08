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
    @brief Implementation of an  FSI with fixed point iterations.

    @author Miguel Fernandez
    @author Gilles Fourestey
    @date 10-06-2010

    @contributor Simone Deparis <simone.deparis@epfl.ch>
    @maintainer Simone Deparis <simone.deparis@epfl.ch>
 */


#ifndef FSIFIXEDPOINT_HPP
#define FSIFIXEDPOINT_HPP

#include <lifev/core/algorithm/NonLinearAitken.hpp>
#include <lifev/fsi/solver/FSIOperator.hpp>

namespace LifeV
{

//! FSIFixedPont - Implementation of an  FSI with fixed point iterations.
/*!

    @author Miguel Fernandez
    @author Gilles Fourestey
    @author Paolo Crosetto <paolo.crosetto@epfl.ch>
    @see
    \cite DeparisDiscacciati2006 (Dirichlet--Neumann )
    \cite BadiaNobileVergara2008 (Robin Neumann)
    FSI

    This class implements an FSI that will solve the FSI problem by a
    relaxed fixed point method.
*/

class FSIFixedPoint : public FSIOperator
{
public:

    //! @name Public Types
    //@{
    typedef FSIOperator                     super;

    typedef super::vector_Type              vector_Type;
    typedef super::vectorPtr_Type           vectorPtr_type;

    typedef fluid_Type::matrix_Type         matrix_Type;
    typedef fluid_Type::matrixPtr_Type      matrixPtr_Type;


    //typedef super::fluid_Type               fluid_Type;
    typedef super::solid_Type               solid_Type;

    //     //! OBSOLETE typedefs
    //     //typedef super::fluidBchandler_Type      fluidBchandler_Type;

    //     typedef super::fluid_Type               fluid_type;
    //     typedef super::solid_Type               solid_type;
    //     typedef super::vector_Type              vector_Type;

    //     typedef super::fluidBchandler_Type      bchandler_type;

    //     typedef fluid_Type::matrix_Type         matrix_types;
    //     typedef fluid_Type::matrixPtr_Type      matrix_ptrtype;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    FSIFixedPoint();

    //! Destructor
    ~FSIFixedPoint();

    //! @name Methods
    //@{

    //! solves the Jacobian system
    /**
       The implementation of this method distinguish the various FSI formulations which derive from this class.
       For this reason it must be pure virtual, snd implemented in the child classes.
       \param muk: unknown solution at the k-th nonlinear iteration
       \param res: residual vector (the right hand side of the Jacobian system)
       \param linearRelTol: tolerance for the nonlinear solver
       \todo{replace Real with Real& }
     */
    void solveJac     (vector_Type&       _muk,
                       const vector_Type& _res,
                       const Real       _linearRelTol);


    //! Evaluates the nonlinear residual of the FSI system
    /**
       The implementation of this method also depends on the child classes, though it does not characterize them.
       \param res:  residual vector  to be computed
       \param disp: current unknown solution
       \param iter: nonlinear iteration counter. The part of th rhs related to the time discretization is computed only for iter=0
    */
    void evalResidual (vector_Type&        _res,
                       const vector_Type&  _disp,
                       const UInt           _iter);


    //! sets the space discretization parameters
    /*!
      The FE discretization is set accordingly to what specified in the FSIdata member (order of accuracy for the fluid
      pressure, velocity and for the structure).
     */
    void setupFEspace();

    //! setup of the fluid and solid solver classes
    /**
       This method computes the number of fluxes assigned at the boundaries and calls setupFluidSolid(UInt fluxes)
     */
    void setupFluidSolid();

    //! initializes the GetPot data file
    void setDataFile ( GetPot const& data );

    //! register the product for the factory
    void registerMyProducts( );

    // OBSOLETE
    void setUpBC     ();

    //@}


private:

    //! @name Private Methods
    //@{

    void eval ( const vector_Type& disp, UInt status );

    //@}

    NonLinearAitken<vector_Type> M_nonLinearAitken;

    vectorPtr_Type       M_rhsNew;
    vectorPtr_Type       M_beta;


}; // end class fixedPointBase


inline FSIOperator* createFP()
{
    return new FSIFixedPoint();
}
namespace
{
//static bool registerFP = FSIOperator::FSIFactory_Type::instance().registerProduct( "fixedPoint", &createFP );
}

}   // Namespace LifeV


#endif // FIXEDPOINTBASE_HPP
