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
    @brief NavierStokesPreconditionerOperator - Abstract interface for preconditioners for Navier-Stokes

    @author Davide Forti <davide.forti@epfl.ch>
    @date 08-12-2014

    @maintainer Davide Forti <davide.Forti@epfl.ch>
 */

#ifndef STABILIZATION_HPP
#define STABILIZATION_HPP 1

#include <lifev/core/util/Factory.hpp>
#include <lifev/core/util/FactorySingleton.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/ReferenceFE.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

namespace LifeV
{

class Stabilization
{
public:

    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef MapEpetra  map_Type;

    typedef VectorEpetra  vector_Type;
    typedef boost::shared_ptr<vector_Type>  vectorPtr_Type;

    typedef MatrixEpetra<Real>  matrix_Type;
    typedef boost::shared_ptr<matrix_Type>  matrixPtr_Type;

    typedef FESpace< mesh_Type, map_Type > fespace_Type;
    typedef boost::shared_ptr< fespace_Type > fespacePtr_Type;

    typedef ETFESpace<mesh_Type, map_Type, 3, 3 > ETFESpace_velocity;
    typedef ETFESpace<mesh_Type, map_Type, 3, 1 > ETFESpace_pressure;

    typedef boost::shared_ptr<ETFESpace_velocity > ETFESpacePtr_velocity;
    typedef boost::shared_ptr<ETFESpace_pressure > ETFESpacePtr_pressure;

    Stabilization();

    virtual ~Stabilization();

    //! Build the graphs of each single block
    virtual void buildGraphs() = 0;

    //! @name Interfaces for the NS in fixed domain, fully implicit
    //@{
    
    //! Updates the jacobian matrix in Navier-Stokes simulations in fixed coordinates
    /*!
       @param velocity_previous_newton_step velocity from the previous Newton step
       @param pressure_previous_newton_step pressure from the previous Newton step
       @param velocity_rhs velocity term from approximation time derivative
     
    virtual void apply_matrix(	const vector_Type& velocity_previous_newton_step,
    							const vector_Type& pressure_previous_newton_step,
    							const vector_Type& velocity_rhs ) {};
    */
    
    //! Adds to the residual the contribution coming from the SUPG stabilization
    //  in Navier-Stokes simulations in fixed coordinates
    /*!
       @param residual_velocity velocity component of the residual
       @param residual_pressure pressure component of the residual
       @param velocity_previous_newton_step velocity from the previous Newton step
       @param pressure_previous_newton_step pressure from the previous Newton step
       @param velocity_rhs velocity term from approximation time derivative
     */
    virtual void apply_vector( 	vectorPtr_Type& /*residual_velocity*/,
    							vectorPtr_Type& /*residual_pressure*/,
    							const vector_Type& /*velocity_previous_newton_step*/,
    							const vector_Type& /*pressure_previous_newton_step*/,
    							const vector_Type& /*velocity_rhs*/) {};
    
    //@}
    
    //! @name Interfaces for the NS in fixed domain, semi-implicit
    //@{
    
    //! Updates the system matrix in Navier-Stokes simulations in fixed coordinates
    //  with semi-implicit treatment of the convective term.
    /*!
     @param velocityExtrapolated extrapolation of the fluid velocity
     */
    virtual void apply_matrix(	const vector_Type& /*velocityExtrapolated*/ ) {};

    
    //! Adds to the right hand side the contribution coming from the SUPG stabilization
    //  in Navier-Stokes simulations in fixed coordinates. Usef for NS semi-implicit
    /*!
     @param rhs_velocity velocity component of the right hand side
     @param rhs_pressure pressure component of the right hand side
     @param velocity_extrapolated velocity extrapolated
     @param velocity_rhs velocity term from approximation time derivative
     */
    virtual void apply_vector( 	vectorPtr_Type& /*rhs_velocity*/,
                                vectorPtr_Type& /*rhs_pressure*/,
                                const vector_Type& /*velocity_extrapolated*/,
                                const vector_Type& /*velocity_rhs*/) {};

    //@}
    
    //! @name Interfaces for the NS in fixed domain, VMSLES semi-implicit
    //@{
    
    //! Updates the system matrix in Navier-Stokes simulations in fixed coordinates
    //  with semi-implicit treatment of the convective term.
    /*!
     @param velocityExtrapolated extrapolation of the fluid velocity
     @param pressureExtrapolated extrapolation of the fluid pressure
     @param velocity_rhs velocity term from approximation time derivative
     */
    virtual void apply_matrix(	const vector_Type& /*velocityExtrapolated*/,
                                const vector_Type& /*pressureExtrapolated*/,
                                const vector_Type& /*velocity_rhs*/) {};
    //@}
    
    //! @name Interfaces for the NS in moving domain, fully implicit
    //@{
    
    //! Updates the jacobian matrix in Navier-Stokes simulations in ALE coordinates
    /*!
     * @param convective_velocity_previous_newton_step convective velocity from the previous Newton step
     * @param velocity_previous_newton_step velocity from the previous Newton step
     * @param pressure_previous_newton_step pressure from the previous Newton step
     * @param velocity_rhs velocity term from approximation time derivative
     */
    virtual void apply_matrix(	const vector_Type& /*convective_velocity_previous_newton_step*/,
    							const vector_Type& /*velocity_previous_newton_step*/,
    							const vector_Type& /*pressure_previous_newton_step*/,
    							const vector_Type& /*velocity_rhs*/ ) {};

    //! Adds to the residual the contribution coming from the SUPG stabilization
    //  in Navier-Stokes simulations in ALE coordinates
    /*!
     * @param residual_velocity velocity component of the residual
     * @param residual_pressure pressure component of the residual
     * @param convective_velocity_previous_newton_step convective velocity from the previous Newton step
     * @param velocity_previous_newton_step velocity from the previous Newton step
     * @param pressure_previous_newton_step pressure from the previous Newton step
     * @param velocity_rhs velocity term from approximation time derivative
     */
    virtual void apply_vector(  vectorPtr_Type& /*residual_velocity*/,
    						    vectorPtr_Type& /*residual_pressure*/,
    						    const vector_Type& /*convective_velocity_previous_newton_step*/,
    						    const vector_Type& /*velocity_previous_newton_step*/,
    						    const vector_Type& /*pressure_previous_newton_step*/,
    						    const vector_Type& /*velocity_rhs*/) {};

    //@}
    
    //! Set the constant C_I for the supg
    virtual void setConstant (const int & value) = 0;;

    //! Set the fluid density
    virtual void setDensity (const Real & density) = 0;

    //! Set the bdf order
    virtual void setBDForder (const Real & bdfOrder) = 0;

    //! Set the bdf order
    virtual void setAlpha (const Real & alpha) = 0;

    //! Set the fluid dynamic viscosity
    virtual void setViscosity (const Real & viscosity) = 0;

    //! Set the Epetra communicator
    virtual void setCommunicator (boost::shared_ptr<Epetra_Comm> comm) = 0;

    //! Set the time step size
    virtual void setTimeStep  (const Real & timestep) = 0;

    virtual void setVelocitySpace(fespacePtr_Type velocityFESpace) = 0;

    virtual void setPressureSpace(fespacePtr_Type pressureFESpace) = 0;

    virtual void setETvelocitySpace(const ETFESpacePtr_velocity & velocityEta_fespace) = 0;

    virtual void setETpressureSpace(const ETFESpacePtr_pressure & pressureEta_fespace) = 0;

    //! @name Getters
    //@{

    virtual matrixPtr_Type const& block_00() const = 0;

    virtual matrixPtr_Type const& block_01() const = 0;

    virtual matrixPtr_Type const& block_10() const = 0;

    virtual matrixPtr_Type const& block_11() const = 0;

private:

};

inline Stabilization::Stabilization()
{
}

inline Stabilization::~Stabilization()
{
}

typedef FactorySingleton<Factory<Stabilization, std::string> >  StabilizationFactory;

} // namespace LifeV


#endif // STABILIZATION_HPP
