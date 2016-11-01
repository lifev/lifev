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
    @brief Stabilization - Abstract interface of stabilizations for Navier-Stokes

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
#include <lifev/core/filter/ExporterHDF5.hpp>

#include <lifev/navier_stokes_blocks/solver/FastAssemblerNS.hpp>

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

    //! @name Methods
    //@{

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
    
    //! Updates the system matrix in Navier-Stokes simulations in fixed coordinates
    //  with semi-implicit treatment of the convective term.
    /*!
     @param velocityExtrapolated extrapolation of the fluid velocity
     */
    virtual void apply_matrix(	const vector_Type& /*velocityExtrapolated*/ ) {};

    //! Updates the system matrix in Navier-Stokes simulations in ALE coordinates
    //  with semi-implicit treatment of the convective term.
    /*!
         @param velocityExtrapolated extrapolation of the fluid velocity
         @param velocityALE ALE velocity
     */
    virtual void apply_matrix(	const vector_Type& /*velocityExtrapolated*/, const vector_Type& /*velocityALE*/ ) {};

    
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
    
    //! @name Setters
    //@{

    //! Set the constant C_I for the stabilization
    /*!
     * @param value order of FE used for velocity
     */
    virtual void setConstant (const int & value) = 0;;

    //! Build the graphs of each single block
    virtual void buildGraphs() {};

    //! Set the fluid density
    /*!
     * @param value order of FE used for velocity
     */
    virtual void setDensity (const Real & density) = 0;

    //! Set the bdf order
    /*!
     * @param bdfOrder order of BDF scheme used
     */
    virtual void setBDForder (const Real & bdfOrder) = 0;

    //! Set the alpha coefficient of the BDF scheme used
    /*!
     * @param alpha coefficient BDF in front of u^{n+1}
     */
    virtual void setAlpha (const Real & alpha) = 0;

    //! Set the fluid dynamic viscosity
    /*!
     * @param viscosity value of fluid viscosity
     */
    virtual void setViscosity (const Real & viscosity) = 0;

    //! Set the Epetra communicator
    /*!
     * @param comm Epetra communicator
     */
    virtual void setCommunicator (boost::shared_ptr<Epetra_Comm> comm) = 0;

    //! Set the time step size
    /*!
     * @param timestep time step size
     */
    virtual void setTimeStep  (const Real & timestep) = 0;

    //! Set FE space for velocity
    /*!
     * @param velocityFESpace FE space for velocity
     */
    virtual void setVelocitySpace(fespacePtr_Type velocityFESpace) = 0;

    //! Set Expression Template FE space for velocity
    /*!
     * @param pressureFESpace FE space for pressure
     */
    virtual void setPressureSpace(fespacePtr_Type pressureFESpace) = 0;

    //! Set Expression Template FE space for velocity
    /*!
     * @param velocityEta_fespace Expression Template FE space for velocity
     */
    virtual void setETvelocitySpace(const ETFESpacePtr_velocity & velocityEta_fespace) = 0;

    //! Set Expression Template FE space for velocity
    /*!
     * @param pressureEta_fespace Expression Template FE space for pressure
     */
    virtual void setETpressureSpace(const ETFESpacePtr_pressure & pressureEta_fespace) = 0;

    //! Set if using the graph
    /*!
     * @param useGraph true it uses the graph, false it does not use the graph
     */
    virtual void setUseGraph (const bool& /*useGraph*/) {};
    
    //! Set if using dynamic fine scale model
    virtual void setUseODEfineScale ( const bool& /*M_useODEfineScale*/ ) {};

    //! Set if the user wants to export the fine scale component
    virtual void setExportFineScaleVelocity ( ExporterHDF5<mesh_Type> & /*exporter*/, const int& /*numElementsTotal*/ ) {};

    //! Set if using the fast assembler
    virtual void setFastAssembler ( boost::shared_ptr<FastAssemblerNS>&  /*fast_assembler*/) {};

    //@}

    //! @name Getters
    //@{

    //! Get block00 of the stabilization matrix
    virtual matrixPtr_Type const& block_00() const = 0;

    //! Get block01 of the stabilization matrix
    virtual matrixPtr_Type const& block_01() const = 0;

    //! Get block10 of the stabilization matrix
    virtual matrixPtr_Type const& block_10() const = 0;

    //! Get block11 of the stabilization matrix
    virtual matrixPtr_Type const& block_11() const = 0;

    //! Get name of stabilization being used
    virtual std::string label () = 0;

    //! Updates the fine scale component
    virtual void updateODEfineScale ( const vectorPtr_Type& /*velocity*/, const vectorPtr_Type& /*pressure*/ ) {};

    //! Updates the fine scale component
    virtual void updateODEfineScale ( const vectorPtr_Type& /*velocity*/, const vectorPtr_Type& /*pressure*/, const vectorPtr_Type& /*vel_extrap*/ ) {};

    //@}

private:

    //! @name Private methods
    //@{

    //! Setup of the fine scale component
    virtual void setupODEfineScale () {};

    //! Compute the fine component of the velocity and pressure
    virtual void computeFineScales ( const vectorPtr_Type& /*velocity*/, const vectorPtr_Type& /*pressure*/, const vectorPtr_Type& /*vel_extrap*/ ) {};

    //! Compute the fine component of the velocity and pressure for visualization
    virtual void computeFineScalesForVisualization ( const vectorPtr_Type& /*velocity*/, const vectorPtr_Type& /*pressure*/, const vectorPtr_Type& /*vel_extrap*/ ) {};

    //@}

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
