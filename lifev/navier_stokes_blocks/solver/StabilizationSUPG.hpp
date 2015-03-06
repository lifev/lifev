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
   @brief SUPG stabilization.
   @author Davide Forti <davide.forti@epfl.ch>
   @maintainer Davide Forti <davide.forti@epfl.ch>
   @date 03-02-2015

   This file contains an ETA implementation of SUPG, fully implicit for the moment.

   For reference, see Paper Y. Bazilevs, V.M. Calo, J.A. Cottrell, T.J.R. Hughes, A. Reali, and G. Scovazzi. <i>Variational multiscale
   residual-based turbulence modeling for large eddy simulation of incompressible flows</i>. Comput. Methods Appl. Mech. Engr. 197(1):173–201, 2007.

*/

#ifndef _STABILIZATIONSUPG_HPP_
#define _STABILIZATIONSUPG_HPP_ 1

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
    #include <mpi.h>
    #include <Epetra_MpiComm.h>
#else
    #include <Epetra_SerialComm.h>
#endif

#include <Epetra_FECrsMatrix.h>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/navier_stokes_blocks/solver/Stabilization.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/ReferenceFE.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

// includes for building the matrix graph
#include <Epetra_FECrsGraph.h>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/eta/expression/BuildGraph.hpp>

namespace LifeV
{

class SquareRoot
{
public:
    typedef Real return_Type;

    inline return_Type operator() (const Real& a)
    {
        return std::sqrt(a);
    }

    SquareRoot() {}
    SquareRoot (const SquareRoot&) {}
    ~SquareRoot() {}
};

class StabilizationSUPG : public Stabilization
{
public:

    //@name Public Types
    //@{
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

    typedef Epetra_FECrsGraph graph_Type;
    typedef boost::shared_ptr<Epetra_FECrsGraph> graphPtr_Type;

    //@}

    //! @name Constructor and Destructor
    //@{
    //! Default Constructor
    StabilizationSUPG();

    //! ~Destructor
    virtual ~StabilizationSUPG(){}

    //@}

    //! Build the graphs of each single block
    void buildGraphs();

    //! Updates the jacobian matrix
    /*!
       @param velocity_previous_newton_step velocity from the previous Newton step
       @param pressure_previous_newton_step pressure from the previous Newton step
       @param velocity_rhs velocity term from approximation time derivative
     */
    void apply_matrix(	const vector_Type& velocity_previous_newton_step,
	  	  	  	  		const vector_Type& pressure_previous_newton_step,
	  	  	  	  		const vector_Type& velocity_rhs );

    //! Adds to the residual the contribution coming from the SUPG stabilization
    /*!
       @param residual_velocity velocity component of the residual
       @param residual_pressure pressure component of the residual
       @param velocity_previous_newton_step velocity from the previous Newton step
       @param pressure_previous_newton_step pressure from the previous Newton step
       @param velocity_rhs velocity term from approximation time derivative
     */
    void apply_vector( vectorPtr_Type& residual_velocity,
			  	   	   vectorPtr_Type& residual_pressure,
			  	   	   const vector_Type& velocity_previous_newton_step,
			  	   	   const vector_Type& pressure_previous_newton_step,
			  	   	   const vector_Type& velocity_rhs);

    void setVelocitySpace(fespacePtr_Type velocityFESpace){ M_uFESpace = velocityFESpace;}

    void setPressureSpace(fespacePtr_Type pressureFESpace){ M_pFESpace = pressureFESpace;}

    //! Set the constant C_I for the supg
    void setConstant (const int & value);

    //! Set the fluid density
    void setDensity (const Real & density) { M_density = density;}

    //! Set the bdf order
    void setBDForder (const Real & bdfOrder) { M_bdfOrder = bdfOrder;}

    //! Set the bdf order
    void setAlpha (const Real & alpha) { M_alpha = alpha;}

    //! Set the fluid dynamic viscosity
    void setViscosity (const Real & viscosity) { M_viscosity = viscosity;}

    //! Set the Epetra communicator
    void setCommunicator (boost::shared_ptr<Epetra_Comm> comm) { M_comm = comm;}

    //! Set the time step size
    void setTimeStep  (const Real & timestep)  { M_timestep = timestep;}

    void setETvelocitySpace(const ETFESpacePtr_velocity & velocityEta_fespace){ M_fespaceUETA = velocityEta_fespace;}

    void setETpressureSpace(const ETFESpacePtr_pressure & pressureEta_fespace){ M_fespacePETA = pressureEta_fespace;}

    //! @name Getters
    //@{

    matrixPtr_Type const& block_00() const
    {
    	return M_block_00;
    }

    matrixPtr_Type const& block_01() const
    {
    	return M_block_01;
    }

    matrixPtr_Type const& block_10() const
    {
    	return M_block_10;
    }

    matrixPtr_Type const& block_11() const
    {
    	return M_block_11;
    }

    //@}

private:

    //! @name Private Attributes
    //@{

    //! finite element spaces for velocity and pressure
    fespacePtr_Type M_uFESpace;
    fespacePtr_Type M_pFESpace;

    ETFESpacePtr_velocity M_fespaceUETA;
    ETFESpacePtr_pressure M_fespacePETA;

    //! Epetra communicator
    boost::shared_ptr<Epetra_Comm> M_comm;

    //! fluid dynamic viscosity @f$\nu@f$
    Real         M_viscosity;

    //! fluid density @f$\nu@f$
    Real         M_density;

    //! stabilization parameters for the momentum and continuity equations
    Real         M_timestep;
    bool         M_flag_timestep;
    Real         M_bdfOrder;
    Real         M_alpha;

    Real         M_C_I;

    // Graphs
    graphPtr_Type M_graph_block00;
    graphPtr_Type M_graph_block01;
    graphPtr_Type M_graph_block10;
    graphPtr_Type M_graph_block11;

    // Matrices
    matrixPtr_Type M_block_00;
    matrixPtr_Type M_block_01;
    matrixPtr_Type M_block_10;
    matrixPtr_Type M_block_11;

    //@}
}; // class StabilizationSUPG

//! Factory create function
inline StabilizationSUPG * createStabilizationSUPG()
{
    return new StabilizationSUPG();
}
namespace
{
static bool S_registerStabilizationSUPG = StabilizationFactory::instance().registerProduct ( "SUPG", &createStabilizationSUPG );
}

} // namespace LifeV

#endif /* _STABILIZATIONSUPG_HPP_ */
