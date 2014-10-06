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
   @brief VMS stabilization.
   @author Davide Forti <davide.forti@epfl.ch>
   @maintainer Davide Forti <davide.forti@epfl.ch>
   @date 04-02-2014

   This file contains an ETA implementation of SUPG.

   TODO add more comments, in particular describe the formulation adopted. Paper Y. Bazilevs, V.M. Calo, J.A. Cottrell, T.J.R. Hughes, A. Reali, and G. Scovazzi. <i>Variational multiscale
   residual-based turbulence modeling for large eddy simulation of incompressible flows</i>. Comput. Methods Appl. Mech. Engr. 197(1):173–201, 2007.

   TODO check the parallel execution

*/

#ifndef _SUPGSTABILIZATION_HPP_
#define _SUPGSTABILIZATION_HPP_

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

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/ReferenceFE.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>


// MACRO TO DEFINE TAU_M
#define TAU_M 	       value(1)/( eval(squareroot,TAU_M_DEN) )
#define TAU_M_DEN      TAU_M_DEN_DT + TAU_M_DEN_VEL + TAU_M_DEN_VISC
#define TAU_M_DEN_DT   value(M_density*M_density)*value(4)/value(M_timestep * M_timestep)
#define TAU_M_DEN_VEL  value(M_density*M_density)*dot(value(M_fespaceUETA, velocityExtrapolated), value(M_fespaceUETA, velocityExtrapolated))/(h_K*h_K)
#define TAU_M_DEN_VISC value(M_C_I)*value(M_viscosity*M_viscosity)/(h_K*h_K*h_K*h_K)

// MACRO TO DEFINE TAU_C
#define TAU_C (h_K*h_K)/(TAU_M)

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

class Maximum
{
public:
    typedef Real return_Type;

    inline return_Type operator() (const Real& a, const Real& b)
    {
        return std::max(a,b);
    }

    Maximum() {}
    Maximum(const Maximum&) {}
    ~Maximum() {}
};

class FlagTime
{
public:
    typedef Real return_Type;

    inline return_Type operator() (const Real& a, const Real& b )
    {
        if ( M_flag == false )
	    return a;
         else
            return b;
    }

    bool M_flag;

    FlagTime(const bool& flag){ M_flag = flag; }
    FlagTime(const FlagTime&){}
    ~FlagTime() {}
};

template<typename MeshType, typename MapType, UInt SpaceDim>
class StabilizationSUPG
{
public:

    //@name Public Types
    //@{
    typedef MeshType mesh_Type;
    typedef MapType  map_Type;

    typedef FESpace< mesh_Type, map_Type > fespace_Type;
    typedef boost::shared_ptr< fespace_Type > fespacePtr_Type;

    typedef ETFESpace<mesh_Type, map_Type, SpaceDim, SpaceDim > ETFESpace_velocity;
    typedef ETFESpace<mesh_Type, map_Type, SpaceDim, 1 >        ETFESpace_pressure;

    typedef boost::shared_ptr<ETFESpace_velocity > ETFESpacePtr_velocity;
    typedef boost::shared_ptr<ETFESpace_pressure > ETFESpacePtr_pressure;

    //@}

    //! @name Constructor and Destructor
    //@{
    //! Default Constructor
    StabilizationSUPG(FESpace<mesh_Type, MapEpetra>&  velocityFESpace, FESpace<mesh_Type, MapEpetra>&  pressureFESpace);

    //! ~Destructor
    ~StabilizationSUPG() {};
    //@}

    template <typename MatrixBlockType, typename VectorType>
    void applySUPG_Matrix_semi_implicit( MatrixBlockType& matrix, const VectorType& velocityExtrapolation, const Real& alpha );

    template <typename VectorType, typename VectorBlockType >
    void applySUPG_RHS_semi_implicit( VectorBlockType& rhs, const VectorType& velocityExtrapolated, const VectorType& velocityRhs);

    //! Set the constant C_I for the supg
    void setConstant (const int & value);

    //! Set the fluid density
    void setDensity (const Real & density) { M_density = density;}

    //! Set the fluid dynamic viscosity
    void setViscosity (const Real & viscosity) { M_viscosity = viscosity;}

    //! Set the Epetra communicator
    void setCommunicator (boost::shared_ptr<Epetra_Comm> comm) { M_comm = comm;}

    //! Set the time step size
    void setTimeStep  (const Real & timestep)  { M_timestep = timestep;}

    void setETvelocitySpace(const ETFESpacePtr_velocity & velocityEta_fespace){ M_fespaceUETA = velocityEta_fespace;}

    void setETpressureSpace(const ETFESpacePtr_pressure & pressureEta_fespace){ M_fespacePETA = pressureEta_fespace;}

private:

    //! @name Private Attributes
    //@{

    //! finite element spaces for velocity and pressure
    fespace_Type& M_uFESpace;
    fespace_Type& M_pFESpace;

    ETFESpacePtr_velocity M_fespaceUETA;
    ETFESpacePtr_pressure M_fespacePETA;

    //! Epetra communicator
    boost::shared_ptr<Epetra_Comm> M_comm;

    //! fluid dynamic viscosity @f$\nu@f$
    Real         M_viscosity;

    //! fluid density @f$\nu@f$
    Real         M_density;

    //! stabilization parameters for the momentum and continuity equations
    Real         M_tauM;
    Real         M_tauC;
    Real         M_timestep;
    bool         M_flag_timestep;

    Real         M_C_I;

    //@}
}; // class StabilizationVMS

//=============================================================================
// Constructor
//=============================================================================

template<typename MeshType, typename MapType, UInt SpaceDim>
StabilizationSUPG<MeshType, MapType, SpaceDim>::StabilizationSUPG(FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
																  FESpace<mesh_Type, MapEpetra>&  pressureFESpace):
M_uFESpace (velocityFESpace),
M_pFESpace (pressureFESpace)
{
}

//=============================================================================
// Methods
//=============================================================================

template<typename MeshType, typename MapType, UInt SpaceDim>
void StabilizationSUPG<MeshType, MapType, SpaceDim>::setConstant(const int & value)
{
	if ( value == 1 )
		M_C_I = 30;
	else if ( value == 2 )
		M_C_I = 60;
	else
		ASSERT(0!=0, "Please implement a suitable value for M_C_I for your velocity FE order");
}

// applySUPG_semi_implicit
template<typename MeshType, typename MapType, UInt SpaceDim>
template <typename MatrixBlockType, typename VectorType>
void StabilizationSUPG<MeshType, MapType, SpaceDim>::applySUPG_Matrix_semi_implicit( MatrixBlockType& matrix,
                                                                            const VectorType& velocityExtrapolated,
                                                                            const Real& alpha )
{
	boost::shared_ptr<SquareRoot> squareroot(new SquareRoot());

	Real alfa = alpha*M_timestep;

    MatrixSmall<3, 3> Eye;
    Eye *= 0.0;
    Eye[0][0] = 1;
    Eye[1][1] = 1;
    Eye[2][2] = 1;
    
	using namespace ExpressionAssembly;

	integrate(
			elements(M_uFESpace.mesh()),
			M_uFESpace.qr(),
			M_fespaceUETA, // test  w -> phi_i
			M_fespaceUETA, // trial u -> phi_j
			value(M_density*M_density)*TAU_M*value(alfa/M_timestep) * dot( value(M_fespaceUETA, velocityExtrapolated)*grad(phi_i), phi_j )+
			value(M_density*M_density)*TAU_M*dot( value(M_fespaceUETA, velocityExtrapolated)*grad(phi_i), value(M_fespaceUETA, velocityExtrapolated)*grad(phi_j) )
			+ TAU_C*div(phi_i)*div(phi_j)
            - value(M_density*M_viscosity)*TAU_M*dot( value(M_fespaceUETA, velocityExtrapolated)*grad(phi_i), laplacian(phi_j))
	) >> matrix->block(0,0);

	integrate(
			elements(M_uFESpace.mesh()),
			M_uFESpace.qr(),
			M_fespaceUETA, // test  w -> phi_i
			M_fespacePETA, // trial p -> phi_j
			TAU_M*value(M_density)*dot( value(M_fespaceUETA, velocityExtrapolated)*grad(phi_i), grad(phi_j) )
	) >> matrix->block(0,1);

	integrate(
			elements(M_uFESpace.mesh()),
			M_pFESpace.qr(),
			M_fespacePETA, // test  q -> phi_i
			M_fespaceUETA, // trial u -> phi_j
			TAU_M*value(M_density*alfa/M_timestep)*dot( grad(phi_i), phi_j ) +
			TAU_M*value(M_density)*dot( grad(phi_i), value(M_fespaceUETA, velocityExtrapolated)*grad(phi_j) ) -
            value(M_viscosity)*TAU_M*dot(grad(phi_i), laplacian(phi_j))
    ) >> matrix->block(1,0);

	integrate(
			elements(M_uFESpace.mesh()),
			M_pFESpace.qr(),
			M_fespacePETA, // test   q -> phi_i
			M_fespacePETA, // trial  p -> phi_j
			TAU_M*dot(  grad(phi_j), grad(phi_i) )
	) >> matrix->block(1,1);

	matrix->globalAssemble();

} // applyVMS_semi_implicit(...)

// RHS, semi-implicit, VMS-LES based on Residual Momentum
template<typename MeshType, typename MapType, UInt SpaceDim>
template <typename VectorType, typename VectorBlockType >
void StabilizationSUPG<MeshType, MapType, SpaceDim>::applySUPG_RHS_semi_implicit( VectorBlockType& rhs,
																			const VectorType& velocityExtrapolated,
                                                                            const VectorType& velocityRhs)
{

	boost::shared_ptr<SquareRoot> squareroot(new SquareRoot());

	using namespace ExpressionAssembly;

	integrate(
			elements(M_uFESpace.mesh()),
			M_uFESpace.qr(),
			M_fespaceUETA,
			TAU_M*value(M_density*M_density)*dot( value(M_fespaceUETA, velocityExtrapolated)*grad(phi_i), value(M_fespaceUETA, velocityRhs))
	)
	>> rhs->block(0);

	integrate(
			elements(M_uFESpace.mesh()),
			M_pFESpace.qr(),
			M_fespacePETA,
			TAU_M*value(M_density)*dot( grad(phi_i), value(M_fespaceUETA, velocityRhs))
	)
	>> rhs->block(1);

	rhs->globalAssemble();

} // applyRHS_semi_implicit(...)

} // namespace LifeV

#endif /* _SUPGSTABILIZATION_HPP_ */
