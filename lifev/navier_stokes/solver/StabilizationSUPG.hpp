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
#define TAU_M 		   value(1)/( eval(squareroot,TAU_M_DEN) )
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

    template <typename MatrixBlockType, typename VectorBlockType, typename VectorType>
    void applySUPG_Matrix_semi_implicit( MatrixBlockType& matrix,
    		const VectorBlockType& velocityExtrapolation,
    		const VectorType& velocityPreviousStep,
    		const Real& alpha );

    template <typename VectorType, typename VectorBlockType >
    void applySUPG_RHS_semi_implicit( VectorType& rhs,
    		const VectorType& velocityRhs,
    		const VectorType& velocityExtrapolation,
    		const VectorType& velocityPreviousStep,
    		const VectorType& pressurePreviousStep  );

    //! Set the fluid density
    void setDensity (const Real & density) { M_density = density;}

    //! Set the fluid dynamic viscosity
    void setViscosity (const Real & viscosity) { M_viscosity = viscosity;}

    //! Set the Epetra communicator
    void setCommunicator (boost::shared_ptr<Epetra_Comm> comm) { M_comm = comm;}

    //! Set the FE space for velocity
    void setVelocitySpace(fespace_Type & uFESpace) { M_uFESpace = uFESpace; };

    //! Set the FE space for pressure
    void setPressureSpace(fespace_Type & pFESpace) { M_pFESpace = pFESpace; }

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

    Real         M_C_I = 30;

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

// applyVMS_semi_implicit
template<typename MeshType, typename MapType, UInt SpaceDim>
template <typename MatrixBlockType, typename VectorBlockType, typename VectorType>
void StabilizationSUPG<MeshType, MapType, SpaceDim>::applySUPG_Matrix_semi_implicit( MatrixBlockType& matrix,
                                                                            const VectorBlockType& velocityExtrapolated,
                                                                            const VectorType& velocityRhs,
                                                                            const Real& alpha )
{

	boost::shared_ptr<SquareRoot> squareroot(new SquareRoot());

	using namespace ExpressionAssembly;

	integrate(
			elements(M_fespaceUETA->mesh()),
			M_uFESpace.qr(),
			M_fespaceUETA, // trial u -> phi_i
			M_fespaceUETA, // test  w -> phi_j
			value(M_density*M_density)*TAU_M*value(alpha/M_timestep) * dot( value(M_fespaceUETA, velocityExtrapolated)*grad(phi_j) , phi_i )
			+ value(M_density*M_density) * TAU_M* dot( value(M_fespaceUETA, velocityExtrapolated)*grad(phi_j), value(M_fespaceUETA, velocityExtrapolated)*grad(phi_i))
			+ TAU_C * div(phi_j) * div(phi_i)
	) >> matrix->block(0,0);

	integrate(
			elements(M_fespaceUETA->mesh()),
			M_uFESpace.qr(),
			M_fespacePETA, // trial p -> phi_i
			M_fespaceUETA, // test  w -> phi_j
			TAU_M*value(M_density)*dot( value(M_fespaceUETA, velocityExtrapolated)*grad(phi_j), grad(phi_i) )
		) >> matrix->block(0,1);

	/*
    checkFESpaces();

    etaUspacePtr_Type ETuFESpace( new etaUspace_Type( M_uFESpace->mesh(), &(M_uFESpace->refFE()), M_comm ) );
    etaPspacePtr_Type ETpFESpace( new etaPspace_Type( M_pFESpace->mesh(), &(M_pFESpace->refFE()), M_comm ) );

    boost::shared_ptr<MatrixBlockType> stabMatrix( new MatrixBlockType( ETuFESpace->map() | ETpFESpace->map() ) );
    *stabMatrix     *= 0;

    boost::shared_ptr<SquareRoot> squareroot(new SquareRoot());
    boost::shared_ptr<Maximum> maximum(new Maximum());
    boost::shared_ptr<FlagTime> flagTime(new FlagTime(M_flag_timestep));

    using namespace ExpressionAssembly;

    // Stabilization terms: DIV/DIV + SUPG + VMS
    integrate(
                    elements(ETuFESpace->mesh()), M_uFESpace->qr(), ETuFESpace, ETuFESpace,
                    RES_CONTINUITY * DIVDIV_TEST + dot(RES_MOMENTUM_1, SUPG_TEST) + dot(RES_MOMENTUM_1, VMS_TEST)
                    + dot(outerProduct( value( -1.0 ) * RES_MOMENTUM_1, RES_MOMENTUM_STEP_N ), LES_TEST) // Turbulence model, semi_implicit VMS-LES, Momentum residual
    )
    >> stabMatrix->block(0,0);

    // Stabilization terms: SUPG + VMS
    integrate(
                    elements(ETuFESpace->mesh()), M_uFESpace->qr(), ETuFESpace, ETpFESpace,
                    dot(RES_MOMENTUM_2, SUPG_TEST) + dot(RES_MOMENTUM_2, VMS_TEST)
                    + dot(outerProduct( value( -1.0 ) * RES_MOMENTUM_2, RES_MOMENTUM_STEP_N ), LES_TEST) // Turbulence model, semi_implicit VMS-LES, Momentum residual
    )
    >> stabMatrix->block(0,1);


    // Stabilization terms: PSPG
    integrate(
                    elements(ETuFESpace->mesh()), M_uFESpace->qr(), ETpFESpace, ETuFESpace,
                    dot(RES_MOMENTUM_1, PSPG_TEST)
    )
    >> stabMatrix->block(1,0);

    // Stabilization terms: PSPG
    integrate(
                    elements(ETuFESpace->mesh()), M_uFESpace->qr(), ETpFESpace, ETpFESpace,
                    dot(RES_MOMENTUM_2, PSPG_TEST)
    )
    >> stabMatrix->block(1,1);

    stabMatrix->globalAssemble();
    matrix += *stabMatrix;
	*/
} // applyVMS_semi_implicit(...)

// RHS, semi-implicit, VMS-LES based on Residual Momentum
template<typename MeshType, typename MapType, UInt SpaceDim>
template <typename VectorType, typename VectorBlockType >
void StabilizationSUPG<MeshType, MapType, SpaceDim>::applySUPG_RHS_semi_implicit( VectorType& rhs,
                                                                            const VectorType& velocityRhs,
                                                                            const VectorType& velocityExtrapolation,
                                                                            const VectorType& velocityPreviousStep,
                                                                            const VectorType& pressurePreviousStep  )
{
    /*
    checkFESpaces();

    etaUspacePtr_Type ETuFESpace( new etaUspace_Type( M_uFESpace->mesh(), &(M_uFESpace->refFE()), M_comm ) );
    etaPspacePtr_Type ETpFESpace( new etaPspace_Type( M_pFESpace->mesh(), &(M_pFESpace->refFE()), M_comm ) );

    VectorBlockType NSRhs( ETuFESpace->map() | ETpFESpace->map(), Repeated );
    NSRhs *= 0.0;

    using namespace ExpressionAssembly;

    boost::shared_ptr<SquareRoot> squareroot(new SquareRoot());
    boost::shared_ptr<Maximum> maximum(new Maximum());
    boost::shared_ptr<FlagTime> flagTime(new FlagTime(M_flag_timestep));

    // Consistency terms: SUPG + VMS-LES
    // Changed by Davide: use DVDT_COMPATIBILITY instead of DVDT_COMPATIBILITY_N
    integrate(
                    elements(ETuFESpace->mesh()), M_uFESpace->qr(), ETuFESpace,
                    dot(DVDT_COMPATIBILITY, SUPG_TEST) + dot(DVDT_COMPATIBILITY, VMS_TEST)
                    + dot(outerProduct( value( -1.0 ) * DVDT_COMPATIBILITY, RES_MOMENTUM_STEP_N ), LES_TEST) // Turbulence model, semi_implicit VMS-LES, Momentum residual
    )
    >> NSRhs.block(0);

    // Consistency term: PSPG
    integrate(
                    elements(ETuFESpace->mesh()), M_pFESpace->qr(), ETpFESpace,
                    dot(DVDT_COMPATIBILITY, PSPG_TEST)
    )
    >> NSRhs.block(1);


    VectorBlockType NSRhsUnique( NSRhs, Unique );
    rhs += NSRhsUnique;
	*/

} // applyRHS_semi_implicit(...)

} // namespace LifeV

#endif /* _SUPGSTABILIZATION_HPP_ */
