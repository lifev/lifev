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
   @author Toni Lassila
   @author Luca Dede'
   @contributor Toni Lassila <toni.lassila@epfl.ch>
   @contributor Luca Dede' <luca.dede@epfl.ch>
   @contributor Davide Forti <davide.forti@epfl.ch>
   @maintainer Toni Lassila <toni.lassila@epfl.ch>
   @date 06-12-2013

   This file contains an ETA implementation of variational multiscale stabilization (VMS) with "turbulence model" (LES)
   for the incompressible Navier-Stokes equations. The standard VMS-SUPG stabilization is obtained by deactivating the LES terms.
   Tested with P1/P1 only and BDF2 method. "Turbulence model" not tested.

 */

#ifndef _VMSSTABILIZATION_HPP_
#define _VMSSTABILIZATION_HPP_

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

// Preprocessor macros for the test and trial functions including residuals and stabilization parameters computed elementwise    
#define C_I                     value(30) // for P1 elements; for P2: C_I = 60
#define DT2_VEL                 dot(value(ETuFESpace, velocityExtrapolation),value(ETuFESpace, velocityExtrapolation))/(h_K*h_K)
//#define iDT_MIN                 eval( maximum, value((M_timestep)/(0.025*0.025)), eval(squareroot,DT2_VEL) ) // set a minimum time step to avoid stability issues with the time step is too small, hardcoded
#define iDT_MIN                 eval( flagTime, value((1.0)/(M_timestep)), eval(squareroot,DT2_VEL) ) // set a minimum time step for the ramp and modified time step after initial ramp
#define TAU_M_DEN_DT            (M_density*M_density)*value(4)/(M_timestep)*iDT_MIN   // (M_density*M_density)*value(4)/(M_timestep*M_timestep) 
#define TAU_M_DEN_VEL           (M_density*M_density) * DT2_VEL
#define TAU_M_DEN_VISC          C_I*(M_viscosity*M_viscosity)/(h_K*h_K*h_K*h_K)
#define TAU_M_DEN               TAU_M_DEN_DT + TAU_M_DEN_VEL + TAU_M_DEN_VISC
#define TAU_M                   value(1)/( eval(squareroot,TAU_M_DEN) ) // stabilization parameter tau_M (with density and dynamic viscosity)
#define TAU_C                   (h_K*h_K)/(TAU_M)               // stabilization parameter tau_C
#define SUPG_TEST               TAU_M * M_density * ( grad(phi_i) * value(ETuFESpace, velocityExtrapolation) ) 
#define VMS_TEST                TAU_M * M_density * ( transpose(grad(phi_i)) * value(ETuFESpace, velocityExtrapolation) )
#define PSPG_TEST               TAU_M * grad(phi_i)
#define DIVDIV_TEST             TAU_C * div(phi_i)
#define LES_TEST                TAU_M * TAU_M * M_density * grad(phi_i)
#define RES_MOMENTUM_1          M_density * ( value(ETuFESpace, velocityExtrapolation) * grad(phi_j) ) + ( M_density / M_timestep ) * phi_j  // Momentum residual, components depending on velocity
#define RES_MOMENTUM_2          grad(phi_j) // Momentum residual, component depending on pressure
#define RES_CONTINUITY          div(phi_j)  // Continuity residual
#define RES_MOMENTUM_EXPLICIT   value(ETuFESpace, momentumResidual) // Momentum residual for explicit VMS-LES with RHS
#define RES_MOMENTUM_1_STEP_N   M_density * ( grad(ETuFESpace, velocityPreviousStep) * value(ETuFESpace, velocityPreviousStep) ) + ( M_density / M_timestep ) * ( value(ETuFESpace, velocityExtrapolation) - value(ETuFESpace, velocityPreviousStep) )
#define RES_MOMENTUM_2_STEP_N   grad(ETpFESpace,pressurePreviousStep)
#define RES_MOMENTUM_STEP_N     RES_MOMENTUM_1_STEP_N + RES_MOMENTUM_2_STEP_N // Momentum residual for explicit VMS-LES for Momentum residual, previous time step
#define DVDT_COMPATIBILITY_N    ( M_density / M_timestep ) * value(ETuFESpace, velocityPreviousStep)   // With velocity at the previous time-step.
#define DVDT_COMPATIBILITY_EX   ( M_density / M_timestep ) * value(ETuFESpace, velocityExtrapolation)  // Extrapolated velocity is used.


namespace LifeV
{

//! StabilizationVMS Class
/*!
 * @brief Variational multiscale stabilization for Navier-Stokes
 * @author T. Lassila
 * @author L. Dede'
 *
 * Implementation of variational multiscale stabilization method for the Navier-Stokes equations slved with BDF methods. <br>
 *
 *  Remarks about the implementation:
 *  <ol>
 *  <li> Version of VMS proposed by Bazilevs et al that replaces fine-scale problems by elementwise approximate inversion of the local fine-scale Green's operator.
 *  <li> Does not include the terms @f$2 \mu \textrm{div}(\varepsilon(\mathbf{u}))@f$ in the elementwise residual, should only be used with @f$P_1/P_1@f$
 *  velocity-pressure spaces.
 *  <li> Does not include the forcing term @f$f@f$ in the momentum equation.
 *  TODO: Add the second order terms in the residual. Add forcing term in the momentum equation and residual.
 *  <li> Stabilization terms are for the momentum equation without density rescaling (density and dynamic viscosity are considered in the momentum equation).
 *  <li> Stabilization parameters @f$\tau_M@f$ and @f$\tau_C@f$ computed elementwise (in the quadrature nodes) of the entire fluid domain. Scaling relations for the parameters are<br>
 *  @f$\tau_M \sim \left( \rho^2 \frac{4}{(\Delta t)^2} + \frac{\rho^2 \|\mathbf{u}\|^2}{h_K^2} + C_I \frac{\mu^2}{h_K^4} \right)^{-1/2}@f$;<br>
 *  @f$\tau_C \sim \left( \frac{\tau_M}{h_K^2} \right)^{-1}@f$,<br>
 *  where @f$C_I \approx 30@f$ is given by the inverse constant of the chosen FE (for P1). 
 *  <li> Includes "turbulence model" (fine-fine term of VMS) explicitly in the rhs, can be toggled on/off by the user.
 *  <li> The VMS-LES fine terms of the "turbulence" model are also implemented in semi-implicit form. 
 *  </ol>
 * For more details, see: Y. Bazilevs, V.M. Calo, J.A. Cottrell, T.J.R. Hughes, A. Reali, and G. Scovazzi. <i>Variational multiscale
 * residual-based turbulence modeling for large eddy simulation of incompressible flows</i>. Comput. Methods Appl. Mech. Engr. 197(1):173–201, 2007.
 * For notation, see: S. Quinodoz, <i>Numerical Simulation of Orbitally Shaken Reactors.</i> PhD thesis EPFL 5477. EPFL, 2012. 
 *
 */

class ShowValue
{
public:
    typedef Real return_Type;

    return_Type operator() (const Real& a)
    {
      std::cout.precision(15);
      std::cout << "value is: " << a <<" \n";
      return 1.0;
    }

    return_Type operator() (const VectorSmall<3>& a)
    {
        std::cout << "value is: " << a[0]  <<" \n" << a[1]  <<" \n" << a[2]  <<" \n";
        return 1.0;
    }

    return_Type operator() (const MatrixSmall<3,3>& a)
    {
        std::cout << "value is\n";
        std::cout << a[0][0]  <<" \t" << a[0][1]  <<" \t" << a[0][2]  <<" \n";
        std::cout << a[1][0]  <<" \t" << a[1][1]  <<" \t" << a[1][2]  <<" \n";
        std::cout << a[2][0]  <<" \t" << a[2][1]  <<" \t" << a[2][2]  <<" \n";

        return 1.0;
    }



    ShowValue() {}
    ShowValue (const ShowValue&) {}
    ~ShowValue() {}
};

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
class StabilizationVMS
{
public:

    //@name Public Types
    //@{
    typedef MeshType mesh_Type;
    typedef MapType  map_Type;

    typedef FESpace< mesh_Type, map_Type > fespace_Type;
    typedef boost::shared_ptr< fespace_Type > fespacePtr_Type;

    typedef ETFESpace< mesh_Type, map_Type, SpaceDim, SpaceDim > etaUspace_Type;
    typedef ETFESpace< mesh_Type, map_Type, SpaceDim, 1 > etaPspace_Type;

    typedef boost::shared_ptr<ETFESpace< mesh_Type, map_Type, SpaceDim, SpaceDim > > etaUspacePtr_Type;
    typedef boost::shared_ptr<ETFESpace< mesh_Type, map_Type, SpaceDim, 1 > > etaPspacePtr_Type;

    //@}

    //! @name Constructor and Destructor
    //@{
    //! Default Constructor
    StabilizationVMS( );

    //! ~Destructor
    virtual ~StabilizationVMS() {};
    //@}

    //! @name Methods
    //@{
    //! Compute the VMS stabilization terms and add them into the Navier-Stokes matrix
    /*!
     *  Add the following stabilization terms to the Navier-Stokes system:
     *  <ol>
     *  <li> Block(1,1):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M [\rho \mathbf{u}^* \cdot (\nabla + \nabla^T) ] \mathbf{v} \cdot [\rho \frac{1}{\Delta t} \mathbf{u} + \rho (\mathbf{u}^* \cdot \nabla) \mathbf{u}] + \tau_C (\nabla \cdot \mathbf{v}) (\nabla \cdot \mathbf{u}) \: d\Omega@f$;
     *  <li> Block(1,2):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M [\rho \mathbf{u}^* \cdot (\nabla + \nabla^T) ] \mathbf{v} \cdot \nabla p \: d\Omega@f$.
     *  <li> Block(2,1):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M \nabla q \cdot [\rho \frac{1}{\Delta t} \mathbf{u} + \rho (\mathbf{u}^* \cdot \nabla) \mathbf{u}] \: d\Omega@f$.
     *  <li> Block(2,2):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M \nabla q \cdot \nabla p \: d\Omega@f$.
     *  </ol>
     *
     *  @param matrix MatrixBlockType where the stabilization terms are added
     *  @param velocityExtrapolation VectorBlockType velocity field @f$\mathbf{u}^*@f$ for the linearization of the stabilization
     */
    template<typename MatrixBlockType, typename VectorBlockType>
    void applyVMS( MatrixBlockType& matrix,
                   const VectorBlockType& velocityExtrapolation );

    //! Compute the VMS stabilization terms and add them into the Navier-Stokes matrix. The VMS-LES terms are treated as semi-implicit/
    /*!
     *  Add the following stabilization terms to the Navier-Stokes system:
     *  <ol>
     *  <li> Block(1,1):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M [\rho \mathbf{u}^* \cdot (\nabla + \nabla^T) ] \mathbf{v} \cdot [\rho \frac{1}{\Delta t} \mathbf{u} + \rho (\mathbf{u}^* \cdot \nabla) \mathbf{u}] + \tau_C (\nabla \cdot \mathbf{v}) (\nabla \cdot \mathbf{u}) \: d\Omega@f$;
     *  <li> Block(1,2):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M [\rho \mathbf{u}^* \cdot (\nabla + \nabla^T) ] \mathbf{v} \cdot \nabla p \: d\Omega@f$.
     *  <li> Block(2,1):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M \nabla q \cdot [\rho \frac{1}{\Delta t} \mathbf{u} + \rho (\mathbf{u}^* \cdot \nabla) \mathbf{u}] \: d\Omega@f$.
     *  <li> Block(2,2):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M \nabla q \cdot \nabla p \: d\Omega@f$.
     *  </ol>
     *  The "turbulence model", VMS-LES in semi-implicit is added in:
     *  <ol>
     *  <li> Block(1,1):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \nabla \textbf{v} : [\rho \tau_M \mathbf{r}_M(\mathbf{u}^*;\mathbf{u}^n,p^n) \otimes \tau_M [- \rho \frac{1}{\Delta t} \mathbf{u} - \rho (\mathbf{u}^* \cdot \nabla) \mathbf{u}] \: d\Omega@f$;
     *  <li> Block(1,2):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \nabla \textbf{v} : [\rho \tau_M \mathbf{r}_M(\mathbf{u}^*;\mathbf{u}^n,p^n) \otimes \tau_M [- \nabla p] \: d\Omega@f$;
     *  </ol>
     *  where: @f$\mathbf{r}_M(\mathbf{u}^*;\mathbf{u}^n,p^n) := \rho \frac{1}{\Delta t} (\mathbf{u}^*-\mathbf{u}^n) + \rho (\mathbf{u}^n \cdot \nabla) \mathbf{u}^n + \nabla p^n@f$ (force term f set equal to zero)
     *
     *  @param matrix MatrixBlockType where the stabilization terms are added
     *  @param velocityExtrapolation VectorBlockType velocity field @f$\mathbf{u}^*@f$ for the linearization of the stabilization
     *  @param velocityPreviousStep  VectorType velocity field @f$\mathbf{u}^n@f$ computed at the previous time step
     *  @param pressurePreviousStep  VectorType pressure field @f$p^n@f$ computed at the previous time step
     */
    template<typename MatrixBlockType, typename VectorBlockType, typename VectorType>
    void applyVMS_semi_implicit( MatrixBlockType& matrix,
                                 const VectorBlockType& velocityExtrapolation,
                                 const VectorType& velocityPreviousStep,
                                 const VectorType& pressurePreviousStep );

    //! Compute the VMS consistency terms and add them into the rhs. This yields the <b>VMS-SUPG</b> stabilization. Force term f assumed to be equal to zero. 
    /*!
     *  Add the following consistency terms to the rhs (consistency terms based on velocity at previous time step):
     *  <ol>
     *  <li> Block(1,1):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M [\rho \mathbf{u}^* \cdot (\nabla + \nabla^T) ] \mathbf{v} \cdot [\rho \frac{1}{\Delta t}\mathbf{u}^n + \mathbf{f}]\: d\Omega@f$;
     *  <li> Block(1,2):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M \nabla q \cdot [\rho \frac{1}{\Delta t}\mathbf{u}^n + \mathbf{f}]\: d\Omega@f$.
     *  </ol>
     *
     *  @param rhs    VectorBlockType the rhs @f$\sum_{i=1}^m \frac{\alpha_j}{\Delta t}\mathbf{u}^{n-1} + \mathbf{f}@f$ of the N-S equations before imposition of the essential b.c.
     *  @param velocityExtrapolation VectorType velocity field @f$\mathbf{u}^*@f$ for the linearization of the stabilization
     *  @param velocityPreviousStep  VectorType velocity field @f$\mathbf{u}^n@f$ computed at the previous time step
     *
     *  This is the implementation <b>without</b> the "turbulence model" (fine-fine terms).
     */
    template <typename VectorType, typename VectorBlockType >
    void applyRHS( VectorType& rhs,
                   const VectorType& velocityExtrapolation,
                   const VectorType& velocityPreviousStep );
    
    //! Compute the VMS consistency and turbulence modelling terms and add them into the rhs; turbulence model set as explicit using the extrapolated momentum residual and extrapolated velocity (<b>VMS-LES explicit</b>, based on Momentum residual).
    /*!
     *  Add the following consistency terms to the rhs (consistency terms based on velocity at previous time step):
     *  <ol>
     *  <li> Block(1,1):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M [\rho \mathbf{u}^* \cdot (\nabla + \nabla^T) ] \mathbf{v} \cdot [\rho \frac{1}{\Delta t}\mathbf{u}^n + \mathbf{f}]\: d\Omega@f$;
     *  <li> Block(1,2):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M \nabla q \cdot [\rho \frac{1}{\Delta t}\mathbf{u}^{n} + \mathbf{f}]\: d\Omega@f$.
     *  </ol>
     *  and the "turbulence model" in
     *  Block(1,1): @f$\sum_{K \in \mathcal{T}_h} \int_{K} \nabla \textbf{v} : [\rho \tau_M \mathbf{r}_M(\mathbf{u}^*;\mathbf{u}^n,p^n) \otimes \tau_M \mathbf{r}_M(\mathbf{u}^*;\mathbf{u}^n,p^n)] \: d\Omega@f$,
     *  where: @f$\mathbf{r}_M(\mathbf{u}^*;\mathbf{u}^n,p^n) := \rho \frac{1}{\Delta t} (\mathbf{u}^*-\mathbf{u}^n) + \rho (\mathbf{u}^n \cdot \nabla) \mathbf{u}^n + \nabla p^n@f$ (force term f set equal to zero)
     *
     *  @param rhs    VectorBlockType the rhs @f$\sum_{j=1}^m \frac{\alpha_j}{\Delta t}\mathbf{u}^{n+1-i} + \mathbf{f}@f$ of the N-S equations before imposition of the essential b.c.
     *  @param velocityExtrapolation VectorType velocity field @f$\mathbf{u}^*@f$ for the linearization of the stabilization and consistency terms
     *  @param velocityPreviousStep  VectorType velocity field @f$\mathbf{u}^n@f$ computed at the previous time step
     *  @param pressurePreviousStep  VectorType pressure field @f$p^n@f$ computed at the previous time step
     *
     *  This is the implementation <b>with</b> the "turbulence model" (fine-fine terms). 
     */
    template <typename VectorType, typename VectorBlockType >
    void applyRHS( VectorType& rhs,
                   const VectorType& velocityExtrapolation,
                   const VectorType& velocityPreviousStep,
                   const VectorType& pressurePreviousStep );
    

    //! Compute the VMS consistency and turbulence modelling terms and add them into the rhs; turbulence model set as semi-implicit using the extrapolated momentum residual and extrapolated velocity (<b>VMS-LES semi-implicit</b>, based on Momentum residual).
    /*!
     *  Add the following consistency terms to the rhs (consistency terms based on velocity at previous time step):
     *  <ol>
     *  <li> Block(1,1):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M [\rho \mathbf{u}^* \cdot (\nabla + \nabla^T) ] \mathbf{v} \cdot [\rho \frac{1}{\Delta t}\mathbf{u}^n + \mathbf{f}]\: d\Omega@f$;
     *  <li> Block(1,2):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M \nabla q \cdot [\rho \frac{1}{\Delta t}\mathbf{u}^n + \mathbf{f}]\: d\Omega@f$.
     *  </ol>
     *  and the "turbulence model" in
     *  Block(1,1): @f$\sum_{K \in \mathcal{T}_h} \int_{K} \nabla \textbf{v} : [\rho \tau_M \mathbf{r}_M(\mathbf{u}^*;\mathbf{u}^n,p^n) \otimes \tau_M \rho \frac{1}{\Delta t} (- \mathbf{u}^n) ] \: d\Omega@f$,
     *  where: @f$\mathbf{r}_M(\mathbf{u}^*;\mathbf{u}^n,p^n) := \rho \frac{1}{\Delta t} (\mathbf{u}^*-\mathbf{u}^n) + \rho (\mathbf{u}^n \cdot \nabla) \mathbf{u}^n + \nabla p^n@f$ (force term f set equal to zero)
     *
     *  @param rhs    VectorBlockType the rhs @f$\sum_{j=1}^m \frac{\alpha_j}{\Delta t}\mathbf{u}^{n+1-i} + \mathbf{f}@f$ of the N-S equations before imposition of the essential b.c.
     *  @param velocityExtrapolation VectorBlockType velocity field @f$\mathbf{u}^*@f$ for the linearization of the stabilization and consistency terms
     *  @param velocityPreviousStep  VectorType velocity field @f$\mathbf{u}^n@f$ computed at the previous time step
     *  @param pressurePreviousStep  VectorType pressure field @f$p^n@f$ computed at the previous time step
     *
     *  This is the implementation <b>with</b> the "turbulence model" (fine-fine terms) in semi-implicit form. 
     */
    template <typename VectorType, typename VectorBlockType >
    void applyRHS_semi_implicit( VectorType& rhs,
	                         const VectorType& velocityExtrapolation,
                                 const VectorType& velocityPreviousStep,
                                 const VectorType& pressurePreviousStep );


    //! Compute the VMS consistency and turbulence modelling terms and add them into the rhs; turbulence model set as explicit using the full <b>residual</b> computed as RHS (<b>VMS-LES explicit, based on RHS</b>). 
    /*!
     *  Add the following consistency terms to the rhs (consistency terms based on extrapolated velocity):
     *  <ol>
     *  <li> Block(1,1):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M [\rho \mathbf{u}^* \cdot (\nabla + \nabla^T) ] \mathbf{v} \cdot [\rho \frac{1}{\Delta t}\mathbf{u}^* + \mathbf{f}]\: d\Omega@f$;
     *  <li> Block(1,2):  @f$\sum_{K \in \mathcal{T}_h} \int_{K} \tau_M \nabla q \cdot [\rho \frac{1}{\Delta t}\mathbf{u}^* + \mathbf{f}]\: d\Omega@f$.
     *  </ol>
     *  and the "turbulence model" in
     *  Block(1,1): @f$\sum_{K \in \mathcal{T}_h} \int_{K} \nabla \textbf{v} : [\rho \tau_M \textrm{Res}(\mathbf{u},p)|_{K} \otimes \tau_M \textrm{Res}(\mathbf{u},p)|_{K}] \: d\Omega@f$
     *
     *  @param rhs    VectorBlockType the rhs @f$\sum_{i=1}^m \frac{\alpha_j}{\Delta t}\mathbf{u}^{n-1} + \mathbf{f}@f$ of the N-S equations before imposition of the essential b.c.
     *  @param velocityExtrapolation VectorType velocity field @f$\mathbf{u}^*@f$ for the linearization of the stabilization and consistency terms
     *  @param momentumResidual VectorType elementwise residual of the momentum equations @f$\textrm{Res}(\mathbf{u},p)|_{K}@f$ at the previous time step (explicit) or previous Newton iteration (semi-implicit)
     *
     *  This is the implementation <b>with</b> the "turbulence model" (fine-fine terms). The user must provide also the momentum residual.
     */
    template <typename VectorType, typename VectorBlockType >
    void applyRHS_residual( VectorType& rhs,
                            const VectorType& velocityExtrapolation,
                            const VectorType& momentumResidual);


    //! Display class informations
    void showMe(std::ostream & output = std::cout) const;
    //@}

    //! @name Set Methods
    //@{
    //! Set the stabilization parameter @f$\tau_M@f$ for the momentum equations
    void setTauM (const Real & tauM) { M_tauM  = tauM;}
    //! Set the stabilization parameter @f$\tau_C@f$ for the continuity equations
    void setTauC  (const Real & tauC)  { M_tauC = tauC;}
    //! Set the fluid density
    void setDensity (const Real & density) { M_density = density;}
    //! Set the fluid dynamic viscosity
    void setViscosity (const Real & viscosity) { M_viscosity = viscosity;}
    //! Set the Epetra communicator
    void setCommunicator (boost::shared_ptr< Epetra_Comm > comm) { M_comm = comm;}
    //! Set the FE space for velocity
    void setVelocitySpace(fespacePtr_Type & uFESpace);
    //! Set the FE space for pressure
    void setPressureSpace(fespacePtr_Type & pFESpace) { M_pFESpace = pFESpace; }
    //! Set the time step size
    void setTimeStep  (const Real & timestep)  { M_timestep = timestep;}
    //! Set the flag for correction of stabilization parameters (end of ramp); boolean true/false
    void setFlagTimeStep  (const bool & flag_timestep)  { M_flag_timestep = flag_timestep;}
    //@}
private:

    //! @name Private Constructor
    //@{
    //! Copy Constructor
    StabilizationVMS(const StabilizationVMS< mesh_Type,map_Type,SpaceDim > & original);
    //@}

    //! @name Private Methods
    //@{
    void checkFESpaces( );
    //@}

    //! @name Private Attributes
    //@{

    //! finite element spaces for velocity and pressure
    fespacePtr_Type M_uFESpace;
    fespacePtr_Type M_pFESpace;

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

    //@}
}; // class StabilizationVMS

//=============================================================================
// Constructor
//=============================================================================

template<typename MeshType, typename MapType, UInt SpaceDim>
StabilizationVMS<MeshType, MapType, SpaceDim>::StabilizationVMS():
                                              M_viscosity( 0.035 ),
                                              M_density( 1 ),
                                              M_tauM ( 0. ),
                                              M_tauC ( 0. ),
                                              M_uFESpace( ),
                                              M_pFESpace( ),
					      M_timestep( 0.01 ),
                                              M_flag_timestep( false ),
                                              M_comm( )
                                              {}

//=============================================================================
// Methods
//=============================================================================
// applyVMS
template<typename MeshType, typename MapType, UInt SpaceDim>
template <typename MatrixBlockType, typename VectorBlockType>
void StabilizationVMS<MeshType, MapType, SpaceDim>::applyVMS( MatrixBlockType& matrix,
                                                              const VectorBlockType& velocityExtrapolation )
{
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
    )
    >> stabMatrix->block(0,0);

    // Stabilization terms: SUPG + VMS
    integrate(
                    elements(ETuFESpace->mesh()), M_uFESpace->qr(), ETuFESpace, ETpFESpace,
                    dot(RES_MOMENTUM_2, SUPG_TEST) + dot(RES_MOMENTUM_2, VMS_TEST)
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

} // applyVMS(...)


// applyVMS_semi_implicit
template<typename MeshType, typename MapType, UInt SpaceDim>
template <typename MatrixBlockType, typename VectorBlockType, typename VectorType>
void StabilizationVMS<MeshType, MapType, SpaceDim>::applyVMS_semi_implicit( MatrixBlockType& matrix,
                                                                            const VectorBlockType& velocityExtrapolation, 
                                                                            const VectorType& velocityPreviousStep,
                							    const VectorType& pressurePreviousStep )
{
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

} // applyVMS_semi_implicit(...)

// RHS, VMS-SUPG
template<typename MeshType, typename MapType, UInt SpaceDim>
template < typename VectorType, typename VectorBlockType >
void StabilizationVMS<MeshType, MapType, SpaceDim>::applyRHS( VectorType& rhs,
                                                              const VectorType& velocityExtrapolation,
                                                              const VectorType& velocityPreviousStep )
{
    checkFESpaces();
    boost::shared_ptr<SquareRoot> squareroot(new SquareRoot());
    boost::shared_ptr<Maximum> maximum(new Maximum());
    boost::shared_ptr<FlagTime> flagTime(new FlagTime(M_flag_timestep));
    
    etaUspacePtr_Type ETuFESpace( new etaUspace_Type( M_uFESpace->mesh(), &(M_uFESpace->refFE()), M_comm ) );
    etaPspacePtr_Type ETpFESpace( new etaPspace_Type( M_pFESpace->mesh(), &(M_pFESpace->refFE()), M_comm ) );

    VectorBlockType NSRhs( ETuFESpace->map() | ETpFESpace->map(), Repeated );
    NSRhs *= 0.0;

    using namespace ExpressionAssembly;

    // Consistency terms: SUPG + VMS
    integrate(
                    elements(ETuFESpace->mesh()), M_uFESpace->qr(), ETuFESpace,
                    dot(DVDT_COMPATIBILITY_N, SUPG_TEST) + dot(DVDT_COMPATIBILITY_N, VMS_TEST)
    )
    >> NSRhs.block(0);

    // Consistency term: PSPG
    integrate(
                    elements(ETuFESpace->mesh()), M_pFESpace->qr(), ETpFESpace,
                    dot(DVDT_COMPATIBILITY_N, PSPG_TEST)
    )
    >> NSRhs.block(1);


    VectorBlockType NSRhsUnique( NSRhs, Unique );
    rhs += NSRhsUnique;

} // applyRHS(...)


// RHS, VMS-LES based on Residual Momentum
template<typename MeshType, typename MapType, UInt SpaceDim>
template <typename VectorType, typename VectorBlockType >
void StabilizationVMS<MeshType, MapType, SpaceDim>::applyRHS( VectorType& rhs,
                                                              const VectorType& velocityExtrapolation,
                                                              const VectorType& velocityPreviousStep,
							      const VectorType& pressurePreviousStep  )
{
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
    integrate(
                    elements(ETuFESpace->mesh()), M_uFESpace->qr(), ETuFESpace, 
                    dot(DVDT_COMPATIBILITY_N, SUPG_TEST) + dot(DVDT_COMPATIBILITY_N, VMS_TEST)
                    + dot(outerProduct( RES_MOMENTUM_STEP_N, RES_MOMENTUM_STEP_N ), LES_TEST) // Turbulence model, explicit VMS-LES, Momentum residual
    )
    >> NSRhs.block(0);

    // Consistency term: PSPG
    integrate(
                    elements(ETuFESpace->mesh()), M_pFESpace->qr(), ETpFESpace,
                    dot(DVDT_COMPATIBILITY_N, PSPG_TEST)
    )
    >> NSRhs.block(1);


    VectorBlockType NSRhsUnique( NSRhs, Unique );
    rhs += NSRhsUnique;
    

} // applyRHS(...)


// RHS, semi-implicit, VMS-LES based on Residual Momentum
template<typename MeshType, typename MapType, UInt SpaceDim>
template <typename VectorType, typename VectorBlockType >
void StabilizationVMS<MeshType, MapType, SpaceDim>::applyRHS_semi_implicit( VectorType& rhs,
                                                                            const VectorType& velocityExtrapolation,
                                                                            const VectorType& velocityPreviousStep,
				                 			    const VectorType& pressurePreviousStep  )
{
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
    integrate(
                    elements(ETuFESpace->mesh()), M_uFESpace->qr(), ETuFESpace, 
                    dot(DVDT_COMPATIBILITY_N, SUPG_TEST) + dot(DVDT_COMPATIBILITY_N, VMS_TEST)
                    + dot(outerProduct( value( -1.0 ) * DVDT_COMPATIBILITY_N, RES_MOMENTUM_STEP_N ), LES_TEST) // Turbulence model, semi_implicit VMS-LES, Momentum residual
    )
    >> NSRhs.block(0);

    // Consistency term: PSPG
    integrate(
                    elements(ETuFESpace->mesh()), M_pFESpace->qr(), ETpFESpace,
                    dot(DVDT_COMPATIBILITY_N, PSPG_TEST)
    )
    >> NSRhs.block(1);


    VectorBlockType NSRhsUnique( NSRhs, Unique );
    rhs += NSRhsUnique;
    

} // applyRHS_semi_implicit(...)


// RHS, VMS-LES based on RHS residual
template<typename MeshType, typename MapType, UInt SpaceDim>
template <typename VectorType, typename VectorBlockType >
void StabilizationVMS<MeshType, MapType, SpaceDim>::applyRHS_residual( VectorType& rhs,
                                                                       const VectorType& velocityExtrapolation,
                                                                       const VectorType& momentumResidual )
{
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
    integrate(
                    elements(ETuFESpace->mesh()), M_uFESpace->qr(), ETuFESpace, 
                    dot(DVDT_COMPATIBILITY_EX, SUPG_TEST) + dot(DVDT_COMPATIBILITY_EX, VMS_TEST)
                    + dot(outerProduct( RES_MOMENTUM_EXPLICIT, RES_MOMENTUM_EXPLICIT ), LES_TEST) // Turbulence model, explicit
    )
    >> NSRhs.block(0);

    // Consistency term: PSPG
    integrate(
                    elements(ETuFESpace->mesh()), M_pFESpace->qr(), ETpFESpace,
                    dot(DVDT_COMPATIBILITY_EX, PSPG_TEST)
    )
    >> NSRhs.block(1);


    VectorBlockType NSRhsUnique( NSRhs, Unique );
    rhs += NSRhsUnique;

} // applyRHS_residual(...)

template<typename MeshType, typename MapType, UInt SpaceDim>
void StabilizationVMS<MeshType, MapType, SpaceDim>::showMe(std::ostream & output) const
{
    output << "StabilizationVMS::showMe() "      << std::endl;
    output << "Fluid Viscosity (dynamic): "      << M_viscosity << std::endl;
    output << "Fluid Density: "                  << M_density << std::endl;
    output << "Stabilization coefficient tau_M: "<< M_tauM << std::endl;
    output << "Stabilization coefficient tau_C: "<< M_tauC << std::endl;
}

//=============================================================================
// Set Methods
//=============================================================================

template<typename MeshType, typename MapType, UInt SpaceDim>
void StabilizationVMS<MeshType, MapType, SpaceDim>::setVelocitySpace(fespacePtr_Type & uFESpace)
{
#ifdef DEBUG
        // For now we only support P1 finite elements for velocity in 2D/3D
        if ( (uFESpace->refFE().type() != FE_TYPE.FE_P1_2D) && (uFESpace->refFE().type() != FE_TYPE.FE_P1_3D) )
          std::cout << std:endl << "WARNING: Variational multiscale stabilization only tested for P1-P1 discetization." << std:endl;
#endif
        M_uFESpace = uFESpace;
    }

//=============================================================================
// Private Methods
//=============================================================================

template<typename MeshType, typename MapType, UInt SpaceDim>
void StabilizationVMS<MeshType, MapType, SpaceDim>::checkFESpaces()
{
    ASSERT(M_uFESpace != 0 && M_uFESpace != 0 && M_comm != 0, "Check that FE spaces and Epetra communicator have been set");
}



} // namespace LifeV

#endif /* _VMSSTABILIZATION_HPP_ */
