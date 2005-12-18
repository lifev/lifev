/* -*- mode: c++ -*-

This file is part of the LifeV library

Author(s): Daniele Antonio Di Pietro <dipietro@unibg.it>
Date: 2-1-2005

Copyright (C) 2005 Università degli Studi di Bergamo

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file NSSolver2FluidsMixed
   \author Daniele A. Di Pietro <dipietro@unibg.it>
   \date 2-1-2005
*/

#ifndef _NSSOLVER2FLUIDSMIXED_HPP_
#define _NSSOLVER2FLUIDSMIXED_HPP_

#include <lifeconfig.h>

#define L_NS2F_UMFPACK 0
#define L_NS2F_PETSC 1

#define L_NS2F_LINEAR_SOLVER_U L_NS2F_PETSC
#define L_NS2F_LINEAR_SOLVER_P L_NS2F_PETSC

#if defined(HAVE_PETSC_H)
#include <life/lifealg/SolverPETSC.hpp>
#endif /* HAVE_PETSC_H */

#if defined(HAVE_UMPFACK_H)
#include <life/lifealg/SolverUMFPACK.hpp>
#endif /* HAVE_PETSC_H */

#if ( ! defined(HAVE_PETSC_H) && ( ( L_NS2F_LINEAR_SOLVER_U == L_NS2F_PETSC ) || ( L_NS2F_LINEAR_SOLVER_P == L_NS2F_PETSC ) ) ) || ( ! defined(HAVE_UMFPACK_H) && ( ( L_NS2F_LINEAR_SOLVER_U == L_NS2F_UMFPACK ) || ( L_NS2F_LINEAR_SOLVER_P == L_NS2F_UMFPACK ) ) )
#define L_NS2F_MISSING_SOLVER
#endif

#define L_STOKES 0
#define L_NAVIER_STOKES 1

#define L_NS2F_PROBLEM L_NAVIER_STOKES

#define L_PRESSURE_MATRIX 0
#define L_YOSIDA 1

#define L_NS2F_SOLVER L_YOSIDA

#define L_DEBUG_MODE 0

#include <fstream>

#include <life/lifecore/chrono.hpp>

#include <life/lifesolver/NavierStokesHandler.hpp>
#include <life/lifesolver/dataNS2Fluids.hpp>

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifearray/pattern.hpp>

#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/elemOper2Fluids.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/assemb.hpp>

#include <life/lifesolver/AFSolvers.hpp>
#include <life/lifealg/PressureMatrixSolver.hpp>

#include <life/lifesolver/LevelSetSolver.hpp>

#include <life/lifearray/boostmatrix.hpp>

namespace ublas = boost::numeric::ublas;

namespace LifeV {
    /*!
      \class NSSolver2FluidsMixed
      \brief 2 fluid Navier-Stokes solver

      \c A 2 fluid Navier-Stokes solver

      @author Daniele Antonio Di Pietro <dipietro@unibg.it>
      @see
    */
    template<typename MeshType>
    class NSSolver2FluidsMixed
        :
        public NavierStokesHandler<MeshType, DataNS2Fluids<MeshType> >
    {
    public:
        /** @name Typedefs
         */
        //@{
        typedef MeshType mesh_type;
        typedef DataNS2Fluids<MeshType> data_type;

        typedef LevelSetSolver<MeshType> lss_type;
        typedef typename lss_type::lsfunction_type lsfunction_type;

        typedef typename NavierStokesHandler<mesh_type, data_type>::Function function_type;
        typedef typename NavierStokesHandler<mesh_type, data_type>::source_type source_type;

        typedef BoostMatrix<ublas::row_major> matrix_type_C;
        typedef BoostMatrix<ublas::row_major> matrix_type_M;
        typedef BoostMatrix<ublas::row_major> matrix_type_D;
        typedef BoostMatrix<ublas::column_major> matrix_type_G;
        typedef DiagonalBoostMatrix matrix_type_M_L;

#if ! defined(L_NS2F_MISSING_SOLVER)

#if L_NS2F_LINEAR_SOLVER_U == L_NS2F_PETSC
        typedef SolverPETSC linear_solver_u_type;
#elif L_NS2F_LINEAR_SOLVER_U == L_NS2F_UMFPACK
        typedef SolverUMFPACK linear_solver_u_type;
#endif

#if L_NS2F_LINEAR_SOLVER_P == L_NS2F_PETSC
        typedef SolverPETSC linear_solver_p_type;
#elif L_NS2F_LINEAR_SOLVER_P == L_NS2F_UMFPACK
        typedef SolverUMFPACK linear_solver_p_type;
#endif

#if L_NS2F_SOLVER == L_YOSIDA
        typedef Yosida<matrix_type_M_L,
                       matrix_type_C,
                       matrix_type_D,
                       matrix_type_G,
                       linear_solver_u_type,
                       linear_solver_p_type> solver_type;
#elif L_NS2F_SOLVER == L_PRESSURE_MATRIX
        typedef PressureMatrixSolver<matrix_type_C,
                                     matrix_type_D,
                                     matrix_type_G,
                                     linear_solver_u_type> solver_type;
#endif

#endif /* MISSING_SOLVER */
        //@}

        /** @name Constructors
         */
        //@{
        NSSolver2FluidsMixed(const GetPot& data_file,
                             const RefFE& refFE_u,
                             const RefFE& refFE_p,
                             BCHandler& bc_h_u,
                             const RefFE& refFE_lss,
                             BCHandler& bc_h_lss,
                             const QuadRule& qr,
                             const QuadRule& bd_qr);
        //@}

        /** @name Accessors
         */
        //@{

        const Vector& velocity() const {
            return this->_u;
        }

        const Vector& pressure() const {
            return this->_p;
        }

        const lsfunction_type& lsfunction() {
            return this->_M_lss.lsfunction();
        }

        const mesh_type& mesh() {
            return this->_mesh;
        }

        const Dof& lsDof() {
            return this->_M_lss.dof();
        }

        CurrentFE& fe_ls() {
            return this->_M_lss.fe();
        }
        //@}

        /** @name Methods
         */
        //@{

        /*!
          \Initialize the solver
        */
        void initialize(const function_type&, const function_type&, const function_type&, Real, Real);

        /*!
          \Time advance re-computes the whole problem matrix at every time step.
          \This is necessary because (discontinuous) coefficient are
          \time-dependent being so the location of the interface.
        */
        void timeAdvance(source_type const& source, Real const& time);

        /**
           \Function iterate needs a fake implementation because it is defined
           \as a virtual member of NavierStokesHandler
        */

        void iterate(const Real& /* time */) {}

        /**
           \Call the proper reinitialization function for the level set solver
        */

        void reinitialize(const std::string& method) {
            if(method == "direct")
                this->_M_lss.directReinitialization();
        }

        void setVerbose() {
            this->_M_verbose = true;
        }

        void unsetVerbose() {
            this->_M_verbose = false;
        }
        //@}

    private:
        //! Data file
        const GetPot& _M_data_file;

        //! Pattern for the p-u block (D)
        CSRPatt _M_pattern_D_block;

        //! Pattern for the u-p block (G)
        CSRPatt _M_pattern_G_block;

        //! Pattern for C
        MSRPatt _M_pattern_C;

        //! Pattern for D
        MixedPattern<1, NDIM, CSRPatt> _M_pattern_D;

        //! Pattern for G
        MixedPattern<NDIM, 1, CSRPatt> _M_pattern_G;

        //! Velocity mass matrix
        matrix_type_M _M_M;

        //! Lumped velocity mass matrix
        matrix_type_M_L _M_M_L;

        //! Weighted lumped velocity mass matrix for Yosida solver
        matrix_type_M_L _M_M_L_w;

        //! Matrix C
        matrix_type_C _M_C;

        //! Matrix D
        matrix_type_D _M_D;

        //! Matrix G
        matrix_type_G _M_G;

        //! Elementary matrices
        ElemMat _M_elmat_C;
        ElemMat _M_elmat_M;
        ElemMat _M_elmat_D;
        ElemMat _M_elmat_G;
        ElemMat _M_elmat_P;

        //! Elementary vector
        ElemVec _M_elvec;

        //! The level set solver
        lss_type _M_lss;

        //! Elementary vector for level set solver
        ElemVec _M_elvec_lss;

        //! RHS for the velocity
        Vector _M_rhs_u;

        //! RHS for the pressure
        Vector _M_rhs_p;

        //! Analytical expression for the velocity field to use in
        //! linearized advection term
        const function_type* _M_beta_fct;

        //! Vector of constant pressure
        Vector _M_constant_pressure;

        //! Reference to the level set function
        const lsfunction_type& _M_lsfunction;

#if ! defined(L_NS2F_MISSING_SOLVER)

        //! Saddle point system solver
        solver_type _M_solver;

        //! Linear solver for velocity systems
        linear_solver_u_type _M_solver_u;

        //! Linear solver for pressure systems
        linear_solver_p_type _M_solver_p;

#endif /* MISSING_SOLVER */

        //! Verbose flag
        bool _M_verbose;

        //! The current time
        Real _M_time;

        /**
           \Advance the Navier-Stokes equations one step in time
        */

        void advance_NS(source_type const& source, Real const& time);

        /**
           \Advance the level set equation  one step in time
        */

        inline void advance_LS();
    };

    // Implementations

    template<typename MeshType> NSSolver2FluidsMixed<MeshType>::
    NSSolver2FluidsMixed(const GetPot& data_file,
                         const RefFE& refFE_u,
                         const RefFE& refFE_p,
                         BCHandler& bc_h_u,
                         const RefFE& refFE_lss,
                         BCHandler& bc_h_lss,
                         const QuadRule& qr,
                         const QuadRule& bd_qr)
        :
        NavierStokesHandler<mesh_type, data_type>(data_file, refFE_u, refFE_p, qr, bd_qr, qr, bd_qr, bc_h_u),
        _M_data_file(data_file),
        _M_pattern_D_block(this->pDof(), this->uDof()),
        _M_pattern_G_block(this->uDof(), this->pDof()),
        _M_pattern_C(this->uDof(), NDIM),
        _M_pattern_D(_M_pattern_D_block),
        _M_pattern_G(_M_pattern_G_block),
        _M_M(_M_pattern_C),
        _M_M_L(this->_dim_u * NDIM),
        _M_M_L_w(this->_dim_u * NDIM),
        _M_C(_M_pattern_C),
        _M_D(_M_pattern_D),
        _M_G(_M_pattern_G),
        _M_elmat_C(this->fe_u().nbNode, NDIM, NDIM),
        _M_elmat_M(this->fe_u().nbNode, NDIM, NDIM),
        _M_elmat_D(this->fe_p().nbNode, 1, 0, this->fe_u().nbNode, 0, NDIM),
        _M_elmat_G(this->fe_u().nbNode, NDIM, 0, this->fe_p().nbNode, 0, 1),
        _M_elmat_P(this->fe_u().nbNode, NDIM + 1, NDIM + 1),
        _M_elvec(this->fe_u().nbNode, NDIM),
        _M_lss(this->_mesh, data_file, "levelset", refFE_lss, qr, bd_qr, bc_h_lss, this->_fe_u, this->_dof_u, this->_u),
        _M_elvec_lss(_M_lss.fe().nbNode, 1),
        _M_rhs_u( NDIM * this->_dim_u),
        _M_rhs_p( this->_dim_p ),
        _M_beta_fct(0),
        _M_constant_pressure( this->_dim_p ),
        _M_lsfunction(_M_lss.lsfunction())
#if ! defined(L_NS2F_MISSING_SOLVER)
#if L_NS2F_SOLVER == L_YOSIDA
        ,_M_solver(_M_M_L_w, _M_C, _M_D, _M_G, _M_data_file, "navier-stokes/yosida", _M_solver_u, _M_solver_p)
#elif L_NS2F_SOLVER == L_PRESSURE_MATRIX
        ,_M_solver(_M_C, _M_D, _M_G, _M_data_file, "navier-stokes/pmm")
#endif
#endif /* MISSING_SOLVER */
    {
        if(_M_verbose)
            std::cout << "[NSSolver2FluidsMixed::constructor] Using boost matrix" << std::endl;
#if ! defined(L_NS2F_MISSING_SOLVER)
#if L_NS2F_LINEAR_SOLVER_U == L_NS2F_PETSC
        std::cout << "[NSSolver2FluidsMixed::constructor] Using PETSC linear solver for u" << std::endl;
        _M_solver_u.setOptionsFromGetPot(_M_data_file, "navier-stokes/yosida/solver-u");
#elif L_NS2F_LINEAR_SOLVER_U == L_NS2F_UMFPACK
        std::cout << "[NSSolver2FluidsMixed::constructor] Using UMFPACK linear solver for u" << std::endl;
#endif

#if L_NS2F_SOLVER == L_YOSIDA
        std::cout << "[NSSolver2FluidsMixed::constructor] Using Yosida method to advance in time" << std::endl;
#elif L_NS2F_SOLVER == L_PRESSURE_MATRIX
        std::cout << "[NSSolver2FluidsMixed::constructor] Using PMM method to advance in time" << std::endl;
#endif

#if L_NS2F_LINEAR_SOLVER_P == L_NS2F_PETSC
        std::cout << "[NSSolver2FluidsMixed::constructor] Using PETSC linear solver for p" << std::endl;
        _M_solver_p.setOptionsFromGetPot(_M_data_file, "navier-stokes/yosida/solver-p");
        if ( this->bcHandler().hasOnlyEssential() ) {
            Real constPress = 1. / sqrt( this->_dim_p );
            for(UInt i = 0; i < this->_dim_p; ++i) {
                _M_constant_pressure[ i ] = constPress;
            }
            std::vector<const Vector*> nullSpace(1);
            nullSpace[ 0 ] = &_M_constant_pressure;
            _M_solver_p.setNullSpace(nullSpace);
        }
#elif L_NS2F_LINEAR_SOLVER_P == L_NS2F_UMFPACK
        std::cout << "[NSSolver2FluidsMixed::constructor] Using UMFPACK linear solver for p" << std::endl;
#endif
#else
        std::cerr << "WARNING: Missing linear solver - no systems are solved!"
                  << std::endl;
#endif /* MISSING_SOLVER */
        _M_verbose = _M_data_file("navier-stokes/miscellaneous/verbose", false);
    }

    template<typename MeshType>
    inline void NSSolver2FluidsMixed<MeshType>::advance_NS(source_type const& source, Real const& time) {
        // Set current time
        _M_time = time;

        // Number of components for the velocity
        UInt nbCompU = this->_u.nbcomp();

        // Velocity vector for the linearization of the convection term
        Vector betaVec(this->_u.size());

        if(_M_beta_fct)
            uInterpolate(*_M_beta_fct, betaVec, time);
        else
            betaVec = this->_bdf.bdf_u().extrap();

        // LSH and RHS initialization
        if(_M_verbose)
            std::cout << "[NSSolver2FluidsMixed::advance_NS] Initializing matrices and vectors" << std::endl;
        _M_M.zeros();
        _M_C.zeros();
        _M_D.zeros();
        _M_G.zeros();

        _M_rhs_u = ZeroVector( _M_rhs_u.size() );
        _M_rhs_p = ZeroVector( _M_rhs_p.size() );

        // Timing stuff
        Chrono __chrono;

        // Total time for elementar contribution computation
        Real __cumul1 = 0.;

        // Total time for u-u block assembling
        Real __cumul2 = 0.;

        // Total time for u-p and p-u blocks assembling
        Real __cumul3 = 0.;

        if(_M_verbose)
            std::cout << "[NSSolver2FluidsMixed::advance_NS] Assembling matrices" << std::endl;

        for(UInt iVol = 1; iVol <= this->_mesh.numVolumes(); iVol++) {
            __chrono.start();

            // Do proper element updates
            this->fe_u().updateFirstDeriv( this->_mesh.volumeList( iVol ) );
            this->fe_p().updateFirstDeriv( this->_mesh.volumeList( iVol ) );
            _M_lss.fe().updateJac( this->_mesh.volumeList( iVol ) );

            _M_elmat_C.zero();
            _M_elmat_M.zero();
            _M_elmat_D.zero();
            _M_elmat_G.zero();

            _M_elvec_lss.zero();
            _M_elvec.zero();

            // Extract the vector of local lsfunction dofs
            UInt element_id = _M_lss.fe().currentId();

            for(int node_id = 0; node_id < _M_lss.fe().nbNode; node_id++) {
                UInt iglo = _M_lss.dof().localToGlobal(element_id, node_id + 1) - 1;
                _M_elvec_lss.vec()[ node_id ] = _M_lsfunction[ iglo ];
            }

            // Stiffness strain
            stiff_strain_2f(2. * this->viscosity(this->fluid1), 2. * this->viscosity(this->fluid2), _M_elvec_lss, _M_lss.fe(), _M_elmat_C, this->fe_u());

            // Mass
            Real bdf_coeff = this->_bdf.bdf_u().coeff_der( 0 ) / this->_dt;
            for(UInt iComp = 0; iComp < nbCompU; iComp++) {
                lumped_mass_2f(this->density(this->fluid1), this->density(this->fluid2), _M_elvec_lss, _M_lss.fe(), _M_elmat_M, this->fe_u(), iComp, iComp);
                lumped_mass_2f(bdf_coeff * this->density(this->fluid1), bdf_coeff * this->density(this->fluid2), _M_elvec_lss, _M_lss.fe(),
                               _M_elmat_C, this->fe_u(), iComp, iComp);
            }

#if L_NS2F_PROBLEM == L_NAVIER_STOKES
            // Compute local contributions
            element_id = this->_fe_u.currentId();

            //  Extract the vector of local velocity dofs
            for(UInt node_id = 0; node_id < (UInt)this->fe_u().nbNode; node_id++) {
                for(UInt comp_id = 0; comp_id < nbCompU; comp_id++) {
                    UInt iglo = this->uDof().localToGlobal(element_id, node_id + 1) - 1 + comp_id * this->_dim_u;
                    _M_elvec.vec()[node_id + comp_id * this->fe_u().nbNode] = betaVec(iglo);
                }
            }

            // Advection
            advection_2f(this->density(this->fluid1), this->density(this->fluid2), _M_elvec_lss, _M_lss.fe(), _M_elvec, _M_elmat_C, this->fe_u());
#endif
            __chrono.stop();
            __cumul1 += __chrono.diff();

            // Assemble all contributions
            for(UInt iComp = 0; iComp < nbCompU; iComp++) {
                // Assemble diagonal block
                __chrono.start();
                for(UInt jComp = 0; jComp < nbCompU; jComp++) {
                    assemb_mat(_M_C, _M_elmat_C, this->fe_u(), this->uDof(), iComp, jComp);
                    assemb_mat(_M_M, _M_elmat_M, this->fe_u(), this->uDof(), iComp, jComp);
                }
                __chrono.stop(); __cumul2 += __chrono.diff();

                // Compute off-diagonal blocks
                __chrono.start();
                grad(iComp, 1.0, _M_elmat_G, this->fe_u(), this->fe_p(), iComp, 0);
                __chrono.stop();
                __cumul1 += __chrono.diff();

                // Assemble off-diagonal blocks
                __chrono.start();
                assemb_mat_mixed(_M_G, _M_elmat_G, this->fe_u(), this->fe_p(), this->uDof(), this->pDof(), iComp, 0);
                __chrono.stop();
                __cumul3 += __chrono.diff();

                // Source term
                _M_elvec.zero();
                compute_vec_2f(this->density(this->fluid1), this->density(this->fluid2), _M_elvec_lss, _M_lss.fe(), source, _M_elvec, this->fe_u(), _M_time, iComp);
                assemb_vec(_M_rhs_u, _M_elvec, this->fe_u(), this->uDof(), iComp);
            }

        } // loop over volumes

        if(_M_verbose) {
            std::cout << "[NSSolver2FluidsMixed::advance_NS] Elementary contributions computation : "
                      << __cumul1 << " s" << std::endl;
            std::cout << "[NSSolver2FluidsMixed::advance_NS] Diagonal block assembling: "
                      << __cumul2 << " s" << std::endl;
            std::cout << "[NSSolver2FluidsMixed::advance_NS] Off-diagonal block assembling: "
                      << __cumul3 << " s" << std::endl;
        }

        // Lump mass matrix
        _M_M_L.lumpRowSum(_M_M);
        _M_M_L_w = ( this->_bdf.bdf_u().coeff_der(0) / this->_dt ) *_M_M_L;

        // Add RHS terms stemming from the time derivative
        _M_rhs_u += prod( _M_M_L, this->_bdf.bdf_u().time_der(this->_dt) );

#if L_DEBUG_MODE
        _M_C.spy( "./results/spyCbeforebc" );
        _M_G.spy( "./results/spyGbeforebc" );
#endif

        // Apply boundary conditions
        if(_M_verbose)
            std::cout << "[NSSolver2FluidsMixed::advance_NS] Applying boundary conditions"
                      << std::endl;

        if (!this->bcHandler().bdUpdateDone())
            this->bcHandler().bdUpdate(this->_mesh, this->_feBd_u, this->uDof());

        bcManage(_M_C, _M_G, _M_rhs_u, this->_mesh, this->uDof(), this->bcHandler(), this->_feBd_u, 1.0, _M_time );

        if(_M_verbose)
            std::cout << "[NSSolver2FluidsMixed::advance_NS] Computing discrete divergence operator"
                      << std::endl;
        _M_D = trans( _M_G );

#if L_DEBUG_MODE
        // Export matrices to Matlab format for debugging purposes
        if(_M_verbose)
            std::cout << "[NSSolver2FluidsMixed::advance_NS] Exporting matrices to Matlab format"
                      << std::endl;
        _M_C.spy( "./results/spyC" );
        _M_D.spy( "./results/spyD" );
        _M_G.spy( "./results/spyG" );
        _M_M.spy( "./results/spyM" );
        _M_M_L.spy( "./results/spyML" );
        _M_M_L_w.spy( "./results/spyMLw" );

        std::cout << "[NSSolver2FluidsMixed::advance_NS] Exporting rhs vectors to Matlab format"
                  << std::endl;
        std::ofstream __ofile("./results/rhsu.m");
        __ofile << "bU = [..." << std::endl;
        for(UInt i = 0; i < NDIM * this->_dim_u - 1; i++)
            __ofile << _M_rhs_u( i ) << ";..." << std::endl;
        __ofile << _M_rhs_u( this->_dim_u - 1 ) << "];" << std::endl;
#endif
        // Solve the system
        if(_M_verbose)
            std::cout << "[NSSolver2FluidsMixed::advance_NS] Solving the system" << std::endl;
#if ! defined(L_NS2F_MISSING_SOLVER)
        _M_solver.solve( this->_u, this->_p, _M_rhs_u, _M_rhs_p);
#endif

#if L_NS2F_SOLVER == L_PRESSURE_MATRIX
        std::cout << _M_solver << std::endl;
#endif

        // Update bdf
        this->_bdf.bdf_u().shift_right(this->_u);
    }

    template<typename MeshType>
    inline void NSSolver2FluidsMixed<MeshType>::advance_LS() {
        _M_lss.setVelocity(this->_u);
        _M_lss.timeAdvance();
    }

    template<typename MeshType>
    void NSSolver2FluidsMixed<MeshType>::initialize(const function_type& u0,
                                                    const function_type& p0,
                                                    const function_type& lsfunction0,
                                                    Real t0, Real delta_t)
    {
        UInt nbCompU = this->_u.nbcomp();
        this->_bdf.bdf_u().initialize_unk(u0, this->_mesh, this->refFEu(), this->fe_u(), this->uDof(), t0, delta_t, nbCompU);
        this->_bdf.bdf_p().initialize_unk(p0, this->_mesh, this->refFEp(), this->fe_p(), this->pDof(), t0, delta_t, 1);

        // Initialize _u and _p
        this->_u = *(this->_bdf.bdf_u().unk().begin());
        this->_p = *(this->_bdf.bdf_p().unk().begin());

        // Initialize level set solver
        _M_lss.initialize(lsfunction0, t0, delta_t);

        // Set the velocity field
        _M_lss.setVelocity(this->_u);
    }

    template<typename MeshType>
    void NSSolver2FluidsMixed<MeshType>::timeAdvance(source_type const& source, Real const& time) {
        if(_M_verbose)
            std::cout << "[NSSolver2FluidsMixed::timeAdvance] Advancing Navier-Stokes solver" << std::endl;
        advance_NS(source, time);
        if(_M_verbose)
            std::cout << "[NSSolver2FluidsMixed::timeAdvance] Advancing Level Set solver" << std::endl;
        advance_LS();
    }
}
#endif


