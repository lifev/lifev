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

#define L_NS2F_UMFPACK 0
#define L_NS2F_PETSC 1

#define L_NS2F_LINEAR_SOLVER L_NS2F_PETSC

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

#if L_NS2F_LINEAR_SOLVER == L_NS2F_PETSC
#include <life/lifealg/SolverPETSC.hpp>
#elif L_NS2F_LINEAR_SOLVER == L_NS2F_UMFPACK
#include <life/lifealg/SolverUMFPACK.hpp>
#endif

#include <life/lifesolver/AFSolvers.hpp>


#include <life/lifesolver/LevelSetSolver.hpp>

#include <life/lifearray/boostmatrix.hpp>

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

        typedef BoostMatrix<boost::numeric::ublas::row_major> matrix_type_C;
        typedef BoostMatrix<boost::numeric::ublas::row_major> matrix_type_M;
        typedef BoostMatrix<boost::numeric::ublas::row_major> matrix_type_D;
        typedef BoostMatrix<boost::numeric::ublas::column_major> matrix_type_Dtr;
        typedef DiagonalBoostMatrix matrix_type_M_L;

#if L_NS2F_LINEAR_SOLVER == L_NS2F_PETSC
        typedef SolverPETSC linear_solver_u_type;
        typedef SolverPETSC linear_solver_p_type;
#elif L_NS2F_LINEAR_SOLVER == L_NS2F_UMFPACK
        typedef SolverUMFPACK linear_solver_u_type;
        typedef SolverUMFPACK linear_solver_p_type;
#endif

        typedef Yosida<matrix_type_M_L,
                       matrix_type_C,
                       matrix_type_D,
                       matrix_type_Dtr,
                       linear_solver_u_type,
                       linear_solver_p_type> solver_type;
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
            return _u;
        }

        const Vector& pressure() const {
            return _p;
        }

        const lsfunction_type& lsfunction() {
            return _M_lss.lsfunction();
        }

        const mesh_type& mesh() {
            return _mesh;
        }

        const Dof& lsDof() {
            return _M_lss.dof();
        }
        
        CurrentFE& fe_ls() {
            return _M_lss.fe();
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
                _M_lss.directReinitialization();
        }

        void setVerbose() {
            _M_verbose = true;
        }

        void unsetVerbose() {
            _M_verbose = false;
        }
        //@}

    private:
        //! Data file
        const GetPot& _M_data_file;

        //! Pattern for the p-u block (D)
        CSRPatt _M_pattern_D_block;

        //! Pattern for the u-p block (Dtr)
        CSRPatt _M_pattern_Dtr_block;

        //! Pattern for C
        MSRPatt _M_pattern_C;

        //! Pattern for D
        MixedPattern<1, NDIM, CSRPatt> _M_pattern_D;

        //! Pattern for Dtr
        MixedPattern<NDIM, 1, CSRPatt> _M_pattern_Dtr;

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

        //! Matrix Dtr
        matrix_type_Dtr _M_Dtr;

        //! Elementary matrices
        ElemMat _M_elmat_C;
        ElemMat _M_elmat_M;
        ElemMat _M_elmat_D;
        ElemMat _M_elmat_Dtr;
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
        //!
        Vector _M_constant_pressure;

         //! Reference to the level set function
        const lsfunction_type& _M_lsfunction;

        //! Saddle point system solver
        solver_type _M_solver;

        //! Vector of pressure constants
        Vector _M_pressure_constant;

        //! Linear solver for velocity systems
        linear_solver_u_type _M_solver_u;

        //! Linear solver for pressure systems
        linear_solver_p_type _M_solver_p;

        //! Verbose flag
        bool _M_verbose;

        //! The current time
        Real _M_time;

        //! A flag for steady problems
        bool _M_steady;

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
        _M_pattern_D_block(pDof(), uDof()),
        _M_pattern_Dtr_block(uDof(), pDof()),
        _M_pattern_C(uDof(), NDIM),
        _M_pattern_D(_M_pattern_D_block),
        _M_pattern_Dtr(_M_pattern_Dtr_block),
        _M_M(_M_pattern_C),
        _M_M_L(_dim_u * NDIM),
        _M_M_L_w(_dim_u * NDIM),
        _M_C(_M_pattern_C),
        _M_D(_M_pattern_D),
        _M_Dtr(_M_pattern_Dtr),
        _M_elmat_C(fe_u().nbNode, NDIM, NDIM),
        _M_elmat_M(fe_u().nbNode, NDIM, NDIM),
        _M_elmat_D(fe_p().nbNode, 1, 0, fe_u().nbNode, 0, NDIM),
        _M_elmat_Dtr(fe_u().nbNode, NDIM, 0, fe_p().nbNode, 0, 1),
        _M_elmat_P(fe_u().nbNode, NDIM + 1, NDIM + 1),
        _M_elvec(fe_u().nbNode, NDIM),
        _M_lss(_mesh, data_file, "levelset", refFE_lss, qr, bd_qr, bc_h_lss, _fe_u, _dof_u, _u),
        _M_elvec_lss(_M_lss.fe().nbNode, 1),
        _M_rhs_u( NDIM * _dim_u),
        _M_rhs_p( _dim_p ),
        _M_beta_fct(0),
        _M_constant_pressure( _dim_p ),
        _M_lsfunction(_M_lss.lsfunction()),
        _M_solver(_M_M_L_w, _M_C, _M_D, _M_Dtr, _M_data_file, "navier-stokes/yosida", _M_solver_u, _M_solver_p)
    {
        if(_M_verbose)
            std::cout << "** NS2F ** Using boost matrix" << std::endl;
#if L_NS2F_LINEAR_SOLVER == L_NS2F_PETSC
        std::cout << "** NS2F ** Using PETSC linear solver" << std::endl;
        _M_solver_u.setOptionsFromGetPot(_M_data_file, "navier-stokes/yosida/solver-u");
        _M_solver_p.setOptionsFromGetPot(_M_data_file, "navier-stokes/yosida/solver-p");
        if ( _BCh_u.hasOnlyEssential() ) {
            Real constPress = 1. / sqrt( _dim_p );
            for(UInt i = 0; i < _dim_p; ++i) {
                _M_constant_pressure[ i ] = constPress;
            }
            std::vector<const Vector*> nullSpace(1);
            nullSpace[ 0 ] = &_M_constant_pressure;
            _M_solver_p.setNullSpace(nullSpace);
        }
#elif L_NS2F_LINEAR_SOLVER == L_NS2F_UMFPACK
        std::cout << "** NS2F ** Using UMFPACK linear solver" << std::endl;
#endif
        _M_verbose = _M_data_file("navier-stokes/miscellaneous/verbose", false);
        _M_steady = _M_data_file("navier-stokes/miscellaneous/steady", false);
    }

    template<typename MeshType>
    inline void NSSolver2FluidsMixed<MeshType>::advance_NS(source_type const& source, Real const& time) {
        // Set current time
        _M_time = time;

        // Number of components for the velocity
        UInt nbCompU = _u.nbcomp();

        // Velocity vector for the linearization of the convection term
        Vector betaVec(_u.size());

        if(_M_beta_fct)
            uInterpolate(*_M_beta_fct, betaVec, time);
        else
            if(_M_steady)
                betaVec = _u;
            else
                betaVec = _bdf.bdf_u().extrap();

        // Matrices initialization
        _M_M.zeros();
        _M_C.zeros();
        _M_D.zeros();
        _M_Dtr.zeros();

        // u-RHS vector reinitialization
        _M_rhs_u = ZeroVector( _M_rhs_u.size() );

        // p-RHS vector reinitialization
        _M_rhs_p = ZeroVector( _M_rhs_p.size() );

        // Elementary computation and matrix assembling
        Chrono __chrono;

        // Total time for elementar contribution computation
        Real __cumul1 = 0.;

        // Total time for u-u block assembling
        Real __cumul2 = 0.;

        // Total time for u-p and p-u blocks assembling
        Real __cumul3 = 0.;

        if(_M_verbose)
            std::cout << "** NS2F ** Assembling matrices" << std::endl;

        for(UInt iVol = 1; iVol <= _mesh.numVolumes(); iVol++) {
            __chrono.start();

            fe_u().updateFirstDeriv(_mesh.volumeList(iVol));
            fe_p().updateFirstDeriv(_mesh.volumeList(iVol));
            _M_lss.fe().updateJac(_mesh.volumeList(iVol));

            _M_elmat_C.zero();
            _M_elmat_M.zero();
            _M_elmat_D.zero();
            _M_elmat_Dtr.zero();

            _M_elvec_lss.zero();
            _M_elvec.zero();

            // Extract the vector of local lsfunction dofs
            ElemVec ls_fun_loc();
            UInt element_id = _M_lss.fe().currentId();

            for(int node_id = 0; node_id < _M_lss.fe().nbNode; node_id++) {
                UInt iglo = _M_lss.dof().localToGlobal(element_id, node_id + 1) - 1;
                _M_elvec_lss.vec()[node_id] = _M_lsfunction[iglo];
            }

            element_id = _fe_u.currentId();
            
            // Stiffness strain
            stiff_strain_2f(2. * viscosity(fluid1), 2. * viscosity(fluid2), _M_elvec_lss, _M_lss.fe(), _M_elmat_C, fe_u());


            // Mass
            if(!_M_steady) {
                for(UInt iComp = 0; iComp < nbCompU; iComp++)
                        mass_2f(density(fluid1), density(fluid2), _M_elvec_lss, _M_lss.fe(), _M_elmat_C, fe_u(), iComp, iComp);
            }

            //  Extract the vector of local velocity dofs
            for(UInt node_id = 0; node_id < (UInt)fe_u().nbNode; node_id++) {
                for(UInt comp_id = 0; comp_id < nbCompU; comp_id++) {
                    UInt iglo = uDof().localToGlobal(element_id, node_id + 1) - 1 + comp_id * _dim_u;
                    _M_elvec.vec()[node_id + comp_id * fe_u().nbNode] = betaVec(iglo);
                }
            }

            // Advection
            advection_2f(density(fluid1), density(fluid2), _M_elvec_lss, _M_lss.fe(), _M_elvec, _M_elmat_C, fe_u());

            __chrono.stop();
            __cumul1 += __chrono.diff();

            // Assemble all contributions
            for(UInt iComp = 0; iComp < nbCompU; iComp++) {
                // Assembling u-u block
                __chrono.start();
                for(UInt jComp = 0; jComp < nbCompU; jComp++) {
                    assemb_mat(_M_C, _M_elmat_C, fe_u(), uDof(), iComp, jComp);

                    if(!_M_steady)
                        assemb_mat(_M_M, _M_elmat_M, fe_u(), uDof(), iComp, jComp);
                }
                __chrono.stop(); __cumul2 += __chrono.diff();

                // Off-diagonal blocks
                __chrono.start();
                grad(iComp, 1.0, _M_elmat_Dtr, fe_u(), fe_p(), iComp, 0);
                __chrono.stop();
                __cumul1 += __chrono.diff();

                /**
                   The separate computation of discrete div operator was removed
                   div(iComp, -1.0, _M_elmat_D, fe_p(), fe_u(), 0, iComp);
                */

                // Assembling
                __chrono.start();
                assemb_mat_mixed(_M_Dtr, _M_elmat_Dtr, fe_u(), fe_p(), uDof(), pDof(), iComp, 0);
                __chrono.stop();
                __cumul3 += __chrono.diff();

                /**
                   The separate computation of discrete div operator was removed
                   assemb_mat_mixed(_M_D, _M_elmat_D, fe_p(), fe_u(), pDof(), uDof(), 0, iComp);
                */

                // Source term
                _M_elvec.zero();
                compute_vec_2f(density(fluid1), density(fluid2), _M_elvec_lss, _M_lss.fe(), source, _M_elvec, fe_u(), _M_time, iComp);
                assemb_vec(_M_rhs_u, _M_elvec, fe_u(), uDof(), iComp);
            }

        } // loop over volumes

        // Lump mass matrix
        _M_M_L.lumpRowSum(_M_M);
        _M_M_L_w = ( _bdf.bdf_u().coeff_der(0) / _dt ) *_M_M_L;
        
        // Add unsteady contribution
        _M_C += _M_M_L_w;

        __chrono.stop();

        if(_M_verbose) {
            std::cout << "** NS2F ** Elementary constributions computation : "
                      << __cumul1 << " s" << std::endl;
            std::cout << "** NS2F ** Diagonal block assembling: "
                      << __cumul2 << " s" << std::endl;
            std::cout << "** NS2F ** Off-diagonal block assembling: "
                      << __cumul3 << " s" << std::endl;
        }

        // Add RHS terms stemming from the time derivative
        if(!_M_steady)
            _M_rhs_u += prod( _M_M_L, _bdf.bdf_u().time_der(_dt) );

        // Apply boundary conditions
        if(_M_verbose)
            std::cout << "** NS2F ** Applying boundary conditions"
                      << std::endl;

        if (!_BCh_u.bdUpdateDone())
            _BCh_u.bdUpdate(_mesh, _feBd_u, uDof());

        bcManage(_M_C, _M_Dtr, _M_rhs_u, _mesh, uDof(), _BCh_u, _feBd_u, 1.0, _M_time );

        if(_M_verbose)
            std::cout << "** NS2F ** Computing discrete divergence operator"
                      << std::endl;
        _M_D = trans( _M_Dtr );

        // Export matrices to Matlab format for debugging purposes
//         if(_M_verbose) {
//             std::cout << "** NS2F ** Exporting matrices to Matlab format"
//                       << std::endl;
//             _M_C.spy( "./results/spyC" );
//             _M_D.spy( "./results/spyD" );
//             _M_Dtr.spy( "./results/spyDtr" );
//         }

        // Solve the system
        if(_M_verbose)
            std::cout << "** NS2F ** Solving the system" << std::endl;
        _M_solver.solve(_u, _p, _M_rhs_u, _M_rhs_p);

        // Update bdf
        _bdf.bdf_u().shift_right(_u);
    }

    template<typename MeshType>
    inline void NSSolver2FluidsMixed<MeshType>::advance_LS() {
        _M_lss.setVelocity(_u);
        _M_lss.timeAdvance();
    }

    template<typename MeshType>
    void NSSolver2FluidsMixed<MeshType>::initialize(const function_type& u0,
                                                    const function_type& p0,
                                                    const function_type& lsfunction0,
                                                    Real t0, Real delta_t)
    {
        UInt nbCompU = _u.nbcomp();
        _bdf.bdf_u().initialize_unk(u0, _mesh, refFEu(), fe_u(), uDof(), t0, delta_t, nbCompU);
        _bdf.bdf_p().initialize_unk(p0, _mesh, refFEp(), fe_p(), pDof(), t0, delta_t, 1);

        // Initialize _u and _p
        _u = *(_bdf.bdf_u().unk().begin());
        _p = *(_bdf.bdf_p().unk().begin());

        // Initialize level set solver
        _M_lss.initialize(lsfunction0, t0, delta_t);

        // Set the velocity field
        _M_lss.setVelocity(_u);
    }

    template<typename MeshType>
    void NSSolver2FluidsMixed<MeshType>::timeAdvance(source_type const& source, Real const& time) {
        if(_M_verbose)
            std::cout << "** NS2F ** Advancing Navier-Stokes solver" << std::endl;
        advance_NS(source, time);
        if(_M_verbose)
            std::cout << "** NS2F ** Advancing Level Set solver" << std::endl;
        advance_LS();
    }
}
#endif


