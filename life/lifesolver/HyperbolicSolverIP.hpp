/* -*- mode: c++ -*-

This file is part of the LifeV library

Author(s): Daniele Antonio Di Pietro <dipietro@unibg.it>
Date: 26-1-2005

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
   \file HyperbolicSolverIP.hpp
   \author Daniele Antonio Di Pietro <dipietro@unibg.it>
   \date 1-26-2005
*/

#ifndef _HYPERBOLICSOLVERIP_HPP_
#define _HYPERBOLICSOLVERIP_HPP_

#include <utility>
#include <algorithm>
#include <functional>
#include <vector>
#include <list>
#include <set>

#include <tab.hpp>

#include <dataMesh.hpp>

#include <GetPot.hpp>
#include <bdf.hpp>

#ifdef USE_AZTEC_SOLVER
#include <SolverAztec.hpp>
#else
#include <lifeconfig.h>
#if defined( HAVE_PETSC_H )
#include <SolverPETSC.hpp>
#endif
#endif


#include <bcHandler.hpp>
#include <dof.hpp>
#include <elemOper.hpp>
#include <elemOper2Fluids.hpp>

namespace LifeV {
    /**
       \A type for points and vectors. Points are defined by their position 
       \vectors.
    */
    typedef boost::numeric::ublas::bounded_vector<Real, 3> geo_point_type;


    /*!
      \class HyperbolicSolverIP
      \brief Level set solver class

      \c Hyperbolic solver with IP stabilization

      @author Daniele Antonio Di Pietro <dipietro@unibg.it>
      @see
    */

    template<typename MeshType>
    class HyperbolicSolverIP {
    public:
        /** @name Typedefs
         */
        //@{
        typedef geo_point_type vector_type;

        typedef MeshType mesh_type;
#ifdef USE_AZTEC_SOLVER
        typedef SolverAztec solver_type;
#else
        //typedef SolverPETSC solver_type;
#endif


        typedef MSRPatt pattern_type;

        typedef MSRMatr<Real> matrix_type;

        typedef Vector u_type;
        typedef Vector velocity_type;
        typedef Vector rhs_type;

        typedef Funct function_type;

        /** @name Constructors
         */
        //@{
        HyperbolicSolverIP(mesh_type& mesh,
                           const GetPot& data_file,
                           const std::string& data_section,
                           const RefFE& reffe, 
                           const QuadRule& qr, 
                           const QuadRule& qr_bd, 
                           const BCHandler& bc_h,
                           CurrentFE& fe_velocity,
                           const Dof& dof_velocity,
                           velocity_type& velocity0) 
            :
            _M_mesh(mesh),
            _M_data_file(data_file),
            _M_data_section(data_section),
            _M_reffe(reffe),
            _M_reffe_bd( _M_reffe.boundaryFE() ),
            _M_qr(qr),
            _M_qr_bd(qr_bd),
            _M_fe(_M_reffe, getGeoMap(_M_mesh), _M_qr),
            _M_fe_bd(_M_reffe_bd, getGeoMap(_M_mesh).boundaryMap(), _M_qr_bd),
            _M_fe_velocity(fe_velocity),
            _M_dof(_M_mesh, _M_reffe),
            _M_M_pattern(_M_dof, 1),
            _M_M(_M_M_pattern),
            _M_A_pattern(_M_dof, _M_mesh, 1),
            _M_A_steady(_M_A_pattern),
            _M_A(_M_A_pattern),
            _M_dim( _M_dof.numTotalDof() ),
            _M_u(_M_dim),
            _M_b(_M_dim),
            _M_bc_h(bc_h),
            _M_dof_velocity(dof_velocity),
            _M_dim_velocity( _M_dof_velocity.numTotalDof() ),
            _M_velocity(velocity0),
            _M_bdf_order( data_file((data_section + "/bdf/order").data(), 2) ),
            _M_bdf(_M_bdf_order)
        {
            _M_solver.setOptionsFromGetPot(_M_data_file, (_M_data_section + "/solver").data());
            _M_gamma = _M_data_file((_M_data_section + "/ipstab/gamma").data(), 0.125);
        }
        //@}

        /** @name Accessors
         */
        //@{
        /**
           \Return current time
        */

        inline Real currentTime() {
            return _M_t;
        }

        const Dof& dof(){
            return _M_dof;
        }

        /**
           \Return time at next step
        */

        inline Real nextStepTime() {
            return _M_t + _M_delta_t;
        }
        //@}

        /** @name Mutators
         */
        //@{
        /**
           \Set advection field
        */
        void setVelocity(velocity_type& velocity) const {
            _M_velocity = velocity;
        }
        //@}
        
        /** @name Methods
         */
        //@{

        /**
           \Initialize the solver
        */
        void initialize(const function_type& lsfunction0, Real t0, Real delta_t);

        /**
           \Update the right hand side for time advancement
        */
        void timeAdvance();

        /**
           \Return the level set function
        */
        u_type const & u() const {
            return _M_u;
        }

        /**
           \Update left hand side and solve the system
        */
        void iterate();
       
    protected:
        //! Mesh
        mesh_type& _M_mesh;

        //! Data file
        const GetPot& _M_data_file;

        //! Data section
        const std::string _M_data_section;

        //! The linear solver
        //#if USE_AZTEC_SOLVER
        SolverAztec _M_solver;
        //#else
        //SolverPETSC _M_solver;
        //#endif

        //! Reference FE
        const RefFE& _M_reffe;

        //! Reference boundary FE
        const RefFE& _M_reffe_bd;

        //! Quadrature rule for volume integrals
        const QuadRule& _M_qr;

        //! Quadrature rule for boundary integrals
        const QuadRule& _M_qr_bd;

        //! Current finite element for the unknown
        CurrentFE _M_fe;

        //! Current boundary element for the unknown
        CurrentBdFE _M_fe_bd;

        //! Current finite element for the velocity
        CurrentFE _M_fe_velocity;

        //! Degrees of freedom for the unknown
        Dof _M_dof;

        //! Pattern for the mass matrix
        pattern_type _M_M_pattern;

        //! The mass matrix
        MSRMatr<Real> _M_M;

        //! Pattern for the problem matrix
        pattern_type _M_A_pattern;

        //! The steady part of problem matrix
        matrix_type _M_A_steady;

        //! The full problem matrix (steady + unsteady)
        matrix_type _M_A;

        //!The number of total dofs for the unknown
        UInt _M_dim;

        //! The unknown
        u_type _M_u;

        //! Right hand side vector
        rhs_type _M_b;

        //! The boundary condition handler
        const BCHandler& _M_bc_h;

        //! Degrees of freedom for the advection field
        const Dof& _M_dof_velocity;

        //! Number of dofs for the velocity
        UInt _M_dim_velocity;

        //! The advection field
        velocity_type& _M_velocity;

        //! Order of BDF time discretization scheme
        UInt _M_bdf_order;

        //! The BDF time advance method
        Bdf _M_bdf;

        //! The initial time
        Real _M_t0;

        //! The time step
        Real _M_delta_t;

        //! The current time
        Real _M_t;

        //! Verbose flag
        bool _M_verbose;

        //! The vector of stabilization parameters
        std::vector<Real> _M_h2pK;

        //! The coefficient for stabilization parameters scaling
        Real _M_gamma;

        /** @name Methods
         */
        //@{

        /**
           \Compute the stabilization parameters
        */

        void compute_stabilization_parameters() {
            _M_h2pK.resize( _M_mesh.numVolumes() );

            Real measure;

            for(UInt i = _M_mesh.numBFaces() + 1; i <= _M_mesh.numFaces(); i++){
                _M_fe_bd.updateMeasQuadPt( _M_mesh.faceList(i) );
                measure = _M_fe_bd.measure();

                _M_h2pK[_M_mesh.faceList(i).ad_first() - 1] += measure * 0.5;
                _M_h2pK[_M_mesh.faceList(i).ad_second() - 1] += measure * 0.5;
            }
        }

        /**
           \Compute the mass matrix and the constant part of problem matrix
        */

        void compute_M_A_steady();

        /**
           \Evaluate velocity field on qn_id-th quadrature node of fe_id-th
           \element. It is assumed that the same mesh is used for both the
           \velocity field and the unknown
        */

        inline void evaluate_velocity(UInt fe_id, ElemVec& elvec) {
            UInt dim = _M_dof_velocity.numTotalDof();
            for(UInt j = 0; j < (UInt)_M_fe_velocity.nbNode; j++) {
                for(UInt i_comp = 0; i_comp < NDIM; i_comp++) {
                    UInt jglo = _M_dof_velocity.localToGlobal(fe_id, j + 1) - 1 + i_comp * dim;
                    elvec.vec()[j + i_comp * _M_fe_velocity.nbNode] = _M_velocity[jglo];
                }
            }
        }

        /**
           \Compute the advection part of problem matrix, which is 
           \time-dependent in non-steady multi-fluid problems
        */

        void add_A_unsteady();

        /**
           \Apply boundary conditions. Only Dirichelet (essential) boundary 
           \conditions are considered on the sole inflow boundary nodes. At
           \present time only bc given in functional form are considered.
        */

        void apply_bc();
        //@}
    };

    template<typename MeshType>
    void HyperbolicSolverIP<MeshType>::compute_M_A_steady() {
        ElemMat elmat(_M_fe.nbNode, 1, 1);
        ElemVec elvec(_M_fe.nbNode, NDIM);

        Real coeff = _M_bdf.coeff_der(0) / _M_delta_t;

        for(UInt i = 1; i <= _M_mesh.numVolumes(); i++){
            _M_fe.updateJac( _M_mesh.volumeList(i) );

            elmat.zero();
            mass(coeff, elmat, _M_fe, 0, 0);
            assemb_mat(_M_A_steady, elmat, _M_fe, _M_dof, 0, 0);

            elmat.zero();
            mass(1., elmat, _M_fe, 0, 0);
            assemb_mat(_M_M, elmat, _M_fe, _M_dof, 0, 0);
        }

        CurrentFE fe2(_M_fe);
        Real stab_coeff;

        for(UInt i = _M_mesh.numBFaces() + 1; i <= _M_mesh.numFaces(); i++){
            elmat.zero();

            UInt ad_first = _M_mesh.faceList(i).ad_first();
            UInt ad_second = _M_mesh.faceList(i).ad_second();

            stab_coeff = _M_gamma * ( pow(_M_h2pK[ad_first - 1], 2) * pow(_M_h2pK[ad_second - 1], 2) );

            _M_fe.updateFirstDerivQuadPt( _M_mesh.volumeList(ad_first) );
            fe2.updateFirstDerivQuadPt( _M_mesh.volumeList(ad_second) );

            _M_fe_bd.updateMeasQuadPt( _M_mesh.faceList(i) );

            ipstab_grad(stab_coeff, elmat, _M_fe, fe2, _M_fe_bd, 0, 0, 1);
            assemb_mat(_M_A_steady, elmat, _M_fe, fe2, _M_dof);
        }
    }

    template<typename MeshType>
    void HyperbolicSolverIP<MeshType>::add_A_unsteady() {
        ElemMat elmat(_M_fe.nbNode, 1, 1);
        ElemVec elvec(_M_fe_velocity.nbNode, NDIM);

        for(UInt i = 1; i <= _M_mesh.numVolumes(); i++){
            _M_fe.updateFirstDeriv( _M_mesh.volumeList(i) );
            elmat.zero();

            evaluate_velocity(i, elvec);
            for(UInt i_comp = 0; i_comp < NDIM; i_comp++) {
                grad(i_comp, elvec, _M_fe_velocity, elmat, _M_fe, _M_fe, 0, 0);
            }

            assemb_mat(_M_A, elmat, _M_fe, _M_dof, 0, 0);
        }
    }


    template<typename MeshType>
    void HyperbolicSolverIP<MeshType>::apply_bc() {
        // To be completed
    }
   
    template<typename MeshType>
    void HyperbolicSolverIP<MeshType>::initialize(const function_type& u0, Real t0, Real delta_t) {
        // Set initial time and time step
        _M_t0 = t0;
        _M_delta_t = delta_t;

        // Initialize bdf

        _M_bdf.initialize_unk(u0, _M_mesh, _M_reffe, _M_fe, _M_dof, _M_t0, _M_delta_t, 1);
        _M_u = *_M_bdf.unk().begin();
                
        // Check if mesh has internal faces. If not, build them

        if( !_M_mesh.hasInternalFaces() ) 
            std::cerr << "ERROR: Mesh must store internal faces" << std::endl;
        /*
          {
          UInt num_b_faces = _M_mesh.numBFaces();
          UInt num_i_faces = _M_mesh.numFaces() - _M_mesh.numBFaces();

          buildFaces(_M_mesh, std::cout, std::cerr, num_b_faces, num_i_faces, true, true, false);
          }
        */  
        // Compute the stabilization parameters

        compute_stabilization_parameters();

        // Initialize right hand side vector

        _M_b.resize(_M_dim);

        // Compute the steady part of problem matrix

        compute_M_A_steady();

    }

    template<typename MeshType>
    void HyperbolicSolverIP<MeshType>::iterate() {
    }

    template<typename MeshType>
    void HyperbolicSolverIP<MeshType>::timeAdvance() {
        // Update left hand side
        _M_A = _M_A_steady;
        add_A_unsteady();

        // Update right hand side
        _M_b = ZeroVector( _M_b.size() );
        _M_b += _M_M * _M_bdf.time_der(_M_delta_t);

        // Apply boundary conditions

        apply_bc();

        // Solve the system

        _M_solver.setMatrix(_M_A);
        _M_solver.solve(_M_u, _M_b);

        _M_bdf.shift_right(_M_u);

        // Update current time

        _M_t += _M_delta_t;
    }
}
#endif
