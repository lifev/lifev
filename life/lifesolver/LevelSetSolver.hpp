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
   \file LevelSetHandler.hpp
   \author Daniele Antonio Di Pietro <dipietro@unibg.it>
   \date 1-26-2005
*/

#ifndef _LEVELSETSOLVER_H_
#define _LEVELSETSOLVER_H_

#include <utility>
#include <algorithm>
#include <functional>
#include <vector>
#include <list>
#include <set>

#include <tab.hpp>
#include <LevelSetSolver_utils.hpp>

#include <GetPot.hpp>
#include <SolverAztec.hpp>
#include <bdf.hpp>

#include <bcHandler.hpp>
#include <dof.hpp>
#include <elemMat.hpp>
#include <elemVec.hpp>
#include <elemOper.hpp>

namespace LifeV {
    /*!
      \class LevelSetSolver
      \brief Level set solver class

      \c LevelSetSolver is a solver for the level set equation. It is 
      meant to be used in multi-fluid problems.

      @author Daniele Antonio Di Pietro <dipietro@unibg.it>
      @see
    */

    template<typename MeshType, typename SolverType = SolverAztec>
    class LevelSetSolver {
    public:
        /**
           @name Subclasses
        */
        //@{
        template<int NUMNODES = 3, class PointType = node_type>
        class Face{
        public:
            /** @name Typedefs
             */
            //@{
            typedef PointType point_type;
            typedef PointType vector_type;
            typedef boost::numeric::ublas::bounded_vector<Real, NUMNODES> face_container_type;
            typedef std::vector<point_type> point_container_type;
            typedef UInt id_type;
            //@}

            /** @name Constructors, destructor
             */
            //@{
            Face() {
                _M_points.resize(3);
            }

            Face(ID id) : _M_id(id) {
                _M_points.resize(3);
            }
            //@}

            /** @name Accessors
             */
            //@{

            
            /**
               \Return a const reference to the i-th point
            */
            const point_type& point(const id_type i) const {
                return _M_points[i - 1];
            }
            //@}

            /** @name Mutators
             */
            //@{
            /**
               \Return the i-th point
            */
            point_type& point(const id_type i) {
                return _M_points[i - 1];
            }

            /**
               \Return the id
            */
            id_type& id() {
                return _M_id;
            }

            /**
               \Set i-th point
            */
            void setPoint(const id_type i, const point_type& P) {
                _M_points[i - 1] = P;
            }

            //@}

            /** @name Methods
             */
            //@{

            /**
               \Return the normal versor
            */

            inline vector_type normal() const {
                vector_type v1 = _M_points[1] - _M_points[0];
                vector_type v2 = _M_points[2] - _M_points[0];
                vector_type n = cross_prod(v1, v2);

                Real d = boost::numeric::ublas::norm_2(n);            
                n /= d;
                return n;
            }

            /**
               \Display the face
            */

            void showMe() {
                std::cout << "Face " << _M_id << ":" <<std::endl;
                std::cout << " points" << std::endl
                          << " (1) " << _M_points[0] << std::endl
                          << " (2) " << _M_points[1] << std::endl
                          << " (3) " << _M_points[2] << std::endl;
            }
            //@}

        private:
            point_container_type _M_points;
            id_type _M_id;
        };
        //@}

        /** @name Typedefs
         */
        //@{
        typedef MeshType mesh_type;
        typedef SolverType solver_type;

        typedef Vector lsfunction_type;
        typedef Vector velocity_type;
        typedef Vector rhs_type;

        typedef node_type point_type;
        typedef Face<NDIM> face_type;

        typedef enum {fluid1 = 1, fluid2 = 2} fluid_type;

        typedef Funct function_type;
        //@}

        /** @name Constructors
         */
        //@{
        LevelSetSolver(mesh_type& mesh, 
                       const GetPot& datafile,
                       solver_type& solver,
                       const GeoMap& geomap,
                       const RefFE& reffe, 
                       const QuadRule& qr, 
                       const QuadRule& qr_bd, 
                       const BCHandler& bc_h,
                       CurrentFE& fe_velocity,
                       const Dof& dof_velocity,
                       velocity_type& velocity0) 
            :
            _M_mesh(mesh),
            _M_datafile(datafile),
            _M_solver(solver),
            _M_geomap(geomap),
            _M_geomap_bd( _M_geomap.boundaryMap() ),
            _M_reffe(reffe),
            _M_reffe_bd( _M_reffe.boundaryFE() ),
            _M_qr(qr),
            _M_qr_bd(qr_bd),
            _M_fe(_M_reffe, _M_geomap, _M_qr),
            _M_fe_bd(_M_reffe_bd, _M_geomap_bd, _M_qr_bd),
            _M_fe_velocity(fe_velocity),
            _M_dof(_M_mesh, _M_reffe),
            _M_M_pattern(_M_dof, 1),
            _M_M(_M_M_pattern),
            _M_A_pattern(_M_dof, mesh, 1),
            _M_A_steady(_M_A_pattern),
            _M_A(_M_A_pattern),
            _M_dim( _M_dof.numTotalDof() ),
            _M_lsfunction(_M_dim),
            _M_b(_M_dim),
            _M_bc_h(bc_h),
            _M_dof_velocity(dof_velocity),
            _M_dim_velocity( _M_dof_velocity.numTotalDof() ),
            _M_velocity(velocity0),
            _M_bdf_order( datafile("levelset/bdf/order", 2) ),
            _M_bdf(_M_bdf_order)
        {
            _M_solver.setOptionsFromGetPot(_M_datafile, "levelset/solver");

            _M_t0 = _M_datafile("levelset/bdf/t0", 0);
            _M_delta_t = _M_datafile("levelset/bdf/delta_t", .02);

            _M_gamma = _M_datafile("levelset/ipstab/gamma", 0.125);
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
        void initialize(const function_type& lsfunction0);

        /**
           \Update the right hand side for time advancement
        */
        void timeAdvance();

        /**
           \Return the level set function
        */
        lsfunction_type& lsfunction() {
            return _M_lsfunction;
        }

        /**
           \Update left hand side and solve the system
        */
        void iterate();

        /**
           \Return the mass of i-th fluid
        */
        Real computeMassOfFluid(fluid_type);
        //@}

    protected:
    
        //! Mesh
        mesh_type& _M_mesh;

        //! Data file
        const GetPot& _M_datafile;

        //! The linear solver
        solver_type _M_solver;
        
        //! Geo map
        const GeoMap& _M_geomap;

        //! Boundary geo map

        const GeoMap& _M_geomap_bd;

        //! Reference FE
        const RefFE& _M_reffe;

        //! Reference boundary FE
        const RefFE& _M_reffe_bd;

        //! Quadrature rule for volume integrals
        const QuadRule& _M_qr;

        //! Quadrature rule for boundary integrals
        const QuadRule& _M_qr_bd;

        //! Current finite element for the levelset function
        CurrentFE _M_fe;

        //! Current boundary element for the levelset function
        CurrentBdFE _M_fe_bd;

        //! Current finite element for the velocity
        CurrentFE _M_fe_velocity;

        //! Degrees of freedom for the level set function
        Dof _M_dof;

        //! Pattern for the mass matrix
        MSRPatt _M_M_pattern;

        //! The mass matrix
        MSRMatr<Real> _M_M;

        //! Pattern for the problem matrix
        MSRPatt _M_A_pattern;

        //! The steady part of problem matrix
        MSRMatr<Real> _M_A_steady;

        //! The full problem matrix (steady + unsteady)
        MSRMatr<Real> _M_A;

        //!The number of total dofs for the level set function
        UInt _M_dim;

        //! The level set function
        lsfunction_type _M_lsfunction;

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
           \velocity and the level set function.
        */
        inline void evaluate_velocity(UInt fe_id, ElemVec& elvec) {
            UInt dim = _M_dof_velocity.numTotalDof();

            for(UInt j = 0; j < (UInt)_M_fe_velocity.nbNode; j++) {
                UInt jloc = _M_fe_velocity.patternFirst(j);
                for(UInt i_comp = 0; i_comp < NDIM; i_comp++) {
                    UInt jglo = _M_dof_velocity.localToGlobal(fe_id, jloc + 1) - 1 + i_comp * dim;
                    elvec.vec()[jloc + i_comp * _M_fe_velocity.nbNode] = _M_velocity[jglo];
                }
            }
        }

        /**
           \Compute the advection part of problem matrix, which is 
           \time-dependent in non-steady multi-fluid problems
        */
        void add_A_unsteady();
        //@}
    };

    template<typename MeshType, typename SolverType>
    void LevelSetSolver<MeshType, SolverType>::compute_M_A_steady() {
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

    template<typename MeshType, typename SolverType>
    void LevelSetSolver<MeshType, SolverType>::add_A_unsteady() {
        ElemMat elmat(_M_fe.nbNode, 1, 1);
        ElemVec elvec(_M_fe.nbNode, NDIM);

        for(UInt i = 1; i <= _M_mesh.numVolumes(); i++){
            _M_fe.updateJac( _M_mesh.volumeList(i) );
            elmat.zero();

            evaluate_velocity(i, elvec);

            for(UInt i_comp = 0; i_comp < NDIM; i_comp++) {
                grad(i_comp, elvec, elmat, _M_fe, _M_fe, 0, 0);
            }

            assemb_mat(_M_A, elmat, _M_fe, _M_dof, 0, 0);
        }
    }

    template<typename MeshType, typename SolverType>
    void LevelSetSolver<MeshType, SolverType>::initialize(const function_type& lsfunction0) {
        // Initialize bdf

        _M_bdf.initialize_unk(lsfunction0, _M_mesh, _M_reffe, _M_fe, _M_dof, _M_t0, _M_delta_t, 1);

        // Check if mesh has internal faces. If not, build them

        if( !_M_mesh.hasInternalFaces() ) {
            UInt num_b_faces = _M_mesh.numBFaces();
            UInt num_i_faces = _M_mesh.numFaces() - _M_mesh.numBFaces();

            buildFaces(_M_mesh, std::cout, std::cerr, num_b_faces, num_i_faces, true, true, false);
        }
            
        // Compute the stabilization parameters

        compute_stabilization_parameters();

        // Initialize solution vector

        _M_lsfunction = ZeroVector(_M_dim);

        // Initialize right hand side vector

        _M_b.resize(_M_dim);

        // Compute the steady part of problem matrix

        compute_M_A_steady();

    }

    template<typename MeshType, typename SolverType>
    void LevelSetSolver<MeshType, SolverType>::iterate() {
    }

    template<typename MeshType, typename SolverType>
    void LevelSetSolver<MeshType, SolverType>::timeAdvance() {
        // Update left hand side

        _M_A = _M_A_steady;
        add_A_unsteady();

        // Update right hand side

        _M_b = ZeroVector( _M_b.size() );
        _M_b += _M_M * _M_bdf.time_der(_M_delta_t);

        _M_solver.setMatrix(_M_A);
        _M_solver.solve(_M_lsfunction, _M_b);

        _M_bdf.shift_right(_M_lsfunction);
    }

    template<typename MeshType, typename SolverType>
    Real LevelSetSolver<MeshType, SolverType>::computeMassOfFluid(fluid_type fluid_id) {
        Real mass = 0.;
        Real ls_fun;
        int jg;

        for(UInt i = 1; i <= _M_mesh.numVolumes(); i++) {
            _M_fe.updateJac( _M_mesh.volumeList(i) );

            for(int iq = 0; iq < _M_fe.nbQuadPt; iq++) {
                ls_fun = 0.;
                for(int j = 0; j < _M_fe.nbNode; j++) {
                    jglo = _M_dof.localToGlobal(i, j + 1) - 1;
                    ls_fun += _M_fe.phi(j, iq) * U(jglo);
                }
                if(fluid_id == fluid1)
                    mass += (ls_fun < 0 ? 1. : 0.) * _M_fe.weightDet(iq);
                else
                    mass += (ls_fun >= 0 ? 1. : 0.) * _M_fe.weightDet(iq);
            }
        }
        return mass;
    }
}
#endif
