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
/*!
  \file HyperbolicSolverIP.hpp
  \author Daniele Antonio Di Pietro <dipietro@unibg.it>
  \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>

  \date 1-26-2005
*/

#ifndef _HYPERBOLICSOLVERIP_HPP_
#define _HYPERBOLICSOLVERIP_HPP_

#define L_HSIP_AZTEC 0
#define L_HSIP_UMFPACK 1
#define L_HSIP_PETSC 2

#define L_HSIP_LINEAR_SOLVER L_HSIP_UMFPACK
#define L_HSIP_USE_BOOST_MATRIX 0

#include <utility>
#include <algorithm>
#include <functional>
#include <vector>
#include <list>
#include <set>

#include <boost/numeric/ublas/operation.hpp>

#include <life/lifecore/chrono.hpp>
#include <life/lifecore/GetPot.hpp>

#include <life/lifearray/boostmatrix.hpp>
#include <life/lifearray/tab.hpp>

#include <life/lifemesh/dataMesh.hpp>

#include <life/lifecore/GetPot.hpp>
#include <life/lifefem/bdf.hpp>


#if L_HSIP_LINEAR_SOLVER == L_HSIP_AZTEC
#include <life/lifealg/SolverAztec.hpp>
#elif L_HSIP_LINEAR_SOLVER == L_HSIP_UMFPACK
#include <life/lifealg/SolverUMFPACK.hpp>
#elif L_HSIP_LINEAR_SOLVER == L_HSIP_PETSC
#include <life/lifealg/SolverPETSC.hpp>
#endif

#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/elemOper2Fluids.hpp>

namespace LifeV {
    /*!
      \class HyperbolicSolverIP
      \brief hyperbolic solver class

      Hyperbolic solver with IP stabilization

      \author Daniele Antonio Di Pietro <dipietro@unibg.it>
    */

    template<typename MeshType>
    class HyperbolicSolverIP {
    public:
        /*! @name Typedefs */
        //@{

        typedef MeshType mesh_type;



#if L_HSIP_LINEAR_SOLVER == L_HSIP_AZTEC
        typedef SolverAztec solver_type;
        typedef MSRPatt pattern_type;
        typedef MSRMatr<Real> matrix_type;

#elif L_HSIP_LINEAR_SOLVER == L_HSIP_UMFPACK
        typedef SolverUMFPACK solver_type;
        typedef CSRPatt pattern_type;
#if L_HSIP_USE_BOOST_MATRIX
        typedef BoostMatrix<boost::numeric::ublas::column_major> matrix_type;
#else
        typedef CSRMatr<CSRPatt, Real> matrix_type;
#endif

#elif L_HSIP_LINEAR_SOLVER == L_HSIP_PETSC
        typedef SolverPETSC solver_type;
        typedef CSRPatt pattern_type;
#if L_HSIP_USE_BOOST_MATRIX
        typedef BoostMatrix<boost::numeric::ublas::row_major> matrix_type;
#else
        typedef CSRMatr<CSRPatt, Real> matrix_type;
#endif

#endif /* L_HSIP_LINEAR_SOLVER */

        typedef Vector u_type;
        typedef Vector velocity_type;
        typedef Vector rhs_type;

        typedef Funct function_type;

        /*! @name Constructors */
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
                           velocity_type& velocity0,
                           UInt nbComp = 1)
            :
            _M_mesh(mesh),
            _M_data_file(data_file),
            _M_data_section(data_section),
            _M_reffe(reffe),
            _M_reffe_bd( _M_reffe.boundaryFE() ),
            _M_qr(qr),
            _M_qr_bd(qr_bd),
            _M_fe(_M_reffe, getGeoMap(_M_mesh), _M_qr),
            _M_fe_2(_M_reffe, getGeoMap(_M_mesh), _M_qr),
            _M_fe_bd(_M_reffe_bd, getGeoMap(_M_mesh).boundaryMap(), _M_qr_bd),
            _M_fe_velocity(fe_velocity),
            _M_dof(_M_mesh, _M_reffe),
            _M_M_pattern(_M_dof, nbComp),
            _M_M(_M_M_pattern),
#if L_HSIP_LINEAR_SOLVER == L_HSIP_AZTEC
            _M_A_pattern(_M_dof, _M_mesh, nbComp),
#else
            _M_A_pattern(_M_dof, nbComp, _M_mesh),
#endif
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
            _M_bdf(_M_bdf_order),
            _M_monitored_times(5)
        {
#if L_HSIP_LINEAR_SOLVER == L_HSIP_AZTEC
            _M_solver.setOptionsFromGetPot(_M_data_file,
                                           (_M_data_section +
                                            "/solver_aztec").data());
            std::cout << "** HSIP ** Using AZTEC solver" << std::endl;
#elif L_HSIP_LINEAR_SOLVER == L_HSIP_UMFPACK
            std::cout << "** HSIP ** Using UMFPACK solver" << std::endl;
#elif L_HSIP_LINEAR_SOLVER == L_HSIP_PETSC
             _M_solver.setOptionsFromGetPot(_M_data_file,
                                            (_M_data_section +
                                             "/solver_petsc").data());
            std::cout << "** HSIP ** Using PETSC solver" << std::endl;
#endif

#if L_HSIP_USE_BOOST_MATRIX
            std::cout << "** HSIP ** Using boost matrix" << std::endl;
#else
            std::cout << "** HSIP ** Using MSR/CSRMatr" << std::endl;
#endif

            _M_gamma = _M_data_file( (_M_data_section + "/ipstab/gamma").data(), 0.125 );

            switch( _M_fe.nbNode )
            {
                case 4:
                    _M_fToP = LinearTetra::fToP;
                    break;
                case 10:
                    _M_fToP = QuadraticTetra::fToP;
                    break;
                case 8:
                    _M_fToP = LinearHexa::fToP;
                    break;
                case 20:
                    _M_fToP = QuadraticHexa::fToP;
                    break;
                default:
                    ERROR_MSG( "This refFE is not allowed with IP stabilisation" );
                    break;
            }

            _M_verbose = false;
        }

        //@}

        /*! @name Accessors
         */
        //@{

        /*! \return current time */
        inline Real currentTime() {
            return _M_t;
        }

        /*! \return the dof table */
        const Dof& dof() const {
            return _M_dof;
        }

        /*! \return time at next step */
        inline Real nextStepTime() {
            return _M_t + _M_delta_t;
        }

        /*! \return the current finite element */
        const CurrentFE& fe() const {
            return _M_fe;
        }

        /*! \return the mesh */
        const mesh_type& mesh() const {
            return _M_mesh;
        }

        /*! \return verbose mode */
        bool verboseMode() const {
            return _M_verbose;
        }

        //@}

        /*! @name Mutators */
        //@{

        //! Set verbose mode
        void setVerboseMode() {
            _M_verbose = true;
        }

        //! Set advection field
        void setVelocity(velocity_type& __velocity) {
            _M_velocity = __velocity;
        }

        /*! \return a reference to the current finite element */
        CurrentFE& fe() {
            return _M_fe;
        }

        //@}

        /*! @name Methods */
        //@{

        //! Initialize the solver
        void initialize(const function_type& lsfunction0,
                        Real t0, Real delta_t);

        //! Update the right hand side for time advancement
        void timeAdvance();

        //! Return the current numeric solution
        u_type const & u() const {
            return _M_u;
        }

        //! Update left hand side and solve the system
        void iterate();

        /*!
          Export matrices A and M in a format suitable for use with
          Matlab's spy function
        */
        void spy(std::string __path) {
            _M_A.spy( __path + "spyA" );
            _M_M.spy( __path + "spyM" );
        }
        //@}

        /**
           @name Friend operators
        */
        //@{

        /**
           \Report the result of last call to timeAdvance() method
        */
        template<typename _MeshType>
        friend std::ostream& operator<<(std::ostream&, HyperbolicSolverIP<_MeshType>&);
        //@}

    protected:

        /*! \return a reference the current numeric solution */
        u_type & u() {
            return _M_u;
        }

    private:
        //! Mesh
        mesh_type& _M_mesh;

        //! Data file
        const GetPot& _M_data_file;

        //! Data section
        const std::string _M_data_section;

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

        //! Second current finite element for the unknown (IP terms)
        CurrentFE _M_fe_2;

        //! Current boundary element for the unknown
        CurrentBdFE _M_fe_bd;

        //! Current finite element for the velocity
        CurrentFE _M_fe_velocity;

        //! Degrees of freedom for the unknown
        Dof _M_dof;

        //! Pattern for the mass matrix
        pattern_type _M_M_pattern;

        //! The mass matrix
        matrix_type _M_M;

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
        velocity_type _M_velocity;

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
        std::vector<Real> _M_hpK;

        //! The linear solver
        solver_type _M_solver;

        //! The coefficient for stabilization parameters scaling
        Real _M_gamma;

        //! Time for matrix assembling
        std::vector<Real> _M_monitored_times;

        /*! @name Methods
         */
        //@{

        //! Compute the stabilization parameters
        void compute_stabilization_parameters() {
            _M_hpK.clear();
            _M_hpK.resize( _M_mesh.numVolumes(), 0. );

            for(UInt i = _M_mesh.numBFaces() + 1; i <= _M_mesh.numFaces(); i++){
                _M_fe_bd.updateMeasQuadPt( _M_mesh.faceList(i) );
                Real diameter = pow(_M_fe_bd.measure() * 2., 0.5);

                _M_hpK[_M_mesh.faceList(i).ad_first() - 1] =
                    std::max<Real>(_M_hpK[_M_mesh.faceList(i).ad_first() - 1],
                             diameter);
                _M_hpK[_M_mesh.faceList(i).ad_second() - 1] =
                    std::max<Real>(_M_hpK[_M_mesh.faceList(i).ad_second() - 1],
                             diameter);
            }
        }

        //! Compute the mass matrix and the constant part of problem matrix
        void compute_M_A_steady();

        /*!
          Evaluate velocity field on qn_id-th quadrature node of fe_id-th
          element. It is assumed that the same mesh is used for both the
          velocity field and the unknown.
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

        /*!
          Compute the advection part of problem matrix, which is
          time-dependent in non-steady multi-fluid problems.
        */
        void add_A_unsteady();

        /*!
          Add the stabilization terms. It may be necessary to re-calculate
          stabilization terms at every time step because the stabilization
          parameter depends on the velocity field.
        */
        void add_stabilization();

        /*!
          Apply boundary conditions. Only Dirichlet (essential) boundary
          conditions are considered on the sole inflow boundary nodes. At
          present time only bc given in functional form are considered.
        */

        void apply_bc();
        //@}

        typedef ID ( *FTOP )( ID const localFace, ID const point );
        FTOP _M_fToP;

    };

// Implementations

    template<typename MeshType>
    void HyperbolicSolverIP<MeshType>::compute_M_A_steady() {
        _M_M.zeros();
        _M_A_steady.zeros();

        ElemMat __elmat(_M_fe.nbNode, 1, 1);
        ElemVec __elvec(_M_fe.nbNode, NDIM);

        Real coeff = _M_bdf.coeff_der(0) / _M_delta_t;

        for(UInt i = 1; i <= _M_mesh.numVolumes(); i++){
            _M_fe.updateJac( _M_mesh.volumeList(i) );

            __elmat.zero();
            mass(coeff, __elmat, _M_fe, 0, 0);
            assemb_mat(_M_A_steady, __elmat, _M_fe, _M_dof, 0, 0);

            __elmat.zero();
            mass(1., __elmat, _M_fe, 0, 0);
            assemb_mat(_M_M, __elmat, _M_fe, _M_dof, 0, 0);
        }
    }

    template<typename MeshType>
    void HyperbolicSolverIP<MeshType>::add_A_unsteady() {
        ElemVec __elvec(_M_fe_velocity.nbNode, NDIM);
        ElemMat __elmat(_M_fe.nbNode, 1, 1);

        for(UInt i = 1; i <= _M_mesh.numVolumes(); i++){
            _M_fe.updateFirstDeriv( _M_mesh.volumeList(i) );
            __elmat.zero();

            evaluate_velocity(i, __elvec);
            for(UInt i_comp = 0; i_comp < NDIM; i_comp++) {
                grad(i_comp, __elvec, _M_fe_velocity, __elmat, _M_fe, _M_fe, 0, 0);
            }

            assemb_mat(_M_A, __elmat, _M_fe, _M_dof, 0, 0);
        }
    }

    template<typename MeshType>
    void HyperbolicSolverIP<MeshType>::add_stabilization(){
        ElemMat __elmat11(_M_fe.nbNode, 1, 1);
        ElemMat __elmat12(_M_fe.nbNode, 1, 1);
        ElemMat __elmat21(_M_fe.nbNode, 1, 1);
        ElemMat __elmat22(_M_fe.nbNode, 1, 1);

        //CurrentFE _M_fe_2(_M_fe);
        Real __stab_coeff;

        ElemVec __beta( _M_fe_bd.nbNode, NDIM );

        for(UInt __face_id = _M_mesh.numBFaces() + 1; __face_id <= _M_mesh.numFaces(); __face_id++){
            __elmat11.zero();
            __elmat12.zero();
            __elmat21.zero();
            __elmat22.zero();

            UInt __ad_first = _M_mesh.faceList(__face_id).ad_first();
            UInt __ad_secnd = _M_mesh.faceList(__face_id).ad_second();

            __stab_coeff = .5 * _M_gamma * ( pow(_M_hpK[__ad_first - 1], 2) + pow(_M_hpK[__ad_secnd - 1], 2) );

            _M_fe.updateFirstDerivQuadPt( _M_mesh.volumeList(__ad_first) );
            _M_fe_2.updateFirstDerivQuadPt( _M_mesh.volumeList(__ad_secnd) );

            _M_fe_bd.updateMeasNormalQuadPt( _M_mesh.faceList(__face_id) );

            // Retrieve local velocity dofs

            __beta.zero();

            // Get the position of the face on the first element sharing it

            UInt iFaEl = _M_mesh.faceList(__face_id).pos_first();

            for (int __node_id = 0; __node_id < _M_fe_bd.nbNode; ++__node_id) {
                UInt iloc = _M_fToP(iFaEl, __node_id + 1);

                for ( int __coor_id = 0; __coor_id < _M_fe.nbCoor; ++__coor_id ) {
                    UInt __iglo = _M_dof_velocity.localToGlobal(__ad_first, iloc) - 1 +
                        __coor_id * _M_dof.numTotalDof();
                    __beta.vec()[__coor_id * _M_fe_bd.nbNode + __node_id] = _M_velocity[__iglo];
                }
            }

            ipstab_bagrad(__stab_coeff, __elmat11, _M_fe, _M_fe, __beta, _M_fe_bd);
            ipstab_bagrad(__stab_coeff, __elmat12, _M_fe, _M_fe_2, __beta, _M_fe_bd);
            ipstab_bagrad(__stab_coeff, __elmat21, _M_fe_2, _M_fe, __beta, _M_fe_bd);
            ipstab_bagrad(__stab_coeff, __elmat22, _M_fe_2, _M_fe_2, __beta, _M_fe_bd);

            assemb_mat(_M_A, __elmat11, _M_fe, _M_dof);
            assemb_mat(_M_A, __elmat12, _M_fe, _M_fe_2, _M_dof);
            assemb_mat(_M_A, __elmat21, _M_fe_2, _M_fe, _M_dof);
            assemb_mat(_M_A, __elmat22, _M_fe_2, _M_dof);
        }

    }

    template<typename MeshType>
    void HyperbolicSolverIP<MeshType>::apply_bc() {
        bcManage( _M_A, _M_b, _M_mesh, _M_dof, _M_bc_h, _M_fe_bd, 1.0, _M_t );
    }

    template<typename MeshType>
    void HyperbolicSolverIP<MeshType>::initialize(const function_type& u0, Real t0, Real delta_t) {
        Chrono __chrono;
        __chrono.start();

        // Set initial time, current time and time step
        _M_t0 = t0;
        _M_t = t0;
        _M_delta_t = delta_t;

        // Initialize bdf

        _M_bdf.initialize_unk(u0, _M_mesh, _M_reffe, _M_fe, _M_dof,
                              _M_t0, _M_delta_t, 1);
        _M_u = *_M_bdf.unk().begin();
        _M_bdf.initialize_unk(u0, _M_mesh, _M_reffe, _M_fe, _M_dof,
                              _M_t0-_M_delta_t, _M_delta_t, 1);

        // Check if mesh has internal faces. If not, build them

        if( !_M_mesh.hasInternalFaces() ) {
            std::cerr << "WARNING: mesh did not have internal faces." << std::endl;
            std::cerr << "         The pattern might not have been constructed properly" << std::endl;

            UInt num_b_faces = _M_mesh.numBFaces();
            UInt num_i_faces = _M_mesh.numFaces() - _M_mesh.numBFaces();

            buildFaces(_M_mesh, std::cout, std::cerr, num_b_faces, num_i_faces, true, true, false);
        }

        // Compute stabilization parameters

        compute_stabilization_parameters();

        // Initialize right hand side vector

        _M_b.resize(_M_dim);

        // Compute the steady part of problem matrix

        compute_M_A_steady();
        __chrono.stop();

        _M_monitored_times[0] = __chrono.diff();

    }

    template<typename MeshType>
    void HyperbolicSolverIP<MeshType>::iterate() {
    }

    template<typename MeshType>
    void HyperbolicSolverIP<MeshType>::timeAdvance() {
        Chrono __chrono;

        _M_bdf.shift_right(_M_u);

        if(_M_verbose) std::cout << "** HSIP ** Computing problem matrix" << std::endl;

        _M_A = _M_A_steady;

        __chrono.start();
        add_A_unsteady();
        __chrono.stop();
        _M_monitored_times[1] = __chrono.diff();

        if(_M_verbose) std::cout << "** HSIP ** Adding stabilization" << std::endl;

        __chrono.start();
        add_stabilization();
        __chrono.stop();
        _M_monitored_times[2] = __chrono.diff();

        if(_M_verbose) std::cout << "** HSIP ** Adding RHS unsteady terms" << std::endl;

#if L_HSIP_USE_BOOST_MATRIX
        for_each(_M_b.begin(), _M_b.end(), boost::lambda::_1 = 0.0 );
        axpy_prod( _M_M, _M_bdf.time_der(_M_delta_t), _M_b, false );

        _M_b += prod( _M_M, _M_bdf.time_der( _M_delta_t ) );
#else
        _M_b = ZeroVector( _M_b.size() );
        _M_b += _M_M * _M_bdf.time_der(_M_delta_t);
#endif

        if(_M_verbose) std::cout << "** HSIP ** Applying boundary conditions" << std::endl;

        __chrono.start();
        apply_bc();
        __chrono.stop();
        _M_monitored_times[3] = __chrono.diff();

        __chrono.start();
        if(_M_verbose) std::cout << "** HSIP ** Passing matrix to linear solver" << std::endl;
        _M_solver.setMatrix(_M_A);

        if(_M_verbose) std::cout << "** HSIP ** Solving linear system" << std::endl;
        _M_solver.solve(_M_u, _M_b);

        if(_M_verbose) std::cout << "** HSIP ** Updating bdf object" << std::endl;

        __chrono.stop();
        _M_monitored_times[4] = __chrono.diff();

        _M_t += _M_delta_t;
    }

template<typename _MeshType>
std::ostream& operator<<(std::ostream& __ostr, HyperbolicSolverIP<_MeshType>& __HS) {
    __ostr << "==================================================" << std::endl;
    __ostr << "** HSIP ** Report" << std::endl;
    __ostr << "==================================================" << std::endl;
#if L_HSIP_LINEAR_SOLVER != L_HSIP_UMFPACK
    __ostr << "Convergence                         ";
    __HS._M_solver.converged() ? __ostr << "YES" : __ostr << "NO";
    __ostr << std::endl;
    __ostr << "Number of iterations                "
           << __HS._M_solver.iterations() << std::endl;
    __ostr << "--------------------------------------------------" << std::endl;
#endif
    __ostr << "Timings" << std::endl << std::endl;
    __ostr << "Solver initialization               " <<
        __HS._M_monitored_times[0] << " s" << std::endl;
    __ostr << "Unsteady contribution assembling    " <<
        __HS._M_monitored_times[1] << " s" << std::endl;
    __ostr << "Stabilization terms assembling      " <<
        __HS._M_monitored_times[2] << " s" << std::endl;
    __ostr << "Boundary conditions handling        " <<
        __HS._M_monitored_times[3] << " s" << std::endl;
    __ostr << "System solution                     " <<
        __HS._M_monitored_times[4] << " s" << std::endl;
    __ostr << "==================================================" << std::endl;

    return __ostr;
}

} // namespace LifeV

#endif /* _HYPERBOLICSOLVERIP_HPP_ */
