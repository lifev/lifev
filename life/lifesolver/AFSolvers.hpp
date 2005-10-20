/* -*- mode: c++ -*-

This file is part of the LifeV library

Author(s): Daniele Antonio Di Pietro <dipietro@unibg.it>
Date: 2-22-2005

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
   \file AFSolvers.hpp
   \author Daniele Antonio Di Pietro <dipietro@unibg.it>
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>

   \date 2-22-2005
*/

#ifndef _AFSOLVERS_HPP_
#define _AFSOLVERS_HPP_

#include <life/lifecore/GetPot.hpp>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <life/lifearray/boostmatrix.hpp>

namespace LifeV {
    /*!
      \class AFSolver
      \brief A class storing data common to all solvers based on algebraic
      factorization of N-S system.

      \author Daniele Antonio Di Pietro <dipietro@unibg.it>
      \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
    */
    template<typename __M_L_matrix_type,
             typename __C_matrix_type,
             typename __D_matrix_type,
             typename __D_T_matrix_type>
    class AFSolver {
    public:
        /**
           @name Typedefs
        */
        //@{
        typedef __M_L_matrix_type M_L_matrix_type;
        typedef __C_matrix_type C_matrix_type;
        typedef __D_matrix_type D_matrix_type;
        typedef __D_T_matrix_type D_T_matrix_type;
        //@}

        /**
           @name Constructors
        */
        //@{
        AFSolver(const M_L_matrix_type& __M_L,
                 const C_matrix_type& __C,
                 const D_matrix_type& __D,
                 const D_T_matrix_type& __D_T)
            :
            _M_M_L(__M_L),
            _M_C(__C),
            _M_D(__D),
            _M_D_T(__D_T)
        {
        }
        //@}

//         /**
//            @name Modifiers
//         */
//         //@{
//         /**!
//            \Set lumped mass matrix used for time advancement
//         */
//         void setML(const M_L_matrix_type& __M_L) {
//             _M_M_L = __M_L;
//         }

//         /**!
//            \Set stiffness matrix
//         */
//         void setC(const C_matrix_type& __C) {
//             _M_C = __C;
//         }

//         /**!
//            \Set discrete divergence operator
//         */
//         void setD(const D_matrix_type& __D) {
//             _M_D = __D;
//         }

//         /**!
//            \Set discrete gradient operator
//         */
//         void setDT(const D_T_matrix_type& __D_T) {
//             _M_D_T = __D_T;
//         }
//         //@}

        /**
           @name Members
        */
        //@{
        /**!
           \Solve member
        */
        template<typename __vector_type_u,
                 typename __vector_type_p,
                 typename __vector_type_b_u,
                 typename __vector_type_b_p>
        void solve(__vector_type_u& __u,
                   __vector_type_p& __p,
                   const __vector_type_b_u& __b_u,
                   const __vector_type_b_p& __b_p) const;
        //@}

    protected:
        /**!
           \The lumped mass matrix used for time advancement. If alfa is the
           \first coefficient of the BDF formula, delta_t the time step and
           \M the mass matrix:
           \M_L = alfa / delta_t * lump(M)
        */
        const M_L_matrix_type& _M_M_L;

        //! The stiffness matrix
        const C_matrix_type& _M_C;

        //! The discrete gradient operator
        const D_matrix_type& _M_D;

        //! The discrete divergence operator
        const D_T_matrix_type& _M_D_T;

        /**!
           The inverse of the mass matrix used for time advancement. The matrix
           M_L must be of a type for which the member function invert() is
           available
        */
        //M_L_matrix_type _M_inv_M_L;
    };

    /*!
      \class YosidaSolver
      \brief Base class for all flavours of algebraic solvers based on Yosida method.

      \author Daniele Antonio Di Pietro <dipietro@unibg.it>
    */
    template<typename __M_L_matrix_type,
             typename __C_matrix_type,
             typename __D_matrix_type,
             typename __D_T_matrix_type,
             typename __solver_u_type,
             typename __solver_p_type>
    class YosidaSolver
        :
        public AFSolver<__M_L_matrix_type,
                        __C_matrix_type,
                        __D_matrix_type,
                        __D_T_matrix_type>
    {
    public:
        /**
           @name Typedefs
        */
        //@{
        typedef __solver_u_type solver_u_type;
        typedef __solver_p_type solver_p_type;

        typedef BoostMatrix<boost::numeric::ublas::row_major> matrix_type;

        //@}
        /**!
           @name Constructors
        */
        //@{
        typedef __M_L_matrix_type M_L_matrix_type;
        typedef __C_matrix_type C_matrix_type;
        typedef __D_matrix_type D_matrix_type;
        typedef __D_T_matrix_type D_T_matrix_type;
        //@}
        /**!
           @name Constructors
        */
        //@{
        YosidaSolver(const M_L_matrix_type& __M_L,
                     const C_matrix_type& __C,
                     const D_matrix_type& __D,
                     const D_T_matrix_type& __D_T,
                     const GetPot& __data_file,
                     const std::string& __data_section,
                     solver_u_type& __solver_u,
                     solver_p_type& __solver_p)
            :
            AFSolver<__M_L_matrix_type,
                     __C_matrix_type,
                     __D_matrix_type,
                     __D_T_matrix_type>(__M_L, __C, __D, __D_T),
            _M_data_file(__data_file),
            _M_data_section(__data_section),
            _M_solver_u(__solver_u),
            _M_solver_p(__solver_p)
        {
            // Set options for solver u
            //__solver_u.setOptionsFromGetPot(M_data_file, (_M_data_section + "solver-u").data());

            // Set options for solver p
            //__solver_p.setOptionsFromGetPot(M_data_file, (_M_data_section + "solver-p").data());
        }
        //@}

        /**!
           @name Methods
        */
        //@{
        template<typename __vector_type_u,
                 typename __vector_type_p,
                 typename __vector_type_b_u,
                 typename __vector_type_b_p>
        void solve(__vector_type_u& __u, __vector_type_p& __p,
                   __vector_type_b_u& __b_u, __vector_type_b_p& __b_p);
        //@}

    protected:
        //! Data file
        const GetPot& _M_data_file;

        //! Data section
        const std::string _M_data_section;

        //! Solver for velocity systems
        solver_u_type& _M_solver_u;

        //! Solver for pressure systems
        solver_p_type& _M_solver_p;
    };

    /*!
      \class Yosida
      \brief Yosida solver implementation as described in

      A. Quarteroni, F. Saleri, A. Veneziani
      Factorization methods for the numerical approximation of Navier-Stokes
      equations
      Comput. Methods Appl. Mech. Engrg. 188 (2000) 505--526

      \author Daniele Antonio Di Pietro <dipietro@unibg.it>
    */
    template<typename __M_L_matrix_type,
             typename __C_matrix_type,
             typename __D_matrix_type,
             typename __D_T_matrix_type,
             typename __solver_u_type,
             typename __solver_p_type>
    class Yosida : public
    YosidaSolver<__M_L_matrix_type,
                 __C_matrix_type,
                 __D_matrix_type,
                 __D_T_matrix_type,
                 __solver_u_type,
                 __solver_p_type>
    {
    public:
        /**!
           @name Typedefs
        */
        //@{
        typedef __M_L_matrix_type M_L_matrix_type;
        typedef __C_matrix_type C_matrix_type;
        typedef __D_matrix_type D_matrix_type;
        typedef __D_T_matrix_type D_T_matrix_type;
        typedef __solver_u_type solver_u_type;
        typedef __solver_p_type solver_p_type;

        typedef typename YosidaSolver<__M_L_matrix_type,
                                      __C_matrix_type,
                                      __D_matrix_type,
                                      __D_T_matrix_type,
                                      __solver_u_type,
                                      __solver_p_type>::matrix_type matrix_type;
        //@}
        /**!
           @name Constructors
        */
        //@{
        Yosida(const M_L_matrix_type& __M_L,
               const C_matrix_type& __C,
               const D_matrix_type& __D,
               const D_T_matrix_type& __D_T,
               const GetPot& __data_file,
               const std::string& __data_section,
               solver_u_type& __solver_u,
               solver_p_type& __solver_p)
            :
            YosidaSolver<__M_L_matrix_type,
                         __C_matrix_type,
                         __D_matrix_type,
                         __D_T_matrix_type,
                         __solver_u_type,
                         __solver_p_type>(__M_L, __C, __D, __D_T
                                          , __data_file, __data_section,
                                          __solver_u, __solver_p) {
        }
        //@}
        /**!
           @name Methods
        */
        //@{
        template<typename __vector_type_u,
                 typename __vector_type_p,
                 typename __vector_type_b_u,
                 typename __vector_type_b_p>
        void solve(__vector_type_u& __u,
                   __vector_type_p& __p,
                   const __vector_type_b_u& __b_u,
                   const __vector_type_b_p& __b_p) {
            Chrono chrono;

            // Inverse of (e.g. lumped) mass matrix
            chrono.start();
            std::cout << "[Yosida::solve] invert mass matrix            "
                      << std::flush;
            M_L_matrix_type __H(this->_M_M_L);
            __H.invert();
            chrono.stop();
            std::cout << "in " << chrono.diff() << " s" << std::endl;

            // Approximate Schur complement
            chrono.start();
            std::cout << "[Yosida::solve] calculate Schur product       "
                      << std::flush;
            matrix_type __S( this->_M_D.size1(), this->_M_D_T.size2());
            schurProduct(this->_M_D, __H, this->_M_D_T, __S);
            chrono.stop();
            std::cout << "in " << chrono.diff() << " s" << std::endl;

            std::cout << "[Yosida::solve] exporting S to Matlab format"
                      << std::endl;
            __S.spy( "./results/spyS" );

            // Set matrices for the linear solvers
            this->_M_solver_u.setMatrix(this->_M_C);
            this->_M_solver_p.setMatrix(__S);

            // Intermediate velocity computation
            chrono.start();
            std::cout << "[Yosida::solve] compute intermediate velocity "
                      << std::flush;
            Vector  __u_tilde( __u.size() );
            this->_M_solver_u.solve( __u_tilde, __b_u );
            chrono.stop();
            std::cout << "in " << chrono.diff() << " s" << std::endl;

            // Pressure computation
            chrono.start();
            std::cout << "[Yosida::solve] compute pressure              "
                      << std::flush;
            this->_M_solver_p.solve( __p, prod(this->_M_D, __u_tilde) - __b_p );
            chrono.stop();
            std::cout << "in " << chrono.diff() << " s" << std::endl;

            // End-of-step velocity computation
            chrono.start();
            std::cout << "[Yosida::solve] compute final velocity        "
                      << std::flush;
            this->_M_solver_u.solve( __u, __b_u - prod(this->_M_D_T, __p));
            chrono.stop();
            std::cout << "in " << chrono.diff() << " s" << std::endl;
        }
    };
}

#endif
