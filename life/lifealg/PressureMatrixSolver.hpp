/* -*- mode: c++ -*-

This file is part of the LifeV library

Author(s): Daniele A. Di Pietro <dipietro@unibg.it>
Date: 2-5-2005

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
   \file PressureMatrixSolver.hp
   \author Daniele A. Di Pietro <dipietro@unibg.it>
   \date 2-5-2005
*/

#ifndef _PRESSUREMATRIXSOLVER_HPP_
#define _PRESSUREMATRIXSOLVER_HPP_

#include <life/lifecore/GetPot.hpp>
#include <life/lifearray/pattern.hpp>
#include <life/lifearray/boostmatrix.hpp>

namespace LifeV {
    /*!
      \class PressureMatrixSolver
      \brief A class for solving the linear system arising from the discretization
      of Stokes and Navier-Stokes equations. Matrices are assumed to share the interface
      of boost::numeric::ublas::compressed_matrix format.

      \author Daniele Antonio Di Pietro <dipietro@unibg.it>
    */

    template<typename MatrixTypeC, typename MatrixTypeD, typename MatrixTypeG, typename SolverType>
    class PressureMatrixSolver {
    public:
        /**
           @name Typedefs
        */
        //@{
        typedef Vector vector_type;
        typedef MatrixTypeC matrix_type_C;
        typedef MatrixTypeD matrix_type_D;
        typedef MatrixTypeG matrix_type_G;
        typedef SolverType solver_type;
        //@}

        /**
           @name Constructor
        */
        //@{
        PressureMatrixSolver(const matrix_type_C& C,
                             const matrix_type_D& D,
                             const matrix_type_G& G,
                             const GetPot& data_file,
                             const std::string& data_section)
            :
            _M_C(C),
            _M_D(D),
            _M_G(G),
            _M_data_file(data_file),
            _M_data_section(data_section),
            _M_tol(_M_data_file( (_M_data_section + "cg/tol").data(), 1e-6)),
            _M_maxit(_M_data_file( (_M_data_section + "cg/maxit").data(), 100)),
            _M_nit(0),
            _M_converged(false)
        {
        }
        //@}
        
        /**
           @name Accessors
        */
        //@{
        /**!
           \Return the number of iterations needed in the last call of
           \solve() method
        */
        const UInt numIterations() const {
            return _M_nit;
        }

        /**!
           \Return true if last call to solve() method converged
        */
        const bool converged() const {
            return _M_converged;
        }

        /**!
           \Return residual at last iteration
        */
        const Real residual() const {
            return _M_res;
        }
        //@}

        /**
           @name Methods
        */
        //@{
        void solve(vector_type& u, vector_type& p, vector_type& b1, vector_type& /* b2 */);
        //@}

        /**
           @name Friend operators
        */
        //@{

        /**
           \Report the result of last call to solve() method
        */
        template<typename _MatrixTypeC, typename _MatrixTypeD, typename _MatrixTypeG, typename _SolverType>
        friend std::ostream& operator<<(std::ostream&, PressureMatrixSolver<_MatrixTypeC, _MatrixTypeD, _MatrixTypeG, _SolverType>&);
        //@}

    private:
        //! C block
        const matrix_type_C& _M_C;

        //! D block
        const matrix_type_D& _M_D;

        //! G block
        const matrix_type_G& _M_G;

        //! Data file
        const GetPot& _M_data_file;

        //! Data section
        const std::string _M_data_section;

        //! Tolerance for conjugate gradient method
        Real _M_tol;

        //! Maximum number of iteration for conjugate gradient method
        UInt _M_maxit;

        //! Number of iterations for the last call of solve() method
        UInt _M_nit;

        //! Convergence flag for the last call of solve() method
        bool _M_converged;

        //! Residual after last call of solve() method
        Real _M_res;

        //! Underlying linear solver
        solver_type _M_solver;
    };

    // Implementations
    
    template<typename MatrixTypeC, typename MatrixTypeD, typename MatrixTypeG, typename SolverType>
    void PressureMatrixSolver<MatrixTypeC, MatrixTypeD, MatrixTypeG, SolverType>
    ::
    solve(vector_type& __u, vector_type& __p, vector_type& b1, vector_type& /* b2 */) {
        vector_type r_k( __p.size() );
        vector_type p_k( __p.size() );
        vector_type g( __p.size() );

        vector_type temp1( __u.size() );
        vector_type temp2( __p.size() );

        Real alfa_k;
        Real beta_k;

        _M_solver.setMatrix( _M_C );
        _M_solver.solve( temp1, b1 );

        g = prod( _M_D, temp1 );

        r_k = g;
        p_k = r_k;

        _M_nit = 0;
        _M_res = norm_2(r_k);

        while(_M_res > _M_tol && _M_nit < _M_maxit) {
            _M_solver.solve(temp1, prod( _M_G,  p_k ) );
            temp2 = prod( _M_D, temp1 );

            alfa_k = inner_prod(p_k, r_k) / inner_prod(p_k, temp2);
            __p += alfa_k * p_k;
            
            _M_solver.solve(temp1, _M_G * __p);
            r_k = g - prod( _M_D, temp1 );

            beta_k = inner_prod(temp2, r_k) / inner_prod(temp2, p_k);

            p_k = r_k - beta_k * p_k;

            _M_nit++;
            _M_res = norm_2(r_k);
        }
        _M_solver.solve(__u, b1 -  prod( _M_G, __p ) );
   
        (_M_nit == _M_maxit) ? _M_converged = false : _M_converged = true;
    }

    template<typename _MatrixTypeC, typename _MatrixTypeD, typename _MatrixTypeG, typename SolverType>
    std::ostream& operator<<(std::ostream& __ostr, PressureMatrixSolver<_MatrixTypeC, _MatrixTypeD, _MatrixTypeG, SolverType>& __PMS) {
        __ostr << std::endl;
        __ostr << "==================================================" << std::endl;
        __ostr << "PMS report" << std::endl;
        __ostr << "==================================================" << std::endl;
        __ostr << "Convergence    ";
        __PMS._M_converged ? __ostr << "YES" : __ostr << "NO";
        __ostr << std::endl;
        __ostr << "--------------------------------------------------" << std::endl;
        __ostr << "Final residual " << __PMS._M_res << std::endl;
        __ostr << "==================================================" << std::endl;

        return __ostr;
    }

}

#endif
