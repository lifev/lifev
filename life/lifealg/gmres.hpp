/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2005-04-07

  Copyright (C) 2005 EPFL

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
   \file gmres.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-04-07
 */
#ifndef __GMRES_H
#define __GMRES_H 1

#include <life/lifearray/tab.hpp>

#include <life/lifealg/givens.hpp>
#include <life/lifealg/iteration.hpp>


namespace LifeV
{
/*!
  \fn int gmres(const Matrix_t &A,
  Vector &x,const VectorB &b,
  const Preconditioner &M,
  int m,
  Iteration& outer)
  \brief Generalized Minimum RESidual Restarted

  \latexonly
  \begin{algorithm}
  \begin{algorithmic}
  \STATE{$\epsilon \gets$ tolerance of order $10^{-6}$ for example}
  \STATE{$r^{0} \gets P^{-1} \left( b \,-\, A x^{0} \right)$}
  \STATE{$\beta \gets \| r^0 \|_2$}
  \WHILE{$\| \beta \|_2 > \epsilon$}
  \STATE{$v^1 \gets \displaystyle\frac{r^0}{\beta}$}
  \FOR{$k =1 \ldots m $}
  \STATE{$w \gets P^{-1} v^k$}
  \FOR{$i =1 \ldots k $}
  \STATE{$h_{k,i} \gets (w, v^k)$}
  \STATE{$w \gets w - h_{k,i} v^j$}
  \ENDFOR
  \STATE{$h_{k+1,k} \gets  \| w \|_2$}
  \STATE{$v^{k+1} \gets \displaystyle\frac{w}{h_{k+1,k}}$}
  \ENDFOR
  \STATE{$V^m \gets \left[ v^1 \ldots v^m \right]$}
  \STATE{$H^m \gets (h_{i,j})_{1 \leq i \leq j+1; 1 \leq j \leq m}$}
  \STATE{$y^m \gets \operatorname{argmin}_{y} \| \beta e^1 - H^m y \|_2$}
  \STATE{$x^m \gets x^0 \,+\, V^m y^m$}
  \STATE{$r^m \gets P^{-1} \left( b \,-\, A x^{m} \right)$}
  \STATE{$\beta \gets \| r^m \|_2$}
  \ENDWHILE
  \end{algorithmic}
  \caption{General Minimum Residual(\texttt{m})}
  \end{algorithm}
  \endlatexonly

  @param A matrix to invert
  @param x solution vector
  @param b right hand side vector
  @param M preconditioner
  @param m number of iteration before restart
  @param outer Iteration data structure to control the solver

  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>

  @see Preconditioner
  @see Y. Saad and M. H. Schultz,
  GMRES: A generalized minimal  residual algorithm   for solving nonsymmetric linear systems
  SIAM J. Sci. Stat. Comp. Vol 7 1986
*/
template <
    class Matrix_t,
    class Vector,
    class VectorB,
    class Preconditioner,
    class IterationKrylov
    >
int
gmres( const Matrix_t &A,
       Vector &x,
       const VectorB &b,
       const Preconditioner &M,
       id_type m,
       IterationKrylov& outer )
{
    typedef typename Matrix_t::value_type value_type;
    typedef ublas::matrix<value_type, ublas::column_major> hessenberg_type;
    typedef ublas::vector<Vector> krylov_type;
    typedef ublas::vector<value_type> vector_type;

    hessenberg_type H( m+1, m );
    krylov_type V(m+1);
    for(int i = 0;i < m+1;i++)
    {
        V(i).resize(x.size());
    }
    vector_type s(m+1);
    Vector w(x), r(x), u(x);

    char *text = "gmres inner";

    std::vector<GivensRotation<value_type> > rotations(m+1);

    axpy_prod(A, u, w, true);
    //SAxpy(-1,w,b,w);
    w.assign( b - w );

    M.solve(w, r);

    typename IterationKrylov::value_type beta = std::abs( ublas::norm_2( r ) );
    outer.setInitialResidual( beta );
    outer.reset();
    while (! outer.isFinished( beta ) )
    {
        V( 0 ) = r/beta;
        s = ublas::zero_vector<value_type>( s.size() );
        s(0) = beta;

        int i = 0;
        //Iter inner(outer.normb(), m, outer.tol(),text);
        //IterationBase<value_type> inner(outer.normb(), m, outer.tol());
        iteration_ptrtype inner = iteration_ptrtype( iteration_type::New() );
        inner->setMaximumNumberOfIterations( m );
        inner->setRelativePrecision( outer.relativePrecision() );
        inner->setInitialResidual( outer.initialResidual() );

        try
        {
            do
            {
                int k;
                axpy_prod(A, V( i ), u, true);
                //SMultiplication(A, V(i), u);
                M.solve(u, w);
                for (k = 0; k <= i; k++)
                {
                    H(k, i) = ublas::inner_prod( w, V(k) );
                    w.minus_assign( H( k, i )*V( k ) );
                    //SAxpy(-H(k,i),V(k),w,w);
                }

                H(i+1, i) = ublas::norm_2(w);
                V( i+1 ).assign( w/H(i+1,i) );
                //SScale(1./H(i+1,i), w, V(i+1));

                for (k = 0; k < i; k++)
                {
                    rotations[k].apply(H(k,i), H(k+1,i));
                }

                rotations[i] = GivensRotation<value_type>(H(i,i), H(i+1,i));
                rotations[i].apply(H(i,i), H(i+1,i));
                rotations[i].apply(s(i), s(i+1));

                std::cout << "s(" << i << ")=" << s( i ) << "\n";
                std::cout << "s(" << i+1 << ")=" << s( i+1 ) << "\n";
                ++*inner, ++i;

            }
            while (! inner->isFinished( std::abs( s(i) ) ) );
        }
        catch(std::exception const& e)
        {
            throw;
        }

        //#if 0
        // loop on the columns of the matrix
        for(int j = i-1;j >= 0;j--)
        {
            s(j) /= H(j,j);

            value_type s_j = s(j);
            // loop on the lines of the matrix
            for(int k = j-1;k >= 0;k--)
            {
                s(k) -= H(k,j)*s_j;
            }
        }
        //#endif
#if 0
        // loop on the columns of the matrix
        for(int j = i-1;j >= 0;j--)
        {
            for(int k = i-1;k > j;k--)
            {
                s(j) -= H(j,k)*s(k);
            }
            s(j) /= H(j,j);
        }
#endif
        for(int j = 0;j < i;j++)
        {
            x.plus_assign( s( j )*V( j ) );
            //SAxpy(s(j),V(j), x, x);
        }
        //mult(V(Range(0,x.size()-1),Range(0,i-1)), s, x, x);

        axpy_prod(A, x, w, true);
        //SMultiplication(A, x, w);
        w.assign( b - w );
        ///SAxpy(-1,w,b,w);
        M.solve(w, r);
        beta = std::abs( ublas::norm_2(r) );

        ++outer;
    }

    return 1;//outer.error_code();
}
}
#endif /* __GMRES_H */
