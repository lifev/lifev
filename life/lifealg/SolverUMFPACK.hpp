/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-04

  Copyright (C) 2004 EPFL

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
   \file SolverUMFPACK.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-04
 */
#ifndef __SolverUMFPACK_H
#define __SolverUMFPACK_H 1


extern "C"
{
#include <umfpack.h>
};

#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <sparseArray.hpp>
#include <vecUnknown.hpp>

namespace LifeV
{
/*!
  \class SolverUMFPACK
  \brief Interface for the UMFPACK Solver

  UMFPACK is a direct Solver for (un)symmetric problem \f$ A x = b \f$.

  @author Christophe Prud'homme
*/
class SolverUMFPACK
{
public:


    /** @name Typedefs
     */
    //@{

    typedef double value_type;

    typedef boost::numeric::ublas::compressed_matrix<value_type, boost::numeric::ublas::column_major> matrix_type;

    typedef Vector array_type;



    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
       default constructor

       it sets the umfpack print level to the maximum (ie 6)
     */
    SolverUMFPACK();

    SolverUMFPACK( SolverUMFPACK const & umfpackSolver );

    ~SolverUMFPACK();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{

    //! set matrix from CSRMatr
    void setMatrix( const CSRMatr<CSRPatt, value_type>& m );

    void setMatrix( matrix_type const& m );

    //@}

    /** @name  Methods
     */
    //@{


    //! solve A X = B
    /*!

    \param __X  the solution
    \param __B	the right hand side
    \return the number of iterations
    */
    void solve( array_type& __X, array_type const& __B );

    //! report some info about umfpack
    void reportInfo();


    /**
       report status of umfpack

       \param status status integer returned by umfpack routines
     */
    void reportStatus( int );

    //@}

private:

    void prepareSolve();


private:

    class Pimpl;
    boost::shared_ptr<Pimpl> _M_p;

    bool _M_matrix_reset;
    bool _M_matrix_values_reset;

    void *_M_symbolic;
    void *_M_numeric;

    double* _M_Control;
    double* _M_Info;

};
}
#endif /* __SolverUMFPACK_H */

