/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-25

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
   \file darcySolverBase.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-25
 */
#ifndef __DarcySolverBase_H
#define __DarcySolverBase_H 1

#include <boost/signal.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/singleton.hpp>
#include <life/lifecore/factory.hpp>
#include <life/lifefem/bcHandler.hpp>


namespace LifeV
{
enum darcy_unknown_type
{
    DARCY_PRESSURE_GLOBAL,
    DARCY_PRESSURE,
    DARCY_VELOCITY
};

struct DarcyDefaultSource
{
    double operator()( double, double, double, int )
        {
            return 0;
        }
};
class DarcyDefaultInversePermeability
{
public:
    DarcyDefaultInversePermeability()
        :
        _M_id( IdentityMatrix( 3 ) )
        {}
    Matrix const& operator()( Real const&, Real const&, Real const& )
        {
            return _M_id;
        }
private:
    Matrix _M_id;
};
/*!
  \class DarcySolverBase
  \brief base class for all Darcy solvers

  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
*/
class DarcySolverBase
{
public:


    /** @name Typedefs
     */
    //@{
    typedef boost::function<Matrix const& ( Real const&, Real const&, Real const& )> tensor_type;

    typedef boost::function<double( double, double, double, int )> source_type;
    typedef boost::function<double( double, double, double )> pressure_solution_type;
    typedef boost::function<double( double, double, double, double, UInt )> velocity_solution_type;
    typedef boost::signal<void ( std::string const&, darcy_unknown_type, double, double, double )> error_signal_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    DarcySolverBase()
        :
        _M_source( DarcyDefaultSource() ),
        _M_inv_permeability( DarcyDefaultInversePermeability() )
        {}

    DarcySolverBase( DarcySolverBase const& __ds )
        :
        _M_source( __ds._M_source ),
        _M_inv_permeability( __ds._M_inv_permeability )
        {}

    virtual ~DarcySolverBase()
        {}


    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    source_type sourceTerm() const
        {
            return _M_source;
        }

    /**
       provide the inverse of the permeability tensor at point (x,y,z)
     */
    Matrix const& inversePermeability( Real const& __x, Real const& __y,  Real const& __z ) const
        {
            return _M_inv_permeability( __x, __y, __z );
        }

    //@}

    /** @name  Mutators
     */
    //@{


    //! set the source term for the darcy solver
    void setSourceTerm( source_type __s )
        {
            _M_source = __s;
        }

    //! set the boundary conditions
    virtual void setBC( BCHandler const&  ) = 0;

    //! set the diffusion tensor
    void setInversePermeability( tensor_type const& __d )
        {
            _M_inv_permeability = __d;
        }

    //@}

    /** @name  Methods
     */
    //@{


    /**
       setup the darcy solver
    */
    virtual void setup() = 0;

    /**
       solve the linear system for TP
    */
    virtual void solve() = 0;

    //! set observers for error signal
    template<typename Observer>
    void doOnErrorComputation( Observer __s )
        {
            _M_error_signal.connect( __s );
        }

    //! compute L2 error for pressure wrt analytical solution
    virtual void errorL2( darcy_unknown_type, pressure_solution_type ) = 0;

    //! compute L2 error for velocity wrt analytical solution
    virtual void errorL2( velocity_solution_type  ) = 0;

    //@}



protected:

    //! source function
    source_type _M_source;

    //! diffusion type
    tensor_type _M_inv_permeability;


    //! signal at each error
    error_signal_type _M_error_signal;


};
typedef singleton<factory<DarcySolverBase,std::string> > FactoryDarcy;
}
#endif /* __DarcySolverBase_H */
