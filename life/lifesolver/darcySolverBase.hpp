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

#include <lifeV.hpp>
#include <singleton.hpp>
#include <factory.hpp>
#include <bcHandler.hpp>


namespace LifeV
{
typedef enum darcy_unknown_type
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
struct DarcyDefaultDiffusion
{
    Matrix operator()( double, double, double )
        {
            return ScalarMatrix( 3, 3, 1.0 );
        }
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
    typedef boost::function<Matrix ( double, double, double )> diffusion_type;
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
        _M_diffusion( DarcyDefaultDiffusion() )
        {}

    DarcySolverBase( DarcySolverBase const& __ds )
        :
        _M_source( __ds._M_source ),
        _M_diffusion( __ds._M_diffusion )
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

    diffusion_type diffusion() const
        {
            return _M_diffusion;
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
    void setDiffusion( diffusion_type __d )
        {
            _M_diffusion = __d;
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
    diffusion_type _M_diffusion;


    //! signal at each error
    error_signal_type _M_error_signal;


};
typedef singleton<factory<DarcySolverBase,std::string> > FactoryDarcy;
}
#endif /* __DarcySolverBase_H */
