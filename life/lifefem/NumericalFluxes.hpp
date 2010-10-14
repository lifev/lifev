/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): A. Fumagalli  <alessio.fumagalli@mail.polimi.it>
            M. Kern       <michel.kern@inria.fr>
      Date: 2010-09

 Copyright (C) 2001-2006 EPFL, Politecnico di Milano, INRIA
 Copyright (C) 2006-2010 Politecnico di Milano

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
  @file NumericalFluxes.hpp
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author M. Kern <michel.kern@inria.fr>
  @date 09/2010

  @brief This file contains the numerical fluxes for hyperbolic equations.
*/
#ifndef _numericalFluxes_H_
#define _numericalFluxes_H_

#include <life/lifearray/tab.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

namespace
{

typedef boost::function<LifeV::Vector ( const LifeV::Real&, const LifeV::Real&,
                                        const LifeV::Real&, const LifeV::Real&,
                                        const LifeV::Real& )> vectorFunction;

LifeV::Real functionDotNormal ( const LifeV::Real&            u,
                                const vectorFunction&         function,
                                const LifeV::KN<LifeV::Real>& normal,
                                const LifeV::Real&            t,
                                const LifeV::Real&            x,
                                const LifeV::Real&            y,
                                const LifeV::Real&            z,
                                const LifeV::Real&            plusMinus )
{
    LifeV::Real value(0);
    const LifeV::UInt dimProblem( normal.size() );

    // Compute \sum_{i} physicalFlux(i) * n(i)
    for ( LifeV::UInt nDim(0); nDim < dimProblem; ++nDim )
    {
        value += plusMinus * function(t, x, y, z, u)[nDim] * normal[nDim];
    }

    return value;
}

LifeV::Real absFunctionDotNormal ( const LifeV::Real&            u,
                                   const vectorFunction&         function,
                                   const LifeV::KN<LifeV::Real>& normal,
                                   const LifeV::Real&            t,
                                   const LifeV::Real&            x,
                                   const LifeV::Real&            y,
                                   const LifeV::Real&            z,
                                   const LifeV::Real&            plusMinus )
{
    LifeV::Real value(0);
    const LifeV::UInt dimProblem( normal.size() );

    // Compute \sum_{i} physicalFlux(i) * n(i)
    for ( LifeV::UInt nDim(0); nDim < dimProblem; ++nDim )
    {
        value += function(t, x, y, z, u)[nDim] * normal[nDim];
    }

    return plusMinus * fabs( value );
}

}


// LifeV namespace.
namespace LifeV
{

class AbstractNumericalFlux
{

public:

    typedef boost::function<Vector ( const Real&, const Real&, const Real&,
                                     const Real&, const Real& )>
                                                  vectorFunction;

    typedef boost::function<Real ( const Real& )> scalarFunction;

    AbstractNumericalFlux ( const vectorFunction& physicalFlux,
                            const vectorFunction& firstDerivativePhysicalFlux,
                            const Real&           CFLBrentToll    = 1e-4,
                            const UInt&           CFLBrentMaxIter = 20 );

    virtual ~AbstractNumericalFlux ();

    virtual Real operator() ( const Real&     leftState,
                              const Real&     rightState,
                              const KN<Real>& normal,
                              const Real&     t = 0,
                              const Real&     x = 0,
                              const Real&     y = 0,
                              const Real&     z = 0 ) = 0;

    inline vectorFunction getPhysicalFlux () const
    {
        return M_physicalFlux;
    }

    inline vectorFunction getFirstDerivativePhysicalFlux () const
    {
        return M_firstDerivativePhysicalFlux;
    }

    inline Real getPhysicalFluxDotNormal ( const KN<Real>& normal,
                                           const Real&     t,
                                           const Real&     x,
                                           const Real&     y,
                                           const Real&     z,
                                           const Real&     u )
    {
        return computeFunctionDotNormal ( M_physicalFlux, normal, t, x, y, z, +1 ) ( u );
    }

    inline Real getFirstDerivativePhysicalFluxDotNormal ( const KN<Real>& normal,
                                                          const Real&     t,
                                                          const Real&     x,
                                                          const Real&     y,
                                                          const Real&     z,
                                                          const Real&     u )
    {
        return computeFunctionDotNormal ( M_firstDerivativePhysicalFlux, normal, t, x, y, z, +1 ) ( u );
    }

    Real getNormInfty ( const Real&     leftState,
                        const Real&     rightState,
                        const KN<Real>& normal,
                        const Real&     t = 0,
                        const Real&     x = 0,
                        const Real&     y = 0,
                        const Real&     z = 0 );

protected:

    scalarFunction computeFunctionDotNormal ( const vectorFunction& function,
                                              const KN<Real>&       normal,
                                              const Real&           t,
                                              const Real&           x,
                                              const Real&           y,
                                              const Real&           z,
                                              const Real&           plusMinus );

    scalarFunction computeAbsFunctionDotNormal ( const vectorFunction& function,
                                                 const KN<Real>&       normal,
                                                 const Real&           t,
                                                 const Real&           x,
                                                 const Real&           y,
                                                 const Real&           z,
                                                 const Real&           plusMinus );

    vectorFunction M_physicalFlux;

    vectorFunction M_firstDerivativePhysicalFlux;

    // Tollerance for the brent algorithm for computing the CFL condition
    Real M_CFLBrentToll;

    // Maxiter for the brent algorithm for computing the CFL condition
    UInt M_CFLBrentMaxIter;

}; // AbstractNumericalFlux

// ######################################################################### //

class GodunovNumericalFlux : public AbstractNumericalFlux
{

public:

    typedef AbstractNumericalFlux::vectorFunction vectorFunction;
    typedef AbstractNumericalFlux::scalarFunction scalarFunction;

    GodunovNumericalFlux ( const vectorFunction& physicalFlux,
                           const vectorFunction& firstDerivativePhysicalFlux,
                           const Real&           CFLBrentToll    = 1e-4,
                           const UInt&           CFLBrentMaxIter = 20,
                           const Real&           brentToll       = 1e-4,
                           const UInt&           brentMaxIter    = 20 );

    Real operator() ( const Real&     leftState,
                      const Real&     rightState,
                      const KN<Real>& normal,
                      const Real&     t = 0,
                      const Real&     x = 0,
                      const Real&     y = 0,
                      const Real&     z = 0 );
protected:

    // Tollerance for the brent algorithm
    Real M_brentToll;

    // Maxiter for the brent algorithm
    UInt M_brentMaxIter;

}; // GodunovNumericalFlux

}

#endif //_numericalFluxes_H_
