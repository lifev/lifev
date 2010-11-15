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
#include <life/lifefem/FESpace.hpp>
#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifealg/brent.hpp>


namespace
{

typedef boost::function<LifeV::Vector ( const LifeV::Real&, const LifeV::Real&,
                                        const LifeV::Real&, const LifeV::Real&,
                                        const std::vector<LifeV::Real>& )>
                                      vectorFunction;

LifeV::Real functionDotNormal ( const LifeV::Real&              u,
                                const vectorFunction&           function,
                                const LifeV::KN<LifeV::Real>&   normal,
                                const LifeV::Real&              t,
                                const LifeV::Real&              x,
                                const LifeV::Real&              y,
                                const LifeV::Real&              z,
                                const LifeV::Real&              plusMinus,
                                const std::vector<LifeV::Real>& fieldsValues )
{
    LifeV::Real value(0);
    const LifeV::UInt dimProblem( normal.size() );
    std::vector<LifeV::Real> uAndFields( 1, 0. );

    uAndFields[0] = u;
    uAndFields.insert ( uAndFields.begin() + 1, fieldsValues.begin(), fieldsValues.end() );

    // Compute \sum_{i} physicalFlux(i) * n(i)
    for ( LifeV::UInt nDim(0); nDim < dimProblem; ++nDim )
    {
        value += plusMinus * function( t, x, y, z, uAndFields )[nDim] * normal[nDim];
    }

    return value;
}

LifeV::Real absFunctionDotNormal ( const LifeV::Real&              u,
                                   const vectorFunction&           function,
                                   const LifeV::KN<LifeV::Real>&   normal,
                                   const LifeV::Real&              t,
                                   const LifeV::Real&              x,
                                   const LifeV::Real&              y,
                                   const LifeV::Real&              z,
                                   const LifeV::Real&              plusMinus,
                                   const std::vector<LifeV::Real>& fieldsValues)
{
    LifeV::Real value(0);
    const LifeV::UInt dimProblem( normal.size() );
    std::vector<LifeV::Real> uAndFields( fieldsValues.size() + 1, 0. );

    uAndFields[0] = u;
    uAndFields.insert ( uAndFields.begin() + 1, fieldsValues.begin(), fieldsValues.end() );

    // Compute \sum_{i} physicalFlux(i) * n(i)
    for ( LifeV::UInt nDim(0); nDim < dimProblem; ++nDim )
    {
        value += function( t, x, y, z, uAndFields )[nDim] * normal[nDim];
    }

    return plusMinus * fabs( value );
}

}


// LifeV namespace.
// We suppose that u[0] contains the non-linear value
namespace LifeV
{

template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class AbstractNumericalFlux
{

public:

    // Policies.
    //! @name Policies
    //@{

    typedef boost::function<Vector ( const Real&, const Real&,
                                     const Real&, const Real&,
                                     const std::vector<Real>& )>
                                                  vectorFunction;

    typedef boost::function<Real ( const Real& )> scalarFunction;

    typedef typename SolverType::vector_type      vector_type;
    typedef boost::shared_ptr<vector_type>        vector_ptrtype;

    typedef GetPot                                dataFile;

    //@}

    // Constructors and destructor.
    //! @name Constructors and destructor
    //@{

    AbstractNumericalFlux ( const vectorFunction&            physicalFlux,
                            const vectorFunction&            firstDerivativePhysicalFlux,
                            const FESpace<Mesh, EpetraMap>&  fESpace,
                            const dataFile&                  data,
                            const std::string&               section = "numerical_flux/");


    virtual ~AbstractNumericalFlux ();

    //@}

    // Add one field
    inline void setField ( const vector_ptrtype & field ) { M_fields.push_back( &field ); };

    virtual Real operator() ( const Real&     leftState,
                              const Real&     rightState,
                              const KN<Real>& normal,
                              const UInt&     iElem,
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
                                           const UInt&     iElem,
                                           const Real&     t,
                                           const Real&     x,
                                           const Real&     y,
                                           const Real&     z,
                                           const Real&     u )
    {
        return computeFunctionDotNormal ( M_physicalFlux, normal, iElem, t, x, y, z, +1 ) ( u );
    }

    inline Real getFirstDerivativePhysicalFluxDotNormal ( const KN<Real>& normal,
                                                          const UInt&     iElem,
                                                          const Real&     t,
                                                          const Real&     x,
                                                          const Real&     y,
                                                          const Real&     z,
                                                          const Real&     u )
    {
        return computeFunctionDotNormal ( M_firstDerivativePhysicalFlux, normal, iElem, t, x, y, z, +1 ) ( u );
    }

    Real getNormInfty ( const Real&     leftState,
                        const Real&     rightState,
                        const KN<Real>& normal,
                        const UInt&     iElem,
                        const Real&     t = 0,
                        const Real&     x = 0,
                        const Real&     y = 0,
                        const Real&     z = 0 );

protected:

    scalarFunction computeFunctionDotNormal ( const vectorFunction& function,
                                              const KN<Real>&       normal,
                                              const UInt&           iElem,
                                              const Real&           t,
                                              const Real&           x,
                                              const Real&           y,
                                              const Real&           z,
                                              const Real&           plusMinus );

    scalarFunction computeAbsFunctionDotNormal ( const vectorFunction& function,
                                                 const KN<Real>&       normal,
                                                 const UInt&           iElem,
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

    // Finite element space
    const FESpace<Mesh, EpetraMap>&     M_fESpace;

    // Vector of pointers for the dependences of the permeability to an external field.
    std::vector< const vector_ptrtype* > M_fields;

}; // AbstractNumericalFlux

// Constructor of the abstract class
template < typename Mesh, typename SolverType >
AbstractNumericalFlux<Mesh, SolverType>::
AbstractNumericalFlux ( const vectorFunction&            physicalFlux,
                        const vectorFunction&            firstDerivativePhysicalFlux,
                        const FESpace<Mesh, EpetraMap>&  fESpace,
                        const dataFile&                  data,
                        const std::string&               section ):
    M_physicalFlux                ( physicalFlux ),
    M_firstDerivativePhysicalFlux ( firstDerivativePhysicalFlux ),
    M_fESpace                     ( fESpace ),
    M_fields                      ( std::vector< const vector_ptrtype* >(0) ),
    M_CFLBrentToll                ( data( ( section + "CFL/brent_toll" ).data(), 1e-4 ) ),
    M_CFLBrentMaxIter             ( data( ( section + "CFL/brent_maxIter" ).data(), 20 ) )
{

    CONSTRUCTOR( "AbstractNumericalFlux" );

} // constructor

// Destructor of the abstract class
template < typename Mesh, typename SolverType >
AbstractNumericalFlux<Mesh, SolverType>::
~AbstractNumericalFlux ()
{

    DESTRUCTOR( "AbstractNumericalFlux" );

} // destructor

template < typename Mesh, typename SolverType >
Real
AbstractNumericalFlux<Mesh, SolverType>::
getNormInfty ( const Real& leftState, const Real& rightState, const KN<Real>& normal, const UInt& iElem,
               const Real& t, const Real& x, const Real& y, const Real& z )
{

    std::vector<Real> values ( M_fields.size() * M_fESpace.fieldDim(), 0 );
    UInt totalDofsPresent( M_fESpace.dof().numTotalDof() );
    UInt fieldDim( M_fESpace.fieldDim() );
    scalarFunction absFunctionDotNormalBound;

    for ( UInt i(0); i < M_fields.size(); ++i )
    {
        for ( UInt iComponent(0); iComponent < fieldDim; ++iComponent )
        {
            values[ i*fieldDim + iComponent ] = (*( *(M_fields)[i] ))[ iComponent*totalDofsPresent + M_fESpace.dof().localToGlobal( iElem, 1)  ];
        }

    }

    absFunctionDotNormalBound = boost::bind ( &absFunctionDotNormal, _1,
                                              M_firstDerivativePhysicalFlux,
                                              normal, t, x, y, z, -1, values );

    // Compute the maximum of the absFunctionDotNormal
    const Real maxValue = brent( absFunctionDotNormalBound, leftState, rightState,
                                 M_CFLBrentToll, M_CFLBrentMaxIter );

    return - absFunctionDotNormalBound ( maxValue );

} // getNormInfty

// Create the fluxDotNormal function
template < typename Mesh, typename SolverType >
boost::function<Real ( const Real& )>
AbstractNumericalFlux<Mesh, SolverType>::
computeFunctionDotNormal ( const vectorFunction& function, const KN<Real>& normal, const UInt& iElem,
                           const Real& t, const Real& x, const Real& y, const Real& z, const Real& plusMinus )
{

    std::vector<Real> values ( M_fields.size() * M_fESpace.fieldDim(), 0 );
    UInt totalDofsPresent( M_fESpace.dof().numTotalDof() );
    UInt fieldDim( M_fESpace.fieldDim() );
    scalarFunction functionDotNormalBound;

    for ( UInt i(0); i < M_fields.size(); ++i )
    {
        for ( UInt iComponent(0); iComponent < fieldDim; ++iComponent )
        {
            values[ i*fieldDim + iComponent ] = (*( *(M_fields)[i] ))[ iComponent*totalDofsPresent + M_fESpace.dof().localToGlobal( iElem, 1)  ];
        }

    }

    // Bind the anonymousnamespace::fluxDotNormal function with known quantities.
    functionDotNormalBound = boost::bind ( &functionDotNormal, _1, function,
                                           normal, t, x, y, z, plusMinus, values );

    return functionDotNormalBound;

} // computeFluxDotNormal

// ######################################################################### //

template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class GodunovNumericalFlux : public AbstractNumericalFlux<Mesh, SolverType>
{

public:

    typedef typename AbstractNumericalFlux<Mesh, SolverType>::vectorFunction vectorFunction;
    typedef typename AbstractNumericalFlux<Mesh, SolverType>::scalarFunction scalarFunction;
    typedef typename AbstractNumericalFlux<Mesh, SolverType>::dataFile dataFile;

    GodunovNumericalFlux ( const vectorFunction&           physicalFlux,
                           const vectorFunction&           firstDerivativePhysicalFlux,
                           const FESpace<Mesh, EpetraMap>& fESpace,
                           const dataFile&                 data,
                           const std::string&              section = "numerical_flux/" );

    virtual ~GodunovNumericalFlux ();

    virtual Real operator() ( const Real&     leftState,
                              const Real&     rightState,
                              const KN<Real>& normal,
                              const UInt&     iElem,
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

template < typename Mesh, typename SolverType >
GodunovNumericalFlux<Mesh, SolverType>::
GodunovNumericalFlux ( const vectorFunction&            physicalFlux,
                       const vectorFunction&            firstDerivativePhysicalFlux,
                       const FESpace<Mesh, EpetraMap>&  fESpace,
                       const dataFile&                  data,
                       const std::string&               section ):
    AbstractNumericalFlux<Mesh, SolverType>::AbstractNumericalFlux  ( physicalFlux,
                                                                      firstDerivativePhysicalFlux,
                                                                      fESpace,
                                                                      data,
                                                                      section ),
    M_brentToll                                   ( data( ( section + "godunov/brent_toll" ).data(), 1e-4 ) ),
    M_brentMaxIter                                ( data( ( section + "godunov/brent_maxIter" ).data(), 20) )
{

    CONSTRUCTOR( "GodunovNumericalFlux" );

} // constructor

// Destructor of the Godunov class
template < typename Mesh, typename SolverType >
GodunovNumericalFlux<Mesh, SolverType>::
~GodunovNumericalFlux ()
{

    DESTRUCTOR( "GodunovNumericalFlux" );

} // destructor

template < typename Mesh, typename SolverType >
Real
GodunovNumericalFlux<Mesh, SolverType>::
operator() ( const Real& leftState, const Real& rightState, const KN<Real>& normal,
             const UInt& iElem,
             const Real& t, const Real& x, const Real& y, const Real& z )
{

    // The value of the flux
    Real fluxValue( static_cast<Real>(0) );

    // The argmin or argmax of flux dot normal
    Real minMax( static_cast<Real>(0) );

    // The normal flux function
    scalarFunction normalFlux;

    if ( rightState > leftState )
    {
        // Create the function f \cdot n
        normalFlux = this->computeFunctionDotNormal( this->M_physicalFlux, normal, iElem, t, x, y, z, +1 );

        // Compute the argmin f \cdot n
        minMax = brent( normalFlux, leftState, rightState, M_brentToll, M_brentMaxIter );

        // Compute the flux value
        fluxValue = normalFlux( minMax );
    }
    else
    {
        // Create the function - f \cdot n
        normalFlux = this->computeFunctionDotNormal( this->M_physicalFlux, normal, iElem, t, x, y, z, -1 );

        // Compute the argmin - f \cdot n
        minMax = brent( normalFlux, leftState, rightState, M_brentToll, M_brentMaxIter );

        // Compute the flux value
        fluxValue = - normalFlux( minMax );
    }

    return fluxValue;

} // operator()

}

#endif //_numericalFluxes_H_
