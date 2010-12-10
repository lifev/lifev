//@HEADER
/*
*******************************************************************************
    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University
    This file is part of LifeV.
    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with LifeV. If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************
*/
//@HEADER
/*!
 * @file
 * @brief Numerical fluxes for hyperbolic scalar equations.
 *
 *
 * @date 30-09-2010
 *
 * @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
 * @author Michel Kern       <michel.kern@inria.fr>
 *
 * @contributor
 *
 * @mantainer Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
 *
 */

#ifndef _NUMERICALFLUXES_H_
#define _NUMERICALFLUXES_H_ 1

#include <boost/bind.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifealg/brent.hpp>

namespace
{

using namespace LifeV;

typedef boost::function<Vector ( const Real&, const Real&,
                                 const Real&, const Real&,
                                 const std::vector<Real>& )>
vectorFunction;

// Compute plus or minus the function dot normal vector
Real functionDotNormal ( const Real&              unknown,
                         const vectorFunction&    function,
                         const KN<Real>&          normal,
                         const Real&              t,
                         const Real&              x,
                         const Real&              y,
                         const Real&              z,
                         const Real&              plusMinus,
                         const std::vector<Real>& fieldsValues )
{
    Real valueFunctionDotNormal(0);
    const UInt problemDimension( normal.size() );
    std::vector<Real> unknownAndFields( 1, 0. );

    // Add to the vector unknownAndFields the values of the unknown then the value of the external fields.
    unknownAndFields[0] = unknown;
    unknownAndFields.insert ( unknownAndFields.begin() + 1,
                              fieldsValues.begin(), fieldsValues.end() );

    // Compute  \sum_{i} physicalFlux(i) * n(i)
    for ( UInt nDim(0); nDim < problemDimension; ++nDim )
    {
        valueFunctionDotNormal += plusMinus * function( t, x, y, z, unknownAndFields )[nDim] * normal[nDim];
    }

    return valueFunctionDotNormal;
}

// Compute plus or minus the absolute value of function dot normal vector
Real absFunctionDotNormal ( const Real&              unknown,
                            const vectorFunction&    function,
                            const KN<Real>&          normal,
                            const Real&              t,
                            const Real&              x,
                            const Real&              y,
                            const Real&              z,
                            const Real&              plusMinus,
                            const std::vector<Real>& fieldsValues)
{
    Real valueFunctionDotNormal(0);
    const UInt problemDimension( normal.size() );
    std::vector<Real> unknownAndFields( fieldsValues.size() + 1, 0. );

    // Add to the vector unknownAndFields the values of the unknown then the value of the external fields.
    unknownAndFields[0] = unknown;
    unknownAndFields.insert ( unknownAndFields.begin() + 1,
                              fieldsValues.begin(), fieldsValues.end() );

    // Compute \sum_{i} physicalFlux(i) * n(i)
    for ( LifeV::UInt nDim(0); nDim < problemDimension; ++nDim )
    {
        valueFunctionDotNormal += function( t, x, y, z, unknownAndFields )[nDim] * normal[nDim];
    }

    return plusMinus * fabs( valueFunctionDotNormal );
}

}


// LifeV namespace.
namespace LifeV
{
//! AbstractNumericalFlux Gives a common interface for hyperbolic's flux function.
/*!
  @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author Michel Kern       <michel.kern@inria.fr>

  This class gives a common interface for hyperbolic's flux functions \f$ \mathbf{F}( u ) \f$, in particoular it computes
  <ol>
  <li> the flux function dot product the outward unit normal \f$ \mathbf{n} \f$;</li>
  <li> the first derivative of the flux function, respect to \f$ u \f$, dot product the outward unit normal;</li>
  <li> the maximum value of \f$ \hat{\mathbf{F}} \cdot \mathbf{n} \f$ between to adjacent elements; </li>
  <li> the value of \f$ \Vert \mathbf{F}^\prime \cdot \mathbf{n}_{e, K} \Vert_{L^\infty(a_0, b_0) } \f$ between to adjacent elements.</li>
  </ol>
  The class can handle the dependeces of the flux function \f$ \mathbf{F} \f$ of external vector or scalar fields.
  @note In the implementation of the physical flux \f$ \mathbf{F} \f$ we suppose that the first parameter of the vector is the unknown.
        See the test case for an example.
*/
template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class AbstractNumericalFlux
{

public:

    //! @name Public Types
    //@{

    typedef boost::function<Vector ( const Real&, const Real&,
                                     const Real&, const Real&,
                                     const std::vector<Real>& )>
    vectorFunction_Type;

    typedef boost::function< Real ( const Real& ) > scalarFunction_Type;

    typedef typename SolverType::vector_type        vector_Type;
    typedef boost::shared_ptr< vector_Type >        vectorPtr_Type;

    typedef GetPot                                  dataFile_Type;
    typedef KN< Real >                              normal_Type;

    //@}

    // Constructors & destructor.
    //! @name Constructors and destructor
    //@{

    //! Constructor for the class
    /*!
      @param physicalFlux Physical flux for the problem.
      @param firstDerivativePhysicalFlux First derivative, respect the the unknown, of the physical flux.
      @param fESpace Finite element space of the hyperbolic problem.
      @param data Data for the problem.
      @param section Section for read the data from GetPot file.
    */
    AbstractNumericalFlux ( const vectorFunction_Type&       physicalFlux,
                            const vectorFunction_Type&       firstDerivativePhysicalFlux,
                            const FESpace<Mesh, EpetraMap>&  fESpace,
                            const dataFile_Type&             data,
                            const std::string&               section = "numerical_flux/");

    //! Virtual destructor.
    virtual ~AbstractNumericalFlux ();

    //@}

    //! @name Operators
    //@{

    //! Computes the face contribution of the flux.
    /*!
      Given a face \f$ e \in \partial K \f$ it evaluates
      \f[
      \max \hat{\mathbf{F}} \cdot \mathbf{n}
      \f]
      between to elements sharing the face.
      @param leftState Left value of the unknown respect to the face.
      @param rightValue Right value of the unknown respect to the face.
      @param normal Normal of the face.
      @param iElem The ID of the current element in the mesh.
      @param t Current time.
      @param x Abscissa.
      @param y Ordinate.
      @param z Quota.
      @note We assume left and right side of \f$ e \f$ is given by: the normal direction goes
      from the left side to the right side.
    */
    virtual Real operator() ( const Real&        leftState,
                              const Real&        rightState,
                              const normal_Type& normal,
                              const UInt&        iElem,
                              const Real&        t = 0,
                              const Real&        x = 0,
                              const Real&        y = 0,
                              const Real&        z = 0 ) const = 0;
    //@}


    //! @name Set Methods
    //@{

    //! Add one external field.
    /*!
      Add one extra field for the dependece from \f$ \mathbf{F} \f$.
      @param field The filed to be added.
    */
    inline void setField ( const vectorPtr_Type & field )
    {
        M_fields.push_back( &field );
    }

    //@}

    //! @name Get Methods
    //@{

    //! Return the physical flux.
    /*!
      @return The physical flux.
    */
    inline vectorFunction_Type getPhysicalFlux () const
    {
        return M_physicalFlux;
    }

    //! Return the first derivative, respect to the unknown, of the physical flux.
    /*!
      @return The first derivative of the physical flux.
    */
    inline vectorFunction_Type getFirstDerivativePhysicalFlux () const
    {
        return M_firstDerivativePhysicalFlux;
    }

    //! Evaluate the flux dot normal in a given point of a face.
    /*!
      @param normal The normal of the face.
      @param iElem The ID of the current element in the mesh.
      @param t Current time.
      @param x Abscissa.
      @param y Ordinate.
      @param z Quota.
      @param unknown The value of the unknown.
      @return The value of \f$ \mathbf{ \hat{ F } } \cdot \mathbf{ n }\f$ in the point \f$ (x, y, z, t, u) \f$.
    */
    inline Real getPhysicalFluxDotNormal ( const normal_Type& normal,
                                           const UInt&        iElem,
                                           const Real&        t,
                                           const Real&        x,
                                           const Real&        y,
                                           const Real&        z,
                                           const Real&        unknown ) const
    {
        return computeFunctionDotNormal ( M_physicalFlux, normal, iElem, t, x, y, z, +1 ) ( unknown );
    }

    //! Evaluate the first derivative of the flux dot normal in a given point of a face.
    /*!
      @param normal The normal of the face.
      @param iElem The ID of the current element in the mesh.
      @param t Current time.
      @param x Abscissa.
      @param y Ordinate.
      @param z Quota.
      @param unknown The value of the unknown.
      @return The value of \f$ \mathbf{ \hat{ F^\prime } } \cdot \mathbf{ n }\f$ in the point \f$ (x, y, z, t, u) \f$.
    */
    inline Real getFirstDerivativePhysicalFluxDotNormal ( const normal_Type& normal,
                                                          const UInt&        iElem,
                                                          const Real&        t,
                                                          const Real&        x,
                                                          const Real&        y,
                                                          const Real&        z,
                                                          const Real&        unknown ) const
    {
        return computeFunctionDotNormal ( M_firstDerivativePhysicalFlux, normal, iElem, t, x, y, z, +1 ) ( unknown );
    }

    //! Computes the local infinity norm of the first derivative of the flux dot normal.
    /*!
      Given a face \f$ e \in \partial K \f$ it evaluates
      \f[
      \displaystyle \max_{e \in \partial K^- \cap \partial K^+} \vert \hat{\mathbf{F^\prime}} \cdot \mathbf{n} \vert
      \f]
      between to elements sharing the face.
      @param leftState Left value of the unknown respect to the face.
      @param rightValue Right value of the unknown respect to the face.
      @param normal Normal of the face.
      @param iElem The ID of the current element in the mesh.
      @param t Current time.
      @param x Abscissa.
      @param y Ordinate.
      @param z Quota.
      @note We assume left and right side of \f$ e \f$ is given by: the normal direction goes
      from the left side to the right side.
    */
    Real getNormInfty ( const Real&        leftState,
                        const Real&        rightState,
                        const normal_Type& normal,
                        const UInt&        iElem,
                        const Real&        t = 0,
                        const Real&        x = 0,
                        const Real&        y = 0,
                        const Real&        z = 0 ) const;

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Return a scalar function from a general vector function dot normal in a given point of a face.
    /*!
      @param function Function to be reduced.
      @param normal The normal of the face.
      @param iElem The ID of the current element in the mesh.
      @param t Current time.
      @param x Abscissa.
      @param y Ordinate.
      @param z Quota.
      @param plusMinus
      @return The function \f$ g(u) = \pm \mathbf{f}(x,y,z,t,u) \cdot \mathbf{n} \f$.
    */
    scalarFunction_Type computeFunctionDotNormal ( const vectorFunction_Type& function,
                                                   const normal_Type&         normal,
                                                   const UInt&                iElem,
                                                   const Real&                t,
                                                   const Real&                x,
                                                   const Real&                y,
                                                   const Real&                z,
                                                   const Real&                plusMinus ) const;

    //! Return a scalar function from the absolute value of a general vector function dot normal in a given point of a face.
    /*!
      @param function Function to be reduced.
      @param normal The normal of the face.
      @param iElem The ID of the current element in the mesh.
      @param t Current time.
      @param x Abscissa.
      @param y Ordinate.
      @param z Quota.
      @param plusMinus
      @return The function \f$ g(u) = \pm \vert \mathbf{f}(x,y,z,t,u) \cdot \mathbf{n} \vert \f$.
    */
    scalarFunction_Type computeAbsFunctionDotNormal ( const vectorFunction_Type& function,
                                                      const normal_Type&         normal,
                                                      const UInt&                iElem,
                                                      const Real&                t,
                                                      const Real&                x,
                                                      const Real&                y,
                                                      const Real&                z,
                                                      const Real&                plusMinus ) const;

    //@}

    //! Physical flux function.
    vectorFunction_Type                  M_physicalFlux;

    //! First derivative, respect to unknown, of physical flux function.
    vectorFunction_Type                  M_firstDerivativePhysicalFlux;

    //! Tollerance for the Brent algorithm for computing the CFL condition.
    Real                                 M_CFLBrentToll;

    //! Maximum of iterations for the Brent algorithm for computing the CFL condition.
    UInt                                 M_CFLBrentMaxIter;

    //! Finite element space of the hyperbolic solver.
    const FESpace<Mesh, EpetraMap>&      M_fESpace;

    //! Vector of pointers for the dependences of the permeability to external vector fields.
    std::vector< const vectorPtr_Type* > M_fields;

}; // AbstractNumericalFlux

// ===================================================
// Constructors & Destructor
// ===================================================

// Constructor of the class
template < typename Mesh, typename SolverType >
AbstractNumericalFlux<Mesh, SolverType>::
AbstractNumericalFlux ( const vectorFunction_Type&       physicalFlux,
                        const vectorFunction_Type&       firstDerivativePhysicalFlux,
                        const FESpace<Mesh, EpetraMap>&  fESpace,
                        const dataFile_Type&             data,
                        const std::string&               section ):
        M_physicalFlux                ( physicalFlux ),
        M_firstDerivativePhysicalFlux ( firstDerivativePhysicalFlux ),
        M_CFLBrentToll                ( data( ( section + "CFL/brent_toll" ).data(), 1e-4 ) ),
        M_CFLBrentMaxIter             ( data( ( section + "CFL/brent_maxIter" ).data(), 20 ) ),
        M_fESpace                     ( fESpace ),
        M_fields                      ( std::vector< const vectorPtr_Type* >(0) )
{

} // Constructor

// Destructor of the class
template < typename Mesh, typename SolverType >
AbstractNumericalFlux<Mesh, SolverType>::
~AbstractNumericalFlux ()
{

} // Destructor

// ===================================================
// Methods
// ===================================================

// Return the norm infinity
template< typename Mesh, typename SolverType >
Real
AbstractNumericalFlux<Mesh, SolverType>::
getNormInfty ( const Real& leftState, const Real& rightState, const normal_Type& normal, const UInt& iElem,
               const Real& t, const Real& x, const Real& y, const Real& z ) const
{

    std::vector<Real> values ( M_fields.size() * M_fESpace.fieldDim(), 0 );
    const UInt totalDofsPresent( M_fESpace.dof().numTotalDof() );
    const UInt fieldDim( M_fESpace.fieldDim() );
    scalarFunction_Type absFunctionDotNormalBound;

    // Takes the value of all the external fields in the current element.
    for ( UInt i(0); i < M_fields.size(); ++i )
    {
        // Select if the external field is a scalar or vector field
        for ( UInt iComponent(0); iComponent < fieldDim; ++iComponent )
        {
            values[ i*fieldDim + iComponent ] = (*( *(M_fields)[i] ))[ iComponent*totalDofsPresent + M_fESpace.dof().localToGlobal( iElem, 1)  ];
        }

    }

    // Bind the function which depends on all the parameter to obatin a scalar function, which depend just on the unknown
    absFunctionDotNormalBound = boost::bind ( &absFunctionDotNormal, _1,
                                              M_firstDerivativePhysicalFlux,
                                              normal, t, x, y, z, -1, values );

    // Compute the minumum of minus absFunctionDotNormal
    const Real maxValue = brent( absFunctionDotNormalBound, leftState, rightState,
                                 M_CFLBrentToll, M_CFLBrentMaxIter );

    // Return minus the value
    return - absFunctionDotNormalBound ( maxValue );

} // getNormInfty

// Create the fluxDotNormal function
template < typename Mesh, typename SolverType >
boost::function< Real ( const Real& ) >
AbstractNumericalFlux< Mesh, SolverType >::
computeFunctionDotNormal ( const vectorFunction_Type& function, const normal_Type& normal, const UInt& iElem,
                           const Real& t, const Real& x, const Real& y, const Real& z, const Real& plusMinus ) const
{

    std::vector<Real> values ( M_fields.size() * M_fESpace.fieldDim(), 0 );
    const UInt totalDofsPresent( M_fESpace.dof().numTotalDof() );
    const UInt fieldDim( M_fESpace.fieldDim() );
    scalarFunction_Type functionDotNormalBound;

    // Takes the value of all the external fields in the current element.
    for ( UInt i(0); i < M_fields.size(); ++i )
    {
        // Select if the external field is a scalar or vector field
        for ( UInt iComponent(0); iComponent < fieldDim; ++iComponent )
        {
            values[ i*fieldDim + iComponent ] = (*( *(M_fields)[i] ))[ iComponent*totalDofsPresent + M_fESpace.dof().localToGlobal( iElem, 1)  ];
        }

    }

    // Bind the fluxDotNormal function with known quantities.
    functionDotNormalBound = boost::bind ( &functionDotNormal, _1, function,
                                           normal, t, x, y, z, plusMinus, values );

    return functionDotNormalBound;

} // computeFluxDotNormal

// ######################################################################### //

//! GodunovNumericalFlux Gives an implementation for Godunov solver for hyperbolic's flux function.
/*!
  @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author Michel Kern       <michel.kern@inria.fr>

  This class gives an implementation for Godunov solver for hyperbolic's flux functions \f$ \mathbf{F}( u ) \f$. In particular it implements
  <ol>
  <li> the flux function dot product the outward unit normal \f$ \mathbf{n} \f$;</li>
  <li> the first derivative of the flux function, respect to \f$ u \f$, dot product the outward unit normal;</li>
  <li> the maximum value of \f$ \hat{\mathbf{F}} \cdot \mathbf{n} \f$ between to adjacent elements computed using Godunov solver
  \f[
  \hat{\mathbf{F}} \cdot \mathbf{n} (a,b)=
    \left\{
      \begin{array}{l l}
        \displaystyle \min_{a \leq u \leq b} \mathbf{F} \cdot \mathbf{n} ( u ) & \mathrm{if} \quad a \leq b\,, \\
        \displaystyle \max_{b \leq u \leq a} \mathbf{F} \cdot \mathbf{n} ( u ) & \mathrm{otherwise}\,.
      \end{array}
    \right.
  \f]
  </li>
  <li> the value of \f$ \Vert \mathbf{F}^\prime \cdot \mathbf{n}_{e, K} \Vert_{L^\infty(a_0, b_0) } \f$ between to adjacent elements.</li>
  </ol>
  The class can handle the dependeces of the flux function \f$ \mathbf{F} \f$ of external vector or scalar fields.
  @note In the implementation of the physical flux \f$ \mathbf{F} \f$ we suppose that the first parameter of the vector is the unknown.
        See the test case for an example.
*/
template< typename Mesh,
typename SolverType = LifeV::SolverTrilinos >
class GodunovNumericalFlux : public AbstractNumericalFlux<Mesh, SolverType>
{

public:

    //! @name Public Types
    //@{

    typedef typename AbstractNumericalFlux<Mesh, SolverType>::vectorFunction_Type vectorFunction_Type;
    typedef typename AbstractNumericalFlux<Mesh, SolverType>::scalarFunction_Type scalarFunction_Type;
    typedef typename AbstractNumericalFlux<Mesh, SolverType>::dataFile_Type       dataFile_Type;
    typedef typename AbstractNumericalFlux<Mesh, SolverType>::normal_Type         normal_Type;

    //@}

    // Constructors & destructor.
    //! @name Constructors and destructor
    //@{

    //! Constructor for the class
    /*!
      @param physicalFlux Physical flux for the problem.
      @param firstDerivativePhysicalFlux First derivative, respect the the unknown, of the physical flux.
      @param fESpace Finite element space of the hyperbolic problem.
      @param data Data for the problem.
      @param section Section for read the data from GetPot file.
    */
    GodunovNumericalFlux ( const vectorFunction_Type&      physicalFlux,
                           const vectorFunction_Type&      firstDerivativePhysicalFlux,
                           const FESpace<Mesh, EpetraMap>& fESpace,
                           const dataFile_Type&            data,
                           const std::string&              section = "numerical_flux/" );

    //! Virtual destructor
    virtual ~GodunovNumericalFlux ();

    //@}

    //! @name Operators
    //@{

    //! Computes the face contribution of the flux.
    /*!
      Given a face \f$ e \in \partial K \f$ it evaluates
      \f[
      \max \hat{\mathbf{F}} \cdot \mathbf{n}
      \f]
      between to elements sharing the face, using Godunov flux.
      @param leftState Left value of the unknown respect to the face.
      @param rightValue Right value of the unknown respect to the face.
      @param normal Normal of the face.
      @param iElem The ID of the current element in the mesh.
      @param t Current time.
      @param x Abscissa.
      @param y Ordinate.
      @param z Quota.
      @note We assume left and right side of \f$ e \f$ is given by: the normal direction goes
      from the left side to the right side.
    */
    virtual Real operator() ( const Real&        leftState,
                              const Real&        rightState,
                              const normal_Type& normal,
                              const UInt&        iElem,
                              const Real&        t = 0,
                              const Real&        x = 0,
                              const Real&        y = 0,
                              const Real&        z = 0 ) const;

    //@}

protected:

    //! Tollerance for the Brent algorithm used for Godunov flux.
    Real M_brentToll;

    //! Maximum of iteration for the Brent algorithm used for Godunov flux.
    UInt M_brentMaxIter;

}; // GodunovNumericalFlux

// ===================================================
// Constructors & Destructor
// ===================================================

// Constructor of the class
template < typename Mesh, typename SolverType >
GodunovNumericalFlux<Mesh, SolverType>::
GodunovNumericalFlux ( const vectorFunction_Type&       physicalFlux,
                       const vectorFunction_Type&       firstDerivativePhysicalFlux,
                       const FESpace<Mesh, EpetraMap>&  fESpace,
                       const dataFile_Type&             data,
                       const std::string&               section ):
    AbstractNumericalFlux<Mesh, SolverType>::AbstractNumericalFlux  ( physicalFlux,
                                                                      firstDerivativePhysicalFlux,
                                                                      fESpace,
                                                                      data,
                                                                      section ),
    M_brentToll                                   ( data( ( section + "godunov/brent_toll" ).data(), 1e-4 ) ),
    M_brentMaxIter                                ( data( ( section + "godunov/brent_maxIter" ).data(), 20) )
{

} // Constructor

// Destructor of the class
template < typename Mesh, typename SolverType >
GodunovNumericalFlux<Mesh, SolverType>::
~GodunovNumericalFlux ()
{

} // Destructor

// ===================================================
// Operators
// ===================================================

template < typename Mesh, typename SolverType >
Real
GodunovNumericalFlux<Mesh, SolverType>::
operator() ( const Real& leftState, const Real& rightState, const normal_Type& normal,
             const UInt& iElem, const Real& t, const Real& x, const Real& y, const Real& z ) const
{

    // It will store the value of the flux
    Real fluxValue( static_cast<Real>(0) );

    // It will store the argmin or argmax of flux dot normal
    Real minMax( static_cast<Real>(0) );

    // The normal flux function
    scalarFunction_Type normalFlux;

    // Select where it is the upwind direction
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

} // Namespace LifeV

#endif //_NUMERICALFLUXES_H_
