/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): A. Fumagalli  <alessio.fumagalli@mail.polimi.it>
            M. Kern       <michel.kern@inria.fr>
      Date: 2010-09

 Copyright (C) 2001-2006 EPFL, Politecnico di Milano, INRIA
 Copyright (C) 2006-2010 EPFL, Politecnico di Milano

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
  @file hyperbolicSolver.hpp
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author M. Kern <michel.kern@inria.fr>
  @date 09/2010

  @brief This file contains a solver class for hyperbolic equations
*/
#ifndef _hyperbolicSolver_H_
#define _hyperbolicSolver_H_

#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>

#include <life/lifealg/clapack.hpp>
#include <life/lifesolver/dataHyperbolic.hpp>
#include <life/lifealg/SolverTrilinos.hpp>
//
#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/NumericalFluxes.hpp>
//
#include <life/lifecore/displayer.hpp>

namespace
{

struct HyperbolicDefaultSource
{
    const LifeV::Real operator()( const LifeV::Real&, const LifeV::Real&,
                                  const LifeV::Real&, const LifeV::Real&,
                                  const LifeV::UInt&) const
    {
    	return static_cast<LifeV::Real>( 0 );
    }
};

struct HyperbolicDefaultInitialSolution
{
     LifeV::Real operator()( const LifeV::Real&, const LifeV::Real&,
                             const LifeV::Real&, const LifeV::Real&,
                             const LifeV::UInt& )
    {
    	return static_cast<LifeV::Real>( 0 );
    }
};

}

// LifeV namespace.
namespace LifeV
{
/*!
*/
template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class HyperbolicSolver
{

public:

    // Policies.
    //! @name Policies
    //@{

    typedef boost::function<Real ( const Real&, const Real&, const Real&,
                                   const Real&, const UInt& )>
                                                     Function;

    typedef DataHyperbolic<Mesh>                     data_type;

    typedef BCHandler                                bchandler_raw_type;
    typedef boost::shared_ptr<bchandler_raw_type>    bchandler_type;

    typedef typename SolverType::vector_type         vector_type;
    typedef boost::shared_ptr<vector_type>           vector_ptrtype;

    typedef Epetra_Comm                              comm_type;
    typedef boost::shared_ptr< comm_type >           comm_ptrtype;

    typedef AbstractNumericalFlux<Mesh, SolverType>  flux_type;
    typedef boost::shared_ptr< flux_type >           flux_ptrtype;

    //@}

    // Constructors and destructor.
    //! @name Constructors and destructor
    //@{

    /*!
      Full constructor for the class.
      @param dataFile Data for the problem.
      @param fESpace Discontinuous finite element space.
      @param bcHandler Boundary conditions for the problem.
      @param comm Shared pointer of the Epetra communicator.
    */
    HyperbolicSolver ( const data_type&          dataFile,
                       FESpace<Mesh, EpetraMap>& fESpace,
                       bchandler_raw_type&       bcHandler,
                       comm_ptrtype&             comm );

    /*!
      Constructor for the class without the definition of the boundary handler.
      @param dataFile Data for the problem.
      @param fESpace Discontinuous finite element space.
      @param comm Shared pointer of the Epetra communicator.
    */
    HyperbolicSolver ( const data_type&          dataFile,
                       FESpace<Mesh, EpetraMap>& fESpace,
                       comm_ptrtype&             comm );

    //! Virtual destructor.
    virtual ~HyperbolicSolver ();

    //@}

    // Update & set methos.
    //! @name Update & set methos
    //@{

    /*!
      Set up the linear solver, the preconditioner for the linear system
      and the exporter to save the solution.
    */
    void setup ();

    /*!
      Compute the initial solution as an interpolation on the analytical
      initial solution.
      @param initialSolution The initial solution function.
    */
    inline void setInitialSolution ( const Function& initialSolution )
    {

        // Interpolate the initial solution.
        M_FESpace.interpolate( initialSolution,
                               *M_uOld,
                               M_data.dataTime()->getInitialTime() );

        *M_u = *M_uOld;
    }

    /*!
      Set the boundary conditions.
      @param bcHandler Boundary condition handler for the problem.
    */
    inline void setBC ( bchandler_raw_type& bcHandler )
    {
        M_BCh   = &bcHandler;
        M_setBC = true;
    }

    /*!
      Set the source term, the default setted source term is the zero function.
      @param source Source term for the problem.
    */
    inline void setSourceTerm ( const Function& source )
    {
        M_source = source;
    }


    inline void setNumericalFlux ( const GodunovNumericalFlux<Mesh, SolverType>& flux )
    {

        M_numericalFlux.reset ( new GodunovNumericalFlux<Mesh, SolverType>( flux ) );
    }

    /*!
      Set the solution vector.
      @param solution Constant vector_type reference of the solution.
    */
    inline void setSolution ( const vector_ptrtype& solution )
    {
        // Set both the final step solution and beginning step solution.
        M_u    = solution;
        M_uOld = solution;
    }


    //@}

    // Get methods.
    //! @name Get methods
    //@{

    /*!
      Returns the local hybrid solution vector.
      @return Constant vector_type reference of the hybrid vector.
     */
    inline const vector_ptrtype& solution () const
    {
        return M_u;
    }

    /*!
      Return if the bounday conditions is setted or not.
      @return Constant boolean with value true if the boundary condition is setted,
      false otherwise
    */
    inline const bool BCset () const
    {
        return M_setBC;
    }

    /*!
      Return the boundary conditions handler.
      @return Reference of boundary conditions handler.
    */
    inline bchandler_type& bcHandler ()
    {
        return M_BCh;
    }

    /*!
      Return the Epetra local map.
      @return Constant EpetraMap reference of the problem.
    */
    inline EpetraMap const& getMap () const
    {
        return M_localMap;
    }

    /*!
      Return the Epetra communicator.
      @return Constant of the shared pointer of the communicator of the problem.
    */
    inline const comm_ptrtype& comm () const
    {
        return M_displayer.comm();
    }

    /*!
      Return the displayer.
      @return Constant reference of the displayer of the problem.
    */
    inline Displayer const & getDisplayer() const
    {
        return M_displayer;
    }

    //@}

    // Solve functions.
    //! @name Solve functions
    //@{

    /*!
      Solve one time step of the hyperbolic problem.
     */
    void solveOneStep();

    /*!
      Compute the global CFL condition.
    */
    Real CFL() const;

    //@}

protected:

    /*!
      Return the number of total degrees of freem of the problem.
      @return The number of total degrees of freedom as a constant UInt.
    */
    inline const UInt dim () const
    {
        return M_FESpace.dim();
    }

    // Local solve functions.
    //! @name Local solve functions
    //@{

    void localReconstruct ( const UInt& Elem );

    void localEvolve      ( const UInt& iElem );

    void localAverage     ( const UInt& iElem );

    //@}

    // Parallel stuff.
    //! @name Parallel stuff
    //@{

    //! MPI process identifier.
    UInt                M_me;

    //! Local map.
    EpetraMap           M_localMap;

    //@}

    // Data of the problem.
    //! @name Data of the problem
    //@{

    //! Data for Darcy solvers.
    const data_type&    M_data;

    //! Source function.
    Function            M_source;

    //! Initial solution.
    Function            M_initialSolution;

    boost::shared_ptr< GodunovNumericalFlux<Mesh, SolverType> > M_numericalFlux;

    //! Bondary conditions handler.
    bchandler_raw_type* M_BCh;

    //! Flag if the boundary conditions are setted or not.
    bool                M_setBC;

    //@}

    // Finite element spaces.
    //! @name Finite element spaces
    //@{

    //! Primal finite element space.
    FESpace<Mesh, EpetraMap>& M_FESpace;

    //@}

    // Algebraic stuff.
    //! @name Algebraic stuff
    //@{

    //! Right hand side.
    vector_ptrtype M_rhs;

    //! Primal solution.
    vector_ptrtype M_u;

    vector_ptrtype M_uOld;

    vector_ptrtype M_globalFlux;

    Real           M_CFLrelax;

    //@}

    // Elementary matrices and vectors used for the static condensation.
    //! @name Elementary matrices and vectors used for the static condensation.
    //@{


    std::vector<ElemMat> M_elmatMass;

    ElemVec              M_localFlux;

    //@}

    Displayer M_displayer;

}; // class HyperbolicSolver

//
// IMPLEMENTATION
//

// Complete constructor.
template<typename Mesh, typename SolverType>
HyperbolicSolver<Mesh, SolverType>::
HyperbolicSolver ( const data_type&          dataFile,
                   FESpace<Mesh, EpetraMap>& fESpace,
                   bchandler_raw_type&       bcHandler,
                   comm_ptrtype&             comm ):
    // Parallel stuff.
    M_me              ( comm->MyPID() ),
    M_localMap        ( fESpace.map() ),
    M_displayer       ( comm ),
    // Data of the problem.
    M_data            ( dataFile ),
	M_source          ( HyperbolicDefaultSource() ),
    M_BCh             ( &bcHandler ),
    M_setBC           ( true ),
    // Finite element spaces.
    M_FESpace         ( fESpace ),
    // Algebraic stuff.
    M_rhs             ( new vector_type ( M_localMap ) ),
    M_u    			  ( new vector_type ( M_FESpace.map(), Repeated ) ),
    M_uOld			  ( new vector_type ( M_FESpace.map(), Repeated ) ),
    M_globalFlux      ( new vector_type ( M_FESpace.map(), Repeated ) ),
    // Local matrices and vectors.
    M_localFlux       ( M_FESpace.refFE().nbDof(), 1 ),
    M_elmatMass       ( ),
    M_initialSolution ( HyperbolicDefaultInitialSolution() )
{

    // Interpolate the primal initial value on the dafault initial value function.
    M_FESpace.interpolate( M_initialSolution,
                           *M_u,
                           M_data.dataTime()->getInitialTime() );

    M_elmatMass.reserve( M_FESpace.mesh()->numElements() );

    CONSTRUCTOR( "HyperbolicSolver" );

} // Constructor


// Constructor without boundary condition handler.
template<typename Mesh, typename SolverType>
HyperbolicSolver<Mesh, SolverType>::
HyperbolicSolver ( const data_type&          dataFile,
                   FESpace<Mesh, EpetraMap>& fESpace,
                   comm_ptrtype&             comm ):
    // Parallel stuff.
    M_me              ( comm->MyPID() ),
    M_localMap        ( fESpace.map() ),
    M_displayer       ( comm ),
    // Data of the problem.
    M_data            ( dataFile ),
	M_source          ( HyperbolicDefaultSource() ),
    M_setBC           ( false ),
    // Finite element spaces.
    M_FESpace         ( fESpace ),
    // Algebraic stuff.
    M_rhs             ( new vector_type ( M_localMap ) ),
    M_u    			  ( new vector_type ( M_FESpace.map(), Repeated ) ),
    M_uOld			  ( new vector_type ( M_FESpace.map(), Repeated ) ),
    M_globalFlux      ( new vector_type ( M_FESpace.map(), Repeated ) ),
    // Local matrices and vectors.
    M_localFlux       ( M_FESpace.refFE().nbDof(), 1 ),
    M_elmatMass       ( ),
    M_initialSolution ( HyperbolicDefaultInitialSolution() )
{

    // Interpolate the primal initial value on the dafault initial value function.
    M_FESpace.interpolate( M_initialSolution,
                           *M_u,
                           M_data.dataTime()->getInitialTime() );

    M_elmatMass.reserve( M_FESpace.mesh()->numElements() );

    CONSTRUCTOR( "HyperbolicSolver" );

} // Constructor

// Virtual destructor.
template <typename Mesh, typename SolverType>
HyperbolicSolver<Mesh, SolverType>::
~HyperbolicSolver ()
{

	DESTRUCTOR( "HyperbolicSolver" );

} // Destructor

// Set up the linear solver, the preconditioner and the exporter.
template<typename Mesh, typename SolverType>
void
HyperbolicSolver<Mesh, SolverType>::
setup ()
{

    // Flags for the BLAS and LAPACK routine.
    int INFO[1]      = {0};
    int NB[1]        = { M_FESpace.refFE().nbDof() };

    // Parameter that indicate the Lower storage of matrices.
    char _param_L[1] = {'L'};

    // Total number of elements.
    UInt meshNumberOfElements = M_FESpace.mesh()->numElements();

    // For each element it creates the mass matrix and factorize it using Cholesky.
    for ( UInt iElem(1); iElem <= meshNumberOfElements; ++iElem )
    {

        // Update the element.
        M_FESpace.fe().update( M_FESpace.mesh()->element( iElem ),
                               UPDATE_QUAD_NODES | UPDATE_WDET );

        // Local mass matrix
        ElemMat matElem(M_FESpace.refFE().nbDof(), 1, 1);
        matElem.zero();

        // Compute the mass matrix for the current element
        mass( static_cast<Real>(1.), matElem, M_FESpace.fe(), 0, 0);

        /* Put in M the matrix L and L^T, where L and L^T is the Cholesky factorization of M.
           For more details see http://www.netlib.org/lapack/double/dpotrf.f */
        dpotrf_( _param_L, NB, matElem.mat(), NB, INFO );
        ASSERT_PRE( !INFO[0], "Lapack factorization of M is not achieved." );

        // Save the local mass matrix in the global vector of mass matrices
        M_elmatMass.push_back( matElem );

    }

} // setup


template<typename Mesh,typename SolverType>
void
HyperbolicSolver<Mesh, SolverType>::
localReconstruct ( const UInt& iElem )
{

    ;

} // localReconstruct

template<typename Mesh, typename SolverType>
void
HyperbolicSolver<Mesh, SolverType>::
localEvolve ( const UInt& iElem )
{

    // Flags for the BLAS and LAPACK routine.
    int INFO[1]      = {0};
    int NB[1]        = { M_FESpace.refFE().nbDof() };

    // Parameter that indicate the Lower storage of matrices.
    char _param_L[1] = {'L'};
    char _param_N[1] = {'N'};

    // Paramater that indicate the Transpose of matrices.
    char _param_T[1] = {'T'};

    // Numbers of columns of the right hand side := 1.
    int NBRHS[1]     = {1};

    // Clean the local flux
    M_localFlux.zero();

    // Loop on the faces of the element iElem and compute the local contribution
    for ( UInt iFace(1); iFace <= M_FESpace.mesh()->numLocalFaces(); ++iFace )
    {

        const UInt iGlobalFace( M_FESpace.mesh()->localFaceId( iElem, iFace ) );

        // Take the left element to the face, see regionMesh for the meaning of left element
        UInt leftElement( M_FESpace.mesh()->faceElement( iGlobalFace, 1 ) );

        // Take the right element to the face, see regionMesh for the meaning of right element
        UInt rightElement( M_FESpace.mesh()->faceElement( iGlobalFace, 2 ) );

        // Update the normal vector of the current face in each quadrature point
        M_FESpace.feBd().updateMeasNormalQuadPt( M_FESpace.mesh()->bElement( iGlobalFace ) );

        // Local flux of a face times the integration weight
        ElemVec localFaceFluxWeight ( M_FESpace.refFE().nbDof(), 1 );

        // Solution in the left element
        ElemVec leftValue  ( M_FESpace.refFE().nbDof(), 1 );

        // Solution in the right element
        ElemVec rightValue ( M_FESpace.refFE().nbDof(), 1 );

        // Extract the solution in the current element, now is the leftElement
        extract_vec( *M_uOld,
                     leftValue,
                     M_FESpace.refFE(),
                     M_FESpace.dof(),
                     leftElement , 0 );

        // Check if the current face is a boundary face, that is rightElement == 0
        if ( !rightElement )
        {

            // Clean the value of the right element
            rightValue.zero();

            // Check if the boundary conditions were updated.
            if ( !M_BCh->bdUpdateDone() )
            {
                // Update the boundary conditions handler. We use the finite element of the boundary of the dual variable.
                M_BCh->bdUpdate( *M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );
            }

            // Take the boundary marker for the current boundary face
            const Index_t faceMarker ( M_FESpace.mesh()->bElement( iGlobalFace ).marker() );

            // Take the corrispective boundary function
            const BCBase& bcBase ( M_BCh->GetBCWithFlag( faceMarker ) );

            // Check if the bounday condition is of type Essential, useful for operator splitting strategies
            if ( bcBase.type() == Essential )
            {

                // Loop on all the quadrature points
                for ( UInt ig(0); ig < M_FESpace.feBd().nbQuadPt; ++ig)
                {

                    // Coordinates of the current quadrature point
                    Real x(0), y(0), z(0);
                    x = M_FESpace.feBd().quadPt( ig, static_cast<UInt>(0) );
                    y = M_FESpace.feBd().quadPt( ig, static_cast<UInt>(1) );
                    z = M_FESpace.feBd().quadPt( ig, static_cast<UInt>(2) );

                    // Compute the boundary contribution
                    rightValue[0] = bcBase( M_data.dataTime()->getTime(), x, y, z, 0 );

                    // Compute the outward unit normal of the boundary
                    const KN<Real> normal ( M_FESpace.feBd().normal('.', static_cast<Int>(ig) ) );

                    const  Real localFaceFlux = M_numericalFlux->getFirstDerivativePhysicalFluxDotNormal ( normal,
                                                                                                           iElem,
                                                                                                           M_data.dataTime()->getTime(),
                                                                                                           x, y, z, rightValue[ 0 ] );
                    // Update the local flux of the current face with the quadrature weight
                    localFaceFluxWeight[0] += localFaceFlux * M_FESpace.feBd().weightMeas( ig );
                }

            }
            else
            {
                /* If the boundary flag is not Essential then is automatically an outflow boundary.
                   We impose to localFaceFluxWeight a positive value. */
                localFaceFluxWeight[0] = 1.;
            }

            // It is an outflow face, we use a ghost cell
            if ( localFaceFluxWeight[0] > 1e-4 )
            {
                rightValue = leftValue;
            }

            // Clean the localFaceFluxWeight
            localFaceFluxWeight.zero();

        }
        else
        {
            // The current element is not a boundary element

            // Clean the rightValue
            rightValue.zero();

            // Extract the solution in the right element
            extract_vec( *M_uOld,
                         rightValue,
                         M_FESpace.refFE(),
                         M_FESpace.dof(),
                         rightElement, 0 );

        }

        // Clean the localFaceFluxWeight
        localFaceFluxWeight.zero();

        // Loop on all the quadrature points
        for ( UInt ig(0); ig < M_FESpace.feBd().nbQuadPt; ++ig )
        {

            // Coordinates of the current quadrature point
            Real x(0), y(0), z(0);

            // Set the coordinates of the current quatrature point
            x = M_FESpace.feBd().quadPt( ig, static_cast<UInt>(0) );
            y = M_FESpace.feBd().quadPt( ig, static_cast<UInt>(1) );
            z = M_FESpace.feBd().quadPt( ig, static_cast<UInt>(2) ) ;

            KN<Real> normal ( M_FESpace.feBd().normal( '.', static_cast<Int>(ig) ) );


            // If the normal is orientated inward, we change its sign and swap the left and value of the solution
            if ( iElem == rightElement )
            {
                normal *= -1.;
                std::swap( leftValue, rightValue );
            }

            const Real localFaceFlux = (*M_numericalFlux)( leftValue[ 0 ],
                                                           rightValue[ 0 ],
                                                           normal,
                                                           iElem,
                                                           M_data.dataTime()->getTime(), x, y, z );

            // Update the local flux of the current face with the quadrature weight
            localFaceFluxWeight[0] += localFaceFlux * M_FESpace.feBd().weightMeas( ig );

        }

        /* Put in localFlux the vector L^{-1} * localFlux
           For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
        dtrtrs_( _param_L, _param_N, _param_N, NB, NBRHS, M_elmatMass[ iElem - 1 ].mat(), NB, localFaceFluxWeight, NB, INFO);
        ASSERT_PRE( !INFO[0], "Lapack Computation M_elvecSource = LB^{-1} rhs is not achieved." );

        /* Put in localFlux the vector L^{-T} * localFlux
           For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
        dtrtrs_( _param_L, _param_T, _param_N, NB, NBRHS, M_elmatMass[ iElem - 1 ].mat(), NB, localFaceFluxWeight, NB, INFO);
        ASSERT_PRE( !INFO[0], "Lapack Computation M_elvecSource = LB^{-1} rhs is not achieved." );

        // Add to the local flux the local flux of the current face
        M_localFlux += localFaceFluxWeight;

    }

} // localEvolve

template<typename Mesh, typename SolverType>
void
HyperbolicSolver<Mesh, SolverType>::
localAverage ( const UInt& iElem )
{

    ;

} // localAverage

template<typename Mesh, typename SolverType>
Real
HyperbolicSolver<Mesh, SolverType>::
CFL() const
{
//    return M_data.dataTime()->getTimeStep();

    // Total number of elements in the mesh
    const UInt meshNumberOfElements( M_FESpace.mesh()->numElements() );

    // The local value for the CFL condition, without the time step
    Real localCFL(0.), localCFLOld( - 1. );

    // Loop on all the elements to perform the fluxes
    for ( UInt iElem(1); iElem <= meshNumberOfElements; ++iElem )
    {
        // Update the property of the current element
        M_FESpace.fe().update( M_FESpace.mesh()->element( iElem ),
                               UPDATE_QUAD_NODES | UPDATE_WDET | UPDATE_PHI );

        // Volumetric measure of the current element
        const Real K( M_FESpace.fe().measure() );

        // Loop on the faces of the element iElem and compute the local contribution
        for ( UInt iFace(1); iFace <= M_FESpace.mesh()->numLocalFaces(); ++iFace )
        {

            const UInt iGlobalFace( M_FESpace.mesh()->localFaceId( iElem, iFace ) );

            // Update the normal vector of the current face in each quadrature point
            M_FESpace.feBd().updateMeasNormalQuadPt( M_FESpace.mesh()->bElement( iGlobalFace ) );

            // Take the left element to the face, see regionMesh for the meaning of left element
            UInt leftElement( M_FESpace.mesh()->faceElement( iGlobalFace, 1 ) );

            // Take the right element to the face, see regionMesh for the meaning of right element
            UInt rightElement( M_FESpace.mesh()->faceElement( iGlobalFace, 2 ) );

            // Solution in the left element
            ElemVec leftValue  ( M_FESpace.refFE().nbDof(), 1 );

            // Solution in the right element
            ElemVec rightValue ( M_FESpace.refFE().nbDof(), 1 );

            // Extract the solution in the current element, now is the leftElement
            extract_vec( *M_uOld,
                         leftValue,
                         M_FESpace.refFE(),
                         M_FESpace.dof(),
                         leftElement , 0 );

            if ( rightElement )
            {

                // Extract the solution in the current element, now is the leftElement
                extract_vec( *M_uOld,
                             rightValue,
                             M_FESpace.refFE(),
                             M_FESpace.dof(),
                             rightElement , 0 );
            }
            else
            {
                rightValue = leftValue;
            }

            // Area of the current face
            const Real e( M_FESpace.feBd().measure() );

            // Loop on all the quadrature points
            for ( UInt ig(0); ig < M_FESpace.feBd().nbQuadPt; ++ig )
            {

                // Coordinates of the current quadrature point
                Real x(0), y(0), z(0);

                // Set the coordinates of the current quatrature point
                x = M_FESpace.feBd().quadPt( ig, static_cast<UInt>(0) );
                y = M_FESpace.feBd().quadPt( ig, static_cast<UInt>(1) );
                z = M_FESpace.feBd().quadPt( ig, static_cast<UInt>(2) ) ;

                KN<Real> normal ( M_FESpace.feBd().normal( '.', static_cast<Int>(ig) ) );

                // Compute the local CFL without the time step
                localCFL = e / K * M_numericalFlux->getNormInfty ( leftValue[0], rightValue[0], normal, iElem,
                                                                   M_data.dataTime()->getTime(), x, y, z );

                // Select the maximum between the old CFL condition and the new CFL condition
                if ( localCFL > localCFLOld  )
                {
                    localCFLOld = localCFL;
                }
                else
                {
                    localCFL = localCFLOld;
                }

            }

        }

    }

    // Compute the timeStep according to CLF
    const Real timeStep( M_data.getCFLrelax() / localCFL );

    // Return the correct value of the time step
    return timeStep;

} //CFL

// Solve one time step of the hyperbolic problem.
template<typename Mesh, typename SolverType>
void
HyperbolicSolver<Mesh, SolverType>::
solveOneStep ()
{

    // Total number of elements in the mesh
    const UInt meshNumberOfElements( M_FESpace.mesh()->numElements() );

    // Loop on all the elements to perform the fluxes
    for ( UInt iElem(1); iElem <= meshNumberOfElements; ++iElem )
    {

        // Update the property of the current element
        M_FESpace.fe().update( M_FESpace.mesh()->element( iElem ),
                               UPDATE_QUAD_NODES | UPDATE_WDET | UPDATE_PHI );

        // Reconstruct step of the current element
        localReconstruct( iElem );

        // Evolve step of the current element
        localEvolve( iElem  );

        // Put the total flux of the current element in the global vector of fluxes
        assembleVector( *M_globalFlux,
    	    			M_FESpace.fe().currentLocalId(),
       				    M_localFlux,
       				    M_FESpace.refFE().nbDof(),
       				    M_FESpace.dof(), 0 );

        // Average step of the current element
        localAverage( iElem );

    }

    // Assemble the global hybrid vector.
    M_globalFlux->GlobalAssemble();

    // Update the value of the solution
    (*M_u) = (*M_uOld) - M_data.dataTime()->getTimeStep() * (*M_globalFlux);

    // Clean the vector of fluxes
    M_globalFlux.reset( new vector_type( M_FESpace.map(), Repeated ) );

    // Update the solution at previous time step
    *M_uOld = *M_u;

} // run

} // namespace LifeV


#endif //_hyperbolicSolver_H_


/*

*/
