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
 * @brief Solver class for hyperbolic scalar equations.
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

#ifndef _HYPERBOLICSOLVER_H_
#define _HYPERBOLICSOLVER_H_ 1

#include <Epetra_LAPACK.h>

#include <lifev/core/algorithm/SolverAztecOO.hpp>

#include <lifev/core/fem/AssemblyElemental.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/fem/HyperbolicFluxNumerical.hpp>

#include <lifev/core/solver/HyperbolicData.hpp>

namespace
{

LifeV::Real _ZeroFun ( const LifeV::Real&, const LifeV::Real&,
                       const LifeV::Real&, const LifeV::Real&,
                       const LifeV::UInt& )
{
    return static_cast<LifeV::Real> ( 0. );
}

}

// LifeV namespace.
namespace LifeV
{
//! HyperbolicSolver Implements an hyperbolic solver.
/*!

  @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author Michel Kern       <michel.kern@inria.fr>
  @see For applications related to two-phase flow see \cite Fumagalli2011a

  This class implements an hyperbolic solver.
  <br>
  <br>
  This class solves a general hyperbolic scalar equation in the conservative form: find \f$ u \in V \f$ such that
  \f[
  \left\{
  \begin{array}{l l}
  \displaystyle \frac{\partial u}{\partial t} + \nabla \cdot \mathbf{F} ( u ) = f & \mathrm{ in } \quad \Omega \times [0,T]\,, \\
  u(\mathbf{x}, 0) = u_0( \mathbf{x} ) & \mathrm{ in } \quad \Omega\,,
  \end{array}
  \right.
  \f]
  with prescribed inflow boundary conditions on \f$ \partial \Omega_{in} \f$
  \f[
  u = \bar{u} \quad \mathrm{ on } \quad \partial \Omega_{in}\,.
  \f]
  We have \f$ \partial \Omega = \partial \Omega_{in} \cap \partial \Omega_{out}\f$.
  Since this solver can be used coupled with an elliptic solver, the Neumann and Robin boundary
  are trated as outflow boundary, while the Dirichlet boundary is automatically dedived into
  outflow boundary and inflow boundary. This choice is due to the general form of the flux function \f$ \mathbf{F} \f$.
  <br>
  Introducing a conforming triangolation \f$ \mathcal{T}_h = \cup K \f$ of the domain \f$ \Omega \f$, we write the strong formulation
  in a weak form in each element \f$ K \f$: give a test function \f$ v \f$
  \f[
  \displaystyle \int_K \frac{\partial u}{\partial t} v + \int_K \nabla \cdot \mathbf{F} ( u ) v = \int_K f v \quad \forall v \in V\,,
  \f]
  integrating by part the integral with the divergence we find
  \f[
  \displaystyle \int_K \frac{\partial u}{\partial t} v - \int_K \mathbf{F} ( u ) \cdot \nabla v +
  \int_{\partial K} \mathbf{F}( u ) \cdot \mathbf{n} v = \int_K f v \quad \forall v \in V\,.
  \f]
  The Discontinuous Galerkin formulation of the problem is given approximating \f$ V \f$ with a discontinuous finite element
  space \f$ V_h \f$, while the numerical flux at the boundaries is approximated using a suitable numerical scheme obtaining
  \f$ \hat{\mathbf{F}} \f$.
  <br>
  Summing on all the elements \f$ K \in \mathcal{T}_h \f$, the finite element problem reads: find \f$ u_h \in V_h \f$ such that
  \f[
  \displaystyle \sum_{K \in \mathcal{T}_h} \int_K \frac{\partial u_h}{\partial t} v_h - \int_K \mathbf{F} ( u_h ) \cdot \nabla v_h +
  \int_{\partial K} \hat{\mathbf{F}}( u_h ) \cdot \mathbf{n} v_h =  \sum_{K \in \mathcal{T}_h} \int_K f v_h \quad \forall v_h \in V_h\,.
  \f]
  Each time step \f$ \Delta t^n \f$ and mesh size \f$ h \f$ are coupled via \f$ CFL \f$ condition at step \f$ n \f$, using esplicit Euler
  scheme we find the relation between \f$ \Delta t^n \f$ and \f$ h \f$ satysfying \f$ CFL^n < 1 \f$.
  <br>
  A general form for computing the \f$ CFL \f$ condition is
   \f[
  CFL^n = \displaystyle \sup_{e \in \partial K, K \in \mathcal{T}_h } \Delta t^n \frac{\vert e \vert}{ \vert K \vert }
          \Vert \mathbf{F}^\prime \cdot \mathbf{n}_{e, K} \Vert_{L^\infty(a_0, b_0) }\,,
  \f]
  where
  \f[
  \begin{array}{l}
  \displaystyle a_0 = \inf_{x \in \Omega, t \in (t^{n-1}, t^n), y \in \partial \Omega_{in} } \left\{ u_0, \bar{u}(y, t) \right\}\,, \\
  \displaystyle b_0 = \sup_{x \in \Omega, t \in (t^{n-1}, t^n), y \in \partial \Omega_{in} } \left\{ u_0, \bar{u}(y, t) \right\}\,.
  \end{array}
  \f]
  Local approximation for the flux function \f$ \hat{\mathbf{F}} \cdot \mathbf{n} \f$ and the computation of
  \f$ \mathbf{F}^\prime \cdot \mathbf{n}_{e, K} \f$ are in NumericalFlux.hpp file.
  @note The implementation is given just for lowest order discontinuous finite elements.
  @todo Implement the forcing term \f$ f \f$ and implement high order finite elements.
  @todo When we will pass to Trilinos >= 10.6 use Epetra wrapper for LAPACK functions.
*/
template < typename Mesh,
         typename SolverType = LifeV::SolverAztecOO >
class HyperbolicSolver
{

public:

    //! @name Public Types
    //@{

    typedef boost::function < Real ( const Real&, const Real&, const Real&,
                                     const Real&, const UInt& ) >
    Function_Type;

    typedef HyperbolicData< Mesh >                   data_Type;

    typedef BCHandler                                bchandler_Type;
    typedef boost::shared_ptr< bchandler_Type >      bchandlerPtr_Type;

    typedef typename SolverType::vector_type         vector_Type;
    typedef boost::shared_ptr< vector_Type >         vectorPtr_Type;

    typedef Epetra_Comm                              comm_Type;
    typedef boost::shared_ptr< comm_Type >           commPtr_Type;

    typedef AbstractNumericalFlux<Mesh, SolverType>  flux_Type;
    typedef boost::shared_ptr< flux_Type >           fluxPtr_Type;

    typedef Real                                     ghostData_Type;
    typedef std::vector< ghostData_Type >            ghostDataContainer_Type;
    typedef std::map< UInt, ghostDataContainer_Type > buffer_Type;
    typedef std::map< ID, ghostData_Type >           ghostDataMap_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Full constructor for the class.
    /*!
      @param dataFile Data for the problem.
      @param fESpace Discontinuous finite element space.
      @param bcHandler Boundary conditions for the problem.
      @param comm Shared pointer of the Epetra communicator.
    */
    HyperbolicSolver ( const data_Type&          dataFile,
                       FESpace<Mesh, MapEpetra>& fESpace,
                       MapEpetra&                ghostMap,
                       bchandler_Type&           bcHandler,
                       commPtr_Type&             comm );

    //! Constructor for the class without the definition of the boundary handler.
    /*!
      @param dataFile Data for the problem.
      @param fESpace Discontinuous finite element space.
      @param comm Shared pointer of the Epetra communicator.
    */
    HyperbolicSolver ( const data_Type&          dataFile,
                       FESpace<Mesh, MapEpetra>& fESpace,
                       MapEpetra&                ghostMap,
                       commPtr_Type&             comm );

    //! Virtual destructor.
    virtual ~HyperbolicSolver ();

    //@}

    //! @name Methods
    //@{

    //! Setup the local mass matrices.
    void setup ();

    //! Solve one time step of the hyperbolic problem.
    void solveOneTimeStep();

    //! Compute the global CFL condition.
    Real CFL();

    //@}

    //! @name Set Methos
    //@{

    //! Set the initial solution for the computation.
    /*!
      Compute the initial solution as an interpolation on the analytical
      initial solution.
      @param initialSolution The initial solution function.
    */
    void setInitialSolution ( const Function_Type& initialSolution );

    //! Set the boundary conditions.
    /*!
      @param bcHandler Boundary condition handler for the problem.
    */
    inline void setBoundaryCondition ( bchandler_Type& bcHandler )
    {
        M_BCh.reset ( new bchandler_Type ( bcHandler ) );
        M_setBC = true;
    }

    //! Set the source term.
    /*
      The default setted source term is \f$ f( \mathbf{x}, t ) \equiv 0 \f$.
      @param source Source term for the problem.
    */
    inline void setSourceTerm ( const Function_Type& source )
    {
        M_source = source;
    }

    //! Set the mass function term.
    /*!
      It does not depend on time. The default setted source term is \f$ \phi( \mathbf{x} ) \equiv 1 \f$.
      @param mass Mass term for the problem.
    */
    inline void setMassTerm ( const Function_Type& mass )
    {
        M_mass = mass;
    }

    //! Set the numerical flux.
    /*!
      @param flux The numerical flux class
    */
    inline void setNumericalFlux ( const flux_Type& flux )
    {
        typedef GodunovNumericalFlux<Mesh, SolverType> godunov_Type;
        M_numericalFlux.reset ( new godunov_Type ( dynamic_cast<const godunov_Type&> ( flux ) ) );
    }

    //! Set the solution vector.
    /*!
      @param solution Constant vector_type reference of the solution.
    */
    inline void setSolution ( const vectorPtr_Type& solution )
    {
        // Set both the final step solution and beginning step solution.
        *M_u    = *solution;
        *M_uOld = *solution;
    }

    //@}

    //! @name Get Methods
    //@{

    //! Returns the solution vector.
    /*!
      @return Constant vector_type reference of the solution vector.
    */
    inline const vectorPtr_Type& solution () const
    {
        return M_u;
    }

    //! Return if the bounday conditions is setted or not.
    /*!
      @return Constant boolean with value true if the boundary condition is setted,
      false otherwise
    */
    inline bool isBoundaryConditionSet() const
    {
        return M_setBC;
    }

    //! Return the boundary conditions handler.
    /*!
      @return Reference of boundary conditions handler.
    */
    inline bchandlerPtr_Type& boundaryConditionHandler ()
    {
        return M_BCh;
    }


    //! Return the Epetra local map.
    /*!
      @return Constant MapEpetra reference of the problem.
    */
    inline MapEpetra const& map () const
    {
        return M_localMap;
    }

    //! Return the displayer.
    /*!
      Useful for parallel print in programs.
      @return Constant reference of the displayer of the problem.
    */
    inline Displayer const& getDisplayer() const
    {
        return M_displayer;
    }

    //! Returns displayer.
    /*!
      @return Reference of the displayer of the problem.
    */
    inline Displayer& getDisplayer()
    {
        return M_displayer;
    }

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Reconstruct locally the solution.
    void localReconstruct ( const UInt& Elem );

    //! Compute the local contribute
    void localEvolve      ( const UInt& iElem );

    //! Apply the flux limiters locally.
    void localAverage     ( const UInt& iElem );

    //@}

    //! MPI process identifier.
    const UInt                M_me;

    //! Local map.
    MapEpetra                 M_localMap;

    //! Parallel displayer
    Displayer                 M_displayer;

    //! Data for Darcy solvers.
    const data_Type&          M_data;

    //! Source function.
    Function_Type             M_source;

    //! Mass function, it does not depend on time
    Function_Type             M_mass;

    //! Bondary conditions handler.
    bchandlerPtr_Type         M_BCh;

    //! Flag if the boundary conditions are setted or not.
    bool                      M_setBC;

    //! Function of initial solution.
    Function_Type             M_initialSolution;

    //! Class of type AbstractNumericalFlux for local flux computations.
    fluxPtr_Type              M_numericalFlux;

    //! Finite element space.
    FESpace<Mesh, MapEpetra>& M_FESpace;

    //! Right hand side.
    vectorPtr_Type            M_rhs;

    //! Solution at current time step.
    vectorPtr_Type            M_u;

    //! Solution at previous time step.
    vectorPtr_Type            M_uOld;

    //! Computed numerical flux.
    vectorPtr_Type            M_globalFlux;

    //! Auxiliary vector for local fluxes.
    VectorElemental           M_localFlux;

    //! Vector of all local mass matrices, possibly with mass function.
    std::vector<MatrixElemental>  M_elmatMass;

private:

    //! @name Private Constructors
    //@{

    //! Inhibited copy constructor.
    HyperbolicSolver ( const HyperbolicSolver<Mesh, SolverType>& );

    //@}

    //! @name Private Operators
    //@{

    //! Inhibited assign operator.
    HyperbolicSolver& operator= ( const HyperbolicSolver<Mesh, SolverType>& );

    //@}

}; // class HyperbolicSolver

// ===================================================
// Constructors & Destructor
// ===================================================

// Complete constructor.
template< typename Mesh, typename SolverType >
HyperbolicSolver< Mesh, SolverType >::
HyperbolicSolver ( const data_Type&          dataFile,
                   FESpace<Mesh, MapEpetra>& fESpace,
                   MapEpetra&                ghostMap,
                   bchandler_Type&           bcHandler,
                   commPtr_Type&             comm ) :
    // Parallel stuff.
    M_me              ( comm->MyPID() ),
    M_localMap        ( fESpace.map() ),
    M_displayer       ( comm ),
    // Data of the problem.
    M_data            ( dataFile ),
    M_source          ( NULL ),
    M_mass            ( NULL ),
    M_BCh             ( new bchandler_Type ( bcHandler ) ),
    M_setBC           ( true ),
    M_initialSolution ( NULL ),
    M_numericalFlux   ( ),
    // Finite element spaces.
    M_FESpace         ( fESpace ),
    // Algebraic stuff.
    M_rhs             ( new vector_Type ( M_localMap ) ),
    M_u               ( new vector_Type ( M_FESpace.map(), Repeated ) ),
    M_uOld            ( new vector_Type ( ghostMap, Repeated ) ),
    M_globalFlux      ( new vector_Type ( ghostMap, Repeated ) ),
    // Local matrices and vectors.
    M_localFlux       ( M_FESpace.refFE().nbDof(), 1 ),
    M_elmatMass       ( )
{

    M_elmatMass.reserve ( M_FESpace.mesh()->numElements() );

} // Constructor


// Constructor without boundary condition handler.
template< typename Mesh, typename SolverType >
HyperbolicSolver< Mesh, SolverType >::
HyperbolicSolver ( const data_Type&          dataFile,
                   FESpace<Mesh, MapEpetra>& fESpace,
                   MapEpetra&                ghostMap,
                   commPtr_Type&             comm ) :
    // Parallel stuff.
    M_me              ( comm->MyPID() ),
    M_localMap        ( fESpace.map() ),
    M_displayer       ( comm ),
    // Data of the problem.
    M_data            ( dataFile ),
    M_source          ( NULL ),
    M_mass            ( NULL ),
    M_BCh             ( ),
    M_setBC           ( false ),
    M_initialSolution ( NULL ),
    M_numericalFlux   ( ),
    // Finite element spaces.
    M_FESpace         ( fESpace ),
    // Algebraic stuff.
    M_rhs             ( new vector_Type ( M_localMap ) ),
    M_u               ( new vector_Type ( M_FESpace.map(), Repeated ) ),
    M_uOld            ( new vector_Type ( ghostMap, Repeated ) ),
    M_globalFlux      ( new vector_Type ( ghostMap, Repeated ) ),
    // Local matrices and vectors.
    M_localFlux       ( M_FESpace.refFE().nbDof(), 1 ),
    M_elmatMass       ( )
{

    M_elmatMass.reserve ( M_FESpace.mesh()->numElements() );

} // Constructor

// Virtual destructor.
template< typename Mesh, typename SolverType >
HyperbolicSolver< Mesh, SolverType >::
~HyperbolicSolver ()
{

} // Destructor

// ===================================================
// Methods
// ===================================================

// Setup the local mass matrices.
template<typename Mesh, typename SolverType>
void
HyperbolicSolver<Mesh, SolverType>::
setup ()
{
    // LAPACK wrapper of Epetra
    Epetra_LAPACK lapack;

    // Flags for LAPACK routines.
    Int INFO[1]  = {0};
    Int NB = M_FESpace.refFE().nbDof();

    // Parameter that indicate the Lower storage of matrices.
    char param_L = 'L';

    // Total number of elements.
    UInt meshNumberOfElements = M_FESpace.mesh()->numElements();

    // Vector of interpolation of the mass function
    vector_Type vectorMass ( M_FESpace.map(), Repeated );

    // If the mass function is given take it, otherwise use one as a value.
    if ( M_mass != NULL )
    {
        // Interpolate the mass function on the finite element space.
        M_FESpace.interpolate ( M_mass, vectorMass );
    }
    else
    {
        vectorMass = 1.;
    }

    // For each element it creates the mass matrix and factorize it using Cholesky.
    for ( UInt iElem (0); iElem < meshNumberOfElements; ++iElem )
    {

        // Update the element.
        M_FESpace.fe().update ( M_FESpace.mesh()->element ( iElem),
                                UPDATE_QUAD_NODES | UPDATE_WDET );

        // Local mass matrix
        MatrixElemental matElem (M_FESpace.refFE().nbDof(), 1, 1);
        matElem.zero();

        // Compute the mass matrix for the current element
        VectorElemental massValue  ( M_FESpace.refFE().nbDof(), 1 );
        extract_vec ( vectorMass, massValue, M_FESpace.refFE(), M_FESpace.dof(), iElem, 0 );
        // TODO: this works only for P0
	AssemblyElemental::mass ( massValue[ 0 ], matElem, M_FESpace.fe(), 0, 0);

        /* Put in M the matrix L and L^T, where L and L^T is the Cholesky factorization of M.
           For more details see http://www.netlib.org/lapack/double/dpotrf.f */
        lapack.POTRF ( param_L, NB, matElem.mat(), NB, INFO );
        ASSERT_PRE ( !INFO[0], "Lapack factorization of M is not achieved." );

        // Save the local mass matrix in the global vector of mass matrices
        M_elmatMass.push_back ( matElem );

    }

    //make sure mesh facets are updated
    if (! M_FESpace.mesh()->hasLocalFacets() )
    {
        M_FESpace.mesh()->updateElementFacets();
    }

} // setup

// Solve one time step of the hyperbolic problem.
template< typename Mesh, typename SolverType >
void
HyperbolicSolver< Mesh, SolverType >::
solveOneTimeStep ()
{
    // Total number of elements in the mesh
    const UInt meshNumberOfElements ( M_FESpace.mesh()->numElements() );

    // Loop on all the elements to perform the fluxes
    for ( UInt iElem (0); iElem < meshNumberOfElements; ++iElem )
    {

        // Update the property of the current element
        M_FESpace.fe().update ( M_FESpace.mesh()->element(iElem), UPDATE_QUAD_NODES | UPDATE_WDET);

        // Reconstruct step of the current element
        localReconstruct ( iElem );

        // Evolve step of the current element
        localEvolve ( iElem  );

        // Put the total flux of the current element in the global vector of fluxes
        assembleVector ( *M_globalFlux,
                         M_FESpace.fe().currentLocalId(),
                         M_localFlux,
                         M_FESpace.refFE().nbDof(),
                         M_FESpace.dof(), 0 );

        // Average step of the current element
        localAverage ( iElem );

    }

    // Assemble the global hybrid vector.
    M_globalFlux->globalAssemble();

    // alternative: instead of modifying M_globalFlux.map, we can make a local copy with the correct map
    //    // this is needed since M_uOld.map != M_globalFlux.map
    //    vector_Type fluxCopy ( M_uOld->map() );
    //    fluxCopy = *M_globalFlux;

    // Update the value of the solution
    (*M_u) = (*M_uOld) - M_data.dataTime()->timeStep() * (*M_globalFlux);
    //    *M_u = *M_uOld - M_data.dataTime()->timeStep() * fluxCopy;

    // Clean the vector of fluxes
    M_globalFlux.reset ( new vector_Type ( M_uOld->map(), Repeated ) );

    // Update the solution at previous time step
    *M_uOld = *M_u;

} // solveOneStep

// Compute the global CFL condition.
template< typename Mesh, typename SolverType >
Real
HyperbolicSolver< Mesh, SolverType >::
CFL()
{
    // Total number of elements in the mesh
    const UInt meshNumberOfElements ( M_FESpace.mesh()->numElements() );

    // The local value for the CFL condition, without the time step
    Real localCFL (0.), localCFLOld ( - 1. );

    // Loop on all the elements to perform the fluxes
    for ( UInt iElem (0); iElem < meshNumberOfElements; ++iElem )
    {
        // Update the property of the current element
        M_FESpace.fe().update ( M_FESpace.mesh()->element (iElem), UPDATE_QUAD_NODES | UPDATE_WDET);

        // Volumetric measure of the current element
        const Real K ( M_FESpace.fe().measure() );

        // Loop on the faces of the element iElem and compute the local contribution
        for ( UInt iFace (0); iFace < M_FESpace.mesh()->numLocalFaces(); ++iFace )
        {

            const UInt iGlobalFace ( M_FESpace.mesh()->localFacetId ( iElem, iFace ) );

            // Update the normal vector of the current face in each quadrature point
            M_FESpace.feBd().update ( M_FESpace.mesh()->boundaryFacet ( iGlobalFace ), UPDATE_W_ROOT_DET_METRIC | UPDATE_NORMALS | UPDATE_QUAD_NODES );

            // Take the left element to the face, see regionMesh for the meaning of left element
            const UInt leftElement ( M_FESpace.mesh()->faceElement ( iGlobalFace, 0 ) );

            // Take the right element to the face, see regionMesh for the meaning of right element
            const UInt rightElement ( M_FESpace.mesh()->faceElement ( iGlobalFace, 1 ) );

            // Solution in the left element
            VectorElemental leftValue  ( M_FESpace.refFE().nbDof(), 1 );

            // Solution in the right element
            VectorElemental rightValue ( M_FESpace.refFE().nbDof(), 1 );

            // Extract the solution in the current element, now is the leftElement
            extract_vec ( *M_uOld,
                          leftValue,
                          M_FESpace.refFE(),
                          M_FESpace.dof(),
                          leftElement , 0 );

            if ( !Flag::testOneSet ( M_FESpace.mesh()->face ( iGlobalFace ).flag(),
                                     EntityFlags::PHYSICAL_BOUNDARY | EntityFlags::SUBDOMAIN_INTERFACE ) )
            {
                // Extract the solution in the current element, now is the leftElement
                extract_vec ( *M_uOld,
                              rightValue,
                              M_FESpace.refFE(),
                              M_FESpace.dof(),
                              rightElement , 0 );
            }
            else if ( Flag::testOneSet ( M_FESpace.mesh()->face ( iGlobalFace ).flag(), EntityFlags::SUBDOMAIN_INTERFACE ) )
            {
                // TODO: this works only for P0 elements
                // but extract_vec works only with lids while RightElement is a gid
                rightValue[ 0 ] = (*M_uOld) [ rightElement ];
            }
            else // Flag::testOneSet ( M_FESpace.mesh()->face ( iGlobalFace ).flag(), PHYSICAL_BOUNDARY )
            {
                rightValue = leftValue;
            }

            // Area of the current face
            const Real e ( M_FESpace.feBd().measure() );

            // Loop on all the quadrature points
            for ( UInt ig (0); ig < M_FESpace.feBd().nbQuadPt(); ++ig )
            {

                // Cuurent quadrature point
                KN<Real> quadPoint (3);
                // normal vector
                KN<Real> normal (3);

                for (UInt icoor (0); icoor < 3; ++icoor)
                {
                    quadPoint (icoor) = M_FESpace.feBd().quadPt ( ig, icoor );
                    normal (icoor)    = M_FESpace.feBd().normal ( icoor, ig ) ;
                }

                // Compute the local CFL without the time step
                localCFL = e / K * M_numericalFlux->normInfinity ( leftValue[0],
                                                                   rightValue[0],
                                                                   normal,
                                                                   iElem,
                                                                   M_data.dataTime()->time(),
                                                                   quadPoint (0),
                                                                   quadPoint (1),
                                                                   quadPoint (2) );

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


    // Compute the time step according to CLF for the current process
    Real timeStepLocal[]  = { M_data.getCFLRelaxParameter() / localCFLOld };
    Real timeStepGlobal[] = { 0. };

    // Compute the minimum of the computed time step for all the processes
    M_displayer.comm()->MinAll ( timeStepLocal, timeStepGlobal, 1 );

    // Return the computed value
    return *timeStepGlobal;

} //CFL

// ===================================================
// Set Methods
// ===================================================

// Set the initial solution for the computation.
template< typename Mesh, typename SolverType >
void
HyperbolicSolver< Mesh, SolverType >::
setInitialSolution ( const Function_Type& initialSolution )
{

    // interpolation must be done on a Unique map
    vector_Type uUnique ( M_u->map(), Unique );

    // Interpolate the initial solution.
    M_FESpace.interpolate ( initialSolution,
                            uUnique,
                            M_data.dataTime()->initialTime() );

    // Update the solutions
    *M_uOld = uUnique;
    *M_u    = uUnique;
} // setInitialSolution

// ===================================================
// Protected Methods
// ===================================================

// Reconstruct locally the solution.
template< typename Mesh, typename SolverType >
void
HyperbolicSolver< Mesh, SolverType >::
localReconstruct ( const UInt& /* iElem */ )
{

    ;

} // localReconstruct

// Compute the local contribute
template< typename Mesh, typename SolverType >
void
HyperbolicSolver< Mesh, SolverType >::
localEvolve ( const UInt& iElem )
{

    // LAPACK wrapper of Epetra
    Epetra_LAPACK lapack;

    // Flags for LAPACK routines.
    Int INFO[1]  = { 0 };
    Int NB = M_FESpace.refFE().nbDof();

    // Parameter that indicate the Lower storage of matrices.
    char param_L = 'L';
    char param_N = 'N';

    // Paramater that indicate the Transpose of matrices.
    char param_T = 'T';

    // Numbers of columns of the right hand side := 1.
    Int NBRHS = 1;

    // Clean the local flux
    M_localFlux.zero();

    // Check if the boundary conditions were updated.
    if ( !M_BCh->bcUpdateDone() )
    {
        // Update the boundary conditions handler. We use the finite element of the boundary of the dual variable.
        M_BCh->bcUpdate ( *M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );
    }

    // Loop on the faces of the element iElem and compute the local contribution
    for ( UInt iFace (0); iFace < M_FESpace.mesh()->numLocalFaces(); ++iFace )
    {
        // Id mapping
        const UInt iGlobalFace ( M_FESpace.mesh()->localFacetId ( iElem, iFace ) );

        // Take the left element to the face, see regionMesh for the meaning of left element
        const UInt leftElement ( M_FESpace.mesh()->faceElement ( iGlobalFace, 0 ) );

        // Take the right element to the face, see regionMesh for the meaning of right element
        const UInt rightElement ( M_FESpace.mesh()->faceElement ( iGlobalFace, 1 ) );

        // Update the normal vector of the current face in each quadrature point
        M_FESpace.feBd().update ( M_FESpace.mesh()->boundaryFacet ( iGlobalFace ), UPDATE_W_ROOT_DET_METRIC | UPDATE_NORMALS | UPDATE_QUAD_NODES );

        // Local flux of a face times the integration weight
        VectorElemental localFaceFluxWeight ( M_FESpace.refFE().nbDof(), 1 );

        // Solution in the left element
        VectorElemental leftValue  ( M_FESpace.refFE().nbDof(), 1 );

        // Solution in the right element
        VectorElemental rightValue ( M_FESpace.refFE().nbDof(), 1 );

        // Extract the solution in the current element, now is the leftElement
        extract_vec ( *M_uOld,
                      leftValue,
                      M_FESpace.refFE(),
                      M_FESpace.dof(),
                      leftElement , 0 );

        // Check if the current face is a boundary face, that is rightElement == NotAnId
        if ( !Flag::testOneSet ( M_FESpace.mesh()->face ( iGlobalFace ).flag(), EntityFlags::PHYSICAL_BOUNDARY | EntityFlags::SUBDOMAIN_INTERFACE ) )
        {
            // Extract the solution in the current element, now is the leftElement
            extract_vec ( *M_uOld,
                          rightValue,
                          M_FESpace.refFE(),
                          M_FESpace.dof(),
                          rightElement , 0 );
        }
        else if ( Flag::testOneSet ( M_FESpace.mesh()->face ( iGlobalFace ).flag(), EntityFlags::SUBDOMAIN_INTERFACE ) )
        {
            // TODO: this works only for P0 elements
            //            rightValue[ 0 ] = M_ghostDataMap[ iGlobalFace ];
            rightValue[ 0 ] = (*M_uOld) [ rightElement ];
        }
        else // Flag::testOneSet ( M_FESpace.mesh()->face ( iGlobalFace ).flag(), PHYSICAL_BOUNDARY )
        {

            // Clean the value of the right element
            rightValue.zero();

            // Take the boundary marker for the current boundary face
            const ID faceMarker ( M_FESpace.mesh()->boundaryFacet ( iGlobalFace ).markerID() );

            // Take the corrispective boundary function
            const BCBase& bcBase ( M_BCh->findBCWithFlag ( faceMarker ) );

            // Check if the bounday condition is of type Essential, useful for operator splitting strategies
            if ( bcBase.type() == Essential )
            {

                // Loop on all the quadrature points
                for ( UInt ig (0); ig < M_FESpace.feBd().nbQuadPt(); ++ig)
                {

                    // Current quadrature point
                    KN<Real> quadPoint (3);

                    // normal vector
                    KN<Real> normal (3);

                    for (UInt icoor (0); icoor < 3; ++icoor)
                    {
                        quadPoint (icoor) = M_FESpace.feBd().quadPt ( ig, icoor );
                        normal (icoor)    = M_FESpace.feBd().normal ( icoor, ig ) ;
                    }

                    // Compute the boundary contribution
                    rightValue[0] = bcBase ( M_data.dataTime()->time(), quadPoint (0), quadPoint (1), quadPoint (2), 0 );

                    const Real localFaceFlux = M_numericalFlux->firstDerivativePhysicalFluxDotNormal ( normal,
                                               iElem,
                                               M_data.dataTime()->time(),
                                               quadPoint (0),
                                               quadPoint (1),
                                               quadPoint (2),
                                               rightValue[ 0 ] );
                    // Update the local flux of the current face with the quadrature weight
                    localFaceFluxWeight[0] += localFaceFlux * M_FESpace.feBd().wRootDetMetric ( ig );
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

        // Clean the localFaceFluxWeight
        localFaceFluxWeight.zero();

        // Loop on all the quadrature points
        for ( UInt ig (0); ig < M_FESpace.feBd().nbQuadPt(); ++ig )
        {

            // Current quadrature point
            KN<Real> quadPoint (3);

            // normal vector
            KN<Real> normal (3);

            for (UInt icoor (0); icoor < 3; ++icoor)
            {
                quadPoint (icoor) = M_FESpace.feBd().quadPt ( ig, icoor );
                normal (icoor)    = M_FESpace.feBd().normal ( icoor, ig ) ;
            }

            // If the normal is orientated inward, we change its sign and swap the left and value of the solution
            if ( iElem == rightElement )
            {
                normal *= -1.;
                std::swap ( leftValue, rightValue );
            }

            const Real localFaceFlux = (*M_numericalFlux) ( leftValue[ 0 ],
                                                            rightValue[ 0 ],
                                                            normal,
                                                            iElem,
                                                            M_data.dataTime()->time(),
                                                            quadPoint (0),
                                                            quadPoint (1),
                                                            quadPoint (2) );

            // Update the local flux of the current face with the quadrature weight
            localFaceFluxWeight[0] += localFaceFlux * M_FESpace.feBd().wRootDetMetric ( ig );

        }

        /* Put in localFlux the vector L^{-1} * localFlux
           For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
        lapack.TRTRS ( param_L, param_N, param_N, NB, NBRHS, M_elmatMass[ iElem ].mat(), NB, localFaceFluxWeight, NB, INFO);
        ASSERT_PRE ( !INFO[0], "Lapack Computation M_elvecSource = LB^{-1} rhs is not achieved." );

        /* Put in localFlux the vector L^{-T} * localFlux
           For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
        lapack.TRTRS ( param_L, param_T, param_N, NB, NBRHS, M_elmatMass[ iElem ].mat(), NB, localFaceFluxWeight, NB, INFO);
        ASSERT_PRE ( !INFO[0], "Lapack Computation M_elvecSource = LB^{-1} rhs is not achieved." );

        // Add to the local flux the local flux of the current face
        M_localFlux += localFaceFluxWeight;

    }

} // localEvolve

// Apply the flux limiters locally.
template< typename Mesh, typename SolverType >
void
HyperbolicSolver< Mesh, SolverType >::
localAverage ( const UInt& /* iElem */ )
{

    ;

} // localAverage

} // namespace LifeV

#endif /*_HYPERBOLICSOLVER_H_ */
