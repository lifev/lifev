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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief This file contains an Oseen equation solver class.

    @author Gilles Fourestey <gilles.fourestey@imag.fr>
            Simone Deparis   <simone.deparis@epfl.ch>
    @contributor Zhen Wang <zhen.wang@emory.edu>

    @date 01-06-2007

    This file contains an Oseen equation solver class.
    The resulting linear systems are solved by GMRES on the full
    matrix ( u and p coupled ).

 */


#ifndef OSEEN_H
#define OSEEN_H 1

#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifealg/EpetraPreconditioner.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/AztecOOPreconditioner.hpp>
#include <life/lifealg/EpetraMap.hpp>

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>

#include <life/lifecore/chrono.hpp>

#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/sobolevNorms.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/postProc.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifesolver/nsipterms.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>

#include <boost/shared_ptr.hpp>

#include <list>

namespace LifeV
{
//! @class Oseen
/*!
    @brief This class contains an Oseen equation solver.

    @author Gilles Fourestey <gilles.fourestey@imag.fr>
            Simone Deparis   <simone.deparis@epfl.ch>
    @contributor Zhen Wang <zhen.wang@emory.edu>

 */

template< typename MeshType, typename SolverType = LifeV::SolverTrilinos >
class Oseen
{

public:

    //! @name Public Types
    //@{

    typedef MeshType                              mesh_Type;
    typedef SolverType                            linearSolver_Type;
    typedef DataNavierStokes                      data_Type;

    typedef boost::function<Real ( const Real& t, const Real& x, const Real& y,
                                   const Real& z, const ID& i )> function_Type;

    typedef boost::function<Real ( const Real& t, const Real& x, const Real& y,
                                   const Real& z, const ID& i )> source_Type;

    typedef BCHandler                             bcHandler_Type;
    typedef boost::shared_ptr<bcHandler_Type>     bcHandlerPtr_Type;

    typedef typename linearSolver_Type::matrix_type     matrix_Type;
    typedef boost::shared_ptr<matrix_Type>        matrixPtr_Type;
    typedef typename linearSolver_Type::vector_type     vector_Type;
    typedef boost::shared_ptr<vector_Type>        vectorPtr_Type;

    typedef typename linearSolver_Type::prec_raw_type   preconditioner_Type;
    typedef typename linearSolver_Type::prec_type       preconditionerPtr_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    Oseen();

    //! Constructor
    /*!
        @param dataType DataNavierStokes class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param communicator MPI communicator
        @param lagrangeMultiplier Lagrange multiplier
     */
    Oseen( boost::shared_ptr<data_Type>    dataType,
           FESpace<mesh_Type, EpetraMap>&  velocityFESpace,
           FESpace<mesh_Type, EpetraMap>&  pressureFESpace,
           boost::shared_ptr<Epetra_Comm>& communicator,
           const Int                       lagrangeMultiplier = 0 );

    //! Constructor
    /*!
        @param dataType DataNavierStokes class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param communicator MPI communicator
        @param monolithicMap EpetraMap class
        @param offset
     */
    Oseen( boost::shared_ptr<data_Type>    dataType,
           FESpace<mesh_Type, EpetraMap>&  velocityFESpace,
           FESpace<mesh_Type, EpetraMap>&  pressureFESpace,
           boost::shared_ptr<Epetra_Comm>& communicator,
           const EpetraMap                 monolithicMap,
           const UInt                      offset = 0 );

    //! Constructor
    /*!
        @param dataType DataNavierStokes class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param lagrangeMultipliers (lagrange multipliers for the flux problem with rufaec flag)
        @param communicator MPI communicator
     */
    Oseen( boost::shared_ptr<data_Type>    dataType,
           FESpace<mesh_Type, EpetraMap>&  velocityFESpace,
           FESpace<mesh_Type, EpetraMap>&  pressureFESpace,
           const std::vector<Int>&         lagrangeMultipliers,
           boost::shared_ptr<Epetra_Comm>& communicator );

    //! virtual destructor
    virtual ~Oseen();

    //@}

    //! @name Methods
    //@{

    //! Set up data from GetPot
    /*!
        @param dataFile GetPot object
     */
    virtual void setUp( const GetPot& dataFile );

    //! Initialize with velocityFunction and pressureFunction
    /*!
        @param velocityFunction
        @param pressureFunction
     */
    void initialize( const function_Type& velocityFunction, const function_Type& pressureFunction );

    //! Initialize with velocityInitialGuess and pressureInitialGuess
    /*!
        @param velocityInitialGuess
        @param pressureInitialGuess
     */
    void initialize( const vector_Type& velocityInitialGuess, const vector_Type& pressureInitialGuess );

    //! Initialize with velocityAndPressure
    /*!
        @param velocityAndPressure
     */
    void initialize( const vector_Type& velocityAndPressure );

    //! Build linear system.
    virtual void buildSystem();

    //! Update system
    /*!
        @param alpha
        @param betaVector
        @param sourceVector
     */
    virtual void updateSystem( const Real         alpha,
                               const vector_Type& betaVector,
                               const vector_Type& sourceVector );

    //! Update system
    /*!
        @param alpha
        @param betaVector
        @param sourceVector
        @param matrix
        @param un
     */
    virtual void updateSystem( const Real         alpha,
                               const vector_Type& betaVector,
                               const vector_Type& sourceVector,
                               matrixPtr_Type     matrix,
                               vectorPtr_Type     un );

    //! Update stabilization term
    /*!
        @param matrixFull
     */
    void updateStabilization( matrix_Type& matrixFull );

    //! Update the right hand side
    /*!
        @param rightHandSide right hand side
     */
    virtual void updateRightHandSide( const vector_Type& rightHandSide )
    {
        M_rightHandSideNoBC = rightHandSide;
        M_rightHandSideNoBC.globalAssemble();
    }

    //! Update convective term, boundary condition and solve the linearized ns system
    /*!
        @param bcHandler BC handler
     */
    virtual void iterate( bcHandler_Type& bcHandler );

    //! Reduce the local solution in global vectors
    /*!
        @param velocity
        @param pressure
     */
    void reduceSolution( Vector& velocity, Vector& pressure );

    //! Reduce the residual
    /*!
        @param residual
     */
    void reduceResidual( Vector& residual );

    //! Set a block preconditioner
    /*!
        @blockPrecconditioner Block preconditioner
     */
    void setBlockPreconditioner( matrixPtr_Type blockPreconditioner );

    //! Update and return the coefficient matrix
    /*!
        @param matrixFull The coefficient matrix
     */
    void getFluidMatrix( matrix_Type& matrixFull );

    //! Set up post processing
    void setupPostProc( const entityFlag_Type& flag, const mesh_Type meshPart );

    //! Compute area on a boundary face with given flag
    /*!
        @param  flag
        @return area
     */
    Real area( const entityFlag_Type& flag );

    //! Compute flux on a boundary face with given flag and a given solution
    /*!
        @param  flag
        @param  solution
        @return flux
     */
    Real flux( const entityFlag_Type& flag, const vector_Type& solution );

    //! Compute flux on a boundary face with given flag
    /*!
        @param flag
        @return flux
     */
    Real flux( const entityFlag_Type& flag );

    //! Compute average pressure on a boundary face with given flag and a given solution
    /*!
        @param  flag
        @param  solution
        @return average pressure
     */
    Real pressure( const entityFlag_Type& flag, const vector_Type& solution );

    //! Compute average pressure on a boundary face with given flag
    /*!
        @param flag
        @return average pressure
     */
    Real pressure( const entityFlag_Type& flag );

    //! Get the Lagrange multiplier related to a flux imposed on a given part of the boundary
    /*!
        @param flag      Flag of the boundary face associated with the flux
                         and the Lagrange multiplier we want.
        @param bcHandler BChandler containing the boundary conditions of the problem.
        @return          Lagrange multiplier
     */
    Real lagrangeMultiplier( const entityFlag_Type& flag, bcHandler_Type& bcHandler );

    //! Get the Lagrange multiplier related to a flux imposed on a given part of the boundary
    /*!
        @param flag      Flag of the boundary face associated
                         with the flux and the Lagrange multiplier we want.
        @param bcHandler BChandler containing the boundary conditions of the problem.
        @param solution  Vector containing the solution of the problem
                         (and also the Lagrange multipliers at the end).
        @return          Lagrange multiplier
     */
    Real lagrangeMultiplier( const entityFlag_Type&  flag,
                             bcHandler_Type& bcHandler,
                             const vector_Type& solution );

    //! Reset the preconditioner.
    /*!
        @param reset Reset preconditioner.
     */
    void resetPreconditioner( bool reset = true )
    {
        if ( reset )
            M_linearSolver.precReset();
    }

    //! Reset stabilization matrix at the same time as the preconditioner
    void resetStabilization()
    {
        M_resetStabilization = true;
    }

    //! Update
    void updateUn()
    {
        *M_un = *M_solution;
    }

    //! Update for the monolithic
    void updateUn( const vector_Type& solution )
    {
        *M_un = solution;
    }

    //! Display general information about the content of the class
    /*!
        @param output specify the output format (std::cout by default)
     */
    void showMe( std::ostream& output = std::cout ) const;

    //@}

    //! @name Set Methods
    //@{

    //! Set
    /*!
        @param recomputeMatrix
     */
    void setRecomputeMatrix( const bool& recomputeMatrix )
    {
        M_recomputeMatrix = recomputeMatrix;
    }

    //! set the source term functor
    /*!
        @param source
     */
    void setSourceTerm( source_Type source )
    {
        M_source = source;
    }

    //! Set the tolerance and the maximum number of iterations of the linear solver
    /*!
        @param tolerance Tolerance
        @param maxIteration maximum number of iterations
     */
    void setTolMaxIteration( const Real& tolerance, const Int& maxIteration = -1 );

    //@}

    //! @name Get Methods
    //@{

    //! Return the density of the fluid
    /*!
        @return Density of the fluid
     */
    const Real& density() const
    {
        return M_dataNavierStokes->density();
    }
    //! Return the viscosity of the fluid
    /*!
        @return Viscosity of the fluid
     */
    const Real& viscosity() const
    {
        return M_dataNavierStokes->viscosity();
    }

    //! Return the local solution vector
    /*!
        @return vectorPtr_Type Solution vector
     */
    const vectorPtr_Type& solution() const
    {
        return M_solution;
    }

    //! Return the local residual vector
    /*!
        @return Residual vector
     */
    const vector_Type& residual() const
    {
        return M_residual;
    }

    //! Return velocity FE space
    /*!
        @return velocity FE space
     */
    FESpace<mesh_Type, EpetraMap>& velocityFESpace()
    {
        return M_velocityFESpace;
    }

    const FESpace<mesh_Type, EpetraMap>& velocityFESpace() const
    {
        return M_velocityFESpace;
    }

    //! Return pressure FE space
    /*!
        @return pressure FE space
     */
    FESpace<mesh_Type, EpetraMap>& pressureFESpace()
    {
        return M_pressureFESpace;
    }

    const FESpace<mesh_Type, EpetraMap>& pressureFESpace() const
    {
        return M_pressureFESpace;
    }

    //! Get the source term
    /*!
        @return Source term
     */
    const source_Type& sourceTerm() const
    {
        return M_source;
    }

    //! Returns the post processing structure
    /*!
        @return Post processing
    */
    PostProc<mesh_Type>& postProcessing()
    {
        return *M_postProcessing;
    }

    const PostProc<mesh_Type>& postProcessing() const
    {
        return *M_postProcessing;
    }

    //! Return EpetraMap.
    /*!
        @return EpetraMap
     */
    const EpetraMap& getMap() const
    {
        return M_localMap;
    }

    //! Return Epetra communicator
    /*!
        @return Epetra communicator
     */
    const boost::shared_ptr<Epetra_Comm>& comm() const
    {
        return M_Displayer.comm();
    }

    //! Return displayer
    /*!
        @return
     */
    const Displayer& getDisplayer() const
    {
        return M_Displayer;
    }

    //! Return
    /*!
        @return recomputeMatrix
     */
    const bool& recomputeMatrix() const
    {
        return M_recomputeMatrix;
    }

    //! Return matrix without boundary conditions
    /*!
        @return Matrix without boundary conditions
     */
    matrix_Type& matrixNoBC()
    {
        return *M_matrixNoBC;
    }

    const matrix_Type& matrixNoBC() const
    {
        return *M_matrixNoBC;
    }

    //! Return mass matrix
    /*!
        @return Mass matrix
     */
    matrix_Type& matrixMass()
    {
        return *M_matrixMass;
    }

    const matrix_Type& matrixMass() const
    {
        return *M_matrixMass;
    }

    //@}

    //@{ unused methods

    //! Set up post processing structures
    void postProcessingSetArea();

    //! Set up post processing
    void postProcessingSetNormal();

    //! Set up post processing
    void postProcessingSetPhi();

    //! Return a bool value if using diagonal block preconditioner
    bool getIsDiagonalBlockPreconditioner()
    {
        return M_isDiagonalBlockPreconditioner;
    }

    const bool& getIsDiagonalBlockPreconditioner() const
    {
        return M_isDiagonalBlockPreconditioner;
    }

    //@}


protected:

    //! @name Constructor
    //@{

    //! Empty copy constructor
    Oseen( const Oseen& oseen);

    //@}

    //! @name Private Methods
    //@{

    //! Removes mean of component of vector x
    /*!
        @param x
        @return
     */
    Real removeMean( vector_Type& x );

    //! Apply boundary conditions.
    /*!
        @param matrix
        @param rightHandSide
        @param bcHandler
     */
    void applyBoundaryConditions( matrix_Type&        matrix,
                                  vector_Type&        rightHandSide,
                                  bcHandler_Type& bcHandler );

    //! Echo message.
    /*!
        @param message
     */
    void echo( std::string message );

    //! Return the dim of velocity FE space
    const UInt& dimVelocity() const
    {
        return M_velocityFESpace.dim();
    }

    //! Return the dim of pressure FE space
    const UInt& dimPressure() const
    {
        return M_pressureFESpace.dim();
    }

    //@}

    //@{ deprecated methods

    const UInt&  __attribute__ ((__deprecated__)) dim_u( ) const
    {
        return dimVelocity( );
    }

    const UInt&  __attribute__ ((__deprecated__)) dim_p( ) const
    {
        return dimPressure( );
    }

    //@}

    //private members

    //! data for Navier-Stokes solvers
    boost::shared_ptr<data_Type>   M_dataNavierStokes;

    // FE spaces
    FESpace<mesh_Type, EpetraMap>& M_velocityFESpace;
    FESpace<mesh_Type, EpetraMap>& M_pressureFESpace;

    //! MPI communicator
    Displayer                      M_Displayer;

    EpetraMap                      M_localMap;

    //! mass matrix
    matrixPtr_Type                 M_matrixMass;

    //! mass matrix
    matrixPtr_Type                 M_matrixMassPtr;

    //! Stokes matrix: nu*stiff
    matrixPtr_Type                 M_matrixStokes;

    //! matrix to be solved
//    matrixPtr_Type               M_matrixFull;

    //! matrix without boundary conditions
    matrixPtr_Type                 M_matrixNoBC;

    //! stabilization matrix
    matrixPtr_Type                 M_matrixStabilization;

    //! source term for Navier-Stokes equations
    source_Type                    M_source;

    //! Right hand side for the velocity component
    vector_Type                    M_rightHandSideNoBC;

    //! Global right hand side
    vector_Type                    M_rightHandSideFull;

    //! Global solution
    vectorPtr_Type                 M_solution;

    //! residual
    vector_Type                    M_residual;

    linearSolver_Type              M_linearSolver;

    bool                           M_steady;

    //! Postprocessing class
    boost::shared_ptr<PostProc<mesh_Type> > M_postProcessing;

    //! Stabilization
    bool                           M_stabilization;
    bool                           M_reuseStabilization;
    bool                           M_resetStabilization;
    Int                            M_iterReuseStabilization;

    details::IPStabilization<mesh_Type, Dof> M_ipStabilization;
    Real                           M_gammaBeta;
    Real                           M_gammaDiv;
    Real                           M_gammaPress;

    const function_Type*                M_betaFunction;

    bool                           M_divBetaUv;

    bool                           M_stiffStrain;

    //
    Real                           M_diagonalize;

    UInt                           M_count;

    bool                           M_recomputeMatrix;

    bool                           M_isDiagonalBlockPreconditioner;

    //! Elementary matrices and vectors
    ElemMat                        M_elementMatrixStiff;      // velocity Stokes
    ElemMat                        M_elementMatrixMass;       // velocity mass
    ElemMat                        M_elementMatrixPreconditioner;          // (p,q) bloc for preconditioners
    ElemMat                        M_elementMatrixDivergence;
    ElemMat                        M_elementMatrixGradient;
    ElemVec                        M_elementRightHandSide;           // Elementary right hand side
    matrixPtr_Type                 M_blockPreconditioner;
    ElemVec                        M_wLoc;
    ElemVec                        M_uLoc;
    boost::shared_ptr<vector_Type> M_un;

}; // class Oseen



// ===================================================
// Constructors & Destructor
// ===================================================

template<typename MeshType, typename SolverType>
Oseen<MeshType, SolverType>::
Oseen( boost::shared_ptr<data_Type>    dataType,
       FESpace<mesh_Type, EpetraMap>&  velocityFESpace,
       FESpace<mesh_Type, EpetraMap>&  pressureFESpace,
       boost::shared_ptr<Epetra_Comm>& communicator,
       const Int                       lagrangeMultiplier ):
        M_dataNavierStokes       ( dataType ),
        M_velocityFESpace        ( velocityFESpace ),
        M_pressureFESpace        ( pressureFESpace ),
        M_Displayer              ( communicator ),
        M_localMap               ( M_velocityFESpace.map() + M_pressureFESpace.map() + lagrangeMultiplier),
        M_matrixMass             ( ),
        M_matrixMassPtr          ( ),
        M_matrixStokes           ( ),
        M_matrixNoBC             ( ),
        M_matrixStabilization    ( ),
        M_rightHandSideNoBC      ( M_localMap ),
        M_rightHandSideFull      ( M_localMap ),
        M_solution               ( new vector_Type( M_localMap ) ),
        M_residual               ( M_localMap ),
        M_linearSolver           ( communicator ),
        M_steady                 ( ),
        M_postProcessing         ( new PostProc<mesh_Type>( M_velocityFESpace.mesh(),
                                                            &M_velocityFESpace.feBd(),
                                                            &M_velocityFESpace.dof(),
                                                            &M_pressureFESpace.feBd(),
                                                            &M_pressureFESpace.dof(),
                                                            M_localMap ) ),
        M_stabilization          ( false ),
        M_reuseStabilization     ( false ),
        M_resetStabilization     ( false ),
        M_iterReuseStabilization ( -1 ),
//        M_ipStabilization        ( M_velocityFESpace.mesh(),
//                                   M_velocityFESpace.dof(),
//                                   M_velocityFESpace.refFE(),
//                                   M_velocityFESpace.feBd(),
//                                   M_velocityFESpace.qr(),
//                                   0., 0., 0.,
//                                   M_dataNavierStokes->viscosity() ),
        M_betaFunction           ( 0 ),
        M_divBetaUv              ( false ),
        M_stiffStrain            ( false ),
        M_diagonalize            ( false ),
        M_count                  ( 0 ),
        M_recomputeMatrix        ( false ),
        M_elementMatrixStiff     ( M_velocityFESpace.fe().nbFEDof(), nDimensions, nDimensions ),
        M_elementMatrixMass      ( M_velocityFESpace.fe().nbFEDof(), nDimensions, nDimensions ),
        M_elementMatrixPreconditioner ( M_pressureFESpace.fe().nbFEDof(), 1, 1 ),
        M_elementMatrixDivergence ( M_pressureFESpace.fe().nbFEDof(), 1, 0,
                                    M_velocityFESpace.fe().nbFEDof(), 0, nDimensions ),
        M_elementMatrixGradient  ( M_velocityFESpace.fe().nbFEDof(), nDimensions, 0,
                                   M_pressureFESpace.fe().nbFEDof(), 0, 1 ),
        M_elementRightHandSide   ( M_velocityFESpace.fe().nbFEDof(), nDimensions ),
        M_blockPreconditioner    ( ),
        M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), nDimensions ),
        M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), nDimensions ),
        M_un                     ( new vector_Type(M_localMap) )
{
    M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() );
    M_ipStabilization.setFeSpaceVelocity(M_velocityFESpace);
    M_ipStabilization.setViscosity(M_dataNavierStokes->viscosity() );
}

template<typename MeshType, typename SolverType>
Oseen<MeshType, SolverType>::
Oseen( boost::shared_ptr<data_Type>    dataType,
       FESpace<mesh_Type, EpetraMap>&  velocityFESpace,
       FESpace<mesh_Type, EpetraMap>&  pressureFESpace,
       boost::shared_ptr<Epetra_Comm>& communicator,
       EpetraMap                       monolithicMap,
       UInt                            offset ):
        M_dataNavierStokes       ( dataType ),
        M_velocityFESpace        ( velocityFESpace ),
        M_pressureFESpace        ( pressureFESpace ),
        M_Displayer              ( communicator ),
        M_localMap               ( monolithicMap ),
        M_matrixMass             ( ),
        M_matrixStokes           ( ),
        M_matrixNoBC             ( ),
        M_matrixStabilization    ( ),
        M_rightHandSideNoBC      ( M_localMap ),
        M_rightHandSideFull      ( M_localMap ),
        M_solution               ( ),
        M_residual               ( M_localMap ),
        M_linearSolver           ( communicator ),
        M_postProcessing         ( new PostProc<mesh_Type>(M_velocityFESpace.mesh(),
                                                           &M_velocityFESpace.feBd(),
                                                           &M_velocityFESpace.dof(),
                                                           &M_pressureFESpace.feBd(),
                                                           &M_pressureFESpace.dof(),
                                                           M_localMap ) ),
        M_stabilization          ( false ),
        M_reuseStabilization     ( false ),
        M_resetStabilization     ( false ),
        M_iterReuseStabilization ( -1 ),
//        M_ipStabilization        ( M_velocityFESpace.mesh(),
//                                   M_velocityFESpace.dof(),
//                                   M_velocityFESpace.refFE(),
//                                   M_velocityFESpace.feBd(),
//                                   M_velocityFESpace.qr(),
//                                   0., 0., 0.,
//                                   M_dataNavierStokes->viscosity() ),
        M_betaFunction           ( 0 ),
        M_divBetaUv              ( false ),
        M_stiffStrain            ( false ),
        M_diagonalize            ( false ),
        M_count                  ( 0 ),
        M_recomputeMatrix        ( false ),
        M_elementMatrixStiff     ( M_velocityFESpace.fe().nbFEDof(), nDimensions, nDimensions ),
        M_elementMatrixMass      ( M_velocityFESpace.fe().nbFEDof(), nDimensions, nDimensions ),
        M_elementMatrixPreconditioner                 ( M_pressureFESpace.fe().nbFEDof(), 1, 1 ),
        M_elementMatrixDivergence ( M_pressureFESpace.fe().nbFEDof(), 1, 0,
                                    M_velocityFESpace.fe().nbFEDof(), 0, nDimensions ),
        M_elementMatrixGradient  ( M_velocityFESpace.fe().nbFEDof(), nDimensions, 0,
                                   M_pressureFESpace.fe().nbFEDof(), 0, 1 ),
        M_elementRightHandSide   ( M_velocityFESpace.fe().nbFEDof(), nDimensions ),
        M_blockPreconditioner    ( ),
        M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), nDimensions ),
        M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), nDimensions ),
        M_un                     ( new vector_Type(M_localMap) )
{
    M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() );
    M_ipStabilization.setFeSpaceVelocity(M_velocityFESpace);
    M_ipStabilization.setViscosity(M_dataNavierStokes->viscosity() );
}

template<typename MeshType, typename SolverType>
Oseen<MeshType, SolverType>::
Oseen( boost::shared_ptr<data_Type>    dataType,
       FESpace<mesh_Type, EpetraMap>&  velocityFESpace,
       FESpace<mesh_Type, EpetraMap>&  pressureFESpace,
       const std::vector<Int> &        lagrangeMultipliers,
       boost::shared_ptr<Epetra_Comm>& communicator ):
        M_dataNavierStokes       ( dataType ),
        M_velocityFESpace        ( velocityFESpace ),
        M_pressureFESpace        ( pressureFESpace ),
        M_Displayer              ( communicator ),
        M_localMap               ( M_velocityFESpace.map() + M_pressureFESpace.map() + lagrangeMultipliers ),
        M_matrixMass             ( ),
        M_matrixStokes           ( ),
        M_matrixNoBC             ( ),
        M_matrixStabilization    ( ),
        M_rightHandSideNoBC      ( M_localMap ),
        M_rightHandSideFull      ( M_localMap ),
        M_solution               ( new vector_Type( M_localMap ) ),
        M_residual               ( M_localMap ),
        M_linearSolver           ( ),
        M_postProcessing         ( new PostProc<mesh_Type>(M_velocityFESpace.mesh(),
                                                           &M_velocityFESpace.feBd(),
                                                           &M_velocityFESpace.dof(),
                                                           &M_pressureFESpace.feBd(),
                                                           &M_pressureFESpace.dof(),
                                                           M_localMap ) ),
        M_stabilization          ( false ),
        M_reuseStabilization     ( false ),
        M_resetStabilization     ( false ),
        M_iterReuseStabilization ( -1 ),
//        M_ipStabilization        ( M_velocityFESpace.mesh(),
//                                   M_velocityFESpace.dof(),
//                                   M_velocityFESpace.refFE(),
//                                   M_velocityFESpace.feBd(),
//                                   M_velocityFESpace.qr(),
//                                   0., 0., 0.,
//                                   M_dataNavierStokes->viscosity() ),
        M_betaFunction           ( 0 ),
        M_divBetaUv              ( false ),
        M_stiffStrain            ( false ),
        M_diagonalize            ( false ),
        M_count                  ( 0 ),
        M_recomputeMatrix        ( false ),
        M_elementMatrixStiff     ( M_velocityFESpace.fe().nbFEDof(), nDimensions, nDimensions ),
        M_elementMatrixMass      ( M_velocityFESpace.fe().nbFEDof(), nDimensions, nDimensions ),
        M_elementMatrixPreconditioner ( M_pressureFESpace.fe().nbFEDof(), 1, 1 ),
        M_elementMatrixDivergence ( M_pressureFESpace.fe().nbFEDof(), 1, 0,
                                    M_velocityFESpace.fe().nbFEDof(), 0, nDimensions ),
        M_elementMatrixGradient  ( M_velocityFESpace.fe().nbFEDof(), nDimensions, 0,
                                   M_pressureFESpace.fe().nbFEDof(), 0, 1 ),
        M_elementRightHandSide   ( M_velocityFESpace.fe().nbFEDof(), nDimensions ),
        M_blockPreconditioner    ( ),
        M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), nDimensions ),
        M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), nDimensions ),
        M_un                     ( new vector_Type(M_localMap) )
{
    M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() );
    M_ipStabilization.setFeSpaceVelocity(M_velocityFESpace);
    M_ipStabilization.setViscosity(M_dataNavierStokes->viscosity() );
}

template<typename MeshType, typename SolverType>
Oseen<MeshType, SolverType>::
~Oseen()
{

}


// ===================================================
// Methods
// ===================================================

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::setUp( const GetPot& dataFile )
{

    M_linearSolver.setUpPrec( dataFile, "fluid/prec" );
    M_linearSolver.setDataFromGetPot( dataFile, "fluid/solver" );

    M_steady        = dataFile( "fluid/miscellaneous/steady", 0 );

    M_gammaBeta     = dataFile( "fluid/ipstab/gammaBeta",  0. );
    M_gammaDiv      = dataFile( "fluid/ipstab/gammaDiv",   0. );
    M_gammaPress    = dataFile( "fluid/ipstab/gammaPress", 0. );
    M_reuseStabilization     = dataFile( "fluid/ipstab/reuse", false );
    M_iterReuseStabilization = dataFile( "fluid/ipstab/max_iter_reuse",
                                         static_cast<Int> ( M_linearSolver.MaxIter() * 8./10. ) );

    // Energetic stabilization term
    M_divBetaUv   = dataFile( "fluid/space_discretization/div_beta_u_v",false);
    // Enable grad( u )^T in stress tensor
    M_stiffStrain = dataFile( "fluid/space_discretization/stiff_strain",false);
    M_diagonalize = dataFile( "fluid/space_discretization/diagonalize", 1. );
    M_isDiagonalBlockPreconditioner = dataFile( "fluid/diagonalBlockPrec", false );

    //    M_linearSolver.setAztecooPreconditioner( dataFile, "fluid/solver" );

    M_ipStabilization.setGammaBeta ( M_gammaBeta );
    M_ipStabilization.setGammaDiv  ( M_gammaDiv );
    M_ipStabilization.setGammaPress( M_gammaPress );
}


template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::
initialize( const function_Type& velocityFunction, const function_Type& pressureFunction )
{
    vector_Type velocityInitialGuess( M_velocityFESpace.map() );
    M_velocityFESpace.interpolate( velocityFunction,
                                   velocityInitialGuess,
                                   M_dataNavierStokes->dataTime()->time() );

    vector_Type pressureInitialGuess( M_pressureFESpace.map() );
    M_pressureFESpace.interpolate( pressureFunction,
                                   pressureInitialGuess,
                                   M_dataNavierStokes->dataTime()->time() );

    initialize( velocityInitialGuess, pressureInitialGuess );
}


template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::
initialize( const vector_Type& velocityInitialGuess, const vector_Type& pressureInitialGuess )
{

    *M_solution = velocityInitialGuess;
    *M_un = velocityInitialGuess;
    M_solution->add( pressureInitialGuess, nDimensions * M_velocityFESpace.dof().numTotalDof() );
    M_un->add( pressureInitialGuess, nDimensions * M_velocityFESpace.dof().numTotalDof() );

}

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::
initialize( const vector_Type& velocityAndPressure )
{

    *M_un = velocityAndPressure;
    if ( M_solution.get() )
        *M_solution = velocityAndPressure;

}

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::buildSystem()
{
    M_matrixMass.reset  ( new matrix_Type( M_localMap ) );
    M_matrixStokes.reset( new matrix_Type( M_localMap ) );

    M_Displayer.leaderPrint( "  F-  Computing constant matrices ...          " );

    Chrono chrono;

    Chrono chronoDer;
    Chrono chronoStiff;
    Chrono chronoMass;
    Chrono chronoGrad;

    Chrono chronoStiffAssemble;
    Chrono chronoMassAssemble;
    Chrono chronoGradAssemble;
    Chrono chronoDivAssemble;
    Chrono chronoStab;
    Chrono chronoZero;

    // Number of velocity components
    UInt numVelocityComponent = nDimensions;

    // Elementary computation and matrix assembling
    // Loop on elements

    UInt velocityTotalDof   = M_velocityFESpace.dof().numTotalDof();
//    UInt pressureTotalDof = M_pressureFESpace.dof().numTotalDof();

    if ( M_isDiagonalBlockPreconditioner == true )
    {
        M_blockPreconditioner.reset( new matrix_Type( M_localMap ) );
    }
    chrono.start();

    for ( UInt iVolume = 1; iVolume <= M_velocityFESpace.mesh()->numVolumes(); iVolume++ )
    {
        chronoDer.start();
        // just to provide the id number in the assem_mat_mixed
        M_pressureFESpace.fe().update( M_velocityFESpace.mesh()->volumeList( iVolume ) );
        // just to provide the id number in the assem_mat_mixed
        // M_pressureFESpace.fe().updateFirstDeriv( M_velocityFESpace.mesh()->volumeList( iVolume ) );
        M_velocityFESpace.fe().updateFirstDeriv( M_velocityFESpace.mesh()->volumeList( iVolume ) );

        chronoDer.stop();

        chronoZero.start();
        M_elementMatrixStiff.zero();
        M_elementMatrixMass.zero();
        M_elementMatrixPreconditioner.zero();
        M_elementMatrixDivergence.zero();
        M_elementMatrixGradient.zero();
        chronoZero.stop();

        // stiffness matrix
        chronoStiff.start();
        if ( M_stiffStrain )
            stiff_strain( 2.0*M_dataNavierStokes->viscosity(),
                          M_elementMatrixStiff,
                          M_velocityFESpace.fe() );
        else
            stiff( M_dataNavierStokes->viscosity(),
                   M_elementMatrixStiff,
                   M_velocityFESpace.fe(), 0, 0, nDimensions );
        //stiff_div( 0.5*M_velocityFESpace.fe().diameter(), M_elementMatrixStiff, M_velocityFESpace.fe() );
        chronoStiff.stop();

        // mass matrix
        if ( !M_steady )
        {
            chronoMass.start();
            mass( M_dataNavierStokes->density(),
                  M_elementMatrixMass,
                  M_velocityFESpace.fe(), 0, 0, nDimensions );
            chronoMass.stop();
        }

        for ( UInt iComponent = 0; iComponent < numVelocityComponent; iComponent++ )
        {
            // stiffness matrix
            chronoStiffAssemble.start();
            if ( M_isDiagonalBlockPreconditioner == true )
            {
                assembleMatrix( *M_blockPreconditioner,
                                M_elementMatrixStiff,
                                M_velocityFESpace.fe(),
                                M_velocityFESpace.fe(),
                                M_velocityFESpace.dof(),
                                M_velocityFESpace.dof(),
                                iComponent, iComponent,
                                iComponent * velocityTotalDof, iComponent * velocityTotalDof);
            }
            else
            {
                if ( M_stiffStrain ) // sigma = 0.5 * mu (grad( u ) + grad ( u )^T)
                {
                    for ( UInt jComp = 0; jComp < numVelocityComponent; jComp++ )
                    {
                        assembleMatrix( *M_matrixStokes,
                                        M_elementMatrixStiff,
                                        M_velocityFESpace.fe(),
                                        M_velocityFESpace.fe(),
                                        M_velocityFESpace.dof(),
                                        M_velocityFESpace.dof(),
                                        iComponent, jComp,
                                        iComponent * velocityTotalDof, jComp * velocityTotalDof);

                    }
                }
                else // sigma = mu grad( u )
                {
                    assembleMatrix( *M_matrixStokes,
                                    M_elementMatrixStiff,
                                    M_velocityFESpace.fe(),
                                    M_velocityFESpace.fe(),
                                    M_velocityFESpace.dof(),
                                    M_velocityFESpace.dof(),
                                    iComponent, iComponent,
                                    iComponent * velocityTotalDof, iComponent * velocityTotalDof);
                }
            }
            chronoStiffAssemble.stop();

            // mass matrix
            if ( !M_steady )
            {
                chronoMassAssemble.start();
                assembleMatrix( *M_matrixMass,
                                M_elementMatrixMass,
                                M_velocityFESpace.fe(),
                                M_velocityFESpace.fe(),
                                M_velocityFESpace.dof(),
                                M_velocityFESpace.dof(),
                                iComponent, iComponent,
                                iComponent * velocityTotalDof, iComponent * velocityTotalDof);
                chronoMassAssemble.stop();
            }

            // divergence
            chronoGrad.start();
            grad( iComponent, 1.0,
                  M_elementMatrixGradient,
                  M_velocityFESpace.fe(),
                  M_pressureFESpace.fe(),
                  iComponent, 0 );
            chronoGrad.stop();

            chronoGradAssemble.start();
            assembleMatrix( *M_matrixStokes,
                            M_elementMatrixGradient,
                            M_velocityFESpace.fe(),
                            M_pressureFESpace.fe(),
                            M_velocityFESpace.dof(),
                            M_pressureFESpace.dof(),
                            iComponent, 0,
                            iComponent * velocityTotalDof, numVelocityComponent * velocityTotalDof );
            chronoGradAssemble.stop();

            chronoDivAssemble.start();
            assembleTransposeMatrix( *M_matrixStokes,
                                     -1.,
                                     M_elementMatrixGradient,
                                     M_pressureFESpace.fe(),
                                     M_velocityFESpace.fe(),
                                     M_pressureFESpace.dof(),
                                     M_velocityFESpace.dof(),
                                     0 , iComponent,
                                     numVelocityComponent * velocityTotalDof, iComponent * velocityTotalDof );
            chronoDivAssemble.stop();
        }
    }

    //    for (UInt ii = nDimensions*dimVelocity(); ii < nDimensions*dimVelocity() + dimPressure(); ++ii)
    //  M_matrixStokes->set_mat_inc( ii ,ii, 0. ); not scalable!!!

    if ( M_isDiagonalBlockPreconditioner == true )
    {
        M_blockPreconditioner->globalAssemble();
        *M_matrixStokes += *M_blockPreconditioner;
    }
    comm()->Barrier();

    chrono.stop();
    M_Displayer.leaderPrintMax( "done in " , chrono.diff() );

    M_Displayer.leaderPrint( "  F-  Finalizing the matrices ...              " );

    chrono.start();

    M_matrixStokes->globalAssemble();
    M_matrixMass->globalAssemble();

    chrono.stop();
    M_Displayer.leaderPrintMax( "done in " , chrono.diff() );

    if ( false )
        std::cout << " partial times:  \n"
                  << " Der            " << chronoDer.diffCumul() << " s.\n"
                  << " Zero           " << chronoZero.diffCumul() << " s.\n"
                  << " Stiff          " << chronoStiff.diffCumul() << " s.\n"
                  << " Stiff Assemble " << chronoStiffAssemble.diffCumul() << " s.\n"
                  << " Mass           " << chronoMass.diffCumul() << " s.\n"
                  << " Mass Assemble  " << chronoMassAssemble.diffCumul() << " s.\n"
                  << " Grad           " << chronoGrad.diffCumul() << " s.\n"
                  << " Grad Assemble  " << chronoGradAssemble.diffCumul() << " s.\n"
                  << " Div Assemble   " << chronoDivAssemble.diffCumul() << " s.\n"
                  << std::endl;

}

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::
updateSystem( const Real         alpha,
              const vector_Type& betaVector,
              const vector_Type& sourceVector )
{
    if ( M_matrixNoBC.get() )
        M_matrixNoBC.reset( new matrix_Type( M_localMap, M_matrixNoBC->meanNumEntries() ) );
    else
        M_matrixNoBC.reset( new matrix_Type( M_localMap ) );

    updateSystem( alpha, betaVector, sourceVector, M_matrixNoBC, M_un );
    if ( alpha != 0. )
        M_matrixNoBC->globalAssemble();

}

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::
updateSystem( const Real         alpha,
              const vector_Type& betaVector,
              const vector_Type& sourceVector,
              matrixPtr_Type     matrixNoBC,
              vectorPtr_Type     un )
{
    Chrono chrono;

    // clearing pressure mass matrix in case we need it in removeMean;
    M_matrixMassPtr.reset( );


    M_Displayer.leaderPrint( "  F-  Updating mass term on right hand side... " );

    chrono.start();

    UInt velocityTotalDof   = M_velocityFESpace.dof().numTotalDof();
//    UInt pressureTotalDof = M_pressureFESpace.dof().numTotalDof();

    // Right hand side for the velocity at time

    updateRightHandSide( sourceVector );

    chrono.stop();

    M_Displayer.leaderPrintMax( "done in ", chrono.diff() );


    //    M_updated = false;

    if ( M_recomputeMatrix )
        buildSystem();

    M_Displayer.leaderPrint( "  F-  Copying the matrices ...                 " );

    chrono.start();

    if ( M_isDiagonalBlockPreconditioner == true )
    {
        matrixPtr_Type tempMatrix( M_blockPreconditioner );
        M_blockPreconditioner.reset( new matrix_Type( M_localMap,
                                                      M_blockPreconditioner->meanNumEntries() ) );
        *M_blockPreconditioner += *tempMatrix;
    }


    chrono.stop();
    M_Displayer.leaderPrintMax( "done in " , chrono.diff() );


    UInt numVelocityComponent = nDimensions;

    //! managing the convective term

    Real normInf;
    betaVector.normInf( &normInf );

    if ( normInf != 0. )
    {
        M_Displayer.leaderPrint( "  F-  Sharing convective term ...              " );
        chrono.start();

        // vector with repeated nodes over the processors

        vector_Type betaVectorRepeated( betaVector, Repeated );
        vector_Type unRepeated( *un, Repeated );

        chrono.stop();

        M_Displayer.leaderPrintMax( "done in " , chrono.diff() );
        M_Displayer.leaderPrint( "  F-  Updating the convective terms ...        " );
        chrono.start();

        for ( UInt iVolume = 1; iVolume <= M_velocityFESpace.mesh()->numVolumes(); ++iVolume )
        {
            // just to provide the id number in the assem_mat_mixed
            M_pressureFESpace.fe().updateFirstDeriv( M_velocityFESpace.mesh()->volumeList( iVolume ) );
            //as updateFirstDer
            M_velocityFESpace.fe().updateFirstDeriv( M_velocityFESpace.mesh()->volumeList( iVolume ) );

            M_elementMatrixStiff.zero();

            UInt elementID = M_velocityFESpace.fe().currentLocalId();
            // Non linear term, Semi-implicit approach
            // M_elementRightHandSide contains the velocity values in the nodes
            for ( UInt iNode = 0 ; iNode < M_velocityFESpace.fe().nbFEDof() ; iNode++ )
            {
                UInt iLocal = M_velocityFESpace.fe().patternFirst( iNode );
                for ( UInt iComponent = 0; iComponent < numVelocityComponent; ++iComponent )
                {
                    UInt iGlobal = M_velocityFESpace.dof().localToGlobal( elementID, iLocal + 1 )
                                   + iComponent * dimVelocity();
                    M_elementRightHandSide.vec() [ iLocal + iComponent * M_velocityFESpace.fe().nbFEDof() ]
                    = betaVectorRepeated[iGlobal]; // BASEINDEX + 1

                    M_uLoc.vec() [ iLocal + iComponent * M_velocityFESpace.fe().nbFEDof() ]
                    = unRepeated(iGlobal);
                    M_wLoc.vec() [ iLocal + iComponent * M_velocityFESpace.fe().nbFEDof() ]
                    = unRepeated(iGlobal) - betaVectorRepeated(iGlobal);
                }
            }


            // ALE term: - rho div w u v
            mass_divw( - M_dataNavierStokes->density(),
                       M_wLoc,
                       M_elementMatrixStiff,
                       M_velocityFESpace.fe(), 0, 0, numVelocityComponent );

            // ALE stab implicit: 0.5 rho div u w v
            mass_divw( 0.5*M_dataNavierStokes->density(),
                       M_uLoc,
                       M_elementMatrixStiff,
                       M_velocityFESpace.fe(), 0, 0, numVelocityComponent );

            // Stabilising term: div u^n u v
            if ( M_divBetaUv )
                mass_divw( 0.5*M_dataNavierStokes->density(),
                           M_elementRightHandSide,
                           M_elementMatrixStiff,
                           M_velocityFESpace.fe(), 0, 0, numVelocityComponent );

            // compute local convective terms
            advection( M_dataNavierStokes->density(),
                       M_elementRightHandSide,
                       M_elementMatrixStiff,
                       M_velocityFESpace.fe(), 0, 0, numVelocityComponent );

            // loop on components
            for ( UInt iComponent = 0; iComponent < numVelocityComponent; ++iComponent )
            {
                // compute local convective term and assembling
                // grad( 0, M_elementRightHandSide, M_elementMatrixStiff, M_velocityFESpace.fe(),
                //       M_velocityFESpace.fe(), iComponent, iComponent );
                // grad( 1, M_elementRightHandSide, M_elementMatrixStiff, M_velocityFESpace.fe(),
                //       M_velocityFESpace.fe(), iComponent, iComponent );
                // grad( 2, M_elementRightHandSide, M_elementMatrixStiff, M_velocityFESpace.fe(),
                //       M_velocityFESpace.fe(), iComponent, iComponent );

                assembleMatrix( *matrixNoBC,
                                M_elementMatrixStiff,
                                M_velocityFESpace.fe(),
                                M_velocityFESpace.fe(),
                                M_velocityFESpace.dof(),
                                M_velocityFESpace.dof(),
                                iComponent, iComponent,
                                iComponent*velocityTotalDof, iComponent*velocityTotalDof );
            }
        }

        chrono.stop();
        M_Displayer.leaderPrintMax( "done in " , chrono.diff() );

        if ( M_stabilization &&
                ( M_resetStabilization || !M_reuseStabilization || ( M_matrixStabilization.get() == 0 ) ) )
        {
            M_Displayer.leaderPrint( "  F-  Updating the stabilization terms ...     " );
            chrono.start();
            M_matrixStabilization.reset ( new matrix_Type( M_localMap ) );
            M_ipStabilization.apply( *M_matrixStabilization, betaVectorRepeated, false );
            M_matrixStabilization->globalAssemble();
            M_resetStabilization = false;
            chrono.stop();
            M_Displayer.leaderPrintMax( "done in " , chrono.diff() );
        }

    }
    else
    {
        if ( M_stabilization )
        {
            M_Displayer.leaderPrint( "  F-  Updating the stabilization terms ...     " );
            chrono.start();

            if ( M_resetStabilization || !M_reuseStabilization || ( M_matrixStabilization.get() == 0 ) )
            {
                M_matrixStabilization.reset( new matrix_Type( M_localMap ) );
                M_ipStabilization.apply( *M_matrixStabilization, betaVector, false );
                M_matrixStabilization->globalAssemble();
                M_resetStabilization = false;
            }
            else
            {
                M_Displayer.leaderPrint( "reusing stab." );
            }
            chrono.stop();
            M_Displayer.leaderPrintMax( "done in " , chrono.diff() );
        }
    }

    if ( alpha != 0. )
    {
        *matrixNoBC += (*M_matrixMass) * alpha;
        if ( M_isDiagonalBlockPreconditioner == true )
        {
            matrixNoBC->globalAssemble();
            *M_blockPreconditioner += *matrixNoBC;
            matrix_Type tempMatrix( *matrixNoBC );
            matrixNoBC.reset( new matrix_Type( M_localMap, tempMatrix.meanNumEntries() ) );
            *matrixNoBC += tempMatrix;
            M_blockPreconditioner->globalAssemble();
        }
    }
    *matrixNoBC += *M_matrixStokes;
}


template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::updateStabilization( matrix_Type& matrixFull )
{

    if ( M_stabilization )
    {
        matrixFull += *M_matrixStabilization;
    }

}

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::iterate( bcHandler_Type& bcHandler )
{

    Chrono chrono;

    // matrix and vector assembling communication
    M_Displayer.leaderPrint( "  F-  Updating the boundary conditions ...     " );

    chrono.start();

    M_matrixNoBC->globalAssemble();

    matrixPtr_Type matrixFull( new matrix_Type( M_localMap, M_matrixNoBC->meanNumEntries() ) );

    updateStabilization( *matrixFull );
    getFluidMatrix( *matrixFull );

    vector_Type rightHandSideFull ( M_rightHandSideNoBC );

//     matrixFull.reset( new matrix_Type( *M_matrixNoBC ) );
//     M_rightHandSideFull = M_rightHandSideNoBC;

    chrono.stop();

    M_Displayer.leaderPrintMax( "done in ", chrono.diff() );

    // boundary conditions update
    M_Displayer.leaderPrint("  F-  Applying boundary conditions ...         ");

    chrono.start();
    applyBoundaryConditions( *matrixFull, rightHandSideFull, bcHandler );

    matrixFull->globalAssemble();
    chrono.stop();

    M_Displayer.leaderPrintMax( "done in " , chrono.diff() );

    //matrixFull->spy( "matrixFull" );
    // solving the system
    M_linearSolver.setMatrix( *matrixFull );

    Int numIter = M_linearSolver.solveSystem( rightHandSideFull, *M_solution, matrixFull );

    // if the preconditioner has been rese the stab terms are to be updated
    if ( numIter < 0 || numIter > M_iterReuseStabilization )
    {
        resetStabilization();
    }

    M_residual  = M_rightHandSideNoBC;
    M_residual -= (*M_matrixNoBC) * (*M_solution);

    //M_residual.spy("residual");
} // iterate()


template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::reduceSolution( Vector& velocityVector, Vector& pressureVector )
{
    vector_Type solution( *M_solution, 0 );

    if ( false /*S_verbose*/ )
    {
        for ( UInt iDof = 0; iDof < nDimensions * dimVelocity(); ++iDof )
        {
            // BASEINDEX + 1
            velocityVector[ iDof ] = solution[ iDof + 1 ];
        }

        for ( UInt iDof = 0; iDof < dimPressure(); ++iDof )
        {
            // BASEINDEX + 1
            pressureVector[ iDof ] = solution[ iDof + nDimensions * dimVelocity() + 1 ];
        }
    }

}

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::reduceResidual( Vector& residualVector )
{
    vector_Type residual( M_residual, 0 );

    if ( false /*S_verbose*/ )
    {
        for ( UInt iDof = 0; iDof < nDimensions * dimVelocity(); ++iDof )
        {
            // BASEINDEX + 1
            residualVector[ iDof ] = residual[ iDof + 1 ];
        }

    }
}

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::setBlockPreconditioner( matrixPtr_Type blockPreconditioner )
{
    // blockPreconditioner.reset(new matrix_Type(M_monolithicMap, M_solid->getMatrixPtr()->getMeanNumEntries()));
    *blockPreconditioner += *M_blockPreconditioner;
}

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::getFluidMatrix( matrix_Type& matrixFull )
{
    M_matrixNoBC->globalAssemble();
    matrixFull += *M_matrixNoBC;
}

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::postProcessingSetArea()
{
    M_postProcessing->set_area();
}

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::postProcessingSetNormal()
{
    M_postProcessing->set_normal();
}

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::postProcessingSetPhi()
{
    M_postProcessing->set_phi();
}

template<typename MeshType, typename SolverType>
Real
Oseen<MeshType, SolverType>::flux( const entityFlag_Type& flag )
{
    return flux( flag, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
Oseen<MeshType, SolverType>::flux( const entityFlag_Type& flag,
                                   const vector_Type& solution )
{
    vector_Type velocityAndPressure( solution, Repeated );
    vector_Type velocity( this->M_velocityFESpace.map(), Repeated );
    velocity.subset( velocityAndPressure );

    return M_postProcessing->flux( velocity, flag );
}

template<typename MeshType, typename SolverType>
Real
Oseen<MeshType, SolverType>::area( const entityFlag_Type& flag )
{
    return M_postProcessing->measure( flag );
}

template<typename MeshType, typename SolverType>
Real
Oseen<MeshType, SolverType>::pressure( const entityFlag_Type& flag )
{
    return pressure( flag, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
Oseen<MeshType, SolverType>::pressure(const entityFlag_Type& flag,
                                      const vector_Type& solution)
{
    vector_Type velocityAndPressure( solution, Repeated );
    vector_Type pressure( this->M_pressureFESpace.map(), Repeated );
    pressure.subset( velocityAndPressure,
                     this->M_velocityFESpace.dim()*this->M_velocityFESpace.fieldDim() );

    // third argument is 1, to use the pressure finite element space (see PostProc docs)
    return M_postProcessing->average( pressure, flag, 1 )[0];
}

template<typename MeshType, typename SolverType>
Real
Oseen<MeshType, SolverType>::lagrangeMultiplier( const entityFlag_Type& flag,
                                                 bcHandler_Type& bcHandler )
{
    return lagrangeMultiplier( flag, bcHandler, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
Oseen<MeshType, SolverType>::lagrangeMultiplier( const entityFlag_Type&  flag,
                                                 bcHandler_Type& bcHandler,
                                                 const vector_Type& solution )
{
    // Create a list of Flux bcName_Type ??
    std::vector< bcName_Type > fluxBCVector = bcHandler.findAllBCWithType( Flux );
    bcName_Type fluxbcName_Type = bcHandler.findBCWithFlag( flag ).name();

    // Create a Repeated vector for looking to the lambda
    vector_Type velocityPressureLambda( solution, Repeated );

    // Find the index associated to the correct Lagrange multiplier
    for ( UInt lmIndex = 0; lmIndex < static_cast <UInt> ( fluxBCVector.size() ); ++lmIndex )
        if ( fluxbcName_Type.compare( fluxBCVector[ lmIndex ] ) == 0 )
            return velocityPressureLambda[3 * M_velocityFESpace.dof().numTotalDof()
                                          + M_pressureFESpace.dof().numTotalDof() + 1 + lmIndex];

    // If lmIndex has not been found a warning message is printed
    std::cout << "!!! Warning - Lagrange multiplier for Flux BC not found!" << std::endl;
    return 0;
}

template<typename MeshType, typename SolverType>
Real
Oseen<MeshType, SolverType>::removeMean( vector_Type& x )
{

    Chrono chrono;
    chrono.start();

    const UInt numVelocityComponent ( nDimensions );
    const UInt velocityTotalDof ( M_velocityFESpace.dof().numTotalDof() );


    if ( M_matrixMassPtr.get() == 0 )
        M_matrixMassPtr.reset( new matrix_Type( M_localMap ) );

    for ( UInt iVolume = 1; iVolume <= M_velocityFESpace.mesh()->numVolumes(); iVolume++ )
    {
        chrono.start();
        // just to provide the id number in the assem_mat_mixed
        M_pressureFESpace.fe().update( M_pressureFESpace.mesh()->volumeList( iVolume ) );

        M_elementMatrixPreconditioner.zero();
        // mass
        chrono.start();
        mass( 1, M_elementMatrixPreconditioner, M_pressureFESpace.fe(), 0, 0, nDimensions );
        chrono.stop();

        chrono.start();
        assembleMatrix( *M_matrixMassPtr,
                        M_elementMatrixPreconditioner,
                        M_pressureFESpace.fe(),
                        M_pressureFESpace.fe(),
                        M_pressureFESpace.dof(),
                        M_pressureFESpace.dof(),
                        numVelocityComponent,
                        numVelocityComponent,
                        numVelocityComponent * velocityTotalDof,
                        numVelocityComponent * velocityTotalDof );
        chrono.stop();
    }

    M_matrixMassPtr->GlobalAssemble();

    vector_Type ones( *M_solution );
    ones = 1.0;

    Real mean;
    mean = ones* ( M_matrixMassPtr * x );
    x += ( -mean );

    ASSERT( abs( ones* ( M_matrixMassPtr * x ) ) < 1e-9 , "after removeMean the mean pressure should be zero!");

    return mean;

} // removeMean()

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::applyBoundaryConditions( matrix_Type&       matrix,
                                                      vector_Type&       rightHandSide,
                                                      bcHandler_Type& bcHandler )
{
    // M_rightHandSideFull = M_rightHandSideNoBC;

    // BC manage for the velocity
    if ( !bcHandler.bcUpdateDone() || M_recomputeMatrix )
    {
        bcHandler.bcUpdate( *M_velocityFESpace.mesh(),
                            M_velocityFESpace.feBd(),
                            M_velocityFESpace.dof() );
    }

    // ignoring non-local entries, Otherwise they are summed up lately
    //vector_Type rightHandSideFull( rightHandSide, Repeated, Zero );
    // ignoring non-local entries, Otherwise they are summed up lately
    vector_Type rightHandSideFull( rightHandSide, Unique );

    bcManage( matrix, rightHandSideFull,
              *M_velocityFESpace.mesh(),
              M_velocityFESpace.dof(),
              bcHandler,
              M_velocityFESpace.feBd(),
              1.,
              M_dataNavierStokes->dataTime()->time() );

    rightHandSide = rightHandSideFull;

    if ( bcHandler.hasOnlyEssential() && M_diagonalize )
    {
        matrix.diagonalize( nDimensions*dimVelocity(),
                            M_diagonalize,
                            rightHandSide,
                            0. );
    }

} // applyBoundaryCondition

// ===================================================
// Set Methods
// ===================================================

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::setupPostProc( const entityFlag_Type& flag, const mesh_Type meshPart )
{
    M_postProcessing.reset( new PostProc<mesh_Type>( M_velocityFESpace.mesh(),
                                                     &M_velocityFESpace.feBd(),
                                                     &M_velocityFESpace.dof(),
                                                     &M_pressureFESpace.feBd(),
                                                     &M_pressureFESpace.dof(),
                                                     M_localMap ) );
}

template<typename MeshType, typename SolverType>
void
Oseen<MeshType, SolverType>::setTolMaxIteration( const Real& tolerance, const Int& maxIteration )
{
    M_linearSolver.setTolMaxIteration( tolerance, maxIteration );
}


} // namespace LifeV

#endif // OSEEN_H
