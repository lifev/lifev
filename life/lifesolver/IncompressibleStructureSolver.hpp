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
    @brief File containing the implementation of the file IncompressibleStructureData.hpp

    @author Simone Rossi <simone.rossi@epfl.ch>
    @contributor
    @maintainer Simone Rossi <simone.rossi@epfl.ch>

    @date 04-2011

    @brief This class contains an IncompressibleStructure equation solver.

    This file contains an IncompressibleStructure equation solver class.
    The resulting linear systems are solved by GMRES on the full
    matrix ( u and p coupled ).

 */


#ifndef INCOMPRESSIBLESTRUCTURESOLVER_H
#define INCOMPRESSIBLESTRUCTURESOLVER_H 1

#include <life/lifealg/SolverAztecOO.hpp>
#include <life/lifealg/Preconditioner.hpp>
#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerAztecOO.hpp>
#include <life/lifearray/MapEpetra.hpp>

#include <life/lifearray/MatrixElemental.hpp>
#include <life/lifearray/VectorElemental.hpp>
#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifearray/VectorEpetra.hpp>

#include <life/lifecore/LifeChrono.hpp>

#include <life/lifefem/Assembly.hpp>
#include <life/lifefem/BCManage.hpp>
#include <life/lifefem/AssemblyElemental.hpp>
#include <life/lifefem/SobolevNorms.hpp>
#include <life/lifefem/GeometricMap.hpp>
#include <life/lifefem/PostProcessingBoundary.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifesolver/StabilizationIP.hpp>
#include <life/lifesolver/IncompressibleStructureData.hpp>

#include <boost/shared_ptr.hpp>

#include <list>

namespace LifeV
{
//! @class IncompressibleStructure
/*!

 */

template< typename MeshType, typename SolverType = LifeV::SolverAztecOO >
class IncompressibleStructureSolver
{

public:

    //! @name Public Types
    //@{

    typedef MeshType                                    mesh_Type;
    typedef SolverType                                  linearSolver_Type;
    typedef IncompressibleStructureData                 data_Type;

    typedef boost::function<Real ( const Real& t, const Real& x, const Real& y,
                                   const Real& z, const ID& i )> function_Type;

    typedef boost::function<Real ( const Real& t, const Real& x, const Real& y,
                                   const Real& z, const ID& i )> source_Type;

    typedef BCHandler                                   bcHandler_Type;
    typedef boost::shared_ptr<bcHandler_Type>           bcHandlerPtr_Type;

    typedef typename linearSolver_Type::matrix_type     matrix_Type;
    typedef boost::shared_ptr<matrix_Type>              matrixPtr_Type;
    typedef typename linearSolver_Type::vector_type     vector_Type;
    typedef boost::shared_ptr<vector_Type>              vectorPtr_Type;

    typedef typename linearSolver_Type::prec_raw_type   preconditioner_Type;
    typedef typename linearSolver_Type::prec_type       preconditionerPtr_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    IncompressibleStructureSolver();

    //! Constructor
    /*!
        @param dataType IncompressibleStructureData class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param communicator MPI communicator
        @param lagrangeMultiplier Lagrange multiplier
     */

    IncompressibleStructureSolver( boost::shared_ptr<data_Type>    dataType,
                                   FESpace<mesh_Type, MapEpetra>&  displacementFESpace,
                                   FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
                                   boost::shared_ptr<Epetra_Comm>& communicator,
                                   const Int                       lagrangeMultiplier = 0 );


    //! virtual destructor
    virtual ~IncompressibleStructureSolver();

    //@}

    //! @name Methods
    //@{

    //! Set up data from GetPot
    /*!
        @param dataFile GetPot object
     */
    virtual void setUp( const GetPot& dataFile );



    //! Build linear system.
    virtual void buildSystem();

    //! Update system
    /*!
        @param betaVector
        @param sourceVector
     */
    virtual void updateSystem( const vector_Type& sourceVector );

    //! Update system
    /*!
        @param sourceVector
        @param matrix
        @param un
     */
    virtual void updateSystem( const vector_Type& sourceVector,
                               matrixPtr_Type     matrix);



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

    //! Update and return the coefficient matrix
    /*!
        @param matrixFull The coefficient matrix
     */
    void getMatrix( matrix_Type& matrixFull );









    //! Display general information about the content of the class
    /*!
        @param output specify the output format (std::cout by default)
     */
    void showMe( std::ostream& output = std::cout ) const;

    //@}

    //! @name Set Methods
    //@{


    //! set the source term functor
    /*!
        @param source
     */
//    void setSourceTerm( source_Type source )
//    {
//        M_source = source;
//    }


    //@}

    //! @name Get Methods
    //@{

    //! Return the viscosity of the fluid
    /*!
        @return Viscosity of the fluid
     */
    const Real& viscosity() const
    {
        return M_incompressibleStructureData->viscosity();
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
    FESpace<mesh_Type, MapEpetra>& displacementFESpace()
    {
        return M_displacementFESpace;
    }

    const FESpace<mesh_Type, MapEpetra>& displacementFESpace() const
    {
        return M_displacementFESpace;
    }

    //! Return pressure FE space
    /*!
        @return pressure FE space
     */
    FESpace<mesh_Type, MapEpetra>& pressureFESpace()
    {
        return M_pressureFESpace;
    }

    const FESpace<mesh_Type, MapEpetra>& pressureFESpace() const
    {
        return M_pressureFESpace;
    }

    //! Get the source term
    /*!
        @return Source term
     */
//    const source_Type& sourceTerm() const
//    {
//        return M_source;
//    }



    //! Return MapEpetra.
    /*!
        @return MapEpetra
     */
    const MapEpetra& getMap() const
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


    //! Return
    /*!
        @return recomputeMatrix
     */
//    const bool& recomputeMatrix() const
//    {
//        return M_recomputeMatrix;
//    }

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


    //@}


    //! Return a shared pointer to the preconditioner (of type derived from EpetraPreconditioner)
//    preconditionerPtr_Type& preconditioner(){return M_linearSolver.preconditioner();}

protected:

    //! @name Constructor
    //@{

    //! Empty copy constructor
    IncompressibleStructureSolver( const IncompressibleStructureSolver& incompressibleStructure);

    //@}

    //! @name Private Methods
    //@{


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
//    void echo( std::string message );

    //! Return the dim of velocity FE space
    const UInt& dimDisplacement() const
    {
        return M_displacementFESpace.dim();
    }

    //! Return the dim of pressure FE space
    const UInt& dimPressure() const
    {
        return M_pressureFESpace.dim();
    }

    //@}

    //private members

    //! data for Navier-Stokes solvers
    boost::shared_ptr<data_Type>   M_incompressibleStructureData;

    // FE spaces
    FESpace<mesh_Type, MapEpetra>& M_displacementFESpace;
    FESpace<mesh_Type, MapEpetra>& M_pressureFESpace;

    //! MPI communicator
    Displayer                      M_Displayer;

    MapEpetra                      M_localMap;

    //! mass matrix
    matrixPtr_Type                 M_matrixMass;

    //! mass matrix
    matrixPtr_Type                 M_matrixMassPtr;

    //! Stokes matrix: nu*stiff
    matrixPtr_Type                 M_matrixStokes;

    //! matrix to be solved

    //! matrix without boundary conditions
    matrixPtr_Type                 M_matrixNoBC;

    //! stabilization matrix
//    matrixPtr_Type                 M_matrixStabilization;

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

//    bool                           M_steady;

    //! Postprocessing class
//    boost::shared_ptr<PostProcessingBoundary<mesh_Type> > M_postProcessing;


//    bool                           M_stiffStrain;

    //
//    Real                           M_diagonalize;

//    UInt                           M_count;

    bool                           M_recomputeMatrix;

//    bool                           M_isDiagonalBlockPreconditioner;

    //! Elementary matrices and vectors
    MatrixElemental                M_elementMatrixStiff;      // velocity Stokes
    MatrixElemental                M_elementMatrixMass;       // velocity mass
    MatrixElemental                M_elementMatrixPreconditioner;          // (p,q) bloc for preconditioners
    MatrixElemental                M_elementMatrixDivergence;
    MatrixElemental                M_elementMatrixGradient;
    VectorElemental                M_elementRightHandSide;           // Elementary right hand side
    matrixPtr_Type                 M_blockPreconditioner;
//    VectorElemental                M_wLoc;
//    VectorElemental                M_uLoc;
//    boost::shared_ptr<vector_Type> M_un;

}; // class IncompressibleStructureSolver



// ===================================================
// Constructors & Destructor
// ===================================================

template<typename MeshType, typename SolverType>
IncompressibleStructureSolver<MeshType, SolverType>::
IncompressibleStructureSolver( boost::shared_ptr<data_Type>    dataType,
                               FESpace<mesh_Type, MapEpetra>&  displacementFESpace,
                               FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
                               boost::shared_ptr<Epetra_Comm>& communicator,
                               const Int                       lagrangeMultiplier ):
        M_incompressibleStructureData       ( dataType ),
        M_displacementFESpace        ( displacementFESpace ),
        M_pressureFESpace        ( pressureFESpace ),
        M_Displayer              ( communicator ),
        M_localMap               ( M_displacementFESpace.map() + M_pressureFESpace.map() + lagrangeMultiplier),
        M_matrixMass             ( ),
        M_matrixMassPtr          ( ),
        M_matrixStokes           ( ),
        M_matrixNoBC             ( ),
//        M_matrixStabilization    ( ),
        M_rightHandSideNoBC      ( M_localMap ),
        M_rightHandSideFull      ( M_localMap ),
        M_solution               ( new vector_Type( M_localMap ) ),
        M_residual               ( M_localMap ),
        M_linearSolver           ( communicator ),
//        M_steady                 ( ),
//        M_postProcessing         ( new PostProcessingBoundary<mesh_Type>( M_displacementFESpace.mesh(),
//                                                            &M_displacementFESpace.feBd(),
//                                                            &M_displacementFESpace.dof(),
//                                                            &M_pressureFESpace.feBd(),
//                                                            &M_pressureFESpace.dof(),
//                                                            M_localMap ) ),
//        M_stiffStrain            ( false ),
//        M_diagonalize            ( false ),
//        M_count                  ( 0 ),
        M_recomputeMatrix        ( false ),
        M_elementMatrixStiff     ( M_displacementFESpace.fe().nbFEDof(), nDimensions, nDimensions ),
        M_elementMatrixMass      ( M_displacementFESpace.fe().nbFEDof(), nDimensions, nDimensions ),
        M_elementMatrixPreconditioner ( M_pressureFESpace.fe().nbFEDof(), 1, 1 ),
        M_elementMatrixDivergence ( M_pressureFESpace.fe().nbFEDof(), 1, 0,
                                    M_displacementFESpace.fe().nbFEDof(), 0, nDimensions ),
        M_elementMatrixGradient  ( M_displacementFESpace.fe().nbFEDof(), nDimensions, 0,
                                   M_pressureFESpace.fe().nbFEDof(), 0, 1 ),
        M_elementRightHandSide   ( M_displacementFESpace.fe().nbFEDof(), nDimensions ),
        M_blockPreconditioner    ( )
//        M_wLoc                   ( M_displacementFESpace.fe().nbFEDof(), nDimensions ),
//        M_uLoc                   ( M_displacementFESpace.fe().nbFEDof(), nDimensions )
//        M_un                     ( new vector_Type(M_localMap) )
{

    std::cout << "costruttore 1 " << std::endl << std::endl;
}


template<typename MeshType, typename SolverType>
IncompressibleStructureSolver<MeshType, SolverType>::
~IncompressibleStructureSolver()
{
    std::cout << "distruttore 1 " << std::endl << std::endl;

}


// ===================================================
// Methods
// ===================================================

template<typename MeshType, typename SolverType>
void
IncompressibleStructureSolver<MeshType, SolverType>::setUp( const GetPot& dataFile )
{

    M_linearSolver.setupPreconditioner( dataFile, "fluid/prec" );
    M_linearSolver.setDataFromGetPot( dataFile, "fluid/solver" );

//    M_steady        = dataFile( "fluid/miscellaneous/steady", 0 );

    // Energetic stabilization term
    // Enable grad( u )^T in stress tensor
//    M_stiffStrain = dataFile( "fluid/space_discretization/stiff_strain",false);
//    M_diagonalize = dataFile( "fluid/space_discretization/diagonalize", 1. );
//    M_isDiagonalBlockPreconditioner = dataFile( "fluid/diagonalBlockPrec", false );

}



template<typename MeshType, typename SolverType>
void
IncompressibleStructureSolver<MeshType, SolverType>::buildSystem()
{
    M_matrixMass.reset  ( new matrix_Type( M_localMap ) );
    M_matrixStokes.reset( new matrix_Type( M_localMap ) );

    M_Displayer.leaderPrint( "  F-  Computing constant matrices ...          " );

    LifeChrono chrono;

    LifeChrono chronoDer;
    LifeChrono chronoStiff;
    LifeChrono chronoMass;
    LifeChrono chronoGrad;

    LifeChrono chronoStiffAssemble;
    LifeChrono chronoMassAssemble;
    LifeChrono chronoGradAssemble;
    LifeChrono chronoDivAssemble;
    LifeChrono chronoStab;
    LifeChrono chronoZero;

    // Number of velocity components
    UInt numVelocityComponent = nDimensions;

    // Elementary computation and matrix assembling
    // Loop on elements

    UInt velocityTotalDof   = M_displacementFESpace.dof().numTotalDof();
//    UInt pressureTotalDof = M_pressureFESpace.dof().numTotalDof();

//    if ( M_isDiagonalBlockPreconditioner == true )
//    {
//        M_blockPreconditioner.reset( new matrix_Type( M_localMap ) );
//    }
    chrono.start();

    for ( UInt iVolume = 0; iVolume < M_displacementFESpace.mesh()->numVolumes(); iVolume++ )
    {
        chronoDer.start();
        // just to provide the id number in the assem_mat_mixed
        M_pressureFESpace.fe().update( M_displacementFESpace.mesh()->volumeList( iVolume ) );
        // just to provide the id number in the assem_mat_mixed
        M_displacementFESpace.fe().updateFirstDeriv( M_displacementFESpace.mesh()->volumeList( iVolume ) );

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
//        if ( M_stiffStrain )
//            stiff_strain( 2.0*M_incompressibleStructureData->viscosity(),
//                          M_elementMatrixStiff,
//                          M_displacementFESpace.fe() );
//        else
            stiff( M_incompressibleStructureData->viscosity(),
                   M_elementMatrixStiff,
                   M_displacementFESpace.fe(), 0, 0, nDimensions );
//        chronoStiff.stop();



        for ( UInt iComponent = 0; iComponent < numVelocityComponent; iComponent++ )
        {
            // stiffness matrix
            chronoStiffAssemble.start();
//            if ( M_isDiagonalBlockPreconditioner == true )
//            {
//                assembleMatrix( *M_blockPreconditioner,
//                                M_elementMatrixStiff,
//                                M_displacementFESpace.fe(),
//                                M_displacementFESpace.fe(),
//                                M_displacementFESpace.dof(),
//                                M_displacementFESpace.dof(),
//                                iComponent, iComponent,
//                                iComponent * velocityTotalDof, iComponent * velocityTotalDof);
//            }
//            else
//            {
//                if ( M_stiffStrain ) // sigma = 0.5 * mu (grad( u ) + grad ( u )^T)
//                {
//                    for ( UInt jComp = 0; jComp < numVelocityComponent; jComp++ )
//                    {
//                        assembleMatrix( *M_matrixStokes,
//                                        M_elementMatrixStiff,
//                                        M_displacementFESpace.fe(),
//                                        M_displacementFESpace.fe(),
//                                        M_displacementFESpace.dof(),
//                                        M_displacementFESpace.dof(),
//                                        iComponent, jComp,
//                                        iComponent * velocityTotalDof, jComp * velocityTotalDof);
//
//                    }
//                }
//                else // sigma = mu grad( u )
//                {
                    assembleMatrix( *M_matrixStokes,
                                    M_elementMatrixStiff,
                                    M_displacementFESpace.fe(),
                                    M_displacementFESpace.fe(),
                                    M_displacementFESpace.dof(),
                                    M_displacementFESpace.dof(),
                                    iComponent, iComponent,
                                    iComponent * velocityTotalDof, iComponent * velocityTotalDof);
//                }
//            }
            chronoStiffAssemble.stop();

            // mass matrix
//            if ( !M_steady )
//            {
//                chronoMassAssemble.start();
//                assembleMatrix( *M_matrixMass,
//                                M_elementMatrixMass,
//                                M_displacementFESpace.fe(),
//                                M_displacementFESpace.fe(),
//                                M_displacementFESpace.dof(),
//                                M_displacementFESpace.dof(),
//                                iComponent, iComponent,
//                                iComponent * velocityTotalDof, iComponent * velocityTotalDof);
//                chronoMassAssemble.stop();
//            }

            // divergence
            chronoGrad.start();
            grad( iComponent, 1.0,
                  M_elementMatrixGradient,
                  M_displacementFESpace.fe(),
                  M_pressureFESpace.fe(),
                  iComponent, 0 );
            chronoGrad.stop();

            chronoGradAssemble.start();
            assembleMatrix( *M_matrixStokes,
                            M_elementMatrixGradient,
                            M_displacementFESpace.fe(),
                            M_pressureFESpace.fe(),
                            M_displacementFESpace.dof(),
                            M_pressureFESpace.dof(),
                            iComponent, 0,
                            iComponent * velocityTotalDof, numVelocityComponent * velocityTotalDof );
            chronoGradAssemble.stop();

            chronoDivAssemble.start();
            assembleTransposeMatrix( *M_matrixStokes,
                                     -1.,
                                     M_elementMatrixGradient,
                                     M_pressureFESpace.fe(),
                                     M_displacementFESpace.fe(),
                                     M_pressureFESpace.dof(),
                                     M_displacementFESpace.dof(),
                                     0 , iComponent,
                                     numVelocityComponent * velocityTotalDof, iComponent * velocityTotalDof );
            chronoDivAssemble.stop();
        }
    }


//    if ( M_isDiagonalBlockPreconditioner == true )
//    {
//        M_blockPreconditioner->globalAssemble();
//        *M_matrixStokes += *M_blockPreconditioner;
//    }
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
IncompressibleStructureSolver<MeshType, SolverType>::
updateSystem( const vector_Type& sourceVector )
{
    if ( M_matrixNoBC.get() )
        M_matrixNoBC.reset( new matrix_Type( M_localMap, M_matrixNoBC->meanNumEntries() ) );
    else
        M_matrixNoBC.reset( new matrix_Type( M_localMap ) );

    updateSystem( sourceVector, M_matrixNoBC);

}

template<typename MeshType, typename SolverType>
void
IncompressibleStructureSolver<MeshType, SolverType>::
updateSystem( const vector_Type& sourceVector,
              matrixPtr_Type     matrixNoBC)
{
    LifeChrono chrono;

    // clearing pressure mass matrix in case we need it in removeMean;
    M_matrixMassPtr.reset( );


    M_Displayer.leaderPrint( "  F-  Updating mass term on right hand side... " );

    chrono.start();

    std::string   sv="sv";
    sourceVector.spy(sv);
    updateRightHandSide( sourceVector );


    chrono.stop();

    M_Displayer.leaderPrintMax( "done in ", chrono.diff() );



    if ( M_recomputeMatrix )
        buildSystem();

    M_Displayer.leaderPrint( "  F-  Copying the matrices ...                 " );

    chrono.start();

//    if ( M_isDiagonalBlockPreconditioner == true )
//    {
//        matrixPtr_Type tempMatrix( M_blockPreconditioner );
//        M_blockPreconditioner.reset( new matrix_Type( M_localMap,
//                                                      M_blockPreconditioner->meanNumEntries() ) );
//        *M_blockPreconditioner += *tempMatrix;
//    }


    chrono.stop();
    M_Displayer.leaderPrintMax( "done in " , chrono.diff() );


//    if ( alpha != 0. )
//    {
//        *matrixNoBC += (*M_matrixMass) * alpha;
//        if ( M_isDiagonalBlockPreconditioner == true )
//        {
//            matrixNoBC->globalAssemble();
//            *M_blockPreconditioner += *matrixNoBC;
//            matrix_Type tempMatrix( *matrixNoBC );
//            matrixNoBC.reset( new matrix_Type( M_localMap, tempMatrix.meanNumEntries() ) );
//            *matrixNoBC += tempMatrix;
//            M_blockPreconditioner->globalAssemble();
//        }
//    }
    *matrixNoBC += *M_matrixStokes;

}




template<typename MeshType, typename SolverType>
void
IncompressibleStructureSolver<MeshType, SolverType>::iterate( bcHandler_Type& bcHandler )
{

    LifeChrono chrono;

    // matrix and vector assembling communication
    M_Displayer.leaderPrint( "  F-  Updating the boundary conditions ...     " );

    chrono.start();

    M_matrixNoBC->globalAssemble();

    matrixPtr_Type matrixFull( new matrix_Type( M_localMap, M_matrixNoBC->meanNumEntries() ) );

    //updateStabilization( *matrixFull );
    getMatrix( *matrixFull );

    vector_Type rightHandSideFull ( M_rightHandSideNoBC );

    chrono.stop();

    M_Displayer.leaderPrintMax( "done in ", chrono.diff() );

    // boundary conditions update
    M_Displayer.leaderPrint("  F-  Applying boundary conditions ...         ");

    chrono.start();
    applyBoundaryConditions( *matrixFull, rightHandSideFull, bcHandler );

    matrixFull->globalAssemble();
    chrono.stop();

    M_Displayer.leaderPrintMax( "done in " , chrono.diff() );

    // solving the system
    M_linearSolver.setMatrix( *matrixFull );

    M_linearSolver.solveSystem( rightHandSideFull, *M_solution, matrixFull );


    M_residual  = M_rightHandSideNoBC;
    M_residual -= (*M_matrixNoBC) * (*M_solution);

    //M_residual.spy("residual");
} // iterate()



template<typename MeshType, typename SolverType>
void
IncompressibleStructureSolver<MeshType, SolverType>::getMatrix( matrix_Type& matrixFull )
{
    M_matrixNoBC->globalAssemble();
    matrixFull += *M_matrixNoBC;
}



template<typename MeshType, typename SolverType>
void
IncompressibleStructureSolver<MeshType, SolverType>::applyBoundaryConditions( matrix_Type&       matrix,
                                                      vector_Type&       rightHandSide,
                                                      bcHandler_Type& bcHandler )
{

    // BC manage for the velocity
    if ( !bcHandler.bcUpdateDone() || M_recomputeMatrix )
    {
        bcHandler.bcUpdate( *M_displacementFESpace.mesh(),
                            M_displacementFESpace.feBd(),
                            M_displacementFESpace.dof() );
    }

    // ignoring non-local entries, Otherwise they are summed up lately
    //vector_Type rightHandSideFull( rightHandSide, Repeated, Zero );
    // ignoring non-local entries, Otherwise they are summed up lately
    vector_Type rightHandSideFull( rightHandSide, Unique );



    bcManage( matrix, rightHandSideFull,
              *M_displacementFESpace.mesh(),
              M_displacementFESpace.dof(),
              bcHandler,
              M_displacementFESpace.feBd(),
              1.);

    rightHandSide = rightHandSideFull;

//    if ( bcHandler.hasOnlyEssential() && M_diagonalize )
//    {
//        matrix.diagonalize( nDimensions*dimDisplacement(),
//                            M_diagonalize,
//                            rightHandSide,
//                            0. );
//    }

} // applyBoundaryCondition

// ===================================================
// Set Methods
// ===================================================




} // namespace LifeV

#endif // INCOMPRESSIBLESTRUCTURESOLVER_H
