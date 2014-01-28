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
    @brief File containing the implementation of the file OseenSolver.hpp

    @author Davide Forti <davide.forti@epfl.ch>

    @date 28-01-2014

 */

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

template<typename MeshType, typename SolverType>
OseenSolver<MeshType, SolverType>::
OseenSolver ( boost::shared_ptr<data_Type>    dataType,
              FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
              FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
              boost::shared_ptr<Epetra_Comm>& communicator,
              const Int                       lagrangeMultiplier ) :
    M_oseenData       ( dataType ),
    M_velocityFESpace        ( velocityFESpace ),
    M_pressureFESpace        ( pressureFESpace ),
    M_Displayer              ( communicator ),
    M_localMap               ( M_velocityFESpace.map() + M_pressureFESpace.map() + lagrangeMultiplier),
    M_velocityMatrixMass     ( ),
    M_pressureMatrixMass     ( ),
    M_matrixStokes           ( ),
    M_matrixNoBC             ( ),
    M_matrixStabilization    ( ),
    M_rightHandSideNoBC      ( ),
    M_solution               ( new vector_Type ( M_localMap ) ),
    M_residual               ( new vector_Type (M_localMap ) ),
    M_linearSolver           ( new linearSolver_Type (communicator) ),
    M_steady                 ( ),
    M_postProcessing         ( new PostProcessingBoundary<mesh_Type> ( M_velocityFESpace.mesh(),
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
    //                                   M_oseenData->viscosity() ),
    M_betaFunction           ( 0 ),
    M_divBetaUv              ( false ),
    M_stiffStrain            ( false ),
    M_diagonalize            ( false ),
    M_count                  ( 0 ),
    M_recomputeMatrix        ( false ),
    M_elementMatrixStiff     ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim(), velocityFESpace.fieldDim() ),
    M_elementMatrixMass      ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim(), velocityFESpace.fieldDim() ),
    M_elementMatrixPreconditioner ( M_pressureFESpace.fe().nbFEDof(), 1, 1 ),
    M_elementMatrixDivergence ( M_pressureFESpace.fe().nbFEDof(), 1, 0,
                                M_velocityFESpace.fe().nbFEDof(), 0, velocityFESpace.fieldDim() ),
    M_elementMatrixGradient  ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim(), 0,
                               M_pressureFESpace.fe().nbFEDof(), 0, 1 ),
    M_elementRightHandSide   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
    M_blockPreconditioner    ( ),
    M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
    M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
    M_un                     ( new vector_Type (M_localMap) )
{
    // if(M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ))
    {
        M_ipStabilization.setFeSpaceVelocity (M_velocityFESpace);
        M_ipStabilization.setViscosity (M_oseenData->viscosity() );
    }
}

template<typename MeshType, typename SolverType>
OseenSolver<MeshType, SolverType>::
OseenSolver ( boost::shared_ptr<data_Type>    dataType,
              FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
              FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
              boost::shared_ptr<Epetra_Comm>& communicator,
              MapEpetra                       monolithicMap,
              UInt                            /*offset*/ ) :
    M_oseenData              ( dataType ),
    M_velocityFESpace        ( velocityFESpace ),
    M_pressureFESpace        ( pressureFESpace ),
    M_Displayer              ( communicator ),
    M_localMap               ( monolithicMap ),
    M_velocityMatrixMass     ( ),
    M_matrixStokes           ( ),
    M_matrixNoBC             ( ),
    M_matrixStabilization    ( ),
    M_rightHandSideNoBC      ( ),
    M_solution               ( ),
    M_residual               (  ),
    M_linearSolver           ( ),
    M_postProcessing         ( new PostProcessingBoundary<mesh_Type> (M_velocityFESpace.mesh(),
                                                                      &M_velocityFESpace.feBd(),
                                                                      &M_velocityFESpace.dof(),
                                                                      &M_pressureFESpace.feBd(),
                                                                      &M_pressureFESpace.dof(),
                                                                      M_localMap ) ),
    M_stabilization          ( false ),
    M_reuseStabilization     ( false ),
    M_resetStabilization     ( false ),
    M_iterReuseStabilization ( -1 ),
    M_betaFunction           ( 0 ),
    M_divBetaUv              ( false ),
    M_stiffStrain            ( false ),
    M_diagonalize            ( false ),
    M_count                  ( 0 ),
    M_recomputeMatrix        ( false ),
    M_elementMatrixStiff     ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), M_velocityFESpace.fieldDim() ),
    M_elementMatrixMass      ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), M_velocityFESpace.fieldDim() ),
    M_elementMatrixPreconditioner                 ( M_pressureFESpace.fe().nbFEDof(), 1, 1 ),
    M_elementMatrixDivergence ( M_pressureFESpace.fe().nbFEDof(), 1, 0,
                                M_velocityFESpace.fe().nbFEDof(), 0, M_velocityFESpace.fieldDim() ),
    M_elementMatrixGradient  ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), 0,
                               M_pressureFESpace.fe().nbFEDof(), 0, 1 ),
    M_elementRightHandSide   ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim() ),
    M_blockPreconditioner    ( ),
    M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim() ),
    M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim() ),
    M_un                     ( /*new vector_Type(M_localMap)*/ )
{
    // if(M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ))
    {
        M_ipStabilization.setFeSpaceVelocity (M_velocityFESpace);
        M_ipStabilization.setViscosity (M_oseenData->viscosity() );
    }
}

template<typename MeshType, typename SolverType>
OseenSolver<MeshType, SolverType>::
OseenSolver ( boost::shared_ptr<data_Type>    dataType,
              FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
              FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
              const std::vector<Int>&         lagrangeMultipliers,
              boost::shared_ptr<Epetra_Comm>& communicator ) :
    M_oseenData       ( dataType ),
    M_velocityFESpace        ( velocityFESpace ),
    M_pressureFESpace        ( pressureFESpace ),
    M_Displayer              ( communicator ),
    M_localMap               ( M_velocityFESpace.map() + M_pressureFESpace.map() + lagrangeMultipliers ),
    M_velocityMatrixMass     ( ),
    M_matrixStokes           ( ),
    M_matrixNoBC             ( ),
    M_matrixStabilization    ( ),
    M_rightHandSideNoBC      ( ),
    M_solution               ( new vector_Type ( M_localMap ) ),
    M_residual               (  ),
    M_linearSolver           ( new linearSolver_Type (communicator) ),
    M_postProcessing         ( new PostProcessingBoundary<mesh_Type> (M_velocityFESpace.mesh(),
                                                                      &M_velocityFESpace.feBd(),
                                                                      &M_velocityFESpace.dof(),
                                                                      &M_pressureFESpace.feBd(),
                                                                      &M_pressureFESpace.dof(),
                                                                      M_localMap ) ),
    M_stabilization          ( false ),
    M_reuseStabilization     ( false ),
    M_resetStabilization     ( false ),
    M_iterReuseStabilization ( -1 ),
    M_betaFunction           ( 0 ),
    M_divBetaUv              ( false ),
    M_stiffStrain            ( false ),
    M_diagonalize            ( false ),
    M_count                  ( 0 ),
    M_recomputeMatrix        ( false ),
    M_elementMatrixStiff     ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), M_velocityFESpace.fieldDim() ),
    M_elementMatrixMass      ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), M_velocityFESpace.fieldDim() ),
    M_elementMatrixPreconditioner ( M_pressureFESpace.fe().nbFEDof(), 1, 1 ),
    M_elementMatrixDivergence ( M_pressureFESpace.fe().nbFEDof(), 1, 0,
                                M_velocityFESpace.fe().nbFEDof(), 0, M_velocityFESpace.fieldDim() ),
    M_elementMatrixGradient  ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), 0,
                               M_pressureFESpace.fe().nbFEDof(), 0, 1 ),
    M_elementRightHandSide   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
    M_blockPreconditioner    ( ),
    M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
    M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
    M_un                     ( new vector_Type (M_localMap) )
{
    // if(M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ))
    {
        M_ipStabilization.setFeSpaceVelocity (M_velocityFESpace);
        M_ipStabilization.setViscosity (M_oseenData->viscosity() );
    }
}

template<typename MeshType, typename SolverType>
OseenSolver<MeshType, SolverType>::
~OseenSolver()
{

}


// ===================================================
// Methods
// ===================================================

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::setUp ( const GetPot& dataFile )
{
    if (M_linearSolver.get() )
    {
        M_linearSolver->setupPreconditioner ( dataFile, "fluid/prec" );
        M_linearSolver->setDataFromGetPot ( dataFile, "fluid/solver" );
    }

    M_stabilization = dataFile ( "fluid/ipstab/use",  (&M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ) );
    M_steady        = dataFile ( "fluid/miscellaneous/steady", 0 );
    if (M_stabilization)
    {
        M_gammaBeta     = dataFile ( "fluid/ipstab/gammaBeta",  0. );
        M_gammaDiv      = dataFile ( "fluid/ipstab/gammaDiv",   0. );
        M_gammaPress    = dataFile ( "fluid/ipstab/gammaPress", 0. );
        M_reuseStabilization     = dataFile ( "fluid/ipstab/reuse", false );
        if (M_linearSolver.get() )
            M_iterReuseStabilization = dataFile ( "fluid/ipstab/max_iter_reuse",
                                                  static_cast<Int> ( M_linearSolver->maxNumIterations() * 8. / 10. ) );
    }
    // Energetic stabilization term
    M_divBetaUv   = dataFile ( "fluid/space_discretization/div_beta_u_v", false);
    // Enable grad( u )^T in stress tensor
    M_stiffStrain = dataFile ( "fluid/space_discretization/stiff_strain", false);
    M_diagonalize = dataFile ( "fluid/space_discretization/diagonalize", 1. );
    M_isDiagonalBlockPreconditioner = dataFile ( "fluid/diagonalBlockPrec", false );

    //    M_linearSolver.setAztecooPreconditioner( dataFile, "fluid/solver" );

    M_ipStabilization.setGammaBeta ( M_gammaBeta );
    M_ipStabilization.setGammaDiv  ( M_gammaDiv );
    M_ipStabilization.setGammaPress ( M_gammaPress );
}


template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::
initialize ( const function_Type& velocityFunction, const function_Type& pressureFunction )
{
    vector_Type velocityInitialGuess ( M_velocityFESpace.map() );
    M_velocityFESpace.interpolate ( velocityFunction,
                                    velocityInitialGuess,
                                    M_oseenData->dataTime()->time() );

    vector_Type pressureInitialGuess ( M_pressureFESpace.map() );
    M_pressureFESpace.interpolate ( pressureFunction,
                                    pressureInitialGuess,
                                    M_oseenData->dataTime()->time() );

    initialize ( velocityInitialGuess, pressureInitialGuess );
}


template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::
initialize ( const vector_Type& velocityInitialGuess, const vector_Type& pressureInitialGuess )
{

    *M_solution = velocityInitialGuess;
    *M_un = velocityInitialGuess;
    M_solution->add ( pressureInitialGuess, M_velocityFESpace.fieldDim() * M_velocityFESpace.dof().numTotalDof() );
    M_un->add ( pressureInitialGuess, M_velocityFESpace.fieldDim() * M_velocityFESpace.dof().numTotalDof() );

}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::
initialize ( const vector_Type& velocityAndPressure )
{

    *M_un = velocityAndPressure;
    if ( M_solution.get() )
    {
        *M_solution = velocityAndPressure;
    }

}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::buildSystem()
{
    M_velocityMatrixMass.reset  ( new matrix_Type ( M_localMap ) );
    M_matrixStokes.reset ( new matrix_Type ( M_localMap ) );

    M_Displayer.leaderPrint ( "  F-  Computing constant matrices ...          " );

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
    UInt numVelocityComponent = M_velocityFESpace.fieldDim();

    // Elementary computation and matrix assembling
    // Loop on elements

    UInt velocityTotalDof   = M_velocityFESpace.dof().numTotalDof();
    //    UInt pressureTotalDof = M_pressureFESpace.dof().numTotalDof();

    if ( M_isDiagonalBlockPreconditioner == true )
    {
        M_blockPreconditioner.reset ( new matrix_Type ( M_localMap ) );
    }
    chrono.start();

    for ( UInt iElement = 0; iElement < M_velocityFESpace.mesh()->numElements(); iElement++ )
    {
        chronoDer.start();
        // just to provide the id number in the assem_mat_mixed
        M_pressureFESpace.fe().update ( M_velocityFESpace.mesh()->element ( iElement ) );
        // just to provide the id number in the assem_mat_mixed
        // M_pressureFESpace.fe().updateFirstDeriv( M_velocityFESpace.mesh()->element( iElement ) );
        M_velocityFESpace.fe().updateFirstDeriv ( M_velocityFESpace.mesh()->element ( iElement ) );

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
            stiff_strain ( 2.0 * M_oseenData->viscosity(),
                           M_elementMatrixStiff,
                           M_velocityFESpace.fe() );
        else
            stiff ( M_oseenData->viscosity(),
                    M_elementMatrixStiff,
                    M_velocityFESpace.fe(), 0, 0, M_velocityFESpace.fieldDim() );
        //stiff_div( 0.5*M_velocityFESpace.fe().diameter(), M_elementMatrixStiff, M_velocityFESpace.fe() );
        chronoStiff.stop();

        // mass matrix
        if ( !M_steady )
        {
            chronoMass.start();
            mass ( M_oseenData->density(),
                   M_elementMatrixMass,
                   M_velocityFESpace.fe(), 0, 0, M_velocityFESpace.fieldDim() );
            chronoMass.stop();
        }

        for ( UInt iComponent = 0; iComponent < numVelocityComponent; iComponent++ )
        {
            // stiffness matrix
            chronoStiffAssemble.start();
            if ( M_isDiagonalBlockPreconditioner == true )
            {
                assembleMatrix ( *M_blockPreconditioner,
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
                        assembleMatrix ( *M_matrixStokes,
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
                    assembleMatrix ( *M_matrixStokes,
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
                assembleMatrix ( *M_velocityMatrixMass,
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
            grad ( iComponent, 1.0,
                   M_elementMatrixGradient,
                   M_velocityFESpace.fe(),
                   M_pressureFESpace.fe(),
                   iComponent, 0 );
            chronoGrad.stop();

            chronoGradAssemble.start();
            assembleMatrix ( *M_matrixStokes,
                             M_elementMatrixGradient,
                             M_velocityFESpace.fe(),
                             M_pressureFESpace.fe(),
                             M_velocityFESpace.dof(),
                             M_pressureFESpace.dof(),
                             iComponent, 0,
                             iComponent * velocityTotalDof, numVelocityComponent * velocityTotalDof );
            chronoGradAssemble.stop();

            chronoDivAssemble.start();
            assembleTransposeMatrix ( *M_matrixStokes,
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

    //    for (UInt ii = M_velocityFESpace.fieldDim()*dimVelocity(); ii < M_velocityFESpace.fieldDim()*dimVelocity() + dimPressure(); ++ii)
    //  M_matrixStokes->set_mat_inc( ii ,ii, 0. ); not scalable!!!

    if ( M_isDiagonalBlockPreconditioner == true )
    {
        M_blockPreconditioner->globalAssemble();
        *M_matrixStokes += *M_blockPreconditioner;
    }
    comm()->Barrier();

    chrono.stop();
    M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

    M_Displayer.leaderPrint ( "  F-  Finalizing the matrices ...              " );

    chrono.start();

    M_matrixStokes->globalAssemble();
    M_velocityMatrixMass->globalAssemble();

    chrono.stop();
    M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

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
OseenSolver<MeshType, SolverType>::
updateSystem ( const Real         alpha,
               const vector_Type& betaVector,
               const vector_Type& sourceVector )
{
    if ( M_matrixNoBC.get() )
    {
        M_matrixNoBC.reset ( new matrix_Type ( M_localMap, M_matrixNoBC->meanNumEntries() ) );
    }
    else
    {
        M_matrixNoBC.reset ( new matrix_Type ( M_localMap ) );
    }

    updateSystem ( alpha, betaVector, sourceVector, M_matrixNoBC, *M_un );
    if ( alpha != 0. )
    {
        M_matrixNoBC->globalAssemble();
    }

}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::
updateSystem ( const Real         alpha,
               const vector_Type& betaVector,
               const vector_Type& sourceVector,
               matrixPtr_Type     matrixNoBC,
               const vector_Type&     un )
{
    LifeChrono chrono;

    // clearing pressure mass matrix in case we need it in removeMean;
    M_pressureMatrixMass.reset( );


    M_Displayer.leaderPrint ( "  F-  Updating mass term on right hand side... " );

    chrono.start();

    UInt velocityTotalDof   = M_velocityFESpace.dof().numTotalDof();
    //    UInt pressureTotalDof = M_pressureFESpace.dof().numTotalDof();

    // Right hand side for the velocity at time

    updateRightHandSide ( sourceVector );

    chrono.stop();

    M_Displayer.leaderPrintMax ( "done in ", chrono.diff() );


    //    M_updated = false;

    if ( M_recomputeMatrix )
    {
        buildSystem();
    }

    M_Displayer.leaderPrint ( "  F-  Copying the matrices ...                 " );

    chrono.start();

    if ( M_isDiagonalBlockPreconditioner == true )
    {
        matrixPtr_Type tempMatrix ( M_blockPreconditioner );
        M_blockPreconditioner.reset ( new matrix_Type ( M_localMap,
                                                        M_blockPreconditioner->meanNumEntries() ) );
        *M_blockPreconditioner += *tempMatrix;
    }


    chrono.stop();
    M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );


    UInt numVelocityComponent = M_velocityFESpace.fieldDim();

    //! managing the convective term

    Real normInf;
    betaVector.normInf ( &normInf );

    if ( normInf != 0. )
    {
        M_Displayer.leaderPrint ( "  F-  Sharing convective term ...              " );
        chrono.start();

        // vector with repeated nodes over the processors

        vector_Type betaVectorRepeated ( betaVector, Repeated );
        vector_Type unRepeated ( un, Repeated );

        chrono.stop();

        M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
        M_Displayer.leaderPrint ( "  F-  Updating the convective terms ...        " );
        chrono.start();

        for ( UInt iElement = 0; iElement < M_velocityFESpace.mesh()->numElements(); ++iElement )
        {
            // just to provide the id number in the assem_mat_mixed
            M_pressureFESpace.fe().updateFirstDeriv ( M_velocityFESpace.mesh()->element ( iElement ) );
            //as updateFirstDer
            M_velocityFESpace.fe().updateFirstDeriv ( M_velocityFESpace.mesh()->element ( iElement ) );

            M_elementMatrixStiff.zero();

            UInt elementID = M_velocityFESpace.fe().currentLocalId();
            // Non linear term, Semi-implicit approach
            // M_elementRightHandSide contains the velocity values in the nodes
            for ( UInt iNode = 0 ; iNode < M_velocityFESpace.fe().nbFEDof() ; iNode++ )
            {
                UInt iLocal = M_velocityFESpace.fe().patternFirst ( iNode );
                for ( UInt iComponent = 0; iComponent < numVelocityComponent; ++iComponent )
                {
                    UInt iGlobal = M_velocityFESpace.dof().localToGlobalMap ( elementID, iLocal )
                                   + iComponent * dimVelocity();
                    M_elementRightHandSide.vec() [ iLocal + iComponent * M_velocityFESpace.fe().nbFEDof() ]
                        = betaVectorRepeated[iGlobal];

                    M_uLoc.vec() [ iLocal + iComponent * M_velocityFESpace.fe().nbFEDof() ]
                        = unRepeated (iGlobal);
                    M_wLoc.vec() [ iLocal + iComponent * M_velocityFESpace.fe().nbFEDof() ]
                        = unRepeated (iGlobal) - betaVectorRepeated (iGlobal);
                }
            }

            if (M_oseenData->conservativeFormulation() )
            {
                // ALE term: - rho div(w) u v
                mass_divw ( - M_oseenData->density(),
                            M_wLoc,
                            M_elementMatrixStiff,
                            M_velocityFESpace.fe(), 0, 0, numVelocityComponent );
            }

            // ALE stab implicit: 0.5 rho div u w v
            mass_divw ( 0.5 * M_oseenData->density(),
                        M_uLoc,
                        M_elementMatrixStiff,
                        M_velocityFESpace.fe(), 0, 0, numVelocityComponent );

            // Stabilising term: -rho div(u^n) u v
            if ( M_divBetaUv )
                mass_divw ( -0.5 * M_oseenData->density(),
                            M_uLoc,
                            M_elementMatrixStiff,
                            M_velocityFESpace.fe(), 0, 0, numVelocityComponent );

            // compute local convective terms
            advection ( M_oseenData->density(),
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

                assembleMatrix ( *matrixNoBC,
                                 M_elementMatrixStiff,
                                 M_velocityFESpace.fe(),
                                 M_velocityFESpace.fe(),
                                 M_velocityFESpace.dof(),
                                 M_velocityFESpace.dof(),
                                 iComponent, iComponent,
                                 iComponent * velocityTotalDof, iComponent * velocityTotalDof );
            }
        }

        chrono.stop();
        M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

        if ( M_stabilization &&
                ( M_resetStabilization || !M_reuseStabilization || ( M_matrixStabilization.get() == 0 ) ) )
        {
            M_Displayer.leaderPrint ( "  F-  Updating the stabilization terms ...     " );
            chrono.start();
            M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
            M_ipStabilization.apply ( *M_matrixStabilization, betaVectorRepeated, false );
            M_matrixStabilization->globalAssemble();
            M_resetStabilization = false;
            chrono.stop();
            M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
        }

    }
    else
    {
        if ( M_stabilization )
        {
            M_Displayer.leaderPrint ( "  F-  Updating the stabilization terms ...     " );
            chrono.start();

            if ( M_resetStabilization || !M_reuseStabilization || ( M_matrixStabilization.get() == 0 ) )
            {
                M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
                M_ipStabilization.apply ( *M_matrixStabilization, betaVector, false );
                M_matrixStabilization->globalAssemble();
                M_resetStabilization = false;
                chrono.stop();
                M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
            }
            else
            {
                M_Displayer.leaderPrint ( "reusing\n" );
            }
        }
    }

    if ( alpha != 0. )
    {
        *matrixNoBC += (*M_velocityMatrixMass) * alpha;
        if ( M_isDiagonalBlockPreconditioner == true )
        {
            matrixNoBC->globalAssemble();
            *M_blockPreconditioner += *matrixNoBC;
            matrix_Type tempMatrix ( *matrixNoBC );
            matrixNoBC.reset ( new matrix_Type ( M_localMap, tempMatrix.meanNumEntries() ) );
            *matrixNoBC += tempMatrix;
            M_blockPreconditioner->globalAssemble();
        }
    }
    *matrixNoBC += *M_matrixStokes;
}


template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::updateStabilization ( matrix_Type& matrixFull )
{

    if ( M_stabilization )
    {
        matrixFull += *M_matrixStabilization;
    }

}

template <typename Mesh, typename SolverType>
void OseenSolver<Mesh, SolverType>::updateSourceTerm ( source_Type const& source )
{
    vector_Type rhs ( M_localMap );

    VectorElemental M_elvec (M_velocityFESpace->fe().nbFEDof(), nDimensions);
    UInt nc = nDimensions;

    // loop on volumes: assembling source term
    for ( UInt i = 1; i <= M_velocityFESpace->mesh()->numVolumes(); ++i )
    {

        M_velocityFESpace->fe().updateFirstDerivQuadPt ( M_velocityFESpace->mesh()->volumeList ( i ) );

        M_elvec.zero();

        for ( UInt ic = 0; ic < nc; ++ic )
        {
            compute_vec ( source, M_elvec, M_velocityFESpace->fe(),  M_oseenData->dataTime()->time(), ic ); // compute local vector
            assembleVector ( *rhs, M_elvec, M_velocityFESpace->fe(), M_velocityFESpace->dof(), ic, ic * M_velocityFESpace->getDim() ); // assemble local vector into global one
        }
    }
    M_rightHandSideNoBC += rhs;
}




template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::iterate ( bcHandler_Type& bcHandler )
{

    LifeChrono chrono;

    // matrix and vector assembling communication
    M_Displayer.leaderPrint ( "  F-  Updating the boundary conditions ...     " );

    chrono.start();

    M_matrixNoBC->globalAssemble();

    matrixPtr_Type matrixFull ( new matrix_Type ( M_localMap, M_matrixNoBC->meanNumEntries() ) );

    updateStabilization ( *matrixFull );
    getFluidMatrix ( *matrixFull );

    vector_Type rightHandSideFull ( *M_rightHandSideNoBC );

    //     matrixFull.reset( new matrix_Type( *M_matrixNoBC ) );
    //     M_rightHandSideFull = M_rightHandSideNoBC;

    chrono.stop();

    M_Displayer.leaderPrintMax ( "done in ", chrono.diff() );

    // boundary conditions update
    M_Displayer.leaderPrint ("  F-  Applying boundary conditions ...         ");

    chrono.start();
    applyBoundaryConditions ( *matrixFull, rightHandSideFull, bcHandler );

    matrixFull->globalAssemble();
    chrono.stop();

    M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

    // solving the system
    M_linearSolver->setMatrix ( *matrixFull );

    boost::shared_ptr<MatrixEpetra<Real> > staticCast = boost::static_pointer_cast<MatrixEpetra<Real> > (matrixFull);

    Int numIter = M_linearSolver->solveSystem ( rightHandSideFull, *M_solution, staticCast );

    // if the preconditioner has been rese the stab terms are to be updated
    if ( numIter < 0 || numIter > M_iterReuseStabilization )
    {
        resetStabilization();
    }

    *M_residual  = *M_rightHandSideNoBC;
    *M_residual -= (*M_matrixNoBC) * (*M_solution);

    //M_residual.spy("residual");
} // iterate()


template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::reduceSolution ( Vector& velocityVector, Vector& pressureVector )
{
    vector_Type solution ( *M_solution, 0 );

    if ( false /*S_verbose*/ )
    {
        for ( UInt iDof = 0; iDof < M_velocityFESpace.fieldDim() * dimVelocity(); ++iDof )
        {
            velocityVector[ iDof ] = solution[ iDof ];
        }

        for ( UInt iDof = 0; iDof < dimPressure(); ++iDof )
        {
            pressureVector[ iDof ] = solution[ iDof + M_velocityFESpace.fieldDim() * dimVelocity() ];
        }
    }

}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::reduceResidual ( Vector& residualVector )
{
    vector_Type residual ( *M_residual, 0 );

    if ( false /*S_verbose*/ )
    {
        for ( UInt iDof = 0; iDof < M_velocityFESpace.fieldDim() * dimVelocity(); ++iDof )
        {
            residualVector[ iDof ] = residual[ iDof ];
        }

    }
}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::setBlockPreconditioner ( matrixPtr_Type blockPreconditioner )
{
    // blockPreconditioner.reset(new matrix_Type(M_monolithicMap, M_solid->getMatrixPtr()->getMeanNumEntries()));
    *blockPreconditioner += *M_blockPreconditioner;
}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::getFluidMatrix ( matrix_Type& matrixFull )
{
    M_matrixNoBC->globalAssemble();
    matrixFull += *M_matrixNoBC;
}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::postProcessingSetArea()
{
    M_postProcessing->set_area();
}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::postProcessingSetNormal()
{
    M_postProcessing->set_normal();
}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::postProcessingSetPhi()
{
    M_postProcessing->set_phi();
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::flux ( const markerID_Type& flag )
{
    return flux ( flag, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::flux ( const markerID_Type& flag,
                                          const vector_Type& solution )
{
    vector_Type velocityAndPressure ( solution, Repeated );
    vector_Type velocity ( this->M_velocityFESpace.map(), Repeated );
    velocity.subset ( velocityAndPressure );

    return M_postProcessing->flux ( velocity, flag );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::kineticNormalStress ( const markerID_Type& flag )
{
    return kineticNormalStress ( flag, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::kineticNormalStress ( const markerID_Type& flag,
                                                         const vector_Type& solution )
{
    vector_Type velocityAndPressure ( solution, Repeated );
    vector_Type velocity ( this->M_velocityFESpace.map(), Repeated );
    velocity.subset ( velocityAndPressure );

    return M_postProcessing->kineticNormalStress ( velocity, M_oseenData->density(), flag );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::area ( const markerID_Type& flag )
{
    return M_postProcessing->measure ( flag );
}

template<typename MeshType, typename SolverType>
Vector
OseenSolver<MeshType, SolverType>::normal ( const markerID_Type& flag )
{
    return M_postProcessing->normal ( flag );
}

template<typename MeshType, typename SolverType>
Vector
OseenSolver<MeshType, SolverType>::geometricCenter ( const markerID_Type& flag )
{
    return M_postProcessing->geometricCenter ( flag );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::pressure ( const markerID_Type& flag )
{
    return pressure ( flag, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::pressure (const markerID_Type& flag,
                                             const vector_Type& solution)
{
    vector_Type velocityAndPressure ( solution, Repeated );
    vector_Type pressure ( this->M_pressureFESpace.map(), Repeated );
    pressure.subset ( velocityAndPressure,
                      this->M_velocityFESpace.dim() *this->M_velocityFESpace.fieldDim() );

    // third argument is 1, to use the pressure finite element space (see PostProcessingBoundary docs)
    return M_postProcessing->average ( pressure, flag, 1 ) [0];
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::meanNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return meanNormalStress ( flag, bcHandler, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::meanNormalStress (const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& solution )
{
    if ( bcHandler.findBCWithFlag ( flag ).type() == Flux )
    {
        return -lagrangeMultiplier ( flag, bcHandler, solution );
    }
    else
    {
#ifdef HAVE_LIFEV_DEBUG
        M_Displayer.leaderPrint ( " !!! WARNING - OseenSolver::meanNormalStress( flag, bcHandler, solution) is returning an approximation \n" );
#endif
        return -pressure ( flag, solution ); // TODO: This is an approximation of the mean normal stress as the pressure.
        // A proper method should be coded in the PostprocessingBoundary class
        // to compute the exact mean normal stress.
    }
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::meanTotalNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return meanTotalNormalStress ( flag, bcHandler, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::meanTotalNormalStress (const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& solution )
{
    return meanNormalStress ( flag, bcHandler, solution ) - kineticNormalStress ( flag, solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::lagrangeMultiplier ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return lagrangeMultiplier ( flag, bcHandler, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::lagrangeMultiplier ( const markerID_Type&  flag,
                                                        bcHandler_Type& bcHandler,
                                                        const vector_Type& solution )
{
    // Create a list of Flux bcName_Type ??
    std::vector< bcName_Type > fluxBCVector = bcHandler.findAllBCWithType ( Flux );
    bcName_Type fluxbcName_Type = bcHandler.findBCWithFlag ( flag ).name();

    // Create a Repeated vector for looking to the lambda
    vector_Type velocityPressureLambda ( solution, Repeated );

    // Find the index associated to the correct Lagrange multiplier
    for ( UInt lmIndex = 0; lmIndex < static_cast <UInt> ( fluxBCVector.size() ); ++lmIndex )
        if ( fluxbcName_Type.compare ( fluxBCVector[ lmIndex ] ) == 0 )
            return velocityPressureLambda[M_velocityFESpace.fieldDim() * M_velocityFESpace.dof().numTotalDof()
                                          + M_pressureFESpace.dof().numTotalDof() + lmIndex];

    // If lmIndex has not been found a warning message is printed
    M_Displayer.leaderPrint (  "!!! Warning - Lagrange multiplier for Flux BC not found!\n" );
    return 0;
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::removeMean ( vector_Type& x )
{

    LifeChrono chrono;
    chrono.start();

    const UInt numVelocityComponent ( M_velocityFESpace.fieldDim() );
    const UInt velocityTotalDof ( M_velocityFESpace.dof().numTotalDof() );


    if ( M_pressureMatrixMass.get() == 0 )
    {
        M_pressureMatrixMass.reset ( new matrix_Type ( M_localMap ) );
    }

    for ( UInt iElement = 0; iElement < M_velocityFESpace.mesh()->numElements(); iElement++ )
    {
        chrono.start();
        // just to provide the id number in the assem_mat_mixed
        M_pressureFESpace.fe().update ( M_pressureFESpace.mesh()->element ( iElement ) );

        M_elementMatrixPreconditioner.zero();
        // mass
        chrono.start();
        mass ( 1, M_elementMatrixPreconditioner, M_pressureFESpace.fe(), 0, 0, M_velocityFESpace.fieldDim() );
        chrono.stop();

        chrono.start();
        assembleMatrix ( *M_pressureMatrixMass,
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

    M_pressureMatrixMass->globalAssemble();

    vector_Type ones ( *M_solution );
    ones = 1.0;

    Real mean;
    mean = ones * ( M_pressureMatrixMass * x );
    x += ( -mean );

    ASSERT ( std::fabs ( ones * ( M_pressureMatrixMass * x ) ) < 1e-9 , "after removeMean the mean pressure should be zero!");

    return mean;

} // removeMean()

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::applyBoundaryConditions ( matrix_Type&       matrix,
                                                             vector_Type&       rightHandSide,
                                                             bcHandler_Type& bcHandler )
{
    // M_rightHandSideFull = M_rightHandSideNoBC;

    // BC manage for the velocity
    if ( !bcHandler.bcUpdateDone() || M_recomputeMatrix )
    {
        bcHandler.bcUpdate ( *M_velocityFESpace.mesh(),
                             M_velocityFESpace.feBd(),
                             M_velocityFESpace.dof() );
    }

    // ignoring non-local entries, Otherwise they are summed up lately
    //vector_Type rightHandSideFull( rightHandSide, Repeated, Zero );
    // ignoring non-local entries, Otherwise they are summed up lately
    vector_Type rightHandSideFull ( rightHandSide, Unique );

    bcManage ( matrix, rightHandSideFull,
               *M_velocityFESpace.mesh(),
               M_velocityFESpace.dof(),
               bcHandler,
               M_velocityFESpace.feBd(),
               1.,
               M_oseenData->dataTime()->time() );

    rightHandSide = rightHandSideFull;

    if ( bcHandler.hasOnlyEssential() && M_diagonalize )
    {
        matrix.diagonalize ( M_velocityFESpace.fieldDim() *dimVelocity(),
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
OseenSolver<MeshType, SolverType>::setupPostProc( )
{
    M_postProcessing.reset ( new PostProcessingBoundary<mesh_Type> ( M_velocityFESpace.mesh(),
                                                                     &M_velocityFESpace.feBd(),
                                                                     &M_velocityFESpace.dof(),
                                                                     &M_pressureFESpace.feBd(),
                                                                     &M_pressureFESpace.dof(),
                                                                     M_localMap ) );
}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::setTolMaxIteration ( const Real& tolerance, const Int& maxIteration )
{
    M_linearSolver->setTolerance ( tolerance );
    M_linearSolver->setMaxNumIterations ( maxIteration );
}

} //end namespace LifeV
