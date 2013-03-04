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
 *  @file
 *  @brief File containing the Multiscale Broyden Algorithm
 *
 *  @date 26-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/solver/MultiscaleAlgorithmBroyden.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleAlgorithmBroyden::MultiscaleAlgorithmBroyden() :
    multiscaleAlgorithm_Type     (),
    M_solver                     (),
    M_jacobian                   (),
    M_initializeAsIdentityMatrix ( false ),
    M_iterationsLimitReached     ( false ),
    M_iterationsLimitForReset    ( 1 ),
    M_orthogonalization          ( false ),
    M_orthogonalizationSize      ( 1 ),
    M_orthogonalizationContainer (),
    M_truncate                   ( true )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8014 ) << "MultiscaleAlgorithmBroyden::MultiscaleAlgorithmBroyden() \n";
#endif

    M_type = Broyden;
}

// ===================================================
// Multiscale Algorithm Virtual Methods
// ===================================================
void
MultiscaleAlgorithmBroyden::setupData ( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8014 ) << "MultiscaleAlgorithmBroyden::setupData( fileName ) \n";
#endif

    // Read parameters
    multiscaleParameterListPtr_Type solverParametersList = Teuchos::rcp ( new Teuchos::ParameterList );
    solverParametersList = Teuchos::getParametersFromXmlFile ( fileName );

    // Set main parameters
    setAlgorithmName ( solverParametersList->sublist ( "Multiscale", true, "" ) );
    setAlgorithmParameters ( solverParametersList->sublist ( "Multiscale Algorithm", true, "" ) );

    // Set main parameters
    M_solver.setCommunicator ( M_comm );
    M_solver.setParameters ( solverParametersList->sublist ( "Linear Solver List", true, "" ) );
}

void
MultiscaleAlgorithmBroyden::setupAlgorithm()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8014 ) << "MultiscaleAlgorithmBroyden::setupAlgorithm() \n";
#endif

    multiscaleAlgorithm_Type::setupAlgorithm();

#ifdef HAVE_HDF5
#if ( H5_VERS_MAJOR > 1 || ( H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 8 ) || ( H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 8 && H5_VERS_RELEASE >= 7 ) )
#ifdef BROYDEN_IMPORTJACOBIAN
    // Import Jacobian from previous simulation
    if ( multiscaleProblemStep > 0 )
    {
        importJacobianFromHDF5();
    }
#endif
#endif
#endif
}

void
MultiscaleAlgorithmBroyden::subIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8014 ) << "MultiscaleAlgorithmBroyden::subIterate() \n";
#endif

    multiscaleAlgorithm_Type::subIterate();

    // Verify tolerance
    if ( checkResidual ( 0 ) )
    {
#ifdef HAVE_HDF5
#if ( H5_VERS_MAJOR > 1 || ( H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 8 ) || ( H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 8 && H5_VERS_RELEASE >= 7 ) )
        exportJacobianToHDF5();
#endif
#endif
        return;
    }

    M_multiscale->exportCouplingVariables ( *M_couplingVariables );

    multiscaleVectorPtr_Type delta ( new multiscaleVector_Type ( *M_couplingResiduals, Unique ) );
    *delta = 0.0;

    for ( UInt subIT (1); subIT <= M_subiterationsMaximumNumber; ++subIT )
    {
        // Compute the Jacobian
        if ( subIT == 1 )
        {
            if ( M_jacobian.get() == 0 || M_iterationsLimitReached || M_multiscale->topologyChange() )
            {
                assembleJacobianMatrix();
                M_iterationsLimitReached = false;
            }
        }
        else
        {
            broydenJacobianUpdate ( *delta );
        }

        // Set matrix and RHS
        M_solver.setOperator ( M_jacobian );
        M_solver.setRightHandSide ( M_couplingResiduals );

        // Solve Newton (without changing the sign of the residual)
        M_solver.solve ( delta );

        // Changing the sign of the solution (this is important for the broydenJacobianUpdate method)
        *delta *= -1;

        // Update Coupling Variables using the Broyden Method
        *M_couplingVariables += *delta;

        // Import Coupling Variables inside the coupling blocks
        M_multiscale->importCouplingVariables ( *M_couplingVariables );

        if ( subIT >= M_iterationsLimitForReset )
        {
            M_iterationsLimitReached = true;
        }

        // Verify tolerance
        if ( checkResidual ( subIT ) )
        {
#ifdef HAVE_HDF5
#if ( H5_VERS_MAJOR > 1 || ( H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 8 ) || ( H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 8 && H5_VERS_RELEASE >= 7 ) )
            exportJacobianToHDF5();
#endif
#endif
            return;
        }
    }

    save ( M_subiterationsMaximumNumber, M_couplingResiduals->norm2() );

    multiscaleErrorCheck ( Tolerance, "Broyden algorithm residual: " + number2string ( M_couplingResiduals->norm2() ) +
                           " (required: " + number2string ( M_tolerance ) + ")\n", M_multiscale->communicator() == 0 );
}

void
MultiscaleAlgorithmBroyden::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        multiscaleAlgorithm_Type::showMe();

        std::cout << "Initialize as identity matrix        = " << M_initializeAsIdentityMatrix << std::endl;
        std::cout << "Reset limit                          = " << M_iterationsLimitForReset << std::endl;
        std::cout << "Enable orthogonalization             = " << M_orthogonalization << std::endl;
        std::cout << "Orthogonalization memory size        = " << M_orthogonalizationSize << std::endl;
        std::cout << std::endl << std::endl;
    }
}

// ===================================================
// Set Methods
// ===================================================
void
MultiscaleAlgorithmBroyden::setAlgorithmParameters ( const multiscaleParameterList_Type& parameterList )
{
    multiscaleAlgorithm_Type::setAlgorithmParameters ( parameterList );

    M_initializeAsIdentityMatrix = parameterList.get<bool> ( "Initialize as identity matrix" );
    M_iterationsLimitForReset    = static_cast <UInt> ( M_subiterationsMaximumNumber
                                                        * parameterList.get<Real> ( "Iterations limit for reset" ) );

    M_orthogonalization          = parameterList.get<bool> ( "Orthogonalization" );
    M_orthogonalizationSize      = parameterList.get<UInt> ( "Orthogonalization size" );
}

// ===================================================
// Private Methods
// ===================================================
void
MultiscaleAlgorithmBroyden::assembleJacobianMatrix()
{
    // Compute the Jacobian matrix
    M_jacobian.reset ( new multiscaleMatrix_Type ( M_couplingVariables->map(), 50 ) );

    if ( M_initializeAsIdentityMatrix )
    {
        M_jacobian->insertValueDiagonal ( 1 );
    }
    else
    {
        M_multiscale->exportJacobian ( *M_jacobian );
    }

    M_jacobian->globalAssemble();

    //M_jacobian->spy( "Jacobian" );

    // Reset orthogonalization
    M_orthogonalizationContainer.clear();
}

void
MultiscaleAlgorithmBroyden::broydenJacobianUpdate ( const multiscaleVector_Type& delta )
{
    M_jacobian->openCrsMatrix();

    // Compute the Broyden update
    if ( M_orthogonalization )
    {
        // Orthogonalize the vector
        multiscaleVector_Type orthogonalization ( delta, Unique );
        for ( containerIterator_Type i = M_orthogonalizationContainer.begin(); i != M_orthogonalizationContainer.end() ; ++i )
        {
            orthogonalization -= orthogonalization.dot ( *i ) * *i;
        }
        orthogonalization /= orthogonalization.norm2();

        // Update orthogonalization
        orthogonalizationUpdate ( delta );

        // Update the Jacobian
        M_jacobian->addDyadicProduct ( ( *M_couplingResiduals ) / orthogonalization.dot ( delta ), orthogonalization );
    }
    else
    {
        // Update the Jacobian
        //M_jacobian->addDyadicProduct( ( *M_couplingResiduals + minusCouplingResidual - *M_jacobian * delta ) / delta.dot( delta ), delta );
        M_jacobian->addDyadicProduct ( ( *M_couplingResiduals ) / delta.dot ( delta ), delta );
    }

    M_jacobian->globalAssemble();

    //M_jacobian->spy( "Jacobian" );
}

void
MultiscaleAlgorithmBroyden::orthogonalizationUpdate ( const multiscaleVector_Type& delta )
{
    if ( M_orthogonalizationContainer.size() < M_orthogonalizationSize )
    {
        M_orthogonalizationContainer.push_back ( delta );
    }
    else
    {
        M_orthogonalizationContainer.pop_front();
        M_orthogonalizationContainer.push_back ( delta );
    }
}

#ifdef HAVE_HDF5

void
MultiscaleAlgorithmBroyden::exportJacobianToHDF5()
{
    // If no coupling variable are present, the Jacobian is 0x0
    if ( M_couplingVariables->size() == 0 )
    {
        return;
    }

    if ( M_multiscale->globalData()->dataTime()->timeStepNumber() % multiscaleSaveEachNTimeSteps == 0 || M_multiscale->globalData()->dataTime()->isLastTimeStep() )
        if ( !M_jacobian.get() == 0 )
        {
            if ( M_comm->MyPID() == 0 )
            {
                std::cout << " MS-  Exporting Jacobian matrix ...            " << std::flush;
            }

            LifeChrono exportJacobianChrono;
            exportJacobianChrono.start();

            M_jacobian->exportToHDF5 ( multiscaleProblemFolder + multiscaleProblemPrefix + "_AlgorithmJacobian" + "_" + number2string ( multiscaleProblemStep ), number2string ( M_multiscale->globalData()->dataTime()->timeStepNumber() ), M_truncate );
            M_truncate = false;

            exportJacobianChrono.stop();
            Real jacobianChrono ( exportJacobianChrono.globalDiff ( *M_comm ) );
            if ( M_comm->MyPID() == 0 )
            {
                std::cout << "done in " << jacobianChrono << " s." << std::endl;
            }

            //M_jacobian->spy( multiscaleProblemFolder + multiscaleProblemPrefix + "_AlgorithmJacobianBroydenExported" + "_" + number2string( multiscaleProblemStep ) + "_" + number2string( M_multiscale->globalData()->dataTime()->timeStepNumber() ) );
        }
}

void
MultiscaleAlgorithmBroyden::importJacobianFromHDF5()
{
    // If no coupling variable are present, the Jacobian is 0x0
    if ( M_couplingVariables->size() == 0 )
    {
        return;
    }

    if ( M_comm->MyPID() == 0 )
    {
        std::cout << " MS-  Importing Jacobian matrix                " << std::flush;
    }

    LifeChrono importJacobianChrono;
    importJacobianChrono.start();

    M_jacobian.reset ( new multiscaleMatrix_Type ( M_couplingVariables->map(), 50 ) );
    M_jacobian->importFromHDF5 ( multiscaleProblemFolder + multiscaleProblemPrefix + "_AlgorithmJacobian" + "_" + number2string ( multiscaleProblemStep - 1 ), number2string ( M_multiscale->globalData()->dataTime()->timeStepNumber() ) );

    importJacobianChrono.stop();
    Real jacobianChrono ( importJacobianChrono.globalDiff ( *M_comm ) );
    if ( M_comm->MyPID() == 0 )
    {
        std::cout << "done in " << jacobianChrono << " s. (Time " << M_multiscale->globalData()->dataTime()->time() << ", Iteration " << M_multiscale->globalData()->dataTime()->timeStepNumber() << ")" << std::endl;
    }

    //M_jacobian->spy( multiscaleProblemFolder + multiscaleProblemPrefix + "_AlgorithmJacobianImported" + "_" + number2string( multiscaleProblemStep ) + "_" + number2string( timeInteger ) );
}

#endif

} // Namespace Multiscale
} // Namespace LifeV
