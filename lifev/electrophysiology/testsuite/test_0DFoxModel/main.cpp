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
    @brief 0D test with the Fox model

    @date 04 - 2013
    @author Marie Dupraz <dupraz.marie@gmail.com>

    @contributor
    @mantainer Marie Dupraz <dupraz.marie@gmail.com>
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <fstream>

#include <lifev/electrophysiology/solver/IonicModels/IonicFox.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

using namespace LifeV;

//This is the norm of the precomputed solution
//we check the test against this value
#define SolutionTestNorm -593906.94845120306127

Int main ( Int argc, char** argv )
{
    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    Epetra_MpiComm Comm (MPI_COMM_WORLD);
    if ( Comm.MyPID() == 0 )
    {
        std::cout << "% using MPI" << std::endl;
    }

    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    std::cout << "Importing parameters list...";
    Teuchos::ParameterList parameterList = * ( Teuchos::getParametersFromXmlFile ( "FoxParameters.xml" ) );
    std::cout << " Done!" << std::endl;

    //********************************************//
    // Creates a new model object representing the//
    // model from Fox 2002. The                   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    std::cout << "Building Constructor for Fox Model with parameters ... ";
    IonicFox  model;
    std::cout << " Done!" << std::endl;

    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//
    model.showMe();

    //********************************************//
    // Initialize the solution with the default   //
    // values                                     //
    //********************************************//
    std::cout << "Initializing solution vector...";
    std::vector<Real> unknowns (model.Size(), 0 );
    model.initialize (unknowns);
    std::cout << " Done!" << std::endl;

    //********************************************//
    // Initialize the rhs to 0. The rhs is the    //
    // vector containing the numerical values of  //
    // the time derivatives of the state          //
    // variables, that is, the right hand side of //
    // the differential equation.                 //
    //********************************************//
    std::cout << "Initializing rhs..." ;
    std::vector<Real> rhs (model.Size() + 1, 0);
    std::cout << " Done! "  << std::endl;

    //********************************************//
    // The model needs as external informations   //
    // the contraction velocity and the Calcium   //
    // concentration.                             //
    //********************************************//
    Real Iapp (0.0);

    //********************************************//
    // Simulation starts on t=0 and ends on t=TF. //
    // The timestep is given by dt                //
    //********************************************//
    Real TF     = parameterList.get ( "endTime", 5.0 );
    Real dt     = parameterList.get ( "timeStep", 0.005 );
    Real timeSt = parameterList.get ( "stimuliTime", 1.0 );
    Real stInt  = parameterList.get ( "stimuliInterval", 1000.0 );

    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//
    std::string filename = "output.txt";
    std::ofstream output  ("output.txt");


    //********************************************//
    // Time loop starts.                          //
    //********************************************//
    std::cout << "Time loop starts...\n";


    int iter (0);
    int savedt ( parameterList.get ( "savedt", 1.0) / dt );

    //********************************************//
    // We record the norm of the solution to      //
    // check the failure of the test              //
    //********************************************//
    Real SolutionNorm = unknowns[0];

    for ( Real t = 0; t < TF; )
    {

        //********************************************//
        // Compute the applied current. This is a     //
        // simple switch.                             //
        //********************************************//
        if ( t >= timeSt && t <= timeSt + 1.0 )
        {
            Iapp = 80.;  //0.516289;
            if ( t >= timeSt + 1.0 - dt && t <= timeSt + 1.0 )
            {
                timeSt = timeSt + stInt;
            }
        }
        else
        {
            Iapp = 0.;
        }

        std::cout << "\r " << t << " ms.       " << std::flush;

        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//
        model.setAppliedCurrent (Iapp);
        model.computeRhs ( unknowns, rhs );
        model.addAppliedCurrent (rhs);


        //********************************************//
        // Writes solution on file.                   //
        //********************************************//


        if ( iter % savedt == 0)
        {
            //********************************************//
            // Update the norm of the solution to check   //
            // test failure                               //
            //********************************************//
            SolutionNorm += unknowns[0];
            output << t << ", " << unknowns.at (0) << ", " << unknowns.at (10) << ", " << unknowns.at (11)
                   << ", " << rhs.at (13) << "\n"; //

        }

        //********************************************//
        // Use forward Euler method to advance the    //
        // solution in time.                          //
        //********************************************//

        for (int j (0); j <= 12; ++j)
        {
            unknowns.at (j) = unknowns.at (j)   + dt * rhs.at (j);
        }

        //********************************************//
        // Update the time.                           //
        //********************************************//
        t = t + dt;
    }

    std::cout << "\n...Time loop ends.\n";
    std::cout << "Solution written on file: " << filename << "\n";

    //********************************************//
    // Close exported file.                       //
    //********************************************//
    output.close();


    //! Finalizing Epetra communicator
    MPI_Finalize();
    Real returnValue;
    Real err = std::abs (SolutionTestNorm - SolutionNorm) / std::abs (SolutionTestNorm);
    std::cout << std::setprecision (20) << "\nError: " << err << "\nSolution norm: " << SolutionNorm << "\n";
    if ( err > 1e-12 )

    {
        std::cout << "\nTest Failed!\n";
        returnValue = EXIT_FAILURE; // Norm of solution did not match
    }
    else
    {
        returnValue = EXIT_SUCCESS;
    }
    return ( returnValue );
}


#undef SolutionTestNorm
