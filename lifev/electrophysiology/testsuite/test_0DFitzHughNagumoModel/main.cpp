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
    @brief 0D test with the Fitz-Hugh Nagumo model.

    @date 01-03-2013
    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @contributor
    @mantainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
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

#include <lifev/core/util/LifeChrono.hpp>

#include <lifev/electrophysiology/solver/IonicModels/IonicFitzHughNagumo.hpp>
#include <lifev/core/LifeV.hpp>


using namespace LifeV;


Real EulerExplicit (Real& dt, const Real& TF, IonicFitzHughNagumo model, const Real& I, std::ofstream& output);

//This is the norm of the precomputed solution
//we check the test against this value
#define SolutionTestNorm 456913.59773800277617


Int main ( Int argc, char** argv )
{
    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm>  Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "% using MPI" << std::endl;
    }

    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    std::cout << "Importing parameters list...";
    Teuchos::ParameterList FHNParameterList = * ( Teuchos::getParametersFromXmlFile ( "FitzHughNagumoParameters.xml" ) );
    std::cout << " Done!" << std::endl;


    //********************************************//
    // Creates a new model object representing the//
    // model from Fitz-Hugh Nagumo. The           //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    std::cout << "Building Constructor for Fitz-Hugh Nagumo Model with parameters ... ";
    IonicFitzHughNagumo  model ( FHNParameterList );
    std::cout << " Done!" << std::endl;


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//
    model.showMe();

    //********************************************//
    // Simulation starts on t=0 and ends on t=TF. //
    // The timestep is given by dt                //
    //********************************************//

    Real TF ( FHNParameterList.get ("TF", 300.0) );
    Real dt ( FHNParameterList.get ("dt", 0.01) );

    std::cout << "Time parameters : " << std::endl;
    std::cout << "TF = " << TF << std::endl;
    std::cout << "dt = " << dt << std::endl;

    std::string filename = "test_expeuler.txt";


    std::ofstream output (filename.c_str() );

    //********************************************//
    // Time loop starts.                          //
    //********************************************//
    std::cout << "Time loop starts...\n\n\n";

    LifeChrono chrono;
    chrono.start();

    Real valueToTest;

    valueToTest = EulerExplicit (dt, TF, model, FHNParameterList.get ("Iapp", 2000.0), output);


    chrono.stop();
    std::cout << "\n...Time loop ends.\n";
    std::cout << "\nElapsed time : " << chrono.diff() << std::endl;
    std::cout << "Solution written on file: " << filename << "\n";
    //********************************************//
    // Close exported file.                       //
    //********************************************//
    output.close();

    //! Finalizing Epetra communicator
    MPI_Finalize();

    Real returnValue;
    Real err = std::abs (valueToTest - SolutionTestNorm) / std::abs (SolutionTestNorm);
    std::cout << std::setprecision (20) << "\nError: " << err << "\nSolution norm: " << valueToTest << "\n";
    if ( err > 1e-12 )
    {
        returnValue = EXIT_FAILURE; // Norm of solution did not match
    }
    else
    {
        returnValue = EXIT_SUCCESS;
    }
    return ( returnValue );
}

Real EulerExplicit (Real& dt, const Real& TF, IonicFitzHughNagumo model, const Real& I, std::ofstream& output)
{
    std::vector<Real> unknowns ( model.Size(), 0.0);
    unknowns.at (0) = 1e-8;
    unknowns.at (1) = 0.3;

    std::vector<Real> rhs ( model.Size(), 0.0);
    Real Iapp;

    std::cout << "Computing using Explicit Euler" << std::endl;

    //********************************************//
    // We record the norm of the solution to      //
    // check the failure of the test              //
    //********************************************//
    Real SolutionNorm = unknowns[0];

    for ( Real t = 0; t < TF; )
    {

        //********************************************//
        // Compute Calcium concentration. Here it is  //
        // given as a function of time.               //
        //********************************************//
        if ( t > 2.0 && t < 2.1 )
        {
            Iapp = I;
        }
        else
        {
            Iapp = 0;
        }
        model.setAppliedCurrent (Iapp);
        std::cout << "\r " << t << " ms.       " << std::flush;

        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//
        model.computeRhs ( unknowns, rhs);
        model.addAppliedCurrent (rhs);

        //********************************************//
        // Use forward Euler method to advance the    //
        // solution in time.                          //
        //********************************************//
        unknowns.at (0) = unknowns.at (0)  + dt * rhs.at (0);
        unknowns.at (1) = unknowns.at (1)  + dt * rhs.at (1);


        //********************************************//
        // Writes solution on file.                   //
        //********************************************//
        output << t << " " << unknowns.at (0) << " " << unknowns.at (1) << "\n";

        //********************************************//
        // Update the time.                           //
        //********************************************//
        t = t + dt;

        //********************************************//
        // Update the norm of the solution to check   //
        // test failure                               //
        //********************************************//
        SolutionNorm += unknowns[0];
    }

    return SolutionNorm;
}


#undef SolutionTestNorm








