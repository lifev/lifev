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
    @brief 0D test with the Goldbeter model for simplified Calcium dynamics

    @date 04 - 2014
    @author Simone Rossi <simone.rossi@epfl.ch>

    @contributor
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

// To save the solution
#include <fstream>

// To solve the model
#include <lifev/electrophysiology/solver/IonicModels/IonicGoldbeter.hpp>

// To use LifeV
#include <lifev/core/LifeV.hpp>


using namespace LifeV;

//This is the norm of the precomputed solution
//we check the test against this value
#define SolutionTestNorm 639.4563875712343588

Int main ()
{

    //********************************************//
    // Creates a new model object representing the//
    // model from  Goldebeter model.              //
    //********************************************//
    std::cout << "Building Constructor for Goldbeter Model with parameters ... ";
    IonicGoldbeter  model;
    std::cout << " Done!" << std::endl;


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//
    model.showMe();


    //********************************************//
    // Initialize the solution to 0. The model    //
    // consist of three state variables.          //
    // The method model.Size()                    //
    // returns the number of state variables of   //
    // the model. rStates is the reference to the //
    // the vector states                          //
    //********************************************//
    std::cout << "Initializing solution vector...";
    std::vector<Real> unknowns (model.Size(), 0);
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
    std::vector<Real> rhs (model.Size(), 0);
    std::cout << " Done! "  << std::endl;

    //********************************************//
    // Simulation starts on t=0 and ends on t=TF. //
    // The timestep is given by dt                //
    //********************************************//
    Real TF (10);
    Real dt (0.01);

    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//
    std::string filename = "output.txt";
    std::ofstream output ("output.txt");

    //********************************************//
    // We record the norm of the solution to      //
    // check the failure of the test              //
    //********************************************//
    Real SolutionNorm = unknowns[0];

    //********************************************//
    // Time loop starts.                          //
    //********************************************//
    std::cout << "Time loop starts...\n";
    for ( Real t = 0; t < TF; )
    {

        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//
        model.computeRhs ( unknowns, rhs);

        //********************************************//
        // Use forward Euler method to advance the    //
        // solution in time.                          //
        //********************************************//
        unknowns.at (0) = unknowns.at (0)  + dt * rhs.at (0);
        unknowns.at (1) = unknowns.at (1)  + dt * rhs.at (1);

        //********************************************//
        // Writes solution on file.                   //
        //********************************************//
        output << t << ", " << unknowns.at (0) << ", " << unknowns.at (1) << "\n";

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
    std::cout << "\n...Time loop ends.\n";
    std::cout << "Solution written on file: " << filename << "\n";
    //********************************************//
    // Close exported file.                       //
    //********************************************//
    output.close();

    Real returnValue;
    Real err = std::abs (SolutionTestNorm - SolutionNorm) / std::abs (SolutionTestNorm);
    std::cout << std::setprecision (20) << "\nError: " << err << "\nSolution norm: " << SolutionNorm << "\n";
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

#undef SolutionTestNorm
