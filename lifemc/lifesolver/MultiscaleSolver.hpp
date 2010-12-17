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
 *  @brief File containing the MultiScale Solver
 *
 *  @date 28-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleSolver_H
#define MultiscaleSolver_H 1

#include <lifemc/lifesolver/MultiscaleDefinitions.hpp>

#include <lifemc/lifesolver/MultiscaleAlgorithm.hpp>
#include <lifemc/lifesolver/MultiscaleAlgorithmAitken.hpp>
#include <lifemc/lifesolver/MultiscaleAlgorithmExplicit.hpp>
#include <lifemc/lifesolver/MultiscaleAlgorithmNewton.hpp>

#include <lifemc/lifesolver/MultiscaleModelMultiscale.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleSolver - The MultiScale problem solver
/*!
 *  @author Cristiano Malossi
 *
 *  The MultiscaleSolver class provides a series of functions to create and
 *  solve a general MultiScale problem.
 */
class MultiscaleSolver
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleSolver();

    //! Destructor
    virtual ~MultiscaleSolver() {}

    //@}


    //! @name Methods
    //@{

    //! Set the epetra communicator for the MultiScale problem
    /*!
     * @param comm Epetra communicator
     */
    void setCommunicator( const multiscaleCommPtr_Type& comm );

    //! Setup the problem
    /*!
     * @param FileName Name of the data file.
     * @param problemName the name of the problem (useful for saving data in a specific folder)
     */
    void setupProblem( const std::string& fileName, const std::string& problemName );

    //! Run the time-loop to solve the MultiScale problem
    /*!
     * @return 0: EXIT_SUCCESS, 1: EXIT_FAILURE
     */
    bool solveProblem( const Real& externalResidual = -1 );

    //! Display some information about the MultiScale problem (to be called after SetupProblem)
    void showMe();

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleSolver( const MultiscaleSolver& solver );

    MultiscaleSolver& operator=( const MultiscaleSolver& solver );

    //@}

    // The main model (can be a specific model or a MultiScale model)
    multiscaleModelPtr_Type          M_model;

    // Algorithm for subiterations
    multiscaleAlgorithmPtr_Type      M_algorithm;

    // PhysicalData container
    multiscaleDataPtr_Type           M_globalData;

    // Communicator
    multiscaleCommPtr_Type           M_comm;

    // Displayer tool for MPI processes
    boost::shared_ptr< Displayer >   M_displayer;

    // Chrono performances
    Chrono                           M_chrono;
};

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleSolver_H */
