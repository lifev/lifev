//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief MultiScale Solver
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 28-09-2009
 */

#ifndef MS_Solver_H
#define MS_Solver_H 1

#include <lifemc/lifesolver/MS_Definitions.hpp>

#include <lifemc/lifesolver/MS_Algorithm.hpp>
#include <lifemc/lifesolver/MS_Algorithm_Aitken.hpp>
#include <lifemc/lifesolver/MS_Algorithm_Newton.hpp>

#include <lifemc/lifesolver/MS_Model_MultiScale.hpp>

namespace LifeV {

//! MS_Solver - The MultiScale problem solver
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_Solver class provides a series of functions to create and
 *  solve a general MultiScale problem.
 */
class MS_Solver
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Solver();

    //! Copy constructor
    /*!
     * @param algorithm MS_Solver
     */
    MS_Solver( const MS_Solver& solver );

    //! Destructor
    ~MS_Solver() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param solver MS_solver
     * @return reference to a copy of the class
     */
    MS_Solver& operator=( const MS_Solver& solver );

    //@}


    //! @name Methods
    //@{

    //! Set the epetra communicator for the MultiScale problem
    /*!
     * @param comm Epetra communicator
     */
    void SetCommunicator( const boost::shared_ptr< Epetra_Comm >& comm );

    //! Setup the problem
    /*!
     * @param dataFile Name and path of the data file
     */
    void SetupProblem( const std::string& dataFile );

    //! Run the time-loop to solve the MultiScale problem
    /*!
     * @return 0: EXIT_SUCCESS, 1: EXIT_FAILURE
     */
    bool SolveProblem();

    //! Display some information about the MultiScale problem (to be called after SetupProblem)
    void ShowMe();

    //@}

private:

    //! @name Private Methods
    //@{


    //@}

    // The main multiscale model
    boost::shared_ptr< MS_Model_MultiScale > M_multiscale;

    // Algorithm for subiterations
    Algorithm_ptrType                        M_algorithm;

    // PhysicalData container
    boost::shared_ptr< MS_PhysicalData >     M_dataPhysics;

    // DataTime container
    boost::shared_ptr< DataTime >            M_dataTime;

    // Communicator
    boost::shared_ptr< Epetra_Comm >         M_comm;

    // Displayer tool for MPI processes
    boost::shared_ptr< Displayer >           M_displayer;

    // Chrono performances
    Chrono                                   M_chrono;
};

} // Namespace LifeV

#endif /* MS_Solver_H */
