/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-09-28

 Copyright (C) 2009 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/**
 \file MS_Solver.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-09-28
 */

#ifndef __MS_Solver_H
#define __MS_Solver_H 1

#include "Epetra_ConfigDefs.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <life/lifecore/life.hpp>
#include <life/lifealg/generalizedAitken.hpp>

#include <lifemc/lifesolver/MS_Model_MultiScale.hpp>
#include <lifemc/lifearray/ContainerOfVectors.hpp>

#include <boost/array.hpp>
#include <boost/algorithm/string.hpp>

namespace LifeV {

//! MS_Solver - The MultiScale problem solver
/*!
 *  The MS_Solver class provides a series of functions to create and
 *  solve a general MultiScale problem.
 *
 *  @author Cristiano Malossi
 */
class MS_Solver
{
public:

    typedef singleton< factory< MS_PhysicalModel, modelsTypes > >       FactoryModels;
    typedef singleton< factory< MS_PhysicalCoupling, couplingsTypes > > FactoryCouplings;

    typedef ContainerOfVectors< EpetraVector > VectorType;

    //! @name Constructors, Destructor
    //@{

    //! Constructor
    MS_Solver();

    //! Copy constructor
    /*!
     * \param algorithm - MS_Solver
     */
    MS_Solver( const MS_Solver& solver );

    //! Destructor
    ~MS_Solver() {}

    //@}


    //! @name Methods
    //@{

    //! Operator=
    /*!
     * \param solver - MS_Solver
     */
    MS_Solver& operator=( const MS_Solver& solver );

    //! Set the epetra communicator for the MultiScale problem
    /*!
     * \param comm - Epetra communicator
     */
    void SetCommunicator( const boost::shared_ptr< Epetra_Comm >& comm );

    //! Setup the problem
    /*!
     * \param dataFile - Name and path of data file
     */
    void SetupProblem( const std::string& dataFile );

    //! Run the time-loop to solve the MultiScale problem
    void SolveProblem( void );

    //! Display some information about the MultiScale problem (to be called after SetupProblem)
    void ShowMe( void );

    //@}

private:

    //! @name Private Methods
    //@{


    //@}

    // The main multiscale model
    MS_Model_MultiScale                  M_multiscale;

    // PhysicalData container
    boost::shared_ptr< MS_PhysicalData > M_dataPhysics;

    // DataTime container
    boost::shared_ptr< DataTime >        M_dataTime;

    // Communicator
    boost::shared_ptr< Epetra_Comm >     M_comm;

    // Displayer tool for MPI processes
    boost::shared_ptr< Displayer >       M_displayer;

    // Chrono performances
    Chrono                               M_chrono;

    // Container of coupling variables
    VectorType                           M_couplingVariables;

    // Container of coupling residuals
    VectorType                           M_couplingResiduals;

    // Aitken Method
    generalizedAitken< VectorType >      M_generalizedAitken;
    UInt                                 M_subITMax;
    Real                                 M_tolerance;
};

} // Namespace LifeV

#endif /* __MS_Solver_H */
