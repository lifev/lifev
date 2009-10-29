/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-10-23

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
 \file MS_Algorithm.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-10-23
 */

#ifndef __MS_Algorithm_H
#define __MS_Algorithm_H 1

#include "Epetra_ConfigDefs.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <lifemc/lifesolver/MS_Model_MultiScale.hpp>

namespace LifeV {

//! @name MS_Algorithm global objects
//@{

enum algorithmsTypes
{
    Aitken,
    Newton
};

extern std::map< std::string, algorithmsTypes > algorithmMap;

//@}

//! MS_Algorithm - The MultiScale Algorithm Interface
/*!
 *  The MS_Algorithm class provides a general interface between the
 *  MS_Solver and the specific Algorithm to solve the problem.
 *
 *  @author Cristiano Malossi
 */
class MS_Algorithm
{
public:

    typedef MS_PhysicalModel::VectorType          VectorType;
    typedef MS_PhysicalModel::Vector_ptrType      Vector_ptrType;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Algorithm();

    //! Copy constructor
    /*!
     * \param algorithm - MS_Algorithm
     */
    MS_Algorithm( const MS_Algorithm& algorithm );

    //! Destructor
    virtual ~MS_Algorithm() {}

    //@}


    //! @name Methods
    //@{

    //! Operator=
    /*!
     * \param algorithm - MS_Algorithm
     */
    MS_Algorithm& operator=( const MS_Algorithm& algorithm );

    //@}


    //! @name Set Methods
    //@{

    //! Set the epetra communicator for the model
    /*!
     * \param comm - Epetra communicator
     */
    void SetCommunicator( const boost::shared_ptr< Epetra_Comm >& comm );

    //! Set the main MultiScale Problem
    /*!
     * \param multiscale - Main MultiScale problem
     */
    void SetMultiScaleProblem( const boost::shared_ptr< MS_Model_MultiScale > multiscale );

    //@}


    //! @name Get Methods
    //@{

    //! Get the type of the algorithm
    const algorithmsTypes& GetType() const
    {
        return M_type;
    }

    //! Get the Multiscale problem
    const boost::shared_ptr< MS_Model_MultiScale > GetMultiScaleProblem() const
    {
        return M_multiscale;
    }

    //! Get the coupling variables
    const Vector_ptrType GetCouplingVariables() const
    {
        return M_couplingVariables;
    }

    //! Get the coupling residuals
    const Vector_ptrType GetCouplingResiduals() const
    {
        return M_couplingResiduals;
    }

    //! Get the communicator
    const boost::shared_ptr< Epetra_Comm > GetCommunicator() const
    {
        return M_comm;
    }

    //! Get the subiterations maximum number
    const UInt& GetSubiterationsMaximumNumber() const
    {
        return M_SubiterationsMaximumNumber;
    }

    //! Get the required tolerance
    const Real& GetTolerance() const
    {
        return M_Tolerance;
    }

    //@}


    //! @name MultiScale Algorithm Virtual Methods
    //@{

    //! Setup the data of the algorithm
    virtual void SetupData( const GetPot& DataFile );

    //! Perform sub-iteration on the couplings
    virtual void SubIterate( void ) = 0;

    //! Display some information about the algorithm
    virtual void ShowMe( void );

    //@}

protected:

    //! @name Protected Methods
    //@{

    //@}

    algorithmsTypes                          M_type;

    boost::shared_ptr< MS_Model_MultiScale > M_multiscale;

    Vector_ptrType                           M_couplingVariables;
    Vector_ptrType                           M_couplingResiduals;

    boost::shared_ptr< Epetra_Comm >         M_comm;
    boost::shared_ptr< Displayer >           M_displayer;

    UInt                                     M_SubiterationsMaximumNumber;
    Real                                     M_Tolerance;
};

} // Namespace LifeV

#endif /* __MS_Algorithm_H */
