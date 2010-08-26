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
 *  @brief MultiScale Algorithm
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 23-10-2009
 */

#ifndef MS_Algorithm_H
#define MS_Algorithm_H 1

#include <lifemc/lifesolver/MS_Definitions.hpp>
#include <lifemc/lifesolver/MS_Model_MultiScale.hpp>

namespace LifeV {

//! MS_Algorithm - The MultiScale Algorithm Interface
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_Algorithm class provides a general interface between the
 *  MS_Solver and the specific Algorithm to solve the problem.
 *
 */
class MS_Algorithm
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Algorithm();

    //! Copy constructor
    /*!
     * @param algorithm MS_Algorithm
     */
    MS_Algorithm( const MS_Algorithm& algorithm );

    //! Destructor
    virtual ~MS_Algorithm() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param algorithm MS_Algorithm
     * @return reference to a copy of the class
     */
    MS_Algorithm& operator=( const MS_Algorithm& algorithm );

    //@}


    //! @name MultiScale Algorithm Virtual Methods
    //@{

    //! Setup the data of the algorithm using a data file
    /*!
     * @param FileName Name of the data file.
     */
    virtual void SetupData( const std::string& FileName );

    //! Perform sub-iteration on the coupling variables
    virtual void SubIterate();

    //! Update coupling variables for the next time step.
    virtual void UpdateCouplingVariables();

    //! Display some information about the algorithm
    virtual void ShowMe();

    //@}


    //! @name Methods
    //@{

    //! Initialize coupling variables for the first time step.
    void InitializeCouplingVariables();

    //@}


    //! @name Set Methods
    //@{

    //! Set the epetra communicator for the model
    /*!
     * @param comm Epetra communicator
     */
    void SetCommunicator( const MS_Comm_PtrType& comm );

    //! Set the main MultiScale model
    /*!
     * @param model MultiScale model
     */
    void SetModel( const MS_Model_PtrType model );

    //@}


    //! @name Get Methods
    //@{

    //! Get the type of the algorithm
    /*!
     * @return type of the algorithm
     */
    const algorithmsTypes& GetType() const;

    //! Get the Multiscale problem
    /*!
     * @return shared_ptr to the MultiScale problem
     */
    const boost::shared_ptr< MS_Model_MultiScale > GetMultiScaleProblem() const;

    //! Get the coupling variables
    /*!
     * @return pointer to the coupling variables vector
     */
    const MS_Vector_PtrType GetCouplingVariables() const;

    //! Get the coupling residuals
    /*!
     * @return pointer to the coupling residuals vector
     */
    const MS_Vector_PtrType GetCouplingResiduals() const;

    //! Get the communicator
    /*!
     * @return pointer to the communicator
     */
    const MS_Comm_PtrType GetCommunicator() const;

    //! Get the subiterations maximum number
    /*!
     * @return maximum number of subiterations
     */
    const UInt& GetSubiterationsMaximumNumber() const;

    //! Get the required tolerance
    /*!
     * @return tolerance
     */
    const Real& GetTolerance() const;

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Save on a Matlab file the information about the convergence of the algorithm.
    /*!
     * @param SubiterationsNumber Number of subiterations performed.
     * @param residual Residual.
     */
    void Save( const UInt& SubiterationsNumber, const Real& residual );

    //! Check if the tolerance has been satisfied
    /*!
     * @return true if the tolerance is satisfied
     */
    bool ToleranceSatisfied();

    //@}

    algorithmsTypes                          M_type;

    boost::shared_ptr< MS_Model_MultiScale > M_multiscale;

    MS_Vector_PtrType                        M_couplingVariables;
    MS_Vector_PtrType                        M_couplingResiduals;

    MS_Comm_PtrType                          M_comm;
    boost::shared_ptr< Displayer >           M_displayer;

    UInt                                     M_SubiterationsMaximumNumber;
    Real                                     M_Tolerance;
};

} // Namespace LifeV

#endif /* MS_Algorithm_H */
