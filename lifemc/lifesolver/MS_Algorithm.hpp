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
 *  @brief File containing the MultiScale Algorithm
 *
 *  @date 23-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MS_Algorithm_H
#define MS_Algorithm_H 1

#include <lifemc/lifesolver/MS_Definitions.hpp>
#include <lifemc/lifesolver/MS_Model_MultiScale.hpp>

namespace LifeV
{

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

    //! @name Type definitions
    //@{

    typedef boost::shared_ptr< MS_Model_MultiScale >                   multiscaleModelPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Algorithm();

    //! Destructor
    virtual ~MS_Algorithm() {}

    //@}


    //! @name MultiScale Algorithm Virtual Methods
    //@{

    //! Setup the data of the algorithm using a data file
    /*!
     * @param FileName Name of the data file.
     */
    virtual void setupData( const std::string& fileName );

    //! Perform sub-iteration on the coupling variables
    virtual void subIterate();

    //! Update coupling variables for the next time step.
    virtual void updateCouplingVariables();

    //! Display some information about the algorithm
    virtual void showMe();

    //@}


    //! @name Methods
    //@{

    //! Initialize coupling variables for the first time step.
    void initializeCouplingVariables();

    Real computeResidual() const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the epetra communicator for the model
    /*!
     * @param comm Epetra communicator
     */
    void setCommunicator( const MS_Comm_PtrType& comm );

    //! Set the main MultiScale model
    /*!
     * @param model MultiScale model
     */
    void setModel( const MS_Model_PtrType model );

    //@}


    //! @name Get Methods
    //@{

    //! Get the type of the algorithm
    /*!
     * @return type of the algorithm
     */
    const algorithmsTypes& type() const;

    //! Get the Multiscale problem
    /*!
     * @return shared_ptr to the MultiScale problem
     */
    const multiscaleModelPtr_Type multiScaleProblem() const;

    //! Get the coupling variables
    /*!
     * @return pointer to the coupling variables vector
     */
    const MS_Vector_PtrType couplingVariables() const;

    //! Get the coupling residuals
    /*!
     * @return pointer to the coupling residuals vector
     */
    const MS_Vector_PtrType couplingResiduals() const;

    //! Get the communicator
    /*!
     * @return pointer to the communicator
     */
    const MS_Comm_PtrType communicator() const;

    //! Get the subiterations maximum number
    /*!
     * @return maximum number of subiterations
     */
    const UInt& subiterationsMaximumNumber() const;

    //! Get the required tolerance
    /*!
     * @return tolerance
     */
    const Real& tolerance() const;

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! save on a Matlab file the information about the convergence of the algorithm.
    /*!
     * @param subiterationsNumber Number of subiterations performed.
     * @param computeResidual computeResidual.
     */
    void save( const UInt& subiterationsNumber, const Real& residual );

    //! Check if the tolerance has been satisfied
    /*!
     * @return true if the tolerance is satisfied
     */
    bool toleranceSatisfied();

    //@}

    algorithmsTypes                          M_type;

    multiscaleModelPtr_Type                  M_multiscale;

    MS_Vector_PtrType                        M_couplingVariables;
    MS_Vector_PtrType                        M_couplingResiduals;

    MS_Comm_PtrType                          M_comm;
    boost::shared_ptr< Displayer >           M_displayer;

    UInt                                     M_subiterationsMaximumNumber;
    Real                                     M_tolerance;
};

} // Namespace LifeV

#endif /* MS_Algorithm_H */
