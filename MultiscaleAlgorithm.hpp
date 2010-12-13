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

#ifndef MultiscaleAlgorithm_H
#define MultiscaleAlgorithm_H 1

#include <lifemc/lifesolver/MultiscaleDefinitions.hpp>
#include <lifemc/lifesolver/MultiscaleModelMultiscale.hpp>

namespace LifeV
{
namespace multiscale
{

//! MultiscaleAlgorithm - The MultiScale Algorithm Interface
/*!
 *  @author Cristiano Malossi
 *
 *  The MultiscaleAlgorithm class provides a general interface between the
 *  MS_Solver and the specific Algorithm to solve the problem.
 *
 */
class MultiscaleAlgorithm
{
public:

    //! @name Type definitions
    //@{

    typedef boost::shared_ptr< MultiscaleModelMultiscale >                   multiscaleModelMultiscalePtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleAlgorithm();

    //! Destructor
    virtual ~MultiscaleAlgorithm() {}

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
    virtual void updateCouplingVariables() { M_multiscale->extrapolateCouplingVariables(); }

    //! Display some information about the algorithm
    virtual void showMe();

    //@}


    //! @name Methods
    //@{

    //! Initialize coupling variables for the first time step.
    void initializeCouplingVariables() { M_multiscale->initializeCouplingVariables(); }

    Real computeResidual() const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the epetra communicator for the model
    /*!
     * @param comm Epetra communicator
     */
    void setCommunicator( const multiscaleCommPtr_Type& comm );

    //! Set the main MultiScale model
    /*!
     * @param model MultiScale model
     */
    void setModel( const multiscaleModelPtr_Type model );

    //@}


    //! @name Get Methods
    //@{

    //! Get the type of the algorithm
    /*!
     * @return type of the algorithm
     */
    const algorithms_Type& type() const { return M_type; }

    //! Get the Multiscale problem
    /*!
     * @return shared_ptr to the MultiScale problem
     */
    const multiscaleModelMultiscalePtr_Type multiScaleProblem() const { return M_multiscale; }

    //! Get the coupling variables
    /*!
     * @return pointer to the coupling variables vector
     */
    const multiscaleVectorPtr_Type couplingVariables() const { return M_couplingVariables; }

    //! Get the coupling residuals
    /*!
     * @return pointer to the coupling residuals vector
     */
    const multiscaleVectorPtr_Type couplingResiduals() const { return M_couplingResiduals; }

    //! Get the communicator
    /*!
     * @return pointer to the communicator
     */
    const multiscaleCommPtr_Type communicator() const { return M_comm; }

    //! Get the subiterations maximum number
    /*!
     * @return maximum number of subiterations
     */
    const UInt& subiterationsMaximumNumber() const { return M_subiterationsMaximumNumber; }

    //! Get the required tolerance
    /*!
     * @return tolerance
     */
    const Real& tolerance() const { return M_tolerance; }

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

    algorithms_Type                          M_type;

    multiscaleModelMultiscalePtr_Type        M_multiscale;

    multiscaleVectorPtr_Type                 M_couplingVariables;
    multiscaleVectorPtr_Type                 M_couplingResiduals;

    multiscaleCommPtr_Type                   M_comm;
    boost::shared_ptr< Displayer >           M_displayer;

    UInt                                     M_subiterationsMaximumNumber;
    Real                                     M_tolerance;

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleAlgorithm( const MultiscaleAlgorithm& algorithm );

    MultiscaleAlgorithm& operator=( const MultiscaleAlgorithm& algorithm );

    //@}
};

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleAlgorithm_H */
