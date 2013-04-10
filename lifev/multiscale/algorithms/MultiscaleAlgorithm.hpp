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
 *  @brief File containing the Multiscale Algorithm
 *
 *  @date 23-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleAlgorithm_H
#define MultiscaleAlgorithm_H 1

#include <lifev/multiscale/framework/MultiscaleDefinitions.hpp>
#include <lifev/multiscale/models/MultiscaleModelMultiscale.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Teuchos_XMLParameterListHelpers.hpp>

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace LifeV
{
namespace Multiscale
{

// Forward declaration
class MultiscaleModelMultiscale;

//! MultiscaleAlgorithm - The Multiscale Algorithm Interface
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleAlgorithm class provides a general interface between the
 *  MultiscaleSolver and the specific Algorithm to solve the problem.
 *
 */
class MultiscaleAlgorithm
{
public:

    //! @name Type definitions
    //@{

    typedef MultiscaleModelMultiscale                                        multiscaleModelMultiscale_Type;
    typedef multiscaleModelMultiscale_Type*                                  multiscaleModelMultiscalePtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleAlgorithm();

    //! Destructor
    virtual ~MultiscaleAlgorithm() {}

    //@}


    //! @name Multiscale Algorithm Virtual Methods
    //@{

    //! Setup the data of the algorithm using a data file
    /*!
     * @param fileName Name of the data file.
     */
    virtual void setupData ( const std::string& fileName ) = 0;

    //! Setup coupling variables and other quantities of the algorithm
    virtual void setupAlgorithm();

    //! Perform sub-iteration on the coupling variables
    virtual void subIterate();

    //! Display some information about the algorithm
    virtual void showMe();

    //@}


    //! @name Methods
    //@{

    Real computeResidual() const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the epetra communicator for the model
    /*!
     * @param comm Epetra communicator
     */
    void setCommunicator ( const multiscaleCommPtr_Type& comm )
    {
        M_comm = comm;
    }

    //! Set the main Multiscale model
    /*!
     * @param model Multiscale model
     */
    void setMultiscaleModel ( const multiscaleModelMultiscalePtr_Type model )
    {
        M_multiscale = model;
    }

    //! Set the maximum number of subiterations
    /*!
     * @param subiterationsMaximumNumber maximum number of subiterations
     */
    void setSubiterationsMaximumNumber ( const UInt& subiterationsMaximumNumber )
    {
        M_subiterationsMaximumNumber = subiterationsMaximumNumber;
    }

    //! Set the tolerance
    /*!
     * @param tolerance coupling tolerance
     */
    void setTolerance ( const Real& tolerance )
    {
        M_tolerance = tolerance;
    }

    //! Set the algorithm name
    /*!
     * @param parameterList teuchos list of parameters
     */
    void setAlgorithmName ( const multiscaleParameterList_Type& parameterList );

    //! Set the the main parameters of the algorithm (tolerance, maximum number of subiterations, etc.)
    /*!
     * @param parameterList teuchos list of parameters
     */
    virtual void setAlgorithmParameters ( const multiscaleParameterList_Type& parameterList );

    //@}


    //! @name Get Methods
    //@{

    //! Get the type of the algorithm
    /*!
     * @return type of the algorithm
     */
    const algorithms_Type& type() const
    {
        return M_type;
    }

    //! Get the Multiscale problem
    /*!
     * @return shared_ptr to the Multiscale problem
     */
    const multiscaleModelMultiscalePtr_Type& multiScaleProblem() const
    {
        return M_multiscale;
    }

    //! Get the coupling variables
    /*!
     * @return pointer to the coupling variables vector
     */
    const multiscaleVectorPtr_Type& couplingVariables() const
    {
        return M_couplingVariables;
    }

    //! Get the coupling residuals
    /*!
     * @return pointer to the coupling residuals vector
     */
    const multiscaleVectorPtr_Type& couplingResiduals() const
    {
        return M_couplingResiduals;
    }

    //! Get the communicator
    /*!
     * @return pointer to the communicator
     */
    const multiscaleCommPtr_Type& communicator() const
    {
        return M_comm;
    }

    //! Get the subiterations maximum number
    /*!
     * @return maximum number of subiterations
     */
    const UInt& subiterationsMaximumNumber() const
    {
        return M_subiterationsMaximumNumber;
    }

    //! Get the required tolerance
    /*!
     * @return tolerance
     */
    const Real& tolerance() const
    {
        return M_tolerance;
    }

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! save on a Matlab file the information about the convergence of the algorithm.
    /*!
     * @param subiterationsNumber Number of subiterations performed.
     * @param computeResidual computeResidual.
     */
    void save ( const UInt& subiterationsNumber, const Real& residual ) const;

    //! Update the residual and check if the tolerance has been satisfied
    /*!
     * @param subIT subiteration number (for output purpose)
     * @return true if the tolerance is satisfied
     */
    bool checkResidual ( const UInt& subIT = 0 ) const;

    //@}

    algorithms_Type                          M_type;
    std::string                              M_name;

    multiscaleModelMultiscalePtr_Type        M_multiscale;

    multiscaleVectorPtr_Type                 M_couplingVariables;
    multiscaleVectorPtr_Type                 M_couplingResiduals;

    multiscaleCommPtr_Type                   M_comm;

    UInt                                     M_subiterationsMaximumNumber;
    Real                                     M_tolerance;

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleAlgorithm ( const MultiscaleAlgorithm& algorithm );

    MultiscaleAlgorithm& operator= ( const MultiscaleAlgorithm& algorithm );

    //@}
};

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleAlgorithm_H */
