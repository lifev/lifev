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
 *  @brief File containing the Multiscale Model Multiscale
 *
 *  @date 12-03-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleModelMultiscale_H
#define MultiscaleModelMultiscale_H 1

#include <lifev/multiscale/framework/MultiscaleCommunicatorsManager.hpp>

#include <lifev/multiscale/algorithms/MultiscaleAlgorithm.hpp>

#include <lifev/multiscale/couplings/MultiscaleCoupling.hpp>

#include <lifev/multiscale/models/MultiscaleModel.hpp>

#if defined(LIFEV_HAS_ZERODIMENSIONAL)
#include <lifev/multiscale/models/MultiscaleModelWindkessel0D.hpp>
#include <lifev/multiscale/models/MultiscaleModel0D.hpp>
#endif

#if defined(LIFEV_HAS_ONEDFSI)
#include <lifev/multiscale/models/MultiscaleModelFSI1D.hpp>
#endif

#if defined(LIFEV_HAS_NAVIERSTOKES)
#include <lifev/multiscale/models/MultiscaleModelFluid3D.hpp>
#endif

#if defined(LIFEV_HAS_FSI)
#include <lifev/multiscale/models/MultiscaleModelFSI3D.hpp>
#endif

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleModelMultiscale - Multiscale model
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleModelMultiscale class is an implementation of the MultiscaleModel
 *  for a general multiscale problem.
 */
class MultiscaleModelMultiscale: public virtual multiscaleModel_Type
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleModelMultiscale();

    //! Destructor
    virtual ~MultiscaleModelMultiscale();

    //@}


    //! @name MultiscaleModel Methods
    //@{

    //! Setup the data of the model.
    /*!
     * @param fileName Name of data file.
     */
    void setupData ( const std::string& fileName );

    //! Setup the model.
    void setupModel();

    //! Build the initial model.
    void buildModel();

    //! Update the model.
    void updateModel();

    //! Solve the model.
    void solveModel();

    //! Update the solution.
    void updateSolution();

    //! Save the solution
    void saveSolution();

    //! Display some information about the multiscale problem (call after SetupProblem)
    void showMe();

    //! Return a specific scalar quantity to be used for a comparison with a reference value.
    /*!
     * This method is meant to be used for night checks.
     * @return reference quantity.
     */
    Real checkSolution() const;

    //@}


    //! @name Methods
    //@{

    //! Build the global map for the coupling vectors
    /*!
     * @param couplingMap Global coupling map
     */
    void createCouplingMap ( MapEpetra& couplingMap );

    //! Import the values of the coupling variables
    void importCouplingVariables ( const multiscaleVector_Type& couplingVariables );

    //! Export the values of the coupling variables
    void exportCouplingVariables ( multiscaleVector_Type& couplingVariables );

    //! Compute the values of the interface residuals
    void computeCouplingResiduals();

    //! Export the values of the interface residuals
    void exportCouplingResiduals ( multiscaleVector_Type& couplingResiduals );

    //! Export the Jacobian matrix
    /*!
     * @param jacobian Matrix
     */
    void exportJacobian ( multiscaleMatrix_Type& jacobian );

    //! Check if the topology is changed
    /*!
     * A topology change can be caused by a change in the coupling equations by,
     * for example, the opening/closure of a valve (see MultiscaleCouplingFlowRateValve).
     * @return true if the topology is changed, false otherwise
     */
    bool topologyChange();

    //@}


    //! @name Get Methods
    //@{

    //! Get the number of the coupling variables
    /*!
     * @return number of the coupling variables
     */
    UInt couplingVariablesNumber();

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleModelMultiscale ( const MultiscaleModelMultiscale& model );

    MultiscaleModelMultiscale& operator= ( const MultiscaleModelMultiscale& model );

    //@}

    // Models & Couplings
    MultiscaleCommunicatorsManager     M_commManager;

    // Models & Couplings
    multiscaleModelsContainer_Type     M_modelsList;
    multiscaleCouplingsContainer_Type  M_couplingsList;

    // Algorithm for subiterations
    multiscaleAlgorithmPtr_Type        M_algorithm;
};

//! Factory create function
inline multiscaleModel_Type* createMultiscaleModelMultiscale()
{
    return new MultiscaleModelMultiscale();
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleModelMultiscale_H */
