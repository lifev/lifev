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

#include <lifemc/lifesolver/MultiscaleModel.hpp>
#include <lifemc/lifesolver/MultiscaleModelFluid3D.hpp>
#include <lifemc/lifesolver/MultiscaleModelFSI3D.hpp>
#include <lifemc/lifesolver/MultiscaleModel1D.hpp>
#include <lifemc/lifesolver/MultiscaleModelWindkessel0D.hpp>

#include <lifemc/lifesolver/MultiscaleCoupling.hpp>
#include <lifemc/lifesolver/MultiscaleCouplingBoundaryCondition.hpp>
//#include <lifemc/lifesolver/MultiscaleCouplingFlowRate.hpp>
#include <lifemc/lifesolver/MultiscaleCouplingFlowRateStress.hpp>
#include <lifemc/lifesolver/MultiscaleCouplingStress.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleModelMultiscale - Multiscale model
/*!
 *  @author Cristiano Malossi
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
    void setupData( const std::string& fileName );

    //! Setup the model.
    void setupModel();

    //! Build the initial model.
    void buildModel();

    //! Update the model.
    void updateModel();

    //! Solve the model.
    void solveModel();

    //! Save the solution
    void saveSolution();

    //! Display some information about the multiscale problem (call after SetupProblem)
    void showMe();

    //@}


    //! @name Methods
    //@{

    //! Build the global map for the coupling vectors
    /*!
     * @param couplingMap Global coupling map
     */
    void createCouplingMap( MapEpetra& couplingMap );

    //! Initialize coupling variables for the first time step
    void initializeCouplingVariables();

    //! Extrapolate coupling variables for the next time step
    void extrapolateCouplingVariables();

    //! Import the values of the coupling variables
    void importCouplingVariables( const multiscaleVector_Type& couplingVariables );

    //! Export the values of the coupling variables
    void exportCouplingVariables( multiscaleVector_Type& couplingVariables );

    //! Export the values of the coupling residuals
    void exportCouplingResiduals( multiscaleVector_Type& couplingResiduals );

    //! Export the Jacobian matrix
    /*!
     * @param jacobian Matrix
     */
    void exportJacobian( multiscaleMatrix_Type& jacobian );

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

    MultiscaleModelMultiscale( const MultiscaleModelMultiscale& model );

    MultiscaleModelMultiscale& operator=( const MultiscaleModelMultiscale& model );

    //@}

    // Models & Couplings
    multiscaleModelsVector_Type        M_modelsList;
    multiscaleCouplingsVector_Type     M_couplingsList;
};

//! Factory create function
inline multiscaleModel_Type* createMultiscaleModelMultiscale()
{
    return new MultiscaleModelMultiscale();
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleModelMultiscale_H */
