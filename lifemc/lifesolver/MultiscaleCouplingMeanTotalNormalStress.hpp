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
 *  @brief File containing the Multiscale Coupling Stress
 *
 *  @date 20-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleCouplingMeanTotalNormalStress_H
#define MultiscaleCouplingMeanTotalNormalStress_H 1

#include <lifemc/lifesolver/MultiscaleCoupling.hpp>
#include <lifemc/lifesolver/MultiscaleModelFluid3D.hpp>
#include <lifemc/lifesolver/MultiscaleModelFSI3D.hpp>
#include <lifemc/lifesolver/MultiscaleModel1D.hpp>
#include <lifemc/lifesolver/MultiscaleModelWindkessel0D.hpp>
#include <lifemc/lifesolver/MultiscaleModel0D.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleCouplingMeanTotalNormalStress - Stress coupling condition
/*!
 *  @author Cristiano Malossi
 *
 *  The MultiscaleCouplingMeanTotalNormalStress class is an implementation of the multiscaleCoupling_Type
 *  for applying Stress coupling conditions on the models.
 */
class MultiscaleCouplingMeanTotalNormalStress: public virtual multiscaleCoupling_Type
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleCouplingMeanTotalNormalStress();

    //! Destructor
    virtual ~MultiscaleCouplingMeanTotalNormalStress() {}

    //@}


    //! @name Multiscale PhysicalCoupling Implementation
    //@{

    //! Setup the coupling
    void setupCoupling();

    //! Initialize the values of the coupling variables
    void initializeCouplingVariables();

    //! Update the coupling
    /*!
     * Nothing to do for this coupling.
     */
    void updateCoupling() {};

    //! Update the values of the coupling residuals
    /*!
     * @param couplingResiduals Global vector of variables
     */
    void exportCouplingResiduals( multiscaleVector_Type& couplingResiduals );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleCouplingMeanTotalNormalStress( const MultiscaleCouplingMeanTotalNormalStress& coupling );

    MultiscaleCouplingMeanTotalNormalStress& operator=( const MultiscaleCouplingMeanTotalNormalStress& coupling );

    //@}


    //! @name Private Multiscale PhysicalCoupling Implementation
    //@{

    //! Build the list of models affected by the perturbation of a local coupling variable
    /*!
     * @param localCouplingVariableID local coupling variable (perturbed)
     * @return list of models affected by the perturbation
     */
    multiscaleModelsContainer_Type listOfPerturbedModels( const UInt& localCouplingVariableID );

    //! Insert constant coefficients into the Jacobian matrix
    /*!
     * @param jacobian the Jacobian matrix
     */
    void insertJacobianConstantCoefficients( multiscaleMatrix_Type& jacobian );

    //! Insert the Jacobian coefficient(s) depending on a perturbation of the model, due to a specific variable (the column)
    /*!
     * @param jacobian the Jacobian matrix
     * @param column the column related to the perturbed variable
     * @param ID the global ID of the model which is perturbed by the variable
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     */
    void insertJacobianDeltaCoefficients( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem );

    //@}
};

//! Factory create function
inline multiscaleCoupling_Type* createMultiscaleCouplingMeanTotalNormalStress()
{
    return new MultiscaleCouplingMeanTotalNormalStress();
}

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleCouplingMeanTotalNormalStress_H */
