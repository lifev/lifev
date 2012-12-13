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

#include <lifev/multiscale/solver/MultiscaleCoupling.hpp>
#include <lifev/multiscale/solver/MultiscaleModelFluid3D.hpp>
#include <lifev/multiscale/solver/MultiscaleModelFSI3D.hpp>
#include <lifev/multiscale/solver/MultiscaleModelFSI1D.hpp>
#include <lifev/multiscale/solver/MultiscaleModelWindkessel0D.hpp>
#include <lifev/multiscale/solver/MultiscaleModel0D.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleCouplingMeanTotalNormalStress - Stress coupling condition
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
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

    //! Setup the coupling variables number.
    virtual void setupCouplingVariablesNumber();

    //! Setup the coupling
    virtual void setupCoupling();

    //! Initialize the values of the coupling variables
    virtual void initializeCouplingVariables();

    //! Update the coupling
    /*!
     * Nothing to do for this coupling.
     */
    virtual void updateCoupling() {};

    //! Compute the local coupling residuals vector
    virtual void computeCouplingResiduals();

    //@}

protected:

    //! @name Private Multiscale PhysicalCoupling Implementation
    //@{

    //! Build the list of models affected by the perturbation of a local coupling variable
    /*!
     * @param localCouplingVariableID id of the perturbed local coupling variable
     * @param perturbedModelsList list of models affected by the perturbation
     */
    virtual void exportListOfPerturbedModels( const UInt& localCouplingVariableID, multiscaleModelsContainer_Type& perturbedModelsList );

    //! Insert constant coefficients into the Jacobian matrix
    /*!
     * @param jacobian the Jacobian matrix
     */
    virtual void insertJacobianConstantCoefficients( multiscaleMatrix_Type& jacobian );

    //! Insert the Jacobian coefficient(s) depending on a perturbation of the model, due to a specific variable (the column)
    /*!
     * @param jacobian the Jacobian matrix
     * @param column the column related to the perturbed variable
     * @param ID the global ID of the model which is perturbed by the variable
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     */
    virtual void insertJacobianDeltaCoefficients( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleCouplingMeanTotalNormalStress( const MultiscaleCouplingMeanTotalNormalStress& coupling );

    MultiscaleCouplingMeanTotalNormalStress& operator=( const MultiscaleCouplingMeanTotalNormalStress& coupling );

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
