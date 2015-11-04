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
 *  @brief File containing the multiscale mean normal stress coupling class
 *
 *  @date 11-10-2012
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleCouplingMeanNormalStressArea_H
#define MultiscaleCouplingMeanNormalStressArea_H 1

#include <lifev/multiscale/couplings/MultiscaleCouplingMeanNormalStress.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleCouplingMeanNormalStressArea - Mean normal stress with area coupling condition
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleCouplingMeanNormalStressArea class is an implementation of the multiscaleCoupling_Type
 *  for applying mean normal stress with area coupling conditions to the models interfaces.
 */
class MultiscaleCouplingMeanNormalStressArea: public virtual MultiscaleCouplingMeanNormalStress
{
public:

    //! @name Type definitions
    //@{

    typedef MultiscaleCouplingMeanNormalStress                  super_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleCouplingMeanNormalStressArea();

    //! Destructor
    virtual ~MultiscaleCouplingMeanNormalStressArea() {}

    //@}


    //! @name Multiscale PhysicalCoupling Implementation
    //@{

    //! Setup the coupling variables number.
    void setupCouplingVariablesNumber();

    //! Setup the coupling
    void setupCoupling();

    //! Initialize the values of the coupling variables
    void initializeCouplingVariables();

    //! Compute the local coupling residuals vector
    void computeCouplingResiduals();

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleCouplingMeanNormalStressArea ( const MultiscaleCouplingMeanNormalStressArea& coupling );

    MultiscaleCouplingMeanNormalStressArea& operator= ( const MultiscaleCouplingMeanNormalStressArea& coupling );

    //@}


    //! @name Private Multiscale PhysicalCoupling Implementation
    //@{

    //! Build the list of models affected by the perturbation of a local coupling variable
    /*!
     * @param localCouplingVariableID id of the perturbed local coupling variable
     * @param perturbedModelsList list of models affected by the perturbation
     */
    void exportListOfPerturbedModels ( const UInt& localCouplingVariableID, multiscaleModelsContainer_Type& perturbedModelsList );


    //! Insert constant coefficients into the Jacobian matrix
    /*!
     * @param jacobian the Jacobian matrix
     */
    void insertJacobianConstantCoefficients ( multiscaleMatrix_Type& jacobian );

    //! Insert the Jacobian coefficient(s) depending on a perturbation of the model, due to a specific variable (the column)
    /*!
     * @param jacobian the Jacobian matrix
     * @param column the column related to the perturbed variable
     * @param ID the global ID of the model which is perturbed by the variable
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     */
    void insertJacobianDeltaCoefficients ( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem );

    //@}

};

//! Factory create function
inline multiscaleCoupling_Type* createMultiscaleCouplingMeanNormalStressArea()
{
    return new MultiscaleCouplingMeanNormalStressArea();
}

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleCouplingMeanNormalStressArea_H */
