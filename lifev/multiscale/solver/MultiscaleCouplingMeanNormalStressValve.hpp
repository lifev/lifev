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
 *  @brief File containing the multiscale mean normal stress coupling class with simple valve
 *
 *  @date 05-04-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleCouplingMeanNormalStressValve_H
#define MultiscaleCouplingMeanNormalStressValve_H 1

#include <lifev/multiscale/solver/MultiscaleCouplingMeanNormalStress.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleCouplingMeanNormalStressValve - Mean normal stress coupling condition with simple valve
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleCouplingMeanNormalStressValve class is an implementation of the multiscaleCoupling_Type
 *  for applying flow rate coupling conditions on the models including a valve between the first model and the others.
 *
 *  NOTE: The current implementation works only for two models!
 */
class MultiscaleCouplingMeanNormalStressValve: public virtual MultiscaleCouplingMeanNormalStress
{
public:

    //! @name Type definitions
    //@{

    typedef MultiscaleCouplingMeanNormalStress                  super_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleCouplingMeanNormalStressValve();

    //! Destructor
    virtual ~MultiscaleCouplingMeanNormalStressValve() {}

    //@}


    //! @name Multiscale PhysicalCoupling Implementation
    //@{

    //! Setup the coupling
    void setupCoupling();

    //! Initialize the values of the coupling variables
    void initializeCouplingVariables();

    //! Update the coupling
    /*!
     * Update the position of the valve.
     */
    void updateCoupling();

    //! Compute the local coupling residuals vector
    void computeCouplingResiduals();

    //! Check if the topology is changed
    /*!
     * The opening/closure of the valve change the topology.
     *
     * @return true if the topology is changed, false otherwise
     */
    bool topologyChange()
    {
        return M_topologyChange;
    }

    //@}

private:

    //! @name Private Multiscale PhysicalCoupling Implementation
    //@{

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


    //! @name Unimplemented Methods
    //@{

    MultiscaleCouplingMeanNormalStressValve ( const MultiscaleCouplingMeanNormalStressValve& coupling );

    MultiscaleCouplingMeanNormalStressValve& operator= ( const MultiscaleCouplingMeanNormalStressValve& coupling );

    //@}

    bool                     M_valveIsOpen;
    bool                     M_topologyChange;
};

//! Factory create function
inline multiscaleCoupling_Type* createMultiscaleCouplingMeanNormalStressValve()
{
    return new MultiscaleCouplingMeanNormalStressValve();
}

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleCouplingMeanNormalStressValve_H */
