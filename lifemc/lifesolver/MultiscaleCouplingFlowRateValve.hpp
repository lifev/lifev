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
 *  @brief File containing the Multiscale Coupling FlowRateValve
 *
 *  @date 05-04-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleCouplingFlowRateValve_H
#define MultiscaleCouplingFlowRateValve_H 1

#include <lifemc/lifesolver/MultiscaleCouplingFlowRate.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleCouplingFlowRateValve - FlowRate coupling condition with a valve
/*!
 *  @author Cristiano Malossi
 *
 *  The MultiscaleCouplingFlowRateValve class is an implementation of the multiscaleCoupling_Type
 *  for applying flow rate coupling conditions on the models including a valve between the first model and the others.
 *
 *  NOTE: The current implementation works only for two models!
 */
class MultiscaleCouplingFlowRateValve: public virtual MultiscaleCouplingFlowRate
{
public:

    //! @name Type definitions
    //@{

    typedef MultiscaleCouplingFlowRate                          super_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleCouplingFlowRateValve();

    //! Destructor
    virtual ~MultiscaleCouplingFlowRateValve() {}

    //@}


    //! @name Multiscale PhysicalCoupling Implementation
    //@{

    //! Initialize the values of the coupling variables
    void initializeCouplingVariables();

    //! Setup the coupling
    void setupCoupling();

    //! Update the coupling
    /*!
     * Update the position of the valve.
     */
    void updateCoupling();

    //! Update the values of the coupling residuals
    /*!
     * @param couplingResiduals Global vector of variables
     */
    void exportCouplingResiduals( multiscaleVector_Type& couplingResiduals );

    //! Check if the topology is changed
    /*!
     * The opening/closure of the valve change the topology.
     *
     * @return true if the topology is changed, false otherwise
     */
    bool topologyChange() { return M_topologyChange; }

    //! Display some information about the coupling
    void showMe();

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleCouplingFlowRateValve( const MultiscaleCouplingFlowRateValve& coupling );

    MultiscaleCouplingFlowRateValve& operator=( const MultiscaleCouplingFlowRateValve& coupling );

    //@}


    //! @name Private Multiscale PhysicalCoupling Implementation
    //@{

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

    bool                     M_valveIsOpen;
    bool                     M_topologyChange;
};

//! Factory create function
inline multiscaleCoupling_Type* createMultiscaleCouplingFlowRateValve()
{
    return new MultiscaleCouplingFlowRateValve();
}

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleCouplingFlowRateValve_H */
