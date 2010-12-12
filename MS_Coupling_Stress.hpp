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
 *  @brief File containing the MultiScale Coupling Stress
 *
 *  @date 20-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MS_Coupling_Stress_H
#define MS_Coupling_Stress_H 1

#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>
#include <lifemc/lifesolver/MS_Model_Fluid3D.hpp>
#include <lifemc/lifesolver/MS_Model_FSI3D.hpp>
#include <lifemc/lifesolver/MS_Model_1D.hpp>

namespace LifeV
{

//! MS_Coupling_Stress - Stress coupling condition
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_Coupling_Stress class is an implementation of the MS_PhysicalCoupling
 *  for applying Stress coupling conditions on the models.
 */
class MS_Coupling_Stress: public virtual MS_PhysicalCoupling
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Coupling_Stress();

    //! Destructor
    virtual ~MS_Coupling_Stress() {}

    //@}


    //! @name MultiScale PhysicalCoupling Implementation
    //@{

    //! Setup the data of the coupling
    /*!
     *  @param FileName Name of data file
     */
    void setupData( const std::string& fileName );

    //! Setup the coupling
    void setupCoupling();

    //! Initialize the values of the coupling variables
    void initializeCouplingVariables();

    //! Update the values of the coupling residuals
    void exportCouplingResiduals( MS_Vector_Type& couplingResiduals );

    //! Display some information about the coupling
    void showMe();

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MS_Coupling_Stress( const MS_Coupling_Stress& coupling );

    MS_Coupling_Stress& operator=( const MS_Coupling_Stress& coupling );

    //@}


    //! @name Private MultiScale PhysicalCoupling Implementation
    //@{

    //! Build the list of models affected by the perturbation of a local coupling variable
    /*!
     * @param LocalCouplingVariableID local coupling variable (perturbed)
     * @return list of models affected by the perturbation
     */
    MS_ModelsVector_Type listOfPerturbedModels( const UInt& localCouplingVariableID );

    //! Insert constant coefficients into the Jacobian matrix
    /*!
     * @param jacobian the Jacobian matrix
     */
    void insertJacobianConstantCoefficients( MS_Matrix_Type& jacobian );

    //! Insert the Jacobian coefficient(s) depending on a perturbation of the model, due to a specific variable (the column)
    /*!
     * @param jacobian the Jacobian matrix
     * @param column the column related to the perturbed variable
     * @param ID the global ID of the model which is perturbed by the variable
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     */
    void insertJacobianDeltaCoefficients( MS_Matrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem );

    //! Display some information about the coupling
    /*!
     * @param output specify the output stream
     */
    void displayCouplingValues( std::ostream& output );

    //@}


    //! @name Private Methods
    //@{

    template< class ModelType >
    void imposeStress3D( const UInt& i );

    template< class ModelType >
    void imposeStress1D( const UInt& i );

    Real functionStress( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const UInt& /*id*/);

    //@}

    BCFunctionBase M_baseStress3D;

    OneDimensionalModel_BCFunction M_baseStress1D;

    stress_Type    M_stressType;
};

//! Factory create function
inline MS_PhysicalCoupling* createMultiscaleCouplingStress()
{
    return new MS_Coupling_Stress();
}

// ===================================================
// Template implementation
// ===================================================
template< class ModelType >
inline void
MS_Coupling_Stress::imposeStress3D( const UInt& i )
{
    ModelType *model = MS_DynamicCast< ModelType >( M_models[i] );

    model->bcInterface().addBC( "CouplingStress_Model_" + number2string( model->ID() ) + "_Flag_" + number2string( M_flags[i] ),
                                   M_flags[i], Natural, Normal, M_baseStress3D );
}

template< class ModelType >
inline void
MS_Coupling_Stress::imposeStress1D( const UInt& i )
{
    ModelType *model = MS_DynamicCast< ModelType >( M_models[i] );

    model->bcInterface().setBC( (M_flags[i] == 0) ? OneD_left : OneD_right, OneD_first, OneD_P, M_baseStress1D );
}

} // Namespace LifeV

#endif /* MS_Coupling_Stress_H */
