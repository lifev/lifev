//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief MultiScale Coupling Stress
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 20-10-2009
 */

#ifndef MS_Coupling_Stress_H
#define MS_Coupling_Stress_H 1

#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>
#include <lifemc/lifesolver/MS_Model_Fluid3D.hpp>
#include <lifemc/lifesolver/MS_Model_FSI3D.hpp>
#include <lifemc/lifesolver/MS_Model_1D.hpp>

namespace LifeV {

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

    typedef MS_PhysicalCoupling super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Coupling_Stress();

    //! Copy constructor
    /*!
     * @param Stress MS_Coupling_Stress
     */
    MS_Coupling_Stress( const MS_Coupling_Stress& Stress );

    //! Destructor
    ~MS_Coupling_Stress() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param Stress MS_Coupling_Stress
     * @return reference to a copy of the class
     */
    MS_Coupling_Stress& operator=( const MS_Coupling_Stress& Stress );

    //@}


    //! @name MultiScale PhysicalCoupling Implementation
    //@{

    //! Setup the data of the coupling
    /*!
     *  @param FileName Name of data file
     */
    void SetupData( const std::string& FileName );

    //! Setup the coupling
    void SetupCoupling();

    //! Initialize the values of the coupling variables
    void InitializeCouplingVariables();

    //! Update the values of the coupling residuals
    void ExportCouplingResiduals( MS_Vector_Type& CouplingResiduals );

    //! Display some information about the coupling
    void ShowMe();

    //@}

private:

    //! @name Private MultiScale PhysicalCoupling Implementation
    //@{

    //! Build the list of models affected by the perturbation of a local coupling variable
    /*!
     * @param LocalCouplingVariableID local coupling variable (perturbed)
     * @return list of models affected by the perturbation
     */
    MS_ModelsVector_Type GetListOfPerturbedModels( const UInt& LocalCouplingVariableID );

    //! Insert constant coefficients into the Jacobian matrix
    /*!
     * @param Jacobian the Jacobian matrix
     */
    void InsertJacobianConstantCoefficients( MS_Matrix_Type& Jacobian );

    //! Insert the Jacobian coefficient(s) depending on a perturbation of the model, due to a specific variable (the column)
    /*!
     * @param Jacobian the Jacobian matrix
     * @param Column the column related to the perturbed variable
     * @param ID the global ID of the model which is perturbed by the variable
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     */
    void InsertJacobianDeltaCoefficients( MS_Matrix_Type& Jacobian, const UInt& Column, const UInt& ID, bool& SolveLinearSystem );

    //! Display some information about the coupling
    /*!
     * @param output specify the output stream
     */
    void DisplayCouplingValues( std::ostream& output );

    //@}


    //! @name Private Methods
    //@{

    template< class model >
    inline void ImposeStress3D( const UInt& i );

    template< class model >
    inline void ImposeStress1D( const UInt& i );

    Real FunctionStress( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/);

    //@}

    BCFunctionBase M_baseStress3D;

    OneDimensionalModel_BCFunction M_baseStress1D;

    stressTypes    M_stressType;
};

//! Factory create function
inline MS_PhysicalCoupling* MS_createStress()
{
    return new MS_Coupling_Stress();
}

// ===================================================
// Template implementation
// ===================================================
template< class model >
inline void
MS_Coupling_Stress::ImposeStress3D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    Model->GetBCInterface().addBC( "CouplingStress_Model_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[i] ),
                                   M_flags[i], Natural, Normal, M_baseStress3D );
}

template< class model >
inline void
MS_Coupling_Stress::ImposeStress1D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    Model->GetBCInterface().setBC( (M_flags[i] == 0) ? OneD_left : OneD_right, OneD_first, OneD_P, M_baseStress1D );
}

} // Namespace LifeV

#endif /* MS_Coupling_Stress_H */
