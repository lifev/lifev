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
 *  @brief MultiScale Coupling FluxStress
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 24-08-2009
 */

#ifndef MS_Coupling_FluxStress_H
#define MS_Coupling_FluxStress_H 1

#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>
#include <lifemc/lifesolver/MS_Model_Fluid3D.hpp>
#include <lifemc/lifesolver/MS_Model_FSI3D.hpp>
#include <lifemc/lifesolver/MS_Model_1D.hpp>

namespace LifeV {

//! MS_Coupling_FluxStress - Flux-Stress coupling condition
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_Coupling_FluxStress class is an implementation of the MS_PhysicalCoupling
 *  for applying Flux-Stress coupling conditions on the models.
 *
 *  The coupling equations are:
 *  Q_j = - \sum Q_i
 *  \sigma_i = -P_j
 *  where Q is the flux and P is the pressure (or the total pressure).
 */
class MS_Coupling_FluxStress: public virtual MS_PhysicalCoupling
{
public:

    typedef MS_PhysicalCoupling super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Coupling_FluxStress();

    //! Copy constructor
    /*!
     * @param FluxStress MS_Coupling_FluxStress
     */
    MS_Coupling_FluxStress( const MS_Coupling_FluxStress& FluxStress );

    //! Destructor
    ~MS_Coupling_FluxStress() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param FluxStress MS_Coupling_FluxStress
     * @return reference to a copy of the class
     */
    MS_Coupling_FluxStress& operator=( const MS_Coupling_FluxStress& FluxStress );

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
    /*!
     * @param CouplingResiduals Global vector of variables
     */
    void ExportCouplingResiduals( MS_Vector_Type& CouplingResiduals );

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

    //! Display some information about the coupling
    void ShowMe();

    //@}

private:

    //! @name Private Methods
    //@{

    template< class model >
    inline void ImposeFlux3D( const UInt& i);

    template< class model >
    inline void ImposeStress3D( const UInt& i );

    Real FunctionFlux3D  ( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/);
    Real FunctionStress3D( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/);

    template< class model >
    inline void ImposeDeltaFlux3D( const UInt& i );

    template< class model >
    inline void ImposeDeltaStress3D( const UInt& i );

    Real FunctionDeltaFlux3D  ( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/);
    Real FunctionDeltaStress3D( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/);



    template< class model >
    inline void ImposeFlux1D( const UInt& i);

    template< class model >
    inline void ImposeStress1D( const UInt& i );

    Real FunctionFlux1D  ( const Real& /*t*/ );
    Real FunctionStress1D( const Real& /*t*/ );

    template< class model >
    inline void ImposeDeltaFlux1D( const UInt& i );

    template< class model >
    inline void ImposeDeltaStress1D( const UInt& i );

    Real FunctionDeltaFlux1D  ( const Real& /*t*/ );
    Real FunctionDeltaStress1D( const Real& /*t*/ );

    //@}

    BCFunctionBase M_baseFlux3D;
    BCFunctionBase M_baseStress3D;

    BCFunctionBase M_baseDeltaFlux3D;
    BCFunctionBase M_baseDeltaStress3D;

    OneDimensionalModel_BCFunction M_baseFlux1D;
    OneDimensionalModel_BCFunction M_baseStress1D;

    OneDimensionalModel_BCFunction M_baseDeltaFlux1D;
    OneDimensionalModel_BCFunction M_baseDeltaStress1D;

    stressTypes    M_stressType;
};

//! Factory create function
inline MS_PhysicalCoupling* MS_createFluxStress()
{
    return new MS_Coupling_FluxStress();
}

// ===================================================
// Template implementation
// ===================================================
template< class model >
inline void
MS_Coupling_FluxStress::ImposeFlux3D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    Model->GetBCInterface().addBC( "imposeFlux_Model_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[i] ), M_flags[i], Flux, Full, M_baseFlux3D, 3 );
}

template< class model >
inline void
MS_Coupling_FluxStress::ImposeStress3D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    Model->GetBCInterface().addBC( "imposeStress_Model_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[i] ), M_flags[i], Natural, Normal, M_baseStress3D );
}

template< class model >
inline void
MS_Coupling_FluxStress::ImposeDeltaFlux3D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    Model->GetLinearBCInterface().addBC( "imposeDeltaFlux_Model_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[i] ), M_flags[i], Flux, Full, M_baseDeltaFlux3D, 3 );
}

template< class model >
inline void
MS_Coupling_FluxStress::ImposeDeltaStress3D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    Model->GetLinearBCInterface().addBC( "imposeDeltaStress_Model_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[i] ), M_flags[i], Natural, Normal, M_baseDeltaStress3D );
}

template< class model >
inline void
MS_Coupling_FluxStress::ImposeFlux1D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    Model->GetBCInterface().setBC( (M_flags[i] == 0) ? OneD_left : OneD_right, OneD_first, OneD_Q, M_baseFlux1D );
}

template< class model >
inline void
MS_Coupling_FluxStress::ImposeStress1D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    Model->GetBCInterface().setBC( (M_flags[i] == 0) ? OneD_left : OneD_right, OneD_first, OneD_P, M_baseStress1D );
}

template< class model >
inline void
MS_Coupling_FluxStress::ImposeDeltaFlux1D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    Model->GetLinearBCInterface().setBC( (M_flags[i] == 0) ? OneD_left : OneD_right, OneD_first, OneD_Q, M_baseDeltaFlux1D );
}

template< class model >
inline void
MS_Coupling_FluxStress::ImposeDeltaStress1D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    Model->GetLinearBCInterface().setBC( (M_flags[i] == 0) ? OneD_left : OneD_right, OneD_first, OneD_P, M_baseDeltaStress1D );
}

} // Namespace LifeV

#endif /* MS_Coupling_FluxStress_H */
