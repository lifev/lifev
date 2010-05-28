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
 *  @brief MultiScale Coupling BoundaryCondition
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 02-09-2009
 */

#ifndef MS_Coupling_BoundaryCondition_H
#define MS_Coupling_BoundaryCondition_H 1

#include <lifemc/lifesolver/BCInterface.hpp>
#include <lifemc/lifesolver/BCInterface1D.hpp>

#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>
#include <lifemc/lifesolver/MS_Model_Fluid3D.hpp>
#include <lifemc/lifesolver/MS_Model_FSI3D.hpp>
#include <lifemc/lifesolver/MS_Model_1D.hpp>

namespace LifeV {

//! MS_Coupling_BoundaryCondition - Coupling condition for standard boundary conditions
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_Coupling_BoundaryCondition class is an implementation of the MS_PhysicalCoupling
 *  for applying standard boundary conditions on the models.
 */
class MS_Coupling_BoundaryCondition: public virtual MS_PhysicalCoupling
{
public:

    typedef MS_PhysicalCoupling super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Coupling_BoundaryCondition();

    //! Copy constructor
    /*!
     * @param BoundaryCondition MS_Coupling_BoundaryCondition
     */
    MS_Coupling_BoundaryCondition( const MS_Coupling_BoundaryCondition& BoundaryCondition );

    //! Destructor
    ~MS_Coupling_BoundaryCondition() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param boundaryCondition MS_Coupling_BoundaryCondition
     * @return reference to a copy of the class
     */
    MS_Coupling_BoundaryCondition& operator=( const MS_Coupling_BoundaryCondition& boundaryCondition );

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

    //! Initialize the values of the coupling variables (DO NOTHING)
    void InitializeCouplingVariables() {}

    //! Export the values of the local coupling residuals into a global vector (DO NOTHING)
    /*!
     * @param CouplingResiduals Global vector of variables
     */
    void ExportCouplingResiduals( VectorType& /*CouplingResiduals*/ ) {}

    //! Build the list of models affected by the perturbation of a local coupling variable
    /*!
     * @param LocalCouplingVariableID local coupling variable (perturbed)
     * @return list of models affected by the perturbation
     */
    ModelsVector_Type GetListOfPerturbedModels( const UInt& /*LocalCouplingVariableID*/ );

    //! Insert constant coefficients into the Jacobian matrix (DO NOTHING)
    /*!
     * @param Jacobian the Jacobian matrix
     */
    void InsertJacobianConstantCoefficients( MatrixType& /*Jacobian*/ ) {}

    //! Insert the Jacobian coefficient(s) depending on a perturbation of the model, due to a specific variable (the column) (DO NOTHING)
    /*!
     * @param Jacobian          the Jacobian matrix
     * @param Column            the column related to the perturbed variable
     * @param ID                the global ID of the model which is perturbed by the variable
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     */
    void InsertJacobianDeltaCoefficients( MatrixType& /*Jacobian*/, const UInt& /*Column*/, const UInt& /*ID*/, bool& /*SolveLinearSystem*/ ) {}

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

    //! Apply the boundary condition to the specific 1D model
    template< class model >
    inline void ApplyBoundaryConditions1D( const UInt& i );

    //! Apply the boundary condition to linear version of the specific 1D model
    template< class model >
    inline void ApplyDeltaBoundaryConditions1D( const UInt& i );

    //! Apply the boundary condition to the specific 3D model
    template< class model >
    inline void ApplyBoundaryConditions3D( const UInt& i );

    //! Apply the boundary condition to linear version of the specific 3D model
    template< class model >
    inline void ApplyDeltaBoundaryConditions3D( const UInt& i );

    //@}

    std::string           M_FileName;

    std::vector< BCName > M_list;
    UInt                  M_listSize;
};

//! Factory create function
inline MS_PhysicalCoupling* createBoundaryCondition()
{
    return new MS_Coupling_BoundaryCondition();
}

// ===================================================
// Template implementation
// ===================================================
template< class model >
inline void
MS_Coupling_BoundaryCondition::ApplyBoundaryConditions1D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    for ( UInt j( 0 ); j < M_listSize; ++j )
    {
        Model->GetBCInterface().ReadBC( M_FileName, "boundary_conditions/", M_list[j] );

        Model->GetBCInterface().GetDataContainer().SetSide( (M_flags[i] == 0) ? OneD_left : OneD_right );

        Model->GetBCInterface().InsertBC();
    }
}

template< class model >
inline void
MS_Coupling_BoundaryCondition::ApplyDeltaBoundaryConditions1D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    for ( UInt j( 0 ); j < M_listSize; ++j )
    {
        Model->GetLinearBCInterface().ReadBC( M_FileName, "boundary_conditions/", M_list[j] );

        Model->GetLinearBCInterface().GetDataContainer().SetSide( (M_flags[i] == 0) ? OneD_left : OneD_right );

        Model->GetLinearBCInterface().GetDataContainer().SetBase( make_pair( "function", BCInterface1D_function ) );
        Model->GetLinearBCInterface().GetDataContainer().SetBaseString( "0" );

        Model->GetLinearBCInterface().InsertBC();
    }
}

template< class model >
inline void
MS_Coupling_BoundaryCondition::ApplyBoundaryConditions3D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    for ( UInt j( 0 ); j < M_listSize; ++j )
    {
        Model->GetBCInterface().ReadBC( M_FileName, "boundary_conditions/", M_list[j] );

        Model->GetBCInterface().GetDataContainer().SetName( "BC_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[i] ) + "_" + M_list[j] );
        Model->GetBCInterface().GetDataContainer().SetFlag( M_flags[i] );

        Model->GetBCInterface().InsertBC();
    }
}

template< class model >
inline void
MS_Coupling_BoundaryCondition::ApplyDeltaBoundaryConditions3D( const UInt& i )
{
    model *Model = MS_DynamicCast< model >( M_models[i] );

    for ( UInt j( 0 ); j < M_listSize; ++j )
    {
        Model->GetLinearBCInterface().ReadBC( M_FileName, "boundary_conditions/", M_list[j] );

        Model->GetLinearBCInterface().GetDataContainer().SetName( "BC_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[i] ) + "_" + M_list[j] );
        Model->GetLinearBCInterface().GetDataContainer().SetFlag( M_flags[i] );

        Model->GetLinearBCInterface().GetDataContainer().SetBase( make_pair( "function", BCInterface_function ) );
        Model->GetLinearBCInterface().GetDataContainer().SetBaseString( "0" );

        Model->GetLinearBCInterface().InsertBC();
    }
}

} // Namespace LifeV

#endif /* MS_Coupling_BoundaryCondition_H */
