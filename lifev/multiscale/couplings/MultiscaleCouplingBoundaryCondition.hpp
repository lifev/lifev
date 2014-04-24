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
 *  @brief File containing the Multiscale Coupling BoundaryCondition
 *
 *  @date 02-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleCouplingBoundaryCondition_H
#define MultiscaleCouplingBoundaryCondition_H 1

#if defined(LIFEV_HAS_ZERODIMENSIONAL)
#include <lifev/bc_interface/0D/bc/BCInterface0D.hpp>
#endif

#if defined(LIFEV_HAS_ONEDFSI)
#include <lifev/bc_interface/1D/bc/BCInterface1D.hpp>
#endif

#if defined(LIFEV_HAS_NAVIERSTOKES)
#if defined(LIFEV_HAS_FSI)
#include <lifev/bc_interface/3D/bc/BCInterfaceFSI3D.hpp>
#else
#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#endif
#endif

#include <lifev/multiscale/couplings/MultiscaleCoupling.hpp>

#if defined(LIFEV_HAS_ZERODIMENSIONAL)
#include <lifev/multiscale/models/MultiscaleModelWindkessel0D.hpp>
#include <lifev/multiscale/models/MultiscaleModel0D.hpp>
#endif

#if defined(LIFEV_HAS_ONEDFSI)
#include <lifev/multiscale/models/MultiscaleModelFSI1D.hpp>
#endif

#if defined(LIFEV_HAS_NAVIERSTOKES)
#include <lifev/multiscale/models/MultiscaleModelFluid3D.hpp>
#endif

#if defined(LIFEV_HAS_FSI)
#include <lifev/multiscale/models/MultiscaleModelFSI3D.hpp>
#endif

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleCouplingBoundaryCondition - Coupling condition for standard boundary conditions
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleCouplingBoundaryCondition class is an implementation of the multiscaleCoupling_Type
 *  for applying standard boundary conditions on the models.
 */
class MultiscaleCouplingBoundaryCondition: public virtual multiscaleCoupling_Type
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleCouplingBoundaryCondition();

    //! Destructor
    virtual ~MultiscaleCouplingBoundaryCondition() {}

    //@}


    //! @name Multiscale PhysicalCoupling Implementation
    //@{

    //! Setup the data of the coupling
    /*!
     *  @param fileName Name of data file
     */
    void setupData ( const std::string& fileName );

    //! Setup the coupling variables number.
    void setupCouplingVariablesNumber();

    //! Setup the coupling
    void setupCoupling();

    //! Initialize the values of the coupling variables (DO NOTHING)
    void initializeCouplingVariables() {}

    //! Update the coupling
    /*!
     * Nothing to do for boundary conditions
     */
    void updateCoupling() {};

    //! Compute the local coupling residuals vector
    void computeCouplingResiduals() {}

    //@}

private:

    //! @name Private Multiscale PhysicalCoupling Implementation
    //@{

    //! Build the list of models affected by the perturbation of a local coupling variable (DO NOTHING)
    /*!
     * @param localCouplingVariableID id of the perturbed local coupling variable
     * @param perturbedModelsList list of models affected by the perturbation
     */
    void exportListOfPerturbedModels ( const UInt& /*localCouplingVariableID*/, multiscaleModelsContainer_Type& /*perturbedModelsList*/ ) {}

    //! Insert constant coefficients into the Jacobian matrix (DO NOTHING)
    /*!
     * @param jacobian the Jacobian matrix
     */
    void insertJacobianConstantCoefficients ( multiscaleMatrix_Type& /*jacobian*/ ) {}

    //! Insert the Jacobian coefficient(s) depending on a perturbation of the model, due to a specific variable (the column) (DO NOTHING)
    /*!
     * @param jacobian          the Jacobian matrix
     * @param column            the column related to the perturbed variable
     * @param ID                the global ID of the model which is perturbed by the variable
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     */
    void insertJacobianDeltaCoefficients ( multiscaleMatrix_Type& /*jacobian*/, const UInt& /*column*/, const UInt& /*ID*/, bool& /*solveLinearSystem*/ ) {}

    //@}


    //! @name Private Methods
    //@{

#if defined(LIFEV_HAS_ZERODIMENSIONAL)
    //! Apply the boundary condition to the specific 0D model
    template< class ModelType >
    void applyBoundaryConditions0D ( const UInt& i );
#endif

#if defined(LIFEV_HAS_ONEDFSI)
    //! Apply the boundary condition to the specific 1D model
    template< class ModelType >
    void applyBoundaryConditions1D ( const UInt& i );
#endif

#if defined(LIFEV_HAS_FSI) || defined(LIFEV_HAS_NAVIERSTOKES)
    //! Apply the boundary condition to the specific 3D model
    template< class ModelType >
    void applyBoundaryConditions3D ( const UInt& i );
#endif
    //@}


    //! @name Unimplemented Methods
    //@{

    MultiscaleCouplingBoundaryCondition ( const MultiscaleCouplingBoundaryCondition& coupling );

    MultiscaleCouplingBoundaryCondition& operator= ( const MultiscaleCouplingBoundaryCondition& coupling );

    //@}

    std::string                M_fileName;

    std::vector< bcName_Type > M_list;
    UInt                       M_listSize;
};

//! Factory create function
inline multiscaleCoupling_Type* createMultiscaleCouplingBoundaryCondition()
{
    return new MultiscaleCouplingBoundaryCondition();
}

// ===================================================
// Template implementation
// ===================================================
#if defined(LIFEV_HAS_ZERODIMENSIONAL)
template< class ModelType >
inline void
MultiscaleCouplingBoundaryCondition::applyBoundaryConditions0D ( const UInt& i )
{
    boost::shared_ptr< ModelType > model = multiscaleDynamicCast< ModelType > ( M_models[i] );
    multiscaleID_Type flag ( model->boundaryFlag ( M_boundaryIDs[i] ) );

    for ( UInt j ( 0 ); j < M_listSize; ++j )
    {
        model->bcInterface().readBC ( M_fileName, "boundary_conditions/", M_list[j] );

        model->bcInterface().dataContainer().setFlag ( flag );

        model->bcInterface().insertBC();
    }
}
#endif

#if defined(LIFEV_HAS_ONEDFSI)
template< class ModelType >
inline void
MultiscaleCouplingBoundaryCondition::applyBoundaryConditions1D ( const UInt& i )
{
    boost::shared_ptr< ModelType > model = multiscaleDynamicCast< ModelType > ( M_models[i] );
    multiscaleID_Type flag ( model->boundaryFlag ( M_boundaryIDs[i] ) );

    for ( UInt j ( 0 ); j < M_listSize; ++j )
    {
        model->bcInterface().readBC ( M_fileName, "boundary_conditions/", M_list[j] );

        model->bcInterface().dataContainer().setSide ( ( flag == 0 ) ? OneDFSI::left : OneDFSI::right );

        model->bcInterface().insertBC();
    }
}
#endif

#if defined(LIFEV_HAS_FSI) || defined(LIFEV_HAS_NAVIERSTOKES)
template< class ModelType >
inline void
MultiscaleCouplingBoundaryCondition::applyBoundaryConditions3D ( const UInt& i )
{
    boost::shared_ptr< ModelType > model = multiscaleDynamicCast< ModelType > ( M_models[i] );
    multiscaleID_Type flag ( model->boundaryFlag ( M_boundaryIDs[i] ) );

    for ( UInt j ( 0 ); j < M_listSize; ++j )
    {
        model->bcInterface().readBC ( M_fileName, "boundary_conditions/", M_list[j] );

        model->bcInterface().dataContainer().setName ( "CouplingBC_Model_" + number2string ( model->ID() ) + "_BoundaryID_" + number2string ( M_boundaryIDs[i] ) + "_" + M_list[j] );
        model->bcInterface().dataContainer().setFlag ( flag );

        model->bcInterface().insertBC();
    }
}
#endif

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleCouplingBoundaryCondition_H */
