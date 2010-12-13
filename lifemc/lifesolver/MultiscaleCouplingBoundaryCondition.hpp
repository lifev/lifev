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
 *  @brief File containing the MultiScale Coupling BoundaryCondition
 *
 *  @date 02-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleCouplingBoundaryCondition_H
#define MultiscaleCouplingBoundaryCondition_H 1

#include <lifemc/lifesolver/BCInterface.hpp>
#include <lifemc/lifesolver/BCInterface1D.hpp>

#include <lifemc/lifesolver/MultiscaleCoupling.hpp>
#include <lifemc/lifesolver/MultiscaleModelFluid3D.hpp>
#include <lifemc/lifesolver/MultiscaleModelFSI3D.hpp>
#include <lifemc/lifesolver/MultiscaleModel1D.hpp>

namespace LifeV
{
namespace multiscale
{

//! MultiscaleCouplingBoundaryCondition - Coupling condition for standard boundary conditions
/*!
 *  @author Cristiano Malossi
 *
 *  The MultiscaleCouplingBoundaryCondition class is an implementation of the MS_Coupling_Type
 *  for applying standard boundary conditions on the models.
 */
class MultiscaleCouplingBoundaryCondition: public virtual MS_Coupling_Type
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MultiscaleCouplingBoundaryCondition();

    //! Destructor
    virtual ~MultiscaleCouplingBoundaryCondition() {}

    //@}


    //! @name MultiScale PhysicalCoupling Implementation
    //@{

    //! Setup the data of the coupling
    /*!
     *  @param fileName Name of data file
     */
    void setupData( const std::string& fileName );

    //! Setup the coupling
    void setupCoupling();

    //! Initialize the values of the coupling variables (DO NOTHING)
    void initializeCouplingVariables() {}

    //! Export the values of the local coupling residuals into a global vector (DO NOTHING)
    /*!
     * @param couplingResiduals Global vector of variables
     */
    void exportCouplingResiduals( MS_Vector_Type& /*couplingResiduals*/ ) {}

    //! Display some information about the coupling
    void showMe();

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleCouplingBoundaryCondition( const MultiscaleCouplingBoundaryCondition& coupling );

    MultiscaleCouplingBoundaryCondition& operator=( const MultiscaleCouplingBoundaryCondition& coupling );

    //@}


    //! @name Private MultiScale PhysicalCoupling Implementation
    //@{

    //! Build the list of models affected by the perturbation of a local coupling variable
    /*!
     * @param localCouplingVariableID local coupling variable (perturbed)
     * @return list of models affected by the perturbation
     */
    MS_ModelsVector_Type listOfPerturbedModels( const UInt& /*localCouplingVariableID*/ );

    //! Insert constant coefficients into the Jacobian matrix (DO NOTHING)
    /*!
     * @param jacobian the Jacobian matrix
     */
    void insertJacobianConstantCoefficients( MS_Matrix_Type& /*jacobian*/ ) {}

    //! Insert the Jacobian coefficient(s) depending on a perturbation of the model, due to a specific variable (the column) (DO NOTHING)
    /*!
     * @param jacobian          the Jacobian matrix
     * @param column            the column related to the perturbed variable
     * @param ID                the global ID of the model which is perturbed by the variable
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     */
    void insertJacobianDeltaCoefficients( MS_Matrix_Type& /*jacobian*/, const UInt& /*column*/, const UInt& /*ID*/, bool& /*solveLinearSystem*/ ) {}

    //! Display some information about the coupling
    /*!
     * @param output specify the output stream
     */
    void displayCouplingValues( std::ostream& output );

    //@}


    //! @name Private Methods
    //@{

    //! Apply the boundary condition to the specific 1D model
    template< class ModelType >
    void applyBoundaryConditions1D( const UInt& i );

    //! Apply the boundary condition to the specific 3D model
    template< class ModelType >
    void applyBoundaryConditions3D( const UInt& i );

    //@}

    std::string           M_fileName;

    std::vector< BCName > M_list;
    UInt                  M_listSize;
};

//! Factory create function
inline MS_Coupling_Type* createMultiscaleCouplingBoundaryCondition()
{
    return new MultiscaleCouplingBoundaryCondition();
}

// ===================================================
// Template implementation
// ===================================================
template< class ModelType >
inline void
MultiscaleCouplingBoundaryCondition::applyBoundaryConditions1D( const UInt& i )
{
    ModelType *model = multiscaleDynamicCast< ModelType >( M_models[i] );

    for ( UInt j( 0 ); j < M_listSize; ++j )
    {
        model->bcInterface().readBC( M_fileName, "boundary_conditions/", M_list[j] );

        model->bcInterface().dataContainer().setSide( (M_flags[i] == 0) ? OneD_left : OneD_right );

        model->bcInterface().insertBC();
    }
}

template< class ModelType >
inline void
MultiscaleCouplingBoundaryCondition::applyBoundaryConditions3D( const UInt& i )
{
    ModelType *model = multiscaleDynamicCast< ModelType >( M_models[i] );

    for ( UInt j( 0 ); j < M_listSize; ++j )
    {
        model->bcInterface().readBC( M_fileName, "boundary_conditions/", M_list[j] );

        model->bcInterface().dataContainer().setName( "CouplingBC_Model_" + number2string( model->ID() ) + "_Flag_" + number2string( M_flags[i] ) + "_" + M_list[j] );
        model->bcInterface().dataContainer().setFlag( M_flags[i] );

        model->bcInterface().insertBC();
    }
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleCouplingBoundaryCondition_H */
