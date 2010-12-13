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
 *  @brief File containing the MultiScale Model FSI3D
 *
 *  @date 19-04-2010
 *  @author Paolo Crosetto <paolo.crosetto@epfl.ch>
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MS_MODEL_FSI3D_H
#define MS_MODEL_FSI3D_H 1

// LifeV includes
#include <life/lifesolver/FSIOperator.hpp>
#include <life/lifealg/nonLinRichardson.hpp>

#include <life/lifefilters/ensight.hpp>
#ifdef HAVE_HDF5
#include <life/lifefilters/hdf5exporter.hpp>
#endif

// Mathcard includes
#include <lifemc/lifesolver/MonolithicGE.hpp>
#include <lifemc/lifesolver/MonolithicGI.hpp>

#include <lifemc/lifesolver/BlockMatrix.hpp>
#include <lifemc/lifesolver/BlockMatrixRN.hpp>
#include <lifemc/lifesolver/ComposedDN.hpp>
#include <lifemc/lifesolver/ComposedNN.hpp>
#include <lifemc/lifesolver/ComposedDNND.hpp>

#include <lifemc/lifesolver/BCInterface.hpp>
#include <lifemc/lifesolver/MultiscaleModel.hpp>

namespace LifeV
{
namespace multiscale
{

//! MultiscaleModelFSI3D - MultiScale model for 3D FSI simulations
/*!
 *  @author Paolo Crosetto
 */
class MultiscaleModelFSI3D: public virtual multiscaleModel_Type
{
public:

    //! @name Public Types
    //@{

    typedef FSIOperator                                                                    FSIOperator_Type;
    typedef boost::shared_ptr< FSIOperator_Type>                                           FSIOperatorPtr_Type;

    typedef FSIOperator::data_Type                                                         data_Type;
    typedef FSIOperator::data_PtrType                                                      dataPtr_Type;

    typedef FSIOperator::mesh_type                                                         mesh_Type;

    typedef FSIOperator::fluid_raw_type                                                    fluid_Type;
    typedef FSIOperator::solid_raw_type                                                    solid_Type;

    typedef FSIOperator::vector_type                                                       vector_Type;
    typedef FSIOperator::vector_ptrtype                                                    vectorPtr_Type;

    typedef Exporter< mesh_Type >                                                          IOFile_Type;
    typedef boost::shared_ptr< IOFile_Type >                                               IOFilePtr_Type;

    typedef Ensight< mesh_Type >                                                           ensightIOFile_Type;
#ifdef HAVE_HDF5
    typedef Hdf5exporter< mesh_Type >                                                      hdf5IOFile_Type;
#endif

    typedef BCHandler                                                                      bc_Type;
    typedef boost::shared_ptr< bc_Type >                                                   bcPtr_Type;
    typedef BCInterface< FSIOperator >                                                     bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >                                          bcInterfacePtr_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructor
    explicit MultiscaleModelFSI3D();

    //! Destructor
    virtual ~MultiscaleModelFSI3D() {}

    //@}


    //! @name MultiScale PhysicalModel Virtual Methods
    //@{

    //! Setup the data of the model.
    /*!
     * @param fileName Name of data file.
     */
    void setupData( const std::string& fileName );

    //! Setup the model.
    void setupModel();

    //! Build the initial system (matrix and vectors).
    void buildSystem();

    //! Update the system for (matrix and vectors).
    void updateSystem();

    //! Solve the problem.
    void solveSystem( );

    //! save the solution
    void saveSolution();

    //! Display some information about the model.
    void showMe();

    //@}


    //! @name Methods
    //@{

    //! Setup the linear model
    void setupLinearModel();

    //! Update the linear system matrix and vectors
    void updateLinearModel();

    //! Solve the linear problem
    void solveLinearModel( bool& solveLinearSystem );

    //@}


    //! @name Get Methods (couplings)
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    bcInterface_Type& bcInterface() { return *M_fluidBC; }

    //! Get the density on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return density value
     */
    Real boundaryDensity( const BCFlag& /*flag*/) const { return M_FSIoperator->dataFluid()->density(); }

    //! Get the viscosity on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return viscosity value
     */
    Real boundaryViscosity( const BCFlag& /*flag*/) const { return M_FSIoperator->dataFluid()->viscosity(); }

    //! Get the area on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return area value
     */
    Real boundaryArea( const BCFlag& flag ) const { return M_FSIoperator->fluid().area( flag ); }

    //! Get the flow rate on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return flow rate value
     */
    Real boundaryFlowRate( const BCFlag& flag ) const { return M_FSIoperator->fluid().flux( flag, *M_FSIoperator->solutionPtr() ); }

    //! Get the integral of the pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return pressure value
     */
    Real boundaryPressure( const BCFlag& flag ) const { return M_FSIoperator->fluid().pressure( flag, *M_FSIoperator->solutionPtr() ); }

    //! Get the integral of the dynamic pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return dynamic pressure value
     */
    Real boundaryDynamicPressure( const BCFlag& flag ) const { return 0.5 * boundaryDensity( flag ) * ( boundaryFlowRate( flag ) * boundaryFlowRate( flag ) ) / ( boundaryArea( flag ) * boundaryArea( flag ) ); }

    //! Get the value of the Lagrange multiplier associated to a specific boundary face
    /*!
     * @param flag flag of the boundary face
     * @return Lagrange multiplier value
     */
    Real boundaryLagrangeMultiplier( const BCFlag& flag ) const { return M_FSIoperator->fluid().LagrangeMultiplier(flag, *M_fluidBC->handler(), M_FSIoperator->getSolution() ); }

    //! Get the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param stressType Type of approximation for the stress
     * @return stress value
     */
    Real boundaryStress( const BCFlag& flag, const stress_Type& stressType = StaticPressure ) const;

    //! Get the variation of the flow rate (on a specific boundary face) using the linear model
    /*!
     * @param flag flag of the boundary face on which quantity should be computed
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flow rate
     */
    Real boundaryDeltaFlowRate( const BCFlag& flag, bool& solveLinearSystem );

    //! Get the variation of the pressure (on a specific boundary face) using the linear model
    /*!
     * @param flag flag of the boundary face on which quantity should be computed
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the pressure
     */
    Real boundaryDeltaPressure( const BCFlag& flag, bool& solveLinearSystem );

    //! Get the variation of the total pressure (on a specific boundary face) using the linear model
    /*!
     * @param flag flag of the boundary face on which quantity should be computed
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the dynamic pressure
     */
    Real boundaryDeltaDynamicPressure( const BCFlag& flag, bool& solveLinearSystem );

    //! Get the variation of the Lagrange multiplier associated to a specific boundary face, using the linear model
    /*!
     * @param flag flag of the boundary face
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return Lagrange multiplier value
     */
    Real boundaryDeltaLagrangeMultiplier( const BCFlag& flag, bool& solveLinearSystem );

    //! Get the variation of the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @param stressType Type of approximation for the stress
     * @return variation of the stress
     */
    Real boundaryDeltaStress( const BCFlag& flag, bool& solveLinearSystem, const stress_Type& stressType = StaticPressure );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleModelFSI3D( const MultiscaleModelFSI3D& model );

    MultiscaleModelFSI3D& operator=( const MultiscaleModelFSI3D& model );

    //@}


    //! @name Private Methods
    //@{

    //! Setup the global data of the model.
    /*!
     * In particular, it replaces the default local values with the ones in the global container.
     * If a value is already specified in the data file, do not perform the replacement.
     *
     * @param fileName File name of the specific model.
     */
    void setupGlobalData( const std::string& fileName );

    void setupCommunicator();

    void setupBC( const std::string& fileName );
    //void setupSegregatedBC( const std::string& fileName );
    void updateBC();

    void setupExporter( IOFilePtr_Type& exporter, const GetPot& dataFile, const std::string& label = "" );
    void setupImporter( IOFilePtr_Type& exporter, const GetPot& dataFile, const std::string& label = "" );

    void setExporterFluid( const IOFilePtr_Type& exporter );
    void setExporterSolid( const IOFilePtr_Type& exporter );

    //! Initialize the solution.
    void initializeSolution();

    //! Impose the coupling perturbation on the correct BC inside the BCHandler
    void imposePerturbation();

    //! Reset all the coupling perturbations imposed on the BCHandler
    void resetPerturbation();

    Real bcFunctionDeltaZero( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const UInt& /*id*/ ) { return 0.; }
    Real bcFunctionDeltaOne( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const UInt& /*id*/ ) { return 1.; }

    //@}

    // Operator
    FSIOperatorPtr_Type                    M_FSIoperator;

    // Data
    dataPtr_Type                           M_data;

    // Exporters
    IOFilePtr_Type                         M_exporterFluid;
    IOFilePtr_Type                         M_exporterSolid;

    // Importers
    IOFilePtr_Type                         M_importerFluid;
    IOFilePtr_Type                         M_importerSolid;

    // Solution
    vectorPtr_Type                         M_fluidVelocityPressure;
    vectorPtr_Type                         M_fluidDisplacement;
    vectorPtr_Type                         M_solidVelocity;
    vectorPtr_Type                         M_solidDisplacement;

    // Old Solution
    vectorPtr_Type                         M_fluidVelocityPressure_tn;
    vectorPtr_Type                         M_solidDisplacement_tn;
    vectorPtr_Type                         M_solidDisplacementOld_tn;
    vectorPtr_Type                         M_rhs_tn;
    UInt                                   M_nonLinearRichardsonIteration;

    // Boundary Conditions
    bcInterfacePtr_Type                    M_fluidBC;
    bcInterfacePtr_Type                    M_solidBC;
    bcInterfacePtr_Type                    M_harmonicExtensionBC;

    // Linear Fluid problem
    bcPtr_Type                             M_linearBC;
    vectorPtr_Type                         M_linearRHS;
    vectorPtr_Type                         M_linearSolution;

    // BC Functions for tangent problem
    BCFunctionBase                         M_bcBaseDeltaZero;
    BCFunctionBase                         M_bcBaseDeltaOne;
};

//! Factory create function
inline multiscaleModel_Type* createMultiscaleModelFSI3D()
{
    return new MultiscaleModelFSI3D();
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MS_MODEL_FSI3D_H */
