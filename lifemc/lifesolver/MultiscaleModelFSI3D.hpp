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
 *  @brief File containing the Multiscale Model FSI3D
 *
 *  @date 19-04-2010
 *  @author Paolo Crosetto <paolo.crosetto@epfl.ch>
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleModelFSI3D_H
#define MultiscaleModelFSI3D_H 1

//#define FSI_WITH_EXTERNALPRESSURE

// LifeV includes

#include <life/lifefilters/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <life/lifefilters/ExporterHDF5.hpp>
#endif

#include <life/lifefem/BCInterface3D.hpp>
#include <life/lifealg/NonLinearRichardson.hpp>

#include <life/lifesolver/FSIOperator.hpp>
#include <life/lifesolver/FSIMonolithicGE.hpp>
#include <life/lifesolver/FSIMonolithicGI.hpp>

#include <life/lifesolver/MonolithicBlockMatrix.hpp>
#include <life/lifesolver/MonolithicBlockMatrixRN.hpp>
#include <life/lifesolver/MonolithicBlockComposedDN.hpp>
#include <life/lifesolver/MonolithicBlockComposedNN.hpp>
#include <life/lifesolver/MonolithicBlockComposedDNND.hpp>

// Mathcard includes
#include <lifemc/lifesolver/MultiscaleModel.hpp>
#include <lifemc/lifesolver/MultiscaleInterfaceFluid.hpp>

namespace LifeV
{
namespace Multiscale
{

#ifndef FSI_WITH_EXTERNALPRESSURE
// Forward declaration
class FSI3DCouplingFunction;
#endif

//! MultiscaleModelFSI3D - Multiscale model for 3D FSI simulations
/*!
 *  @author Paolo Crosetto
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 */
class MultiscaleModelFSI3D: public virtual multiscaleModel_Type,
                            public virtual MultiscaleInterfaceFluid
{
public:

    //! @name Public Types
    //@{

    typedef FSIOperator                                FSIOperator_Type;
    typedef boost::shared_ptr< FSIOperator_Type>       FSIOperatorPtr_Type;

    typedef FSIOperator::data_Type                     data_Type;
    typedef FSIOperator::dataPtr_Type                  dataPtr_Type;

    typedef FSIOperator::mesh_Type                     mesh_Type;

    typedef FSIOperator::fluid_Type                    fluid_Type;
    typedef FSIOperator::solid_Type                    solid_Type;

    typedef FSIOperator::vector_Type                   vector_Type;
    typedef FSIOperator::vectorPtr_Type                vectorPtr_Type;

    typedef Exporter< mesh_Type >                      IOFile_Type;
    typedef boost::shared_ptr< IOFile_Type >           IOFilePtr_Type;
    typedef ExporterData< mesh_Type >                  IOData_Type;

    typedef ExporterEnsight< mesh_Type >               ensightIOFile_Type;
#ifdef HAVE_HDF5
    typedef ExporterHDF5< mesh_Type >                  hdf5IOFile_Type;
#endif

    typedef BCHandler                                  bc_Type;
    typedef boost::shared_ptr< bc_Type >               bcPtr_Type;
    typedef BCInterface3D< bc_Type, FSIOperator >      bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >      bcInterfacePtr_Type;

#ifndef FSI_WITH_EXTERNALPRESSURE
    typedef FSI3DCouplingFunction                      couplingFunction_Type;
    typedef boost::shared_ptr< couplingFunction_Type > couplingFunctionPtr_Type;
    typedef std::vector< couplingFunction_Type >       couplingFunctionsContainer_Type;
#endif

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructor
    explicit MultiscaleModelFSI3D();

    //! Destructor
    virtual ~MultiscaleModelFSI3D() {}

    //@}


    //! @name MultiscaleModel Methods
    //@{

    //! Setup the data of the model.
    /*!
     * @param fileName Name of data file.
     */
    void setupData( const std::string& fileName );

    //! Setup the model.
    void setupModel();

    //! Build the initial model.
    void buildModel();

    //! Update the model.
    void updateModel();

    //! Solve the model.
    void solveModel();

    //! Update the solution.
    void updateSolution();

    //! Save the solution
    void saveSolution();

    //! Display some information about the model.
    void showMe();

    //! Return a specific scalar quantity to be used for a comparison with a reference value.
    /*!
     * This method is meant to be used for night checks.
     * @return reference quantity.
     */
    Real checkSolution() const;

    //@}


    //! @name MultiscaleInterfaceFluid Methods
    //@{

    //! Impose the flow rate on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @param function boundary condition function
     */
    void imposeBoundaryFlowRate( const bcFlag_Type& flag, const function_Type& function );

    //! Impose the integral of the normal stress on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @param function boundary condition function
     */
    void imposeBoundaryStress( const bcFlag_Type& flag, const function_Type& function );

    //! Get the flow rate on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return flow rate value
     */
    Real boundaryFlowRate( const bcFlag_Type& flag ) const { return M_FSIoperator->fluid().flux( flag, M_FSIoperator->solution() ); }

    //! Get the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param stressType Type of approximation for the stress
     * @return stress value
     */
    Real boundaryStress( const bcFlag_Type& flag ) const { return -boundaryPressure( flag ); }

    //! Get the variation of the flow rate (on a specific boundary face) using the linear model
    /*!
     * @param flag flag of the boundary face on which quantity should be computed
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flow rate
     */
    Real boundaryDeltaFlowRate( const bcFlag_Type& flag, bool& solveLinearSystem );

    //! Get the variation of the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @param stressType Type of approximation for the stress
     * @return variation of the stress
     */
    Real boundaryDeltaStress( const bcFlag_Type& flag, bool& solveLinearSystem );

    //@}


    //! @name Get Methods
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
    Real boundaryDensity( const bcFlag_Type& /*flag*/) const { return M_FSIoperator->dataFluid()->density(); }

    //! Get the viscosity on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return viscosity value
     */
    Real boundaryViscosity( const bcFlag_Type& /*flag*/) const { return M_FSIoperator->dataFluid()->viscosity(); }

    //! Get the area on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return area value
     */
    Real boundaryArea( const bcFlag_Type& flag ) const { return M_FSIoperator->fluid().area( flag ); }

    //! Get the integral of the pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return pressure value
     */
    Real boundaryPressure( const bcFlag_Type& flag ) const;

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

    //! Initialize the FSI solution.
    void initializeSolution();

    void setupCommunicator();

    void setupBC( const std::string& fileName );
    void updateBC();

    void setupExporter( IOFilePtr_Type& exporter, const GetPot& dataFile, const std::string& label = "" );
    void setupImporter( IOFilePtr_Type& exporter, const GetPot& dataFile, const std::string& label = "" );

    void setExporterFluid( const IOFilePtr_Type& exporter );
    void setExporterSolid( const IOFilePtr_Type& exporter );

    //! Setup the linear model
    void setupLinearModel();

    //! Update the linear system matrix and vectors
    void updateLinearModel();

    //! Solve the linear problem
    void solveLinearModel( bool& solveLinearSystem );

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

#ifndef FSI_WITH_EXTERNALPRESSURE
    // Stress coupling function container
    couplingFunctionsContainer_Type        M_stressCouplingFunction;

    // Scalar external pressure
    Real                                   M_externalPressureScalar;
#endif

    // Vectorial external pressure
    vectorPtr_Type                         M_externalPressureVector;

    // Post processing members TODO NOW SHOULD BE REMOVED
    vectorPtr_Type                         M_fluidVelocityAndPressure;
    vectorPtr_Type                         M_fluidDisplacement;
    vectorPtr_Type                         M_solidVelocity;
    vectorPtr_Type                         M_solidDisplacement;

    // Old Solution (used for subiterations)
    std::vector<vectorPtr_Type>            M_fluidVelocityAndPressure_tn;
    std::vector<vectorPtr_Type>            M_fluidDisplacement_tn;
    std::vector<vectorPtr_Type>            M_solidVelocity_tn;
    std::vector<vectorPtr_Type>            M_solidDisplacement_tn;

    vectorPtr_Type                         M_stateVariable;

    UInt                                   M_nonLinearRichardsonIteration;

    // Boundary Conditions
    bcInterfacePtr_Type                    M_fluidBC;
    bcInterfacePtr_Type                    M_solidBC;
    bcInterfacePtr_Type                    M_harmonicExtensionBC;

    // Linear Fluid problem
    bcPtr_Type                             M_linearBC;
    vectorPtr_Type                         M_linearRHS;
    vectorPtr_Type                         M_linearSolution;
};



#ifndef FSI_WITH_EXTERNALPRESSURE

//! FSI3DCouplingFunction - The FSI3D coupling function
/*!
 *  @author Cristiano Malossi
 *
 *  This simple class provides the implementation for the BC function used by the FSI3D model
 *  in order to apply an offset to the coupling quantity (e.g.: to remove the wall pressure).
 */
class FSI3DCouplingFunction
{
public:

    //! @name Type definitions
    //@{

    typedef MultiscaleInterfaceFluid::function_Type function_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     * @param couplingFunction original coupling function
     * @param delta delta to be applied
     */
    explicit FSI3DCouplingFunction( const function_Type& couplingFunction, const Real& delta ) : M_couplingFunction( couplingFunction ), M_delta( delta ) {}

    //! Destructor
    virtual ~FSI3DCouplingFunction() {}

    //@}


    //! @name Methods
    //@{

    //! Evaluate the coupling quantity
    /*!
     * @return evaluation of the function
     */
    Real function( const Real& t, const Real& x, const Real& y, const Real& z, const UInt& id )
    {
        return M_couplingFunction( t, x, y, z, id ) + M_delta;
    }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    FSI3DCouplingFunction();

    FSI3DCouplingFunction( const MultiscaleCoupling& coupling );

    FSI3DCouplingFunction& operator=( const MultiscaleCoupling& coupling );

    //@}

    function_Type                          M_couplingFunction;
    Real                                   M_delta;
};

#endif

//! Factory create function
inline multiscaleModel_Type* createMultiscaleModelFSI3D()
{
    return new MultiscaleModelFSI3D();
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleModelFSI3D_H */
