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

// If the following macro is defined, the pressure inside the 3-D FSI model
// is initialized to the external pressure. Otherwise, the external pressure
// is applied just as an offset at the interfaces.
//#define FSI_WITH_EXTERNALPRESSURE

// If the following macro is defined, any flow rate BC is prescribed in the
// classical way, i.e., weakly through a Lagrange multiplier. Otherwise, the user
// can chose to prescribe it through an essential BC (a given velocity profile).
//#define FSI_WITHOUT_VELOCITYPROFILE

// This macro now serves just to easily identify the part of the code where the
// methods needed for coupling the area are called. Apart from that, it is useless,
// and in the future can be removed.
#define FSI_WITH_BOUNDARYAREA

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif

#include <lifev/core/algorithm/NonLinearRichardson.hpp>
#include <lifev/bc_interface/fem/BCInterface3D.hpp>

#include <lifev/fsi/solver/FSIMonolithic.hpp>
#include <lifev/fsi/solver/FSIMonolithicGE.hpp>
#include <lifev/fsi/solver/FSIMonolithicGI.hpp>

#include <lifev/fsi/solver/MonolithicBlockMatrix.hpp>
#include <lifev/fsi/solver/MonolithicBlockMatrixRN.hpp>
#include <lifev/fsi/solver/MonolithicBlockComposedDN.hpp>
#include <lifev/fsi/solver/MonolithicBlockComposedNN.hpp>
#include <lifev/fsi/solver/MonolithicBlockComposedDNND.hpp>

#include <lifev/multiscale/solver/MultiscaleModel.hpp>
#include <lifev/multiscale/solver/MultiscaleInterface.hpp>

namespace LifeV
{
namespace Multiscale
{

#ifndef FSI_WITH_EXTERNALPRESSURE
// Forward declaration
class FSI3DBoundaryStressFunction;
#endif

#ifndef FSI_WITHOUT_VELOCITYPROFILE
// Forward declaration
class FSI3DBoundaryFlowRateFunction;
#endif

#ifdef FSI_WITH_BOUNDARYAREA
// Forward declaration
class FSI3DBoundaryAreaFunction;
#endif

//! MultiscaleModelFSI3D - Multiscale model for 3D FSI simulations
/*!
 *  @author Paolo Crosetto
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 */
class MultiscaleModelFSI3D: public virtual multiscaleModel_Type,
                            public virtual MultiscaleInterface
{
public:

    //! @name Public Types
    //@{

    typedef FSIMonolithic                                      FSIOperator_Type;
    typedef boost::shared_ptr< FSIOperator_Type>               FSIOperatorPtr_Type;

    typedef FSIOperator::data_Type                             data_Type;
    typedef FSIOperator::dataPtr_Type                          dataPtr_Type;

    typedef FSIOperator::mesh_Type                             mesh_Type;

    typedef FSIOperator::fluid_Type                            fluid_Type;
    typedef FSIOperator::solid_Type                            solid_Type;

    typedef FSIOperator::vector_Type                           vector_Type;
    typedef FSIOperator::vectorPtr_Type                        vectorPtr_Type;

    typedef Exporter< mesh_Type >                              IOFile_Type;
    typedef boost::shared_ptr< IOFile_Type >                   IOFilePtr_Type;
    typedef ExporterData< mesh_Type >                          IOData_Type;

    typedef ExporterEnsight< mesh_Type >                       ensightIOFile_Type;
#ifdef HAVE_HDF5
    typedef ExporterHDF5< mesh_Type >                          hdf5IOFile_Type;
#endif

    typedef BCHandler                                          bc_Type;
    typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
    typedef BCInterface3D< bc_Type, FSIOperator >              bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;

#ifndef FSI_WITH_EXTERNALPRESSURE
    typedef FSI3DBoundaryStressFunction                        boundaryStressFunction_Type;
    typedef boost::shared_ptr< boundaryStressFunction_Type >   boundaryStressFunctionPtr_Type;
    typedef std::vector< boundaryStressFunctionPtr_Type >      boundaryStressFunctionContainer_Type;
#endif

#ifndef FSI_WITHOUT_VELOCITYPROFILE

    /*! @enum FSI3DBoundaryFlowRate_Type
    */
    enum FSI3DBoundaryFlowRate_Type
    {
        Weak,     /*!< always impose flow rate weakly (DEFAULT) */
        Semiweak, /*!< impose essential when flow rate = 0, otherwise impose flow rate weakly */
        Strong    /*!< always impose essential with user prescribed velocity profile */
    };

    typedef std::vector< FSI3DBoundaryFlowRate_Type >          boundaryFlowRateTypeContainer_Type;

    typedef std::map< multiscaleID_Type, FSI3DBoundaryFlowRate_Type > boundaryFlowRateMap_Type;

    typedef FSI3DBoundaryFlowRateFunction                      boundaryFlowRateFunction_Type;
    typedef boost::shared_ptr< boundaryFlowRateFunction_Type > boundaryFlowRateFunctionPtr_Type;
    typedef std::vector< boundaryFlowRateFunctionPtr_Type >    boundaryFlowRateFunctionsContainer_Type;
    typedef boundaryFlowRateFunctionsContainer_Type::iterator  boundaryFlowRateFunctionsContainerIterator_Type;
#endif

#ifdef FSI_WITH_BOUNDARYAREA
    typedef FSI3DBoundaryAreaFunction                          boundaryAreaFunction_Type;
    typedef boost::shared_ptr< boundaryAreaFunction_Type >     boundaryAreaFunctionPtr_Type;
    typedef std::vector< boundaryAreaFunctionPtr_Type >        boundaryAreaFunctionsContainer_Type;
    typedef boundaryAreaFunctionsContainer_Type::iterator      boundaryAreaFunctionsContainerIterator_Type;
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


    //! @name MultiscaleInterface Methods
    //@{

    //! Impose the flow rate on a specific interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    void imposeBoundaryFlowRate( const multiscaleID_Type& boundaryID, const function_Type& function );

    //! Impose the integral of the mean normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    void imposeBoundaryMeanNormalStress( const multiscaleID_Type& boundaryID, const function_Type& function );

    //! Impose the integral of the mean total normal stress on a specific boundary interface of the model
    /*!
     * Note: mean total normal stress cannot be imposed at the interfaces of this model.
     *
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    void imposeBoundaryMeanTotalNormalStress( const multiscaleID_Type& /*boundaryID*/, const function_Type& /*function*/ )
    {
        multiscaleErrorCheck( ModelInterface, "Invalid interface [MeanTotalNormalStress] for model type [" + enum2String( M_type, multiscaleModelsMap ) +"]", M_comm->MyPID() == 0 );
    }

    //! Impose the area on a specific boundary interface of the model
    /*!
     * Note: area cannot be imposed at the interfaces of this model.
     *
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    void imposeBoundaryArea( const multiscaleID_Type& boundaryID, const function_Type& function );

    //! Get the flow rate on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return flow rate value
     */
    Real boundaryFlowRate( const multiscaleID_Type& boundaryID ) const { return M_FSIoperator->fluid().flux( boundaryFlag( boundaryID ), *M_stateVariable ); }

    //! Get the integral of the mean normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return mean normal stress value
     */
    Real boundaryMeanNormalStress( const multiscaleID_Type& boundaryID ) const;

    //! Get the integral of the mean total normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return mean total normal stress value
     */
    Real boundaryMeanTotalNormalStress( const multiscaleID_Type& boundaryID ) const;

    //! Get the area on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return area value
     */
    Real boundaryArea( const multiscaleID_Type& boundaryID ) const { return M_FSIoperator->fluid().area( boundaryFlag( boundaryID ) ); }

    //! Get the variation of the flow rate (on a specific boundary interface) using the linear model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flow rate
     */
    Real boundaryDeltaFlowRate( const multiscaleID_Type& boundaryID, bool& solveLinearSystem );

    //! Get the variation of the integral of the mean normal stress (on a specific boundary interface) using the linear model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the mean normal stress
     */
    Real boundaryDeltaMeanNormalStress( const multiscaleID_Type& boundaryID, bool& solveLinearSystem );

    //! Get the variation of the integral of the mean total normal stress (on a specific boundary interface) using the linear model
    /*!
     * TODO The integral terms of the derivative of the area have not been coded yet. They are used only by the GI formulation.
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the mean total normal stress
     */
    Real boundaryDeltaMeanTotalNormalStress( const multiscaleID_Type& boundaryID, bool& solveLinearSystem );

    //! Get the variation of the integral of the area (on a specific boundary interface) using the linear model
    /*!
     *  Note: returns always a NaN since this specific method is never used in the current interface equations.
     *
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the area
     */
    Real boundaryDeltaArea( const multiscaleID_Type& /*boundaryID*/, bool& /*solveLinearSystem*/ ) { return NaN; }

    //@}


    //! @name Get Methods
    //@{

    //! Get the FSI3D fluid bcHandler
    /*!
     * @return FSI3D fluid bcHandler
     */
    bc_Type& bcHandlerFluid() { return *M_fluidBC->handler(); }

    //! Get the FSI3D solid bcHandler
    /*!
     * @return FSI3D solid bcHandler
     */
    bc_Type& bcHandlerSolid() { return *M_solidBC->handler(); }

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    bcInterface_Type& bcInterface() { return *M_fluidBC; }

    //! Get the density on a specific boundary face of the model
    /*!
     * @return density value
     */
    Real boundaryDensity() const { return M_FSIoperator->dataFluid()->density(); }

    //! Get the viscosity on a specific boundary face of the model
    /*!
     * @return viscosity value
     */
    Real boundaryViscosity() const { return M_FSIoperator->dataFluid()->viscosity(); }

    //! Get the integral of the pressure (on a specific boundary face)
    /*!
     * @param boundaryID ID of the boundary interface
     * @return pressure value
     */
    Real boundaryPressure( const multiscaleID_Type& boundaryID ) const;

    //! Get the integral of the total pressure (on a specific boundary face)
    /*!
     * @param boundaryID ID of the boundary interface
     * @return total pressure value
     */
    Real boundaryTotalPressure( const multiscaleID_Type& boundaryID ) const;

    //! Get the external wall pressure
    /*!
     * @return external pressure value
     */
    Real externalPressure() const;

    //! Get the FSI3D data container
    /*!
     * @return FSI3D data container
     */
    const dataPtr_Type& data() const { return M_data; }

    //! Get the FSI3D operator
    /*!
     * @return FSI3D operator
     */
    const FSIOperatorPtr_Type& solver() const { return M_FSIoperator; }

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

    void exportFluidSolution();
    void exportSolidSolution();

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
    FSIOperatorPtr_Type                     M_FSIoperator;

    // Data
    dataPtr_Type                            M_data;

    // Exporters
    IOFilePtr_Type                          M_exporterFluid;
    IOFilePtr_Type                          M_exporterSolid;

    // Importers
    IOFilePtr_Type                          M_importerFluid;
    IOFilePtr_Type                          M_importerSolid;

#ifndef FSI_WITH_EXTERNALPRESSURE
    // Stress coupling function container
    boundaryStressFunctionContainer_Type    M_boundaryStressFunctions;

    // Scalar external pressure
    Real                                    M_externalPressureScalar;
#endif

#ifndef FSI_WITHOUT_VELOCITYPROFILE
    // Flow rate coupling function container
    boundaryFlowRateFunctionsContainer_Type M_boundaryFlowRateFunctions;

    // Flow rate type container
    boundaryFlowRateTypeContainer_Type      M_boundaryFlowRateType;
#endif

#ifdef FSI_WITH_BOUNDARYAREA
    // Boundary area coupling function container
    boundaryAreaFunctionsContainer_Type     M_boundaryAreaFunctions;

    // Free flags, available for the couplings
    multiscaleIDContainer_Type              M_boundaryFlagsArea;
    std::vector< bool >                     M_boundaryFlagsAreaPerturbed;
#endif

    // Post processing members
    vectorPtr_Type                         M_fluidVelocity;
    vectorPtr_Type                         M_fluidPressure;
    vectorPtr_Type                         M_fluidDisplacement;
    vectorPtr_Type                         M_solidVelocity;
    vectorPtr_Type                         M_solidDisplacement;

    vectorPtr_Type                         M_stateVariable;

    UInt                                   M_nonLinearRichardsonIteration;

    // Boundary Conditions
    bcInterfacePtr_Type                    M_fluidBC;
    bcInterfacePtr_Type                    M_solidBC;
    bcInterfacePtr_Type                    M_harmonicExtensionBC;

    // Linear problem
    bcPtr_Type                             M_linearFluidBC;
    bcPtr_Type                             M_linearSolidBC;
    vectorPtr_Type                         M_linearRHS;
    vectorPtr_Type                         M_linearSolution;
};



#ifndef FSI_WITH_EXTERNALPRESSURE

//! FSI3DBoundaryStressFunction - The FSI3D coupling function
/*!
 *  @author Cristiano Malossi
 *
 *  This simple class provides the implementation for the BC function used by the FSI3D model
 *  in order to apply an offset to the coupling quantity (e.g.: to remove the wall pressure).
 */
class FSI3DBoundaryStressFunction
{
public:

    //! @name Type definitions
    //@{

    typedef MultiscaleInterface::function_Type function_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit FSI3DBoundaryStressFunction() : M_function(), M_delta() {}

    //! Destructor
    virtual ~FSI3DBoundaryStressFunction() {}

    //@}


    //! @name Methods
    //@{

    //! Evaluate the coupling quantity
    /*!
     * @return evaluation of the function
     */
    Real function( const Real& t, const Real& x, const Real& y, const Real& z, const UInt& id )
    {
        return M_function( t, x, y, z, id ) + M_delta;
    }

    //@}


    //! @name Set methods
    //@{

    //! Set the offset to be applied to the boundary condition
    /*!
     * @param delta offset to be applied to the boundary condition
     */
    void setDelta( const Real& delta ) { M_delta = delta; }

    //! Set the area function
    /*!
     * @param function area function
     */
    void setFunction( const function_Type& function ) { M_function = function; }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    FSI3DBoundaryStressFunction( const FSI3DBoundaryStressFunction& boundaryFunction );

    FSI3DBoundaryStressFunction& operator=( const FSI3DBoundaryStressFunction& boundaryFunction );

    //@}

    function_Type                          M_function;
    Real                                   M_delta;
};

#endif



#ifndef FSI_WITHOUT_VELOCITYPROFILE

//! FSI3DBoundaryFlowRateFunction - The FSI3D coupling function
/*!
 *  @author Cristiano Malossi
 *  @author Toni Lassila
 *
 *  This simple class provides the implementation for the BC function used by the FSI3D model
 *  in order to convert the given flow rate to a velocity profile.
 */
class FSI3DBoundaryFlowRateFunction
{
public:

    //! @name Type definitions
    //@{

    typedef MultiscaleInterface::function_Type function_Type;

    typedef MultiscaleModelFSI3D::FSI3DBoundaryFlowRate_Type FSI3DBoundaryFlowRate_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit FSI3DBoundaryFlowRateFunction() :
        M_FSI3D                (),
        M_function             (),
        M_fluidFlag            (),
        M_boundaryFlowRateType (),
        M_n                    (),
        M_area                 (),
        M_flowRateIsZero       ()
        {}

    //! Destructor
    virtual ~FSI3DBoundaryFlowRateFunction() { /* M_FSI3D is deleted outside */ }

    //@}


    //! @name Methods
    //@{

    //! Update function parameters
    /*!
     *  Instead of evaluating all the parameters at each call, we update them once at each time step.
     */
    void updateParameters()
    {
        // Updating the area
        M_area = M_FSI3D->solver()->fluid().area( M_fluidFlag );

        // Updating the BC
        switch ( M_boundaryFlowRateType )
        {
            case MultiscaleModelFSI3D::Strong:
            {
                // Update the approximate surface normal direction for this b.c.
                useNormalDirectionFlow( M_FSI3D->bcHandlerFluid().findBCWithFlag( M_fluidFlag ) );

                break;
            }
            case MultiscaleModelFSI3D::Semiweak:
            {
                // Check if the flow rate is (almost) zero and impose b.c. accordingly:
                // we use a cutoff to determine when the coupled flow rate is sufficiently small
                M_flowRateIsZero = ( std::abs( M_function( M_FSI3D->data()->dataFluid()->dataTime()->time(), 0.0, 0.0, 0.0, 0 ) ) < 1e-8 );

                BCFunctionBase base;
                if ( M_flowRateIsZero )
                {
                    // Impose essential Dirichlet
                    M_FSI3D->bcHandlerFluid().modifyBC( M_fluidFlag, Essential );

                    /* In case of essential b.c. we must take care of the Lagrange multiplier, otherwise the global
                       system will be singular. For now we circumvent the problem by setting the corresponding diagonal
                       element to zero, i.e., the Lagrange multiplier is retained but becomes trivial. */
                }
                else
                {
                    // Impose weakly the flow rate
                    M_FSI3D->bcHandlerFluid().modifyBC( M_fluidFlag, Flux );
                }

                break;
            }
            case MultiscaleModelFSI3D::Weak:

                std::cerr << "!!! ERROR: in case Weak option has been selected, this class should not be instantiated at all !!!" << std::endl;
                break;

            default:

                break;
        }
    }

    //! Evaluate the coupling quantity
    /*!
     *  For now, we impose a flat profile in the normal direction.
     *  TODO: In the future, this method can be extended in order to prescribe some "predefined" velocity profile
     *  such as Poiseuille, Womersley, etc...
     *
     *  @return evaluation of the function
     */
    Real function( const Real& t, const Real& x, const Real& y, const Real& z, const UInt& id )
    {
        /* CAUTION: Here lies a possible bug. The flat profile is assumed to extend to the boundary of the
         * surface patch, but if the user specific explicitly a no-slip condition on the ring then an error
         * of order h will be made in the flow rate coupling equations. The proper implementation is to
         * make no assumptions but to integrate the effective flow profile across the surface patch and use
         * that value for the scaling. That way it is always mesh-independent. */

        if ( M_boundaryFlowRateType == MultiscaleModelFSI3D::Semiweak && M_flowRateIsZero )
            return 0;

        Real flowRate = M_function( t, x, y, z, id );
        Real meanVelocity = flowRate / M_area;

        return meanVelocity * M_n[id];
    }

    //@}


    //! @name Set Methods
    //@{

    //! Set the FSI3D model
    /*!
     * @param modelFSI3D a pointer to the FSI3D model
     */
    void setModel( MultiscaleModelFSI3D* modelFSI3D ) { M_FSI3D = modelFSI3D; }

    //! Set the fluid flag of the boundary
    /*!
     * @param flag flag of the fluid boundary
     */
    void setFluidFlag( const multiscaleID_Type& flag ) { M_fluidFlag = flag; }

    //! Set the outgoing normal of the fluid boundary
    /*!
     * @param normal outgoing normal of the fluid boundary
     */
    void setNormal( const boost::array< Real, 3 >& normal ) { M_n = normal; }

    //! Set the area function
    /*!
     * @param function area function
     */
    void setFunction( const function_Type& function ) { M_function = function; }

    //! Set the boundary flow rate type
    /*!
     * @param boundaryFlowRateType boundary flow rate type
     */
    void setBoundaryFlowRateType( const FSI3DBoundaryFlowRate_Type& boundaryFlowRateType ) { M_boundaryFlowRateType = boundaryFlowRateType; }

    //@}


    //! @name Get methods
    //@{

    //! Get the fluid flag of the boundary
    /*!
     * @return flag of the fluid boundary
     */
    const multiscaleID_Type& fluidFlag() const { return M_fluidFlag; }

    //! Get the boundary flow rate type
    /*!
     * @return boundary flow rate type
     */
    const FSI3DBoundaryFlowRate_Type& boundaryFlowRateType() const { return M_boundaryFlowRateType; }

    //@}


private:

    //! @name Unimplemented Methods
    //@{

    FSI3DBoundaryFlowRateFunction( const FSI3DBoundaryFlowRateFunction& boundaryFunction );

    FSI3DBoundaryFlowRateFunction& operator=( const FSI3DBoundaryFlowRateFunction& boundaryFunction );

    //@}

    //! @name Private Methods
    //@{

    //! Impose a flow in the outgoing normal direction.
    void useNormalDirectionFlow( const BCBase& boundaryID )
    {
        // Use the PostProcessingBoundary utility class to extract the surface normals
        PostProcessingBoundary<MultiscaleModelFSI3D::mesh_Type> normalExtraction( M_FSI3D->solver()->uFESpacePtr()->mesh(),
                                                                                  &(M_FSI3D->solver()->uFESpacePtr()->feBd()),
                                                                                  &(M_FSI3D->solver()->uFESpacePtr()->dof()),
                                                                                  M_FSI3D->solver()->uFESpacePtr()->map() );

        Vector approxNormal = normalExtraction.normal( boundaryID.flag() );

        // Take the first surface normal direction
        M_n[0] = approxNormal[0];
        M_n[1] = approxNormal[1];
        M_n[2] = approxNormal[2];
    }

    //@}

    MultiscaleModelFSI3D*                          M_FSI3D;
    function_Type                                  M_function;
    multiscaleID_Type                              M_fluidFlag;
    FSI3DBoundaryFlowRate_Type                     M_boundaryFlowRateType;
    boost::array< Real, 3 >                        M_n;
    Real                                           M_area;
    bool                                           M_flowRateIsZero;
};

#endif



#ifdef FSI_WITH_BOUNDARYAREA

//! FSI3DBoundaryAreaFunction - The FSI3D area function
/*!
 *  @author Cristiano Malossi
 *
 *  This simple class provides the implementation for the BC function used by the FSI3D model
 *  in order to apply the area of the fluid vessel at the boundaries.
 */
class FSI3DBoundaryAreaFunction
{
public:

    //! @name Type definitions
    //@{

    typedef MultiscaleInterface::function_Type function_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit FSI3DBoundaryAreaFunction() :
        M_FSI3D           (),
        M_fluidFlag       (),
        M_referenceArea   (),
        M_geometricCenter (),
        M_n               (),
        M_t1              (),
        M_t2              (),
        M_function        ()
        {}

    //! Destructor
    virtual ~FSI3DBoundaryAreaFunction() { /* M_FSI3D is deleted outside */ }

    //@}


    //! @name Methods
    //@{

    //! Evaluate the displacement of a point
    /*!
     * @return displacement of a point
     */
    Real function( const Real& t, const Real& x, const Real& y, const Real& z, const UInt& id )
    {
        return displacement( std::sqrt( M_function( t, x, y, z, id ) / M_referenceArea ) - 1, x , y, z, id );
    }

    //! Evaluate the displacement of a point when solving the tangent problem
    /*!
     * @return displacement of a point
     */
    Real functionLinear( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const UInt& id )
    {
        return displacement( std::sqrt( 1. / M_referenceArea ) - 1, x , y, z, id );
    }

    //! Setup main quantities
    void setup()
    {
        // Compute reference area
        M_referenceArea = M_FSI3D->solver()->fluid().area( M_fluidFlag );

        // Compute normal
        Vector normal = M_FSI3D->solver()->fluid().normal( M_fluidFlag );

        M_n[0] = normal[0];
        M_n[1] = normal[1];
        M_n[2] = normal[2];

        // Compute tangent
        setupTangent();

        // Compute geometric center
        Vector geometricCenter = M_FSI3D->solver()->fluid().geometricCenter( M_fluidFlag );

        M_geometricCenter[0] = geometricCenter[0];
        M_geometricCenter[1] = geometricCenter[1];
        M_geometricCenter[2] = geometricCenter[2];

        // ShowMe
        showMe();
    }

    //! ShowMe
    void showMe() const
    {
        if ( M_FSI3D->communicator()->MyPID() == 0 )
        {
            std::cout << "Fluid flag          = " << M_fluidFlag << std::endl
                      << "Reference area      = " << M_referenceArea << std::endl
                      << "Geometric center    = " << M_geometricCenter[0] << " " << M_geometricCenter[1] << " " << M_geometricCenter[2] << std::endl
                      << "Normal              = " << M_n[0] << " " << M_n[1] << " " << M_n[2] << std::endl
                      << "Tangent 1           = " << M_t1[0] << " " << M_t1[1] << " " << M_t1[2] << std::endl
                      << "Tangent 2           = " << M_t2[0] << " " << M_t2[1] << " " << M_t2[2] << std::endl << std::endl;
        }
    }

    //@}


    //! @name Set methods
    //@{

    //! Set the FSI3D model
    /*!
     * @param modelFSI3D a pointer to the FSI3D model
     */
    void setModel( const MultiscaleModelFSI3D* modelFSI3D ) { M_FSI3D = modelFSI3D; }

    //! Set the fluid flag of the boundary
    /*!
     * @param flag flag of the fluid boundary
     */
    void setFluidFlag( const multiscaleID_Type& flag ) { M_fluidFlag = flag; }

    //! Set the reference area the fluid boundary
    /*!
     * @param referenceArea reference area of the fluid boundary
     */
    void setReferenceArea( const Real& referenceArea ) { M_referenceArea = referenceArea; }

    //! Set the geometric center of the fluid boundary
    /*!
     * @param geometricCenter the x-y-z coordinate of the geometric center of the fluid boundary
     */
    void setGeometricCenter( const boost::array< Real, 3 >& geometricCenter ) { M_geometricCenter = geometricCenter; }

    //! Set the outgoing normal of the fluid boundary
    /*!
     * @param normal outgoing normal of the fluid boundary
     */
    void setNormal( const boost::array< Real, 3 >& normal ) { M_n = normal; setupTangent(); }

    //! Set the area function
    /*!
     * @param function area function
     */
    void setFunction( const function_Type& function ) { M_function = function; }

    //@}

    //! @name Get methods
    //@{

    //! Get the fluid flag of the boundary
    /*!
     * @return flag of the fluid boundary
     */
    const multiscaleID_Type& fluidFlag() const { return M_fluidFlag; }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    FSI3DBoundaryAreaFunction( const FSI3DBoundaryAreaFunction& boundaryFunction );

    FSI3DBoundaryAreaFunction& operator=( const FSI3DBoundaryAreaFunction& boundaryFunction );

    //@}


    //! @name Private Methods
    //@{

    //! Setup tangent vectors
    void setupTangent()
    {
        // Normalization
        Real module = std::sqrt( M_n[0] * M_n[0] + M_n[1] * M_n[1] + M_n[2] * M_n[2] );
        M_n[0] = M_n[0] / module;
        M_n[1] = M_n[1] / module;
        M_n[2] = M_n[2] / module;

        // Compute t1
        M_t1[0] = M_n[1] * M_n[1]  - M_n[0] * M_n[2];
        M_t1[1] = M_n[2] * M_n[2]  - M_n[0] * M_n[1];
        M_t1[2] = M_n[0] * M_n[0]  - M_n[1] * M_n[2];

        // Compute t2
        M_t2[0] = M_n[1] * M_t1[2] - M_n[2] * M_t1[1];
        M_t2[1] = M_n[2] * M_t1[0] - M_n[0] * M_t1[2];
        M_t2[2] = M_n[0] * M_t1[1] - M_n[1] * M_t1[0];
    }

    //! Evaluate the displacement of a point given the scale factor
    /*!
     * @param scale factor
     * @param x x-coordinate of the point
     * @param y y-coordinate of the point
     * @param z z-coordinate of the point
     * @param id id of the component
     * @return displacement of a point
     */
    Real displacement( const Real& scaleFactor, const Real& x, const Real& y, const Real& z, const UInt& id )
    {
        // Compute the RHS
        boost::array< Real, 3 > rhs;
        rhs[0] = 0;
        rhs[1] = scaleFactor * ( ( x - M_geometricCenter[0] ) * M_t1[0] + ( y - M_geometricCenter[1] ) * M_t1[1] + ( z - M_geometricCenter[2] ) * M_t1[2] );
        rhs[2] = scaleFactor * ( ( x - M_geometricCenter[0] ) * M_t2[0] + ( y - M_geometricCenter[1] ) * M_t2[1] + ( z - M_geometricCenter[2] ) * M_t2[2] );

        // Compute the displacement
        Real determinant = M_n[0] * ( M_t1[1] * M_t2[2] - M_t1[2] * M_t2[1] )
                         + M_n[1] * ( M_t1[2] * M_t2[0] - M_t1[0] * M_t2[2] )
                         + M_n[2] * ( M_t1[0] * M_t2[1] - M_t1[1] * M_t2[0] );
        switch ( id )
        {
        case 0:
            return -( rhs[0] * ( M_t1[2] * M_t2[1] - M_t1[1] * M_t2[2] )
                    + rhs[1] * ( M_n[1]  * M_t2[2] - M_n[2]  * M_t2[1] )
                    + rhs[2] * ( M_n[2]  * M_t1[1] - M_n[1]  * M_t1[2] ) ) / determinant;

        case 1:
            return  ( rhs[0] * ( M_t1[2] * M_t2[0] - M_t1[0] * M_t2[2] )
                    + rhs[1] * ( M_n[0]  * M_t2[2] - M_n[2]  * M_t2[0] )
                    + rhs[2] * ( M_n[2]  * M_t1[0] - M_n[0]  * M_t1[2] ) ) / determinant;

        case 2:
            return -( rhs[0] * ( M_t1[1] * M_t2[0] - M_t1[0] * M_t2[1] )
                    + rhs[1] * ( M_n[0]  * M_t2[1] - M_n[1]  * M_t2[0] )
                    + rhs[2] * ( M_n[1]  * M_t1[0] - M_n[0]  * M_t1[1] ) ) / determinant;

        default:

            return 0.;
        }
    }

    //@}


    const MultiscaleModelFSI3D*       M_FSI3D;             // A pointer to the model class
    multiscaleID_Type                 M_fluidFlag;         // Boundary flag
    Real                              M_referenceArea;     // Reference area
    boost::array< Real, 3 >           M_geometricCenter;   // Geometric center
    boost::array< Real, 3 >           M_n;                 // Normal direction
    boost::array< Real, 3 >           M_t1;                // Tangential direction 1
    boost::array< Real, 3 >           M_t2;                // Tangential direction 2
    function_Type                     M_function;          // Function scale factor
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
