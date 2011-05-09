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
 *  @brief File containing the Multiscale Physical Coupling
 *
 *  @date 02-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleCoupling_H
#define MultiscaleCoupling_H 1

#include <lifemc/lifesolver/MultiscaleDefinitions.hpp>
#include <lifemc/lifesolver/MultiscaleData.hpp>
#include <lifemc/lifesolver/MultiscaleModel.hpp>

namespace LifeV
{
namespace Multiscale
{

// Forward declaration
class MultiscaleCouplingFunction;

//! MultiscaleCoupling - The Multiscale Physical Coupling
/*!
 *  @author Cristiano Malossi
 *
 *  The MultiscaleCoupling class provides a general interface between the
 *  MS_Algorithm and all the coupling conditions.
 */
class MultiscaleCoupling
{
public:

    //! @name Type definitions
    //@{

    typedef MultiscaleCouplingFunction                           couplingFunction_Type;
    typedef boost::shared_ptr < couplingFunction_Type >          couplingFunctionPtr_Type;
    typedef std::vector< couplingFunction_Type >                 couplingFunctionsContainer_Type;

    typedef std::vector< multiscaleVectorPtr_Type >              couplingVariablesContainer_Type;

    typedef std::vector< Real >                                  timeContainer_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleCoupling();

    //! Destructor
    virtual ~MultiscaleCoupling() {}

    //@}


    //! @name Multiscale PhysicalCoupling Virtual Methods
    //@{

    //! Setup the data of the coupling.
    /*!
     * @param FileName Name of data file
     */
    virtual void setupData( const std::string& fileName );

    //! Setup the coupling
    virtual void setupCoupling() = 0;

    //! Initialize the values of the coupling variables
    virtual void initializeCouplingVariables() = 0;

    //! Update the coupling
    /*!
     * This method is the analogous of the "updateModel" for the models.
     * It is alternative to initializeCouplingVariables and is called from the second timestep.
     * It is reserved for the update of:
     * <ol>
     *     <li> objects that are not constant with respect to the time but should not be updated during subiterations.
     * </ol>
     */
    virtual void updateCoupling() = 0;

    //! Export the values of the local coupling residuals into a global vector
    /*!
     * @param couplingResiduals Global vector of variables
     */
    virtual void exportCouplingResiduals( multiscaleVector_Type& couplingResiduals ) = 0;

    //! Check if the topology is changed
    /*!
     * A topology change can be caused by a change in the coupling equations by,
     * for example, the opening/closure of a valve (see MultiscaleCouplingFlowRateValve).
     *
     * @return true if the topology is changed, false otherwise
     */
    virtual bool topologyChange() { return false; }

    //! Display some information about the coupling
    virtual void showMe();

    //@}


    //! @name Methods
    //@{

    //! Build the global map for the coupling vectors
    /*!
     * @param couplingMap Global coupling map
     */
    void createCouplingMap( MapEpetra& couplingMap );

    //! Import the values of the coupling variables
    /*!
     * @param couplingVariables Global vector of coupling variables
     */
    void importCouplingVariables( const multiscaleVector_Type& couplingVariables ) { importCouplingVector( couplingVariables, localCouplingVariables( 0 ) ); }

    //! Export the values of the coupling variables
    /*!
     * @param couplingVariables Global vector of coupling variables
     */
    void exportCouplingVariables( multiscaleVector_Type& couplingVariables ) { exportCouplingVector( localCouplingVariables( 0 ), couplingVariables ); }

    //! Extrapolate the values of the coupling variables for the next time step
    void extrapolateCouplingVariables();

    //! Lagrange interpolation/extrapolation of the coupling variables at selected time.
    /*!
     * @param t interpolation time
     * @param interpolatedCouplingVariables variables interpolated/extrapolated at time t
     */
    void interpolateCouplingVariables( const Real& t, multiscaleVector_Type& interpolatedCouplingVariables ) const;

    //! Find if a perturbation is imposed on the coupling.
    /*!
     * @return true if a perturbation is imposed
     */
    bool isPerturbed() const { return M_perturbedCoupling == -1 ? false : true; }

    //! Export the Jacobian matrix
    /*!
     * @param Jacobian Jacobian Matrix
     */
    void exportJacobian( multiscaleMatrix_Type& jacobian );

    //! save the coupling variables information on a file
    void saveSolution();

    //! Clear the list of pointers to the models.
    /*!
     *  This method has to be called before the automatic destructor, in order
     *  to disconnect the coupling classes from the model classes.
     */
    void clearModelsList() { M_models.clear(); }

    //@}


    //! @name Set Methods
    //@{

    //! Set the global ID of the coupling
    /*!
     * @param id Coupling ID
     */
    void setID( const UInt& id ) { M_ID = id; }

    //! Add a pointer to one of the models to couple
    /*!
     * @param model shared_ptr of the model
     */
    void addModel( const multiscaleModelPtr_Type& model ) { M_models.push_back( model ); }

    //! Add a flag of one of the models to couple
    /*!
     * @param flag flag of the model
     */
    void addFlag( const bcFlag_Type& flag ) { M_flags.push_back( flag ); }

    //! Add a flag of one of the models to couple
    /*!
     * @param flagID get from the model the flag with this flagID
     */
    void addFlagID( const UInt& flagID );

    //! Setup the global data of the coupling.
    /*!
     * In particular, it can be used to replace the local values specified in
     * the model data file, with the ones in the global container.
     *
     * @param globalData Global data container.
     */
    void setGlobalData( const multiscaleDataPtr_Type& globalData ) { M_globalData = globalData; }

    //! Set the epetra communicator for the coupling
    /*!
     * @param comm Epetra communicator
     */
    void setCommunicator( const multiscaleCommPtr_Type& comm ) { M_comm = comm; }

    //@}


    //! @name Get Methods
    //@{

    //! Get the global ID of the coupling
    /*!
     * @return global ID of the coupling
     */
    const UInt& ID() const { return M_ID; }

    //! Get the type of the coupling
    /*!
     * @return type of the coupling
     */
    const couplings_Type& type() const { return M_type; }

    //! Get the name of the coupling
    /*!
     * @return name of the coupling
     */
    const std::string& couplingName() const { return M_couplingName; }

    //! Get the number of models connected by the coupling
    UInt modelsNumber() const { return static_cast< UInt > ( M_models.size() ); }

    //! Get the model local ID through global ID
    /*!
     * @param ID global ID of the model
     * @return local ID of the model
     */
    UInt modelGlobalToLocalID( const UInt& ID ) const;

    //! Get the model connected by the coupling through local ID
    /*!
     * @param LocalID local ID of the model
     * @return Pointer to the model
     */
    multiscaleModelPtr_Type model( const UInt& localID ) const { return M_models[localID]; }

    //! Get the model connected by the coupling through local ID
    /*!
     * @param LocalID local ID of the model
     * @return Coupling flag of the model
     */
    const bcFlag_Type& flag( const UInt& localID ) const { return M_flags[localID]; }

    //! Get the number of the coupling variables
    /*!
     * @return number of the coupling variables
     */
    const UInt& couplingVariablesNumber() const { return M_couplingIndex.first; }

    //! Get the container of the local coupling variables
    /*!
     * @return container of the local coupling variables
     */
    const couplingVariablesContainer_Type& couplingVariables() const { return M_localCouplingVariables; }

    //! Get the perturbed coupling.
    /*!
     * If it is unperturbed it returns -1.
     * @return the localID of the perturbed coupling
     */
    const Int& perturbedCoupling() const { return M_perturbedCoupling; }

    //! Get the local residual.
    /*!
     * @return the local residual of the coupling
     */
    const multiscaleVector_Type& residual() const { return *M_localCouplingResiduals; }

    //! Get the time interpolation order.
    /*!
     * @return the value of the time interpolation order.
     */
    const UInt& timeInterpolationOrder() const { return M_timeInterpolationOrder; }


    //@}

protected:

    //! @name Protected MultiscaleCoupling Virtual Methods
    //@{

    //! Build the list of models affected by the perturbation of a local coupling variable
    /*!
     * @param LocalCouplingVariableID local coupling variable (perturbed)
     * @return list of models affected by the perturbation
     */
    virtual multiscaleModelsVector_Type listOfPerturbedModels( const UInt& localCouplingVariableID ) = 0;

    //! Insert constant coefficients into the Jacobian matrix
    /*!
     * @param Jacobian the Jacobian matrix
     */
    virtual void insertJacobianConstantCoefficients( multiscaleMatrix_Type& jacobian ) = 0;

    //! Insert the Jacobian coefficient(s) depending on a perturbation of the model, due to a specific variable (the column)
    /*!
     * @param Jacobian the Jacobian matrix
     * @param Column the column related to the perturbed variable
     * @param ID the global ID of the model which is perturbed by the variable
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     */
    virtual void insertJacobianDeltaCoefficients( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& linearSystemSolved ) = 0;

    //! Display some information about the coupling
    /*!
     * @param output specify the output stream
     */
    virtual void displayCouplingValues( std::ostream& output ) = 0;

    //@}


    //! @name Protected Methods
    //@{

    //! Access by reference to a specific local coupling variable
    /*!
     * This method is used to simplify the access to a specific local coupling variables vector.
     * Note that the returned value is not const!
     * @param id id of the local coupling variables vector
     * @return reference to the local coupling variables vector
     */
    const multiscaleVector_Type& localCouplingVariables( const UInt& id ) const { return *M_localCouplingVariables[id]; }
          multiscaleVector_Type& localCouplingVariables( const UInt& id )       { return *M_localCouplingVariables[id]; }

    //! Create the local vectors of the coupling
    void createLocalVectors();

    //! Reset the history of the couplings
    /*!
     *  This method is used when the topology change.
     */
    void resetCouplingHistory();

    //! Import the content of the Global Vector in the Local vector
    /*!
     * @param globalVector the global vector
     * @param localVector the local vector
     */
    void importCouplingVector( const multiscaleVector_Type& globalVector, multiscaleVector_Type& localVector );

    //! Export the content of the Local Vector in the Global vector
    /*!
     *
     * @param localVector the local vector
     * @param globalVector the global vector
     */
    void exportCouplingVector( const multiscaleVector_Type& localVector, multiscaleVector_Type& globalVector );

    //! Display and error message for the specific model
    /*!
     * @param model shared_ptr to the specific model
     */
    void switchErrorMessage( const multiscaleModelPtr_Type& model );

    //@}

    static UInt                          M_couplingsNumber;

    UInt                                 M_ID;
    couplings_Type                       M_type;

    multiscaleModelsVector_Type          M_models;
    std::string                          M_couplingName;
    std::vector< bcFlag_Type >           M_flags;

    multiscaleDataPtr_Type               M_globalData;

    std::pair< UInt, UInt >              M_couplingIndex;

    couplingFunctionsContainer_Type      M_localCouplingFunctions;
    couplingVariablesContainer_Type      M_localCouplingVariables;
    multiscaleVectorPtr_Type             M_localCouplingResiduals;

    UInt                                 M_timeInterpolationOrder;

    Int                                  M_perturbedCoupling;

    multiscaleCommPtr_Type               M_comm;

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleCoupling( const MultiscaleCoupling& coupling );

    MultiscaleCoupling& operator=( const MultiscaleCoupling& coupling );

    //@}
};

//! MultiscaleCouplingFunction - The multiscale function for the couplings
/*!
 *  @author Cristiano Malossi
 *
 *  This simple class provides the implementation for the BC function used by the couplings.
 */
class MultiscaleCouplingFunction
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleCouplingFunction() : M_coupling(), M_id() {}

    //! Constructor
    /*!
     * @param coupling pointer to the coupling
     * @param id id of the coupling variable
     */
    explicit MultiscaleCouplingFunction( const multiscaleCoupling_Type* coupling, const UInt& id) : M_coupling( coupling ), M_id( id ) {}

    //! Destructor
    virtual ~MultiscaleCouplingFunction() { /* M_coupling is deleted outside */ }

    //@}


    //! @name Methods
    //@{

    //! Evaluate the coupling quantity
    /*!
      @return evaluation of the function
    */
    Real function( const Real& t, const Real&, const Real&, const Real&, const UInt& )
    {
        multiscaleVector_Type interpolatedCouplingVariables( *M_coupling->couplingVariables()[0] );
        M_coupling->interpolateCouplingVariables( t, interpolatedCouplingVariables );

        return interpolatedCouplingVariables[M_id];
    }

    //@}

private:

    const multiscaleCoupling_Type*         M_coupling;
    UInt                                   M_id;
};

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleCoupling_H */
