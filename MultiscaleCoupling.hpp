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
 *  @brief File containing the MultiScale Physical Coupling
 *
 *  @date 02-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MS_PhysicalCoupling_H
#define MS_PhysicalCoupling_H 1

#include <lifemc/lifesolver/MS_Definitions.hpp>
#include <lifemc/lifesolver/MS_PhysicalData.hpp>
#include <lifemc/lifesolver/MS_PhysicalModel.hpp>

namespace LifeV
{

//! MS_PhysicalCoupling - The MultiScale Physical Coupling
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_PhysicalCoupling class provides a general interface between the
 *  MS_Algorithm and all the coupling conditions.
 */
class MS_PhysicalCoupling
{
public:

    //! @name Type definitions
    //@{

    typedef std::vector< MS_Vector_PtrType >                     couplingVariablesContainer_Type;
    typedef std::vector< Real >                                  timeContainer_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_PhysicalCoupling();

    //! Destructor
    virtual ~MS_PhysicalCoupling() {}

    //@}


    //! @name MultiScale PhysicalCoupling Virtual Methods
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

    //! Export the values of the local coupling residuals into a global vector
    /*!
     * @param couplingResiduals Global vector of variables
     */
    virtual void exportCouplingResiduals( MS_Vector_Type& couplingResiduals ) = 0;

    //! Display some information about the coupling
    virtual void showMe();

    //@}


    //! @name Methods
    //@{

    //! Build the global map for the coupling vectors
    /*!
     * @param couplingMap Global coupling map
     */
    void createCouplingMap( EpetraMap& couplingMap );

    //! Import the values of the coupling variables
    /*!
     * @param couplingVariables Global vector of coupling variables
     */
    void importCouplingVariables( const MS_Vector_Type& couplingVariables );

    //! Export the values of the coupling variables
    /*!
     * @param couplingVariables Global vector of coupling variables
     */
    void exportCouplingVariables( MS_Vector_Type& couplingVariables );

    //! Extrapolate the values of the coupling variables for the next time step
    void extrapolateCouplingVariables();

    //! Find if a perturbation is imposed on the coupling.
    /*!
     * @return true if a perturbation is imposed
     */
    bool isPerturbed() const;

    //! Export the Jacobian matrix
    /*!
     * @param Jacobian Jacobian Matrix
     */
    void exportJacobian( MS_Matrix_Type& jacobian );

    //! Export the values of the Jacobian product
    /*!
     * @param deltaCouplingVariables variation of the coupling variables
     * @param jacobianProduct the product of the Jacobian by the varuatuib if tge coupling variables
     */
    void exportJacobianProduct( const MS_Vector_Type& deltaCouplingVariables, MS_Vector_Type& jacobianProduct );

    //! save the coupling variables information on a file
    void saveSolution();

    //! Clear the list of pointers to the models.
    /*!
     *  This method has to be called before the automatic destructor, in order
     *  to disconnect the coupling classes from the model classes.
     */
    void clearModelsList();

    //@}


    //! @name Set Methods
    //@{

    //! Set the global ID of the coupling
    /*!
     * @param id Coupling ID
     */
    void setID( const UInt& id );

    //! Add a pointer to one of the models to couple
    /*!
     * @param model shared_ptr of the model
     */
    void addModel( const MS_Model_PtrType& model );

    //! Add a flag of one of the models to couple
    /*!
     * @param flag flag of the model
     */
    void addFlag( const BCFlag& flag );

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
    void setGlobalData( const MS_GlobalDataContainer_PtrType& globalData );

    //! Set the epetra communicator for the coupling
    /*!
     * @param comm Epetra communicator
     */
    void setCommunicator( const MS_Comm_PtrType& comm );

    //@}


    //! @name Get Methods
    //@{

    //! Get the global ID of the coupling
    /*!
     * @return global ID of the coupling
     */
    const UInt& ID() const;

    //! Get the type of the coupling
    /*!
     * @return type of the coupling
     */
    const couplings_Type& type() const;

    //! Get the name of the coupling
    /*!
     * @return name of the coupling
     */
    const std::string& couplingName() const;

    //! Get the number of models connected by the coupling
    UInt modelsNumber() const;

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
    MS_Model_PtrType model( const UInt& localID ) const;

    //! Get the model connected by the coupling through local ID
    /*!
     * @param LocalID local ID of the model
     * @return Coupling flag of the model
     */
    const BCFlag& flag( const UInt& localID ) const;

    //! Get the number of the coupling variables
    /*!
     * @return number of the coupling variables
     */
    const UInt& couplingVariablesNumber() const;

    //! Get the perturbed coupling.
    /*!
     * If it is unperturbed it returns -1.
     * @return the localID of the perturbed coupling
     */
    const Int& perturbedCoupling() const;

    //! Get the local residual.
    /*!
     * @return the local residual of the coupling
     */
    const MS_Vector_Type& residual() const;

    //! Get the time interpolation order.
    /*!
     * @return the value of the time interpolation order.
     */
    const UInt& timeInterpolationOrder() const;


    //@}

protected:

    //! @name Protected MultiScale PhysicalCoupling Virtual Methods
    //@{

    //! Build the list of models affected by the perturbation of a local coupling variable
    /*!
     * @param LocalCouplingVariableID local coupling variable (perturbed)
     * @return list of models affected by the perturbation
     */
    virtual MS_ModelsVector_Type listOfPerturbedModels( const UInt& localCouplingVariableID ) = 0;

    //! Insert constant coefficients into the Jacobian matrix
    /*!
     * @param Jacobian the Jacobian matrix
     */
    virtual void insertJacobianConstantCoefficients( MS_Matrix_Type& jacobian ) = 0;

    //! Insert the Jacobian coefficient(s) depending on a perturbation of the model, due to a specific variable (the column)
    /*!
     * @param Jacobian the Jacobian matrix
     * @param Column the column related to the perturbed variable
     * @param ID the global ID of the model which is perturbed by the variable
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     */
    virtual void insertJacobianDeltaCoefficients( MS_Matrix_Type& jacobian, const UInt& column, const UInt& ID, bool& linearSystemSolved ) = 0;

    //! Display some information about the coupling
    /*!
     * @param output specify the output stream
     */
    virtual void displayCouplingValues( std::ostream& output ) = 0;

    //@}


    //! @name Protected Methods
    //@{

    //! Create the local vectors of the coupling
    void createLocalVectors();

    //! Import the content of the Global Vector in the Local vector
    /*!
     * @param globalVector the global vector
     * @param localVector the local vector
     */
    void importCouplingVector( const MS_Vector_Type& globalVector, MS_Vector_Type& localVector );

    //! Export the content of the Local Vector in the Global vector
    /*!
     *
     * @param localVector the local vector
     * @param globalVector the global vector
     */
    void exportCouplingVector( const MS_Vector_Type& localVector, MS_Vector_Type& globalVector );

    //! Lagrange interpolation/extrapolation of the coupling variables at selected time.
    /*!
     * @param timeContainer vector of times
     * @param t interpolation time
     * @param interpolatedCouplingVariables variables interpolated/extrapolated at time t
     */
    void interpolateCouplingVariables( const timeContainer_Type& timeContainer, const Real& t,
                                       MS_Vector_Type& interpolatedCouplingVariables );

    //! Display and error message for the specific model
    /*!
     * @param model shared_ptr to the specific model
     */
    void switchErrorMessage( const MS_Model_PtrType& model );

    //@}

    static UInt                          M_couplingsNumber;

    UInt                                 M_ID;
    couplings_Type                       M_type;

    MS_ModelsVector_Type                 M_models;
    std::string                          M_couplingName;
    std::vector< BCFlag >                M_flags;

    MS_GlobalDataContainer_PtrType       M_globalData;

    std::pair< UInt, UInt >              M_couplingIndex;

    couplingVariablesContainer_Type      M_localCouplingVariables;
    MS_Vector_PtrType                    M_localCouplingResiduals;

    UInt                                 M_timeInterpolationOrder;

    Int                                  M_perturbedCoupling;

    MS_Comm_PtrType                      M_comm;
    boost::shared_ptr< Displayer >       M_displayer;

private:

    //! @name Unimplemented Methods
    //@{

    MS_PhysicalCoupling( const MS_PhysicalCoupling& coupling );

    MS_PhysicalCoupling& operator=( const MS_PhysicalCoupling& coupling );

    //@}
};

} // Namespace LifeV

#endif /* MS_PhysicalCoupling_H */
