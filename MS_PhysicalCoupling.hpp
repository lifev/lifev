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
 *  @brief MultiScale Physical Coupling
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 02-09-2009
 */

#ifndef MS_PhysicalCoupling_H
#define MS_PhysicalCoupling_H 1

#include <lifemc/lifesolver/MS_Definitions.hpp>
#include <lifemc/lifesolver/MS_PhysicalData.hpp>
#include <lifemc/lifesolver/MS_PhysicalModel.hpp>

namespace LifeV {

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

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_PhysicalCoupling();

    //! Copy constructor
    /*!
     * @param coupling MS_PhysicalCoupling
     */
    MS_PhysicalCoupling( const MS_PhysicalCoupling& coupling );

    //! Destructor
    virtual ~MS_PhysicalCoupling() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param coupling MS_PhysicalCoupling
     * @return reference to a copy of the class
     */
    MS_PhysicalCoupling& operator=( const MS_PhysicalCoupling& coupling );

    //@}


    //! @name MultiScale PhysicalCoupling Virtual Methods
    //@{

    //! Setup the data of the coupling
    virtual void SetupData() = 0;

    //! Setup the coupling
    virtual void SetupCoupling() = 0;

    //! Initialize the values of the coupling variables
    virtual void InitializeCouplingVariables() = 0;

    //! Export the values of the local coupling residuals into a global vector
    /*!
     * @param CouplingResiduals Global vector of variables
     */
    virtual void ExportCouplingResiduals( VectorType& CouplingResiduals ) = 0;

    //! Build the list of models affected by the perturbation of a local coupling variable
    /*!
     * @param LocalCouplingVariableID local coupling variable (perturbed)
     * @return list of models affected by the perturbation
     */
    virtual ModelsVector_Type GetListOfPerturbedModels( const UInt& LocalCouplingVariableID ) = 0;

    //! Insert constant coefficients into the Jacobian matrix
    /*!
     * @param Jacobian the Jacobian matrix
     */
    virtual void InsertJacobianConstantCoefficients( MatrixType& Jacobian ) = 0;

    //! Insert the Jacobian coefficient(s) depending on a perturbation of the model, due to a specific variable (the column)
    /*!
     * @param Jacobian the Jacobian matrix
     * @param Column the column related to the perturbed variable
     * @param ID the global ID of the model which is perturbed by the variable
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     */
    virtual void InsertJacobianDeltaCoefficients( MatrixType& Jacobian, const UInt& Column, const UInt& ID, bool& LinearSystemSolved ) = 0;

    //! Display some information about the coupling
    /*!
     * @param output specify the output stream
     */
    virtual void DisplayCouplingValues( std::ostream& output ) = 0;

    //! Display some information about the coupling
    virtual void ShowMe();

    //@}


    //! @name Methods
    //@{

    //! Build the global map for the coupling vectors
    /*!
     * @param couplingMap Global coupling map
     */
    void CreateCouplingMap( EpetraMap& couplingMap );

    //! Import the values of the coupling variables
    /*!
     * @param CouplingVariables Global vector of coupling variables
     */
    void ImportCouplingVariables( const VectorType& CouplingVariables );

    //! Export the values of the coupling variables
    /*!
     * @param CouplingVariables Global vector of coupling variables
     */
    void ExportCouplingVariables( VectorType& CouplingVariables );

    //! Export the Jacobian matrix
    /*!
     * @param Jacobian Jacobian Matrix
     */
    void ExportJacobian( MatrixType& Jacobian );

    //! Export the values of the Jacobian product
    /*!
     * @param deltaCouplingVariables variation of the coupling variables
     * @param JacobianProduct the product of the Jacobian by the varuatuib if tge coupling variables
     */
    void ExportJacobianProduct( const VectorType& deltaCouplingVariables, VectorType& JacobianProduct );

    //! Save the coupling variables information on a file
    void SaveSolution();

    //! Clear the list of pointers to the models.
    /*!
     *  This method has to be called before the automatic destructor, in order
     *  to disconnect the coupling classes from the model classes.
     */
    void ClearModelsList();

    //@}


    //! @name Set Methods
    //@{

    //! Set the global ID of the coupling
    /*!
     * @param id Coupling ID
     */
    void SetID( const UInt& id );

    //! Set the data file to load information of the coupling condition
    /*!
     * @param dataFile Name and path of data file
     */
    void SetDataFile( const std::string& dataFile );

    //! Add a pointer to one of the models to couple
    /*!
     * @param model shared_ptr of the model
     */
    void AddModel( const Model_ptrType& model );

    //! Set global data for physical quantities and time
    /*!
     * @param dataPhysics Data container for physical quantities
     * @param dataTime Data container for time parameters
     */
    void SetData( const boost::shared_ptr< MS_PhysicalData >& dataPhysics,
                  const boost::shared_ptr< DataTime >& dataTime );

    //! Add a flag of one of the models to couple
    /*!
     * @param flag flag of the model
     */
    void AddFlag( const BCFlag& flag );

    //! Add a flag of one of the models to couple
    /*!
     * @param flagID get from the model the flag with this flagID
     */
    void AddFlagID( const UInt& flagID );

    //! Set the epetra communicator for the coupling
    /*!
     * @param comm Epetra communicator
     */
    void SetCommunicator( const boost::shared_ptr< Epetra_Comm >& comm );

    //@}


    //! @name Get Methods
    //@{

    //! Get the global ID of the coupling
    /*!
     * @return global ID of the coupling
     */
    const UInt& GetID() const;

    //! Get the type of the coupling
    /*!
     * @return type of the coupling
     */
    const couplingsTypes& GetType() const;

    //! Get the name of the coupling
    /*!
     * @return name of the coupling
     */
    const std::string& GetCouplingName() const;

    //! Get the number of models connected by the coupling
    UInt GetModelsNumber() const;

    //! Get the model local ID through global ID
    /*!
     * @param ID global ID of the model
     * @return local ID of the model
     */
    UInt GetModelLocalID( const UInt& ID ) const;

    //! Get the model connected by the coupling through local ID
    /*!
     * @param ID local ID of the model
     * @return Pointer to the model
     */
    Model_ptrType GetModel( const UInt& LocalID ) const;

    //! Get the number of the coupling variables
    /*!
     * @return number of the coupling variables
     */
    const UInt& GetCouplingVariablesNumber() const;

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Create the local vectors of the coupling
    void CreateLocalVectors();

    //! Import the content of the Global Vector in the Local vector
    /*!
     * @param globalVector the global vector
     * @param localVector the local vector
     */
    void ImportCouplingVector( const VectorType& globalVector, VectorType& localVector );

    //! Export the content of the Local Vector in the Global vector
    /*!
     *
     * @param localVector the local vector
     * @param globalVector the global vector
     */
    void ExportCouplingVector( const VectorType& localVector, VectorType& globalVector );

    //! Display and error message for the specific model
    /*!
     * @param model shared_ptr to the specific model
     */
    void switchErrorMessage( const Model_ptrType& model );

    //@}

    static UInt                          M_couplingsNumber;

    UInt                                 M_ID;
    couplingsTypes                       M_type;

    GetPot                               M_dataFile;
    ModelsVector_Type                    M_models;
    std::string                          M_couplingName;
    std::vector< BCFlag >                M_flags;

    std::pair< UInt, UInt >              M_couplingIndex;

    Vector_ptrType                       M_LocalCouplingVariables;
    Vector_ptrType                       M_LocalCouplingResiduals;
    Vector_ptrType                       M_LocalDeltaCouplingVariables;

    boost::shared_ptr< MS_PhysicalData > M_dataPhysics;
    boost::shared_ptr< DataTime >        M_dataTime;

    boost::shared_ptr< Epetra_Comm >     M_comm;
    boost::shared_ptr< Displayer >       M_displayer;
};

} // Namespace LifeV

#endif /* MS_PhysicalCoupling_H */
