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

    typedef std::vector< MS_Vector_PtrType >                     CouplingVariablesContainer_Type;
    typedef std::vector< Real >                                  TimeContainer_Type;

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
    virtual void SetupData( const std::string& FileName );

    //! Setup the coupling
    virtual void SetupCoupling() = 0;

    //! Initialize the values of the coupling variables
    virtual void InitializeCouplingVariables() = 0;

    //! Export the values of the local coupling residuals into a global vector
    /*!
     * @param CouplingResiduals Global vector of variables
     */
    virtual void ExportCouplingResiduals( MS_Vector_Type& CouplingResiduals ) = 0;

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
    void ImportCouplingVariables( const MS_Vector_Type& CouplingVariables );

    //! Export the values of the coupling variables
    /*!
     * @param CouplingVariables Global vector of coupling variables
     */
    void ExportCouplingVariables( MS_Vector_Type& CouplingVariables );

    //! Extrapolate the values of the coupling variables for the next time step
    void ExtrapolateCouplingVariables();

    //! Find if a perturbation is imposed on the coupling.
    /*!
     * @return true if a perturbation is imposed
     */
    bool IsPerturbed() const;

    //! Export the Jacobian matrix
    /*!
     * @param Jacobian Jacobian Matrix
     */
    void ExportJacobian( MS_Matrix_Type& Jacobian );

    //! Export the values of the Jacobian product
    /*!
     * @param deltaCouplingVariables variation of the coupling variables
     * @param JacobianProduct the product of the Jacobian by the varuatuib if tge coupling variables
     */
    void ExportJacobianProduct( const MS_Vector_Type& deltaCouplingVariables, MS_Vector_Type& JacobianProduct );

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

    //! Add a pointer to one of the models to couple
    /*!
     * @param model shared_ptr of the model
     */
    void AddModel( const MS_Model_PtrType& model );

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

    //! Setup the global data of the coupling.
    /*!
     * In particular, it can be used to replace the local values specified in
     * the model data file, with the ones in the global container.
     *
     * @param globalData Global data container.
     */
    void SetGlobalData( const MS_GlobalDataContainer_PtrType& globalData );

    //! Set the epetra communicator for the coupling
    /*!
     * @param comm Epetra communicator
     */
    void SetCommunicator( const MS_Comm_PtrType& comm );

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
     * @param LocalID local ID of the model
     * @return Pointer to the model
     */
    MS_Model_PtrType GetModel( const UInt& LocalID ) const;

    //! Get the model connected by the coupling through local ID
    /*!
     * @param LocalID local ID of the model
     * @return Coupling flag of the model
     */
    const BCFlag& GetFlag( const UInt& LocalID ) const;

    //! Get the number of the coupling variables
    /*!
     * @return number of the coupling variables
     */
    const UInt& GetCouplingVariablesNumber() const;

    //! Get the perturbed coupling.
    /*!
     * If it is unperturbed it returns -1.
     * @return the localID of the perturbed coupling
     */
    const Int& GetPerturbedCoupling() const;

    //! Get the local residual.
    /*!
     * @return the local residual of the coupling
     */
    const MS_Vector_Type& GetResidual() const;

    //! Get the time interpolation order.
    /*!
     * @return the value of the time interpolation order.
     */
    const UInt& GetTimeInterpolationOrder() const;


    //@}

protected:

    //! @name Protected MultiScale PhysicalCoupling Virtual Methods
    //@{

    //! Build the list of models affected by the perturbation of a local coupling variable
    /*!
     * @param LocalCouplingVariableID local coupling variable (perturbed)
     * @return list of models affected by the perturbation
     */
    virtual MS_ModelsVector_Type GetListOfPerturbedModels( const UInt& LocalCouplingVariableID ) = 0;

    //! Insert constant coefficients into the Jacobian matrix
    /*!
     * @param Jacobian the Jacobian matrix
     */
    virtual void InsertJacobianConstantCoefficients( MS_Matrix_Type& Jacobian ) = 0;

    //! Insert the Jacobian coefficient(s) depending on a perturbation of the model, due to a specific variable (the column)
    /*!
     * @param Jacobian the Jacobian matrix
     * @param Column the column related to the perturbed variable
     * @param ID the global ID of the model which is perturbed by the variable
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     */
    virtual void InsertJacobianDeltaCoefficients( MS_Matrix_Type& Jacobian, const UInt& Column, const UInt& ID, bool& LinearSystemSolved ) = 0;

    //! Display some information about the coupling
    /*!
     * @param output specify the output stream
     */
    virtual void DisplayCouplingValues( std::ostream& output ) = 0;

    //@}


    //! @name Protected Methods
    //@{

    //! Create the local vectors of the coupling
    void CreateLocalVectors();

    //! Import the content of the Global Vector in the Local vector
    /*!
     * @param globalVector the global vector
     * @param localVector the local vector
     */
    void ImportCouplingVector( const MS_Vector_Type& globalVector, MS_Vector_Type& localVector );

    //! Export the content of the Local Vector in the Global vector
    /*!
     *
     * @param localVector the local vector
     * @param globalVector the global vector
     */
    void ExportCouplingVector( const MS_Vector_Type& localVector, MS_Vector_Type& globalVector );

    //! Lagrange interpolation/extrapolation of the coupling variables at selected time.
    /*!
     * @param timeContainer vector of times
     * @param t interpolation time
     * @param interpolatedCouplingVariables variables interpolated/extrapolated at time t
     */
    void InterpolateCouplingVariables( const TimeContainer_Type& timeContainer,
                                       const Real& t,
                                       MS_Vector_Type& interpolatedCouplingVariables );

    //! Display and error message for the specific model
    /*!
     * @param model shared_ptr to the specific model
     */
    void switchErrorMessage( const MS_Model_PtrType& model );

    //@}

    static UInt                          M_couplingsNumber;

    UInt                                 M_ID;
    couplingsTypes                       M_type;

    MS_ModelsVector_Type                 M_models;
    std::string                          M_couplingName;
    std::vector< BCFlag >                M_flags;

    MS_GlobalDataContainer_PtrType       M_globalData;

    std::pair< UInt, UInt >              M_couplingIndex;

    CouplingVariablesContainer_Type      M_LocalCouplingVariables;
    MS_Vector_PtrType                    M_LocalCouplingResiduals;

    UInt                                 M_timeInterpolationOrder;

    Int                                  M_perturbedCoupling;

    MS_Comm_PtrType                      M_comm;
    boost::shared_ptr< Displayer >       M_displayer;
};

} // Namespace LifeV

#endif /* MS_PhysicalCoupling_H */
