/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-09-02

 Copyright (C) 2009 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/**
 \file MS_PhysicalCoupling.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-09-02
 */

#ifndef __MS_PhysicalCoupling_H
#define __MS_PhysicalCoupling_H 1

#include "Epetra_ConfigDefs.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <lifemc/lifesolver/MS_PhysicalModel.hpp>

namespace LifeV {

//! @name MS_PhysicalCoupling global objects
//@{

enum couplingsTypes
{
    BoundaryCondition,
    Stress,
    FluxStress
};

extern std::map< std::string, couplingsTypes > couplingsMap;

//@}

//! MS_PhysicalCoupling - The MultiScale Physical Coupling
/*!
 *  The MS_PhysicalCoupling class provides a general interface between the
 *  MS_Algorithm and all the coupling conditions.
 *
 *  @author Cristiano Malossi
 */
class MS_PhysicalCoupling
{
public:

    typedef MS_PhysicalModel::VectorType          VectorType;
    typedef MS_PhysicalModel::Vector_ptrType      Vector_ptrType;

    typedef boost::shared_ptr< MS_PhysicalModel > PhysicalModel_ptr;
    typedef EntityFlag                            BCFlag;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_PhysicalCoupling();

    //! Copy constructor
    /*!
     * \param coupling - MS_PhysicalCoupling
     */
    MS_PhysicalCoupling( const MS_PhysicalCoupling& coupling );

    //! Destructor
    virtual ~MS_PhysicalCoupling() {}

    //@}


    //! @name Methods
    //@{

    //! Operator=
    /*!
     * \param coupling - MS_PhysicalCoupling
     */
    MS_PhysicalCoupling& operator=( const MS_PhysicalCoupling& coupling );

    //! Build the global map for the coupling vectors
    /*!
     * \param couplingVariablesGlobalNumber - Global number of coupling variables
     * \param MyGlobalElements - Epetra_Map MyGlobalElements
     */
    void BuildCouplingVectorsMap( UInt&             couplingVariablesGlobalNumber,
                                  std::vector<int>& MyGlobalElements );

    //! Import the values of the coupling variables
    void ImportCouplingVariables( const VectorType& CouplingVariables );

    //! Export the values of the coupling variables
    void ExportCouplingVariables( VectorType& CouplingVariables );

    //! Export the values of the Jacobian product
    /*!
     * \param deltaCouplingVariables - variation of the coupling variables
     */
    void ExportJacobianProduct( const VectorType& deltaCouplingVariables, VectorType& JacobianProduct );

    //@}


    //! @name Set Methods
    //@{

    //! Set the global ID of the coupling
    /*!
     * \param id - Coupling ID
     */
    void SetID( const UInt& id )
    {
        M_ID = id;
    }

    //! Set the data file to load information of the coupling condition
    /*!
     * \param dataFile - Name and path of data file
     */
    void SetDataFile( const std::string& dataFile );

    //! Add a pointer to one of the models to couple
    /*!
     * \param model - shared_ptr of the model
     */
    void AddModel( const PhysicalModel_ptr& model )
    {
        M_models.push_back( model );
    }

    //! Add a flag of one of the models to couple
    /*!
     * \param flag - flag of the model
     */
    void AddFlag( const BCFlag& flag )
    {
        M_flags.push_back( flag );
    }

    //! Add a flag of one of the models to couple
    /*!
     * \param flagID - get from the model the flag with this flagID
     */
    void AddFlagID( const UInt& flagID )
    {
        M_flags.push_back( M_models.back()->GetFlag( flagID ) );
    }

    //! Set the epetra communicator for the coupling
    /*!
     * \param comm - Epetra communicator
     */
    void SetCommunicator( const boost::shared_ptr< Epetra_Comm >& comm );

    //@}


    //! @name Get Methods
    //@{

    //! Get the global ID of the coupling
    const UInt& GetID() const
    {
        return M_ID;
    }

    //! Get the type of the coupling
    const couplingsTypes& GetType() const
    {
        return M_type;
    }

    //! Get the name of the coupling
    const std::string& GetCouplingName() const
    {
        return M_couplingName;
    }

    //! Get the number of the coupling variables
    const UInt& GetCouplingVariablesNumber() const
    {
        return M_couplingIndex.first;
    }

    //@}


    //! @name MultiScale PhysicalCoupling Virtual Functions
    //@{

    //! Setup the data of the coupling
    virtual void SetupData( void ) = 0;

    //! Setup the coupling
    virtual void SetupCoupling( void ) = 0;

    //! Initialize the values of the coupling variables
    virtual void InitializeCouplingVariables( void ) = 0;

    //! Export the values of the coupling residuals
    virtual void ExportCouplingResiduals( VectorType& CouplingResiduals ) = 0;

    //! Compute Jacobian product: J(X) * X
    /*!
     * \param deltaCouplingVariables - variation of the coupling variables
     */
    virtual void ComputeJacobianProduct( const VectorType& deltaCouplingVariables ) = 0;

    //! Display some information about the coupling
    virtual void ShowMe( void );

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Get the number of models connected by the coupling
    inline UInt modelsNumber() const
    {
        return static_cast< UInt > ( M_models.size() );
    }

    //! Import the content of the Global Vector in the Local vector
    /*!
     * \param globalVector - the global vector
     * \param localVector - the local vector
     */
    void ImportCouplingVector( const VectorType& globalVector, VectorType& localVector );

    //! Export the content of the Local Vector in the Global vector
    /*!
     * \param localVector - the local vector
     * \param globalVector - the global vector
     */
    void ExportCouplingVector( const VectorType& localVector, VectorType& globalVector );

    //! Display and error message for the specific model
    /*!
     * \param model - shared_ptr to the specific model
     */
    void switchErrorMessage( const PhysicalModel_ptr& model );

    //@}

    static UInt                         M_couplingsNumber;

    UInt                                M_ID;
    couplingsTypes                      M_type;

    GetPot                              M_dataFile;
    std::vector< PhysicalModel_ptr >    M_models;
    std::vector< BCFlag >               M_flags;
    std::string                         M_couplingName;

    boost::shared_ptr< Epetra_Comm >    M_comm;
    boost::shared_ptr< Displayer >      M_displayer;

    std::pair< UInt, UInt >             M_couplingIndex;

    Vector_ptrType                      M_LocalCouplingVariables;
    Vector_ptrType                      M_LocalCouplingResiduals;
    Vector_ptrType                      M_LocalJacobianProduct;
    Vector_ptrType                      M_LocalDeltaCouplingVariables;
};

} // Namespace LifeV

#endif /* __MS_PhysicalCoupling_H */
