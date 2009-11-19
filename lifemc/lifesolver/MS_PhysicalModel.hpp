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
 *  @brief MultiScale Physical Model
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 12-03-2009
 */

#ifndef MS_PhysicalModel_H
#define MS_PhysicalModel_H 1

#include <lifemc/lifesolver/MS_Definitions.hpp>
#include <lifemc/lifesolver/MS_PhysicalData.hpp>
#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>


namespace LifeV {

//! MS_PhysicalModel - The MultiScale Physical Model
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_PhysicalModel class provides a general interface between the
 *  MS_Algorithm and all the other models.
 */
class MS_PhysicalModel
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_PhysicalModel();

    //! Copy constructor
    /*!
     * @param model MS_PhysicalModel
     */
    MS_PhysicalModel( const MS_PhysicalModel& model );

    //! Destructor
    virtual ~MS_PhysicalModel() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param model MS_PhysicalModel
     * @return reference to a copy of the class
     */
    MS_PhysicalModel& operator=( const MS_PhysicalModel& model );

    //@}


    //! @name MultiScale PhysicalModel Virtual Methods
    //@{

    //! Setup the data of the model
    virtual void SetupData() = 0;

    //! Setup the model
    virtual void SetupModel() = 0;

    //! Build the system matrix and vectors
    virtual void BuildSystem() = 0;

    //! Update the system matrix and vectors
    virtual void UpdateSystem() = 0;

    //! Solve the problem
    virtual void SolveSystem() = 0;

    //! Save the solution
    virtual void SaveSolution() = 0;

    //! Display some information about the model
    virtual void ShowMe();

    //@}


    //! @name Set Methods
    //@{

    //! Set the global ID of the model
    /*!
     * @param id Model ID
     */
    void SetID( const UInt& id );

    //! Set the data file to load information of the model
    /*!
     * @param dataFile Name and path of data file
     */
    void SetDataFile( const std::string& dataFile );

    //! Add a pointer to one of the couplings which couple the model
    /*!
     * @param coupling shared_ptr of the coupling
     */
    void AddCoupling( const Coupling_ptrType& coupling );

    //! Rescale, rotate and translate the Model in the 3D space
    /*!
     * @param scale Vector (Sx,Sy,Sz) of scale factors
     * @param rotate Vector (Rx,Ry,Rz) of angles for rotation (degree units)
     * @param translate Vector (Tx,Ty,Tz) of offset for position
     */
    void SetGeometry( const boost::array< Real, nDimensions >& scale,
                      const boost::array< Real, nDimensions >& rotate,
                      const boost::array< Real, nDimensions >& translate );

    //! Set global data for physical quantities and time
    /*!
     * @param dataPhysics Data container for physical quantities
     * @param dataTime Data container for time parameters
     */
    void SetData( const boost::shared_ptr< MS_PhysicalData >& dataPhysics,
                  const boost::shared_ptr< DataTime >& dataTime );

    //! Set the epetra communicator for the model
    /*!
     * @param comm Epetra communicator
     */
    void SetCommunicator( const boost::shared_ptr< Epetra_Comm >& comm );

    //@}


    //! @name Get Methods
    //@{

    //! Get the global ID of the model
    /*!
     * @return global ID of the model
     */
    const UInt& GetID() const;

    //! Get the type of the model
    /*!
     * @return type of the model
     */
    const modelsTypes& GetType() const;

    //! Get one available flag by id
    /*!
     * @param id id of the boundary flag
     * @return flag
     */
    const BCFlag& GetFlag( const UInt& id ) const;

    //! Get the flag vector of the model
    /*!
     * @return vector of flags available for coupling
     */
    const std::vector< BCFlag >& GetFlags() const;

    //! Get the name of the model
    /*!
     * @return name of the model
     */
    const std::string& GetModelName() const;

    //! Get the number of couplings connecting the model
    /*!
     * @return number of couplings connecting the model
     */
    UInt GetCouplingsNumber() const;

    //! Get the coupling local ID through global ID
    /*!
     * @param ID global ID of the coupling
     * @return local ID of the coupling
     */
    UInt GetCouplingLocalID( const UInt& ID ) const;

    //! Get the coupling through local ID
    /*!
     * @param ID local ID of the coupling
     * @return Pointer to the coupling
     */
    Coupling_ptrType GetCoupling( const UInt& LocalID ) const;

    //@}

protected:

    //! @name Protected Methods
    //@{

    //@}

    static UInt                          M_modelsNumber;

    UInt                                 M_ID;
    modelsTypes                          M_type;

    GetPot                               M_dataFile;
    CouplingsVector_Type                 M_couplings;
    std::string                          M_modelName;
    std::vector< BCFlag >                M_flags;

    boost::array< Real, nDimensions >    M_geometryScale;
    boost::array< Real, nDimensions >    M_geometryRotate;
    boost::array< Real, nDimensions >    M_geometryTranslate;

    boost::shared_ptr< MS_PhysicalData > M_dataPhysics;
    boost::shared_ptr< DataTime >        M_dataTime;

    boost::shared_ptr< Epetra_Comm >     M_comm;
    boost::shared_ptr< Displayer >       M_displayer;
};

} // Namespace LifeV

#endif /* MS_PhysicalModel_H */
