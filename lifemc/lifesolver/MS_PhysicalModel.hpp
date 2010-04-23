//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

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
 *  MS_Algorithm and all the other models. Moreover it provides internal methods
 *  for accessing general information and model related couplings.
 */
class MS_PhysicalModel
{
public:

    //! @name Constructors & Destructor
    //@{

    //! The main constructor.
    /*!
     * All the derived classes has to
     * be constructed using an empty constructor which calls (as first operation)
     * this empty constructor.
     */
    MS_PhysicalModel();

    //! Copy constructor
    /*!
     * Should not be used in classical situations.
     * @param model MS_PhysicalModel
     */
    MS_PhysicalModel( const MS_PhysicalModel& model );

    //! Destructor
    virtual ~MS_PhysicalModel() {}

    //@}


    //! @name Operators
    //@{

    //! Copy operator.
    /*!
     * Should not be used in classical situations.
     * @param model MS_PhysicalModel
     * @return reference to a copy of the class
     */
    MS_PhysicalModel& operator=( const MS_PhysicalModel& model );

    //@}


    //! @name MultiScale PhysicalModel Virtual Methods
    //@{

    //! Setup the data of the model.
    /*!
     * In particular it does the following operations:
     * <ol>
     *     <li> read data from files;
     *     <li> set global parameter for the MS simulation (viscosity, time, ...);
     *     <li> perform preliminary operations which don't depend on the couplings.
     * </ol>
     *
     * @param FileName Name of data file.
     */
    virtual void SetupData( const std::string& FileName );

    //! Setup the model.
    /*!
     * In particular it does the following operations:
     * <ol>
     *     <li> initialize the model object;
     *     <li> initialize all the other private objects;
     *     <li> perform preliminary operations which depend on the couplings.
     * </ol>
     */
    virtual void SetupModel() = 0;

    //! Build the initial system (matrix and vectors).
    /*!
     * Matrix and vectors built here are for the first time step. If they are time independent could be reused later.
     * More generally, this method initialize all the objects which are not defined at the first time step.
     *
     * <b>This method should be called only once at the first time step!</b>
     */
    virtual void BuildSystem() = 0;

    //! Update the system for (matrix and vectors)
    /*!
     * Only time dependent objects are updated here.
     */
    virtual void UpdateSystem() = 0;

    //! Solve the problem.
    /*!
     * This method solve the problem.
     */
    virtual void SolveSystem() = 0;

    //! Save the solution.
    /*!
     * This method wrote to file the solution computed during the last call of SolveSystem.
     */
    virtual void SaveSolution() = 0;

    //! Display some information about the model.
    virtual void ShowMe();

    //@}


    //! @name Methods
    //@{

    //! Clear the list of pointers to the couplings.
    /*!
     *  This method has to be called before the automatic destructor, in order
     *  to disconnect the coupling classes from the model classes.
     */
    void ClearCouplingsList();

    //@}


    //! @name Set Methods
    //@{

    //! Set the global ID of the model
    /*!
     * @param id Model ID
     */
    void SetID( const UInt& id );

    //! Add a pointer to one of the couplings which couple the model
    /*!
     * @param coupling shared_ptr of the coupling
     */
    void AddCoupling( const Coupling_ptrType& coupling );

    //! Scale, rotate and translate the Model in the 3D space
    /*!
     * This method apply directly to the mesh before its partitioning.
     * @param scale Vector (Sx,Sy,Sz) of scale factors
     * @param rotate Vector (Rx,Ry,Rz) of angles for rotation (degree units)
     * @param translate Vector (Tx,Ty,Tz) of offset for position
     */
    void SetGeometry( const boost::array< Real, NDIM >& scale,
                      const boost::array< Real, NDIM >& rotate,
                      const boost::array< Real, NDIM >& translate );

    //! Set global data for physical quantities and time
    /*!
     * This method set all the data
     * @param dataPhysics Data container for physical quantities
     */
    void SetGlobalData( const boost::shared_ptr< MS_PhysicalData >& dataPhysics );

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

    static UInt                          M_modelsNumber;       // Total number of models

    UInt                                 M_ID;                 // Global ID of the model
    modelsTypes                          M_type;               // Type of the model (depends on the derived class)

    CouplingsVector_Type                 M_couplings;          // Container for the couplings
    std::string                          M_modelName;          // Name of the model
    std::vector< BCFlag >                M_flags;              // Free flags, available for the couplings

    boost::array< Real, NDIM >           M_geometryScale;      // Global geometrical scale
    boost::array< Real, NDIM >           M_geometryRotate;     // Global geometrical rotation
    boost::array< Real, NDIM >           M_geometryTranslate;  // Global geometrical translation

    boost::shared_ptr< MS_PhysicalData > M_dataPhysics;        // Data container for global physical quantities

    boost::shared_ptr< Epetra_Comm >     M_comm;               // Communicator
    boost::shared_ptr< Displayer >       M_displayer;          // Displayer
};

} // Namespace LifeV

#endif /* MS_PhysicalModel_H */
