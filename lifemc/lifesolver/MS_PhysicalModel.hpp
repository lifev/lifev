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


namespace LifeV
{

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

    //! Destructor
    virtual ~MS_PhysicalModel() {}

    //@}


    //! @name MultiScale PhysicalModel Virtual Methods
    //@{

    //! Setup the data of the model.
    /*!
     * This method is called only once at the beginning of the simulation.
     * It is designed for the following operations:
     * <ol>
     *     <li> read data from files;
     *     <li> set global parameter for the MS simulation (viscosity, time, ...);
     *     <li> perform preliminary operations which do not depend on the couplings.
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

    //! Build the initial system.
    /*!
     * This method is alternative to UpdateSystem and should be called only once at the fist timestep.
     * This method is reserved for the construction of:
     * <ol>
     *     <li> objects that are constant (with respect to the time);
     *     <li> objects that are required for the first time step and should not be updated during subiterations.
     * </ol>
     */
    virtual void BuildSystem() = 0;

    //! Update the system.
    /*!
     * This method is alternative to BuildSystem and should be called from the second timestep.
     * This method is reserved for the update of:
     * <ol>
     *     <li> objects that are not constant with respect to the time but should not be updated during subiterations.
     * </ol>
     */
    virtual void UpdateSystem() = 0;

    //! Solve the System.
    /*!
     * This method is called once for each subiteration (in the case of implicit coupling).
     * It computes the solution at time t_n+1.
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
    void AddCoupling( const MS_Coupling_PtrType& coupling );

    //! Setup the global data of the model.
    /*!
     * In particular, it can be used to replace the local values specified in
     * the model data file, with the ones in the global container.
     *
     * @param globalData Global data container.
     */
    void SetGlobalData( const MS_GlobalDataContainer_PtrType& globalData );

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

    //! Set the epetra communicator for the model
    /*!
     * @param comm Epetra communicator
     */
    void SetCommunicator( const MS_Comm_PtrType& comm );

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
    MS_Coupling_PtrType GetCoupling( const UInt& LocalID ) const;

    //! Get the global data of the model.
    /*!
     * @return Global data container.
     */
    const MS_GlobalDataContainer_PtrType& GetGlobalData() const;

    //@}

protected:

    static UInt                          M_modelsNumber;       // Total number of models

    UInt                                 M_ID;                 // Global ID of the model
    modelsTypes                          M_type;               // Type of the model (depends on the derived class)

    MS_CouplingsVector_Type              M_couplings;          // Container for the couplings
    std::string                          M_modelName;          // Name of the model
    std::vector< BCFlag >                M_flags;              // Free flags, available for the couplings

    MS_GlobalDataContainer_PtrType       M_globalData;         // GlobalDataContainer

    boost::array< Real, NDIM >           M_geometryScale;      // Global geometrical scale
    boost::array< Real, NDIM >           M_geometryRotate;     // Global geometrical rotation
    boost::array< Real, NDIM >           M_geometryTranslate;  // Global geometrical translation

    MS_Comm_PtrType                      M_comm;               // Communicator
    boost::shared_ptr< Displayer >       M_displayer;          // Displayer
};

} // Namespace LifeV

#endif /* MS_PhysicalModel_H */
