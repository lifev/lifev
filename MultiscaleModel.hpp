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
 *  @brief File containing the MultiScale Physical Model
 *
 *  @date 12-03-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MS_PhysicalModel_H
#define MS_PhysicalModel_H 1

#include <lifemc/lifesolver/MultiscaleDefinitions.hpp>
#include <lifemc/lifesolver/MultiscaleData.hpp>
#include <lifemc/lifesolver/MultiscaleCoupling.hpp>


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
     * @param fileName Name of data file.
     */
    virtual void setupData( const std::string& fileName );

    //! Setup the model.
    /*!
     * In particular it does the following operations:
     * <ol>
     *     <li> initialize the model object;
     *     <li> initialize all the other private objects;
     *     <li> perform preliminary operations which depend on the couplings.
     * </ol>
     */
    virtual void setupModel() = 0;

    //! Build the initial system.
    /*!
     * This method is alternative to updateSystem and should be called only once at the fist timestep.
     * This method is reserved for the construction of:
     * <ol>
     *     <li> objects that are constant (with respect to the time);
     *     <li> objects that are required for the first time step and should not be updated during subiterations.
     * </ol>
     */
    virtual void buildSystem() = 0;

    //! Update the system.
    /*!
     * This method is alternative to buildSystem and should be called from the second timestep.
     * This method is reserved for the update of:
     * <ol>
     *     <li> objects that are not constant with respect to the time but should not be updated during subiterations.
     * </ol>
     */
    virtual void updateSystem() = 0;

    //! Solve the System.
    /*!
     * This method is called once for each subiteration (in the case of implicit coupling).
     * It computes the solution at time t_n+1.
     */
    virtual void solveSystem() = 0;

    //! save the solution.
    /*!
     * This method wrote to file the solution computed during the last call of solveSystem.
     */
    virtual void saveSolution() = 0;

    //! Display some information about the model.
    virtual void showMe();

    //@}


    //! @name Methods
    //@{

    //! Clear the list of pointers to the couplings.
    /*!
     *  This method has to be called before the automatic destructor, in order
     *  to disconnect the coupling classes from the model classes.
     */
    void clearCouplingsList();

    //@}


    //! @name Set Methods
    //@{

    //! Set the global ID of the model
    /*!
     * @param id Model ID
     */
    void setID( const UInt& id );

    //! Add a pointer to one of the couplings which couple the model
    /*!
     * @param coupling shared_ptr of the coupling
     */
    void addCoupling( const MS_Coupling_PtrType& coupling );

    //! Setup the global data of the model.
    /*!
     * In particular, it can be used to replace the local values specified in
     * the model data file, with the ones in the global container.
     *
     * @param globalData Global data container.
     */
    void setGlobalData( const MS_GlobalDataContainer_PtrType& globalData );

    //! Scale, rotate and translate the Model in the 3D space
    /*!
     * This method apply directly to the mesh before its partitioning.
     * @param scale Vector (Sx,Sy,Sz) of scale factors
     * @param rotate Vector (Rx,Ry,Rz) of angles for rotation (degree units)
     * @param translate Vector (Tx,Ty,Tz) of offset for position
     */
    void setGeometry( const boost::array< Real, NDIM >& scale,
                      const boost::array< Real, NDIM >& rotate,
                      const boost::array< Real, NDIM >& translate );

    //! Set the epetra communicator for the model
    /*!
     * @param comm Epetra communicator
     */
    void setCommunicator( const MS_Comm_PtrType& comm );

    //@}


    //! @name Get Methods
    //@{

    //! Get the global ID of the model
    /*!
     * @return global ID of the model
     */
    const UInt& ID() const;

    //! Get the type of the model
    /*!
     * @return type of the model
     */
    const models_Type& type() const;

    //! Get one available flag by id
    /*!
     * @param id id of the boundary flag
     * @return flag
     */
    const BCFlag& flag( const UInt& id ) const;

    //! Get the flag vector of the model
    /*!
     * @return vector of flags available for coupling
     */
    const std::vector< BCFlag >& flags() const;

    //! Get the name of the model
    /*!
     * @return name of the model
     */
    const std::string& modelName() const;

    //! Get the number of couplings connecting the model
    /*!
     * @return number of couplings connecting the model
     */
    UInt couplingsNumber() const;

    //! Get the coupling local ID through global ID
    /*!
     * @param ID global ID of the coupling
     * @return local ID of the coupling
     */
    UInt couplingLocalID( const UInt& ID ) const;

    //! Get the coupling through local ID
    /*!
     * @param ID local ID of the coupling
     * @return Pointer to the coupling
     */
    MS_Coupling_PtrType coupling( const UInt& localID ) const;

    //! Get the global data of the model.
    /*!
     * @return Global data container.
     */
    const MS_GlobalDataContainer_PtrType& globalData() const;

    //@}

protected:

    static UInt                          M_modelsNumber;       // Total number of models

    UInt                                 M_ID;                 // Global ID of the model
    models_Type                          M_type;               // Type of the model (depends on the derived class)

    MS_CouplingsVector_Type              M_couplings;          // Container for the couplings
    std::string                          M_modelName;          // Name of the model
    std::vector< BCFlag >                M_flags;              // Free flags, available for the couplings

    MS_GlobalDataContainer_PtrType       M_globalData;         // GlobalDataContainer

    boost::array< Real, NDIM >           M_geometryScale;      // Global geometrical scale
    boost::array< Real, NDIM >           M_geometryRotate;     // Global geometrical rotation
    boost::array< Real, NDIM >           M_geometryTranslate;  // Global geometrical translation

    MS_Comm_PtrType                      M_comm;               // Communicator
    boost::shared_ptr< Displayer >       M_displayer;          // Displayer

private:

    //! @name Unimplemented Methods
    //@{

    MS_PhysicalModel( const MS_PhysicalModel& model );

    MS_PhysicalModel& operator=( const MS_PhysicalModel& model );

    //@}
};

} // Namespace LifeV

#endif /* MS_PhysicalModel_H */
