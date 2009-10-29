/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-03-12

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
 \file MS_PhysicalModel.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-03-12
 */

#ifndef __MS_PhysicalModel_H
#define __MS_PhysicalModel_H 1

#include "Epetra_ConfigDefs.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <life/lifecore/life.hpp>
#include <life/lifecore/util_string.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/displayer.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifemesh/markers.hpp>

#include <life/lifealg/EpetraMap.hpp>
#include <life/lifearray/EpetraVector.hpp>

#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>

#include <lifemc/lifesolver/MS_PhysicalData.hpp>

namespace LifeV {

//! @name MS_PhysicalModel global objects
//@{

enum modelsTypes
{
    MultiScale,
    Fluid3D
};

extern std::map< std::string, modelsTypes > modelsMap;

//@}

//! MS_PhysicalModel - The MultiScale Physical Model
/*!
 *  The MS_PhysicalModel class provides a general interface between the
 *  MS_Algorithm and all the other models.
 *
 *  @author Cristiano Malossi
 */
class MS_PhysicalModel
{
public:

    typedef EpetraVector                       VectorType;
    typedef boost::shared_ptr<VectorType>      Vector_ptrType;

    typedef EntityFlag                         BCFlag;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_PhysicalModel();

    //! Copy constructor
    /*!
     * \param model - MS_PhysicalModel
     */
    MS_PhysicalModel( const MS_PhysicalModel& model );

    //! Destructor
    virtual ~MS_PhysicalModel() {}

    //@}


    //! @name Methods
    //@{

    //! Operator=
    /*!
     * \param model - MS_PhysicalModel
     */
    MS_PhysicalModel& operator=( const MS_PhysicalModel& model );

    //@}


    //! @name Set Methods
    //@{

    //! Set the global ID of the model
    /*!
     * \param id - Model ID
     */
    void SetID( const UInt& id )
    {
        M_ID = id;
    }

    //! Set the data file to load information of the model
    /*!
     * \param dataFile - Name and path of data file
     */
    void SetDataFile( const std::string& dataFile );

    //! Rescale, rotate and translate the Model in the 3D space
    /*!
     * \param scale - Vector (Sx,Sy,Sz) of scale factors
     * \param rotate - Vector (Rx,Ry,Rz) of angles for rotation (degree units)
     * \param translate - Vector (Tx,Ty,Tz) of offset for position
     */
    void SetGeometry( const boost::array< Real, nDimensions >& scale,
                      const boost::array< Real, nDimensions >& rotate,
                      const boost::array< Real, nDimensions >& translate );

    //! Set global data for physical quantities and time
    /*!
     * \param dataPhysics - Data container for physical quantities
     * \param dataTime - Data container for time parameters
     */
    void SetData( const boost::shared_ptr< MS_PhysicalData >& dataPhysics,
                  const boost::shared_ptr< DataTime >& dataTime );

    //! Set the epetra communicator for the model
    /*!
     * \param comm - Epetra communicator
     */
    void SetCommunicator( const boost::shared_ptr< Epetra_Comm >& comm );

    //@}


    //! @name Get Methods
    //@{

    //! Get the global ID of the model
    const UInt& GetID() const
    {
        return M_ID;
    }

    //! Get the type of the model
    const modelsTypes& GetType() const
    {
        return M_type;
    }

    //! Get the flag i-flag of the model
    /*!
     * \param id - id of the boundary flag
     */
    const BCFlag& GetFlag( const UInt& id ) const
    {
        return M_flags[id];
    }

    //! Get the flag vector of the model
    const std::vector< BCFlag >& GetFlags() const
    {
        return M_flags;
    }

    //! Get the name of the model
    const std::string& GetModelName() const
    {
        return M_modelName;
    }

    //@}


    //! @name MultiScale PhysicalModel Virtual Functions
    //@{

    //! Setup the data of the model
    virtual void SetupData( void ) = 0;

    //! Setup the model
    virtual void SetupModel( void ) = 0;

    //! Build the system matrix and vectors
    virtual void BuildSystem( void ) = 0;

    //! Update the system matrix and vectors
    virtual void UpdateSystem( void ) = 0;

    //! Solve the problem
    virtual void SolveSystem( void ) = 0;

    //! Save the solution
    virtual void SaveSolution( void ) = 0;

    //! Display some information about the model
    virtual void ShowMe( void );

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Get the number of flags of the model
    UInt flagsNumber() const
    {
        return static_cast< UInt > ( M_flags.size() );
    }

    //@}

    static UInt                          M_modelsNumber;

    UInt                                 M_ID;
    modelsTypes                          M_type;

    GetPot                               M_dataFile;
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

#endif /* __MS_PhysicalModel_H */
