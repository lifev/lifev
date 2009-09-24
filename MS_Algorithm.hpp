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
   \file MS_Algorithm.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-03-12
 */

#ifndef __MS_Algorithm_H
#define __MS_Algorithm_H 1

#include "Epetra_ConfigDefs.h"
#ifdef EPETRA_MPI
    #include "Epetra_MpiComm.h"
#else
    #include "Epetra_SerialComm.h"
#endif

#include <life/lifecore/life.hpp>
#include <life/lifearray/tab.hpp>

#include <lifemc/lifesolver/MS_PhysicalData.hpp>
#include <lifemc/lifesolver/MS_PhysicalModel.hpp>
#include <lifemc/lifesolver/MS_Model_Fluid3D.hpp>

#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>
#include <lifemc/lifesolver/MS_Coupling_BoundaryCondition.hpp>
#include <lifemc/lifesolver/MS_Coupling_FluxPressure.hpp>

#include <boost/array.hpp>
#include <boost/algorithm/string.hpp>

namespace LifeV {

//! MS_Algorithm - The MultiScale problem generator
/*!
 *  The MS_Algorithm class provides a series of functions to coupling different
 *  models and run a general multiscale problem.
 *
 *  @author Cristiano Malossi
 */
class MS_Algorithm //: public MS_PhysicalModel //(Future idea!)
{
public:

	typedef singleton< factory< MS_PhysicalModel,		modelsTypes > >		FactoryModels;
	typedef singleton< factory< MS_PhysicalCoupling,	couplingsTypes > >	FactoryCouplings;

	typedef boost::shared_ptr< MS_PhysicalModel >							PhysicalModel_ptr;
	typedef boost::shared_ptr< MS_PhysicalCoupling >						PhysicalCoupling_ptr;

	//! @name Constructors, Destructor
    //@{

    //! Constructor
	MS_Algorithm();

	//! Copy constructor
	/*!
	 * \param algorithm - MS_Algorithm
	 */
	MS_Algorithm( const MS_Algorithm& algorithm );

	//! Destructor
    ~MS_Algorithm() {}

    //@}



    //! @name Methods
    //@{

	//! Operator=
	/*!
	 * \param algorithm - MS_Algorithm
	 */
    MS_Algorithm& operator=( const MS_Algorithm& algorithm );

	//! Set the data file to load information of the model
	/*!
	 * \param dataFile - Name and path of data file
	 */
	void SetDataFile( const std::string& dataFile );

	//! Set the epetra communicator for the model
	/*!
	 * \param comm - Epetra communicator
	 */
	void SetCommunicator( const boost::shared_ptr<Epetra_Comm>& comm );

    //! Load data and create the models and the couplings
	void SetupData( void );

    //! Setup the models and the couplings
	void SetupProblem( void );

    //! Run the time-loop to solve the multiscale problem
	void SolveProblem( void );

	//! Display some information about the multiscale problem (call after SetupProblem)
	void ShowMe( void );

    //@}

private:

	//! @name Private Methods
    //@{

	inline void loadModels		( const GetPot& dataFile );
	inline void loadCouplings	( const GetPot& dataFile );
	inline void loadGeometry	( const GetPot& dataFile );

	template <typename number>
	inline std::vector<number> string2numVect( const std::string& string );

    //@}

	// Communicator
	boost::shared_ptr<Epetra_Comm>						M_comm;

	// Displayer tool for MPI processes
	boost::shared_ptr<Displayer>						M_displayer;

    // Name of the GetPot data file
	GetPot		 										M_dataFile;

    // DataTime container
	boost::shared_ptr<DataTime>							M_dataTime;

    // PhysicalData container
	boost::shared_ptr<MS_PhysicalData>					M_dataPhysics;

	// Models & Connections
	std::map<UInt, PhysicalModel_ptr>					M_models;
	std::map<UInt, PhysicalCoupling_ptr>				M_couplings;

	// Models & Connections size
	UInt												M_modelsNumber;
	UInt												M_couplingsNumber;

	// Chrono performances
	Chrono												M_chrono;
};

} // Namespace LifeV

#endif /* __MS_Algorithm_H */
