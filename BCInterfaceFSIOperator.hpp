/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-04-23

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
   \file BCInterfaceFSIOperator.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-04-23
 */

#ifndef __BCInterfaceFSIOperator_H
#define __BCInterfaceFSIOperator_H 1





// ===================================================
//! Include
// ===================================================
#include <life/lifecore/life.hpp>
#include <life/lifefem/bcVector.hpp>

#include <life/lifesolver/exactJacobianBase.hpp>
#include <life/lifesolver/fixedPointBase.hpp>

#include <string>





// ===================================================
//! Namespaces
// ===================================================
using namespace LifeV;
enum FSIList{
				DerFluidLoadToFluid,
				DerFluidLoadToStructure,
				DerHarmonicExtensionVelToFluid,
				DerStructureDispToSolid,
				FluidInterfaceDisp,
				FluidLoadToStructure,
				HarmonicExtensionVelToFluid,
				SolidLoadToStructure,
				StructureDispToHarmonicExtension,
				StructureDispToSolid,
				StructureToFluid
			};





/*!
 * \class BCInterfaceFSIOperator
 * \brief LifeV function parser based on boost::spirit
 *
 *  @author Cristiano Malossi
 *  @see
 */
class BCInterfaceFSIOperator
//     :
//     public LifeV::Application
{
public:

	// ===================================================
	//! Typedef
	// ===================================================



	// ===================================================
	//! Public functions
	// ===================================================

	/** @name Constructors & Destructor
     */
    //@{

    //! Constructor
	BCInterfaceFSIOperator( const std::string& functionString, const boost::shared_ptr<FSIOperator>& oper );

    //! Destructor
    ~BCInterfaceFSIOperator() {}

    //@}



    /** @name Get functions
     */
    //@{

    BCVectorInterface& getBase() { return *M_base; }

    //@}

private:

	// ===================================================
	//! Private variables
	// ===================================================

	std::string											M_baseString;
	boost::shared_ptr<FSIOperator>						M_FSIOperator;
	boost::shared_ptr<BCVectorInterface>				M_base;

	std::map<std::string, FSIList> 						M_mapFSI;



	// ===================================================
	//! Private functions
	// ===================================================

	/** @name Private functions
	*/
	//@{

	//! checkFSIList
	inline void checkFSIList( void );




    //@}
};

#endif /* __BCInterfaceFSIOperator_H */
