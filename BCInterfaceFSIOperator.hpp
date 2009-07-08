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
//#include <life/lifesolver/steklovPoincareBase.hpp>

#include <string>





// ===================================================
//! Namespaces & Enums
// ===================================================
namespace LifeV {

enum FSIMethod
{
	EXACTJACOBIAN,
	FIXEDPOINT,
	MONOLITHIC,
	STEKLOVPOINCARE
};

enum FSIFunction
{
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
 * \brief LifeV bcVector wrapper for BCInterface (FSI problems).
 *
 *  @author Cristiano Malossi
 *  @see
 *
 *  This class allows to use impose interface conditions for FSI problems.
 *
 *
 *
 *  <b>DETAILS:</b>
 *
 *  The constructor of the class takes a string contains the ID of the interface condition to impose,
 *  and the FSIOperator. The list of available conditions is the FSIFunction variable. These are:
 *
 *	- DerFluidLoadToFluid,					(not implemented)
 *	- DerFluidLoadToStructure,
 *	- DerHarmonicExtensionVelToFluid,
 *	- DerStructureDispToSolid,				(not implemented)
 *	- FluidInterfaceDisp,					(not working)
 *	- FluidLoadToStructure,
 *	- HarmonicExtensionVelToFluid,
 *	- SolidLoadToStructure,
 *	- StructureDispToHarmonicExtension,
 *	- StructureDispToSolid, 				(not implemented)
 *	- StructureToFluid
 *
 *	The class automatically recognize which method is used among:
 *
 *	- EXACTJACOBIAN
 *	- FIXEDPOINT
 *	- MONOLITHIC		(not working)
 *	- STEKLOVPOINCARE 	(not working)
 *
 *	To get the base for the boundary condition call the getBase function.
 *
 */
class BCInterfaceFSIOperator
//     :
//     public LifeV::Application
{
public:

	// ===================================================
	//! Public functions
	// ===================================================

	/** @name Constructors & Destructor
     */
    //@{

    //! Constructor
	/*!
	 * functionString - interface condition ID
	 * oper           - fsiOperator
	 */
	BCInterfaceFSIOperator( const std::string& functionString, const boost::shared_ptr<FSIOperator>& oper );

	//! Copy constructor
	/*!
	 * \param fsiOperator - BCInterfaceFSIOperator
	 */
	BCInterfaceFSIOperator( const BCInterfaceFSIOperator& fsiOperator );

	//! Operator =
	/*!
	 * \param fsiOperator - BCInterfaceFSIOperator
	 */
	BCInterfaceFSIOperator& operator=( const BCInterfaceFSIOperator& fsiOperator );

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

	std::map<std::string, FSIMethod> 					M_mapMethod;
	std::map<std::string, FSIFunction> 					M_mapFunction;



	// ===================================================
	//! Private functions
	// ===================================================

	/** @name Private functions
	*/
	//@{

	//! checkMethod
	inline void checkMethod( void );

	template <class method>
	inline void checkFunction( void );

    //@}
};





// ===================================================
//! Template function
// ===================================================
template <class method>
inline void
BCInterfaceFSIOperator::checkFunction( void )
{
	method *operMethod = dynamic_cast<method *>(&*M_FSIOperator);

	switch ( M_mapFunction[M_baseString] )
	{
		case DerFluidLoadToFluid :

#ifdef DEBUG
			Debug( 5022 ) << "BCInterfaceFSIOperator::checkFunction -> DerFluidLoadToFluid" << "\n";
#endif

			break;

		case DerFluidLoadToStructure :

#ifdef DEBUG
			Debug( 5022 ) << "BCInterfaceFSIOperator::checkFunction -> DerFluidLoadToStructure" << "\n";
#endif

			operMethod->setDerFluidLoadToStructure( M_FSIOperator->sigmaSolidRepeated() );

	        M_base = operMethod->bcvDerFluidLoadToStructure();

			break;

		case DerHarmonicExtensionVelToFluid :

#ifdef DEBUG
			Debug( 5022 ) << "BCInterfaceFSIOperator::checkFunction -> DerHarmonicExtensionVelToFluid" << "\n";
#endif

			operMethod->setDerHarmonicExtensionVelToFluid( M_FSIOperator->derVeloFluidMesh() );

	        M_base = operMethod->bcvDerHarmonicExtensionVelToFluid();

			break;

		case DerStructureDispToSolid :

#ifdef DEBUG
			Debug( 5022 ) << "BCInterfaceFSIOperator::checkFunction -> DerStructureDispToSolid" << "\n";
#endif

			break;

		case FluidInterfaceDisp :

#ifdef DEBUG
			Debug( 5022 ) << "BCInterfaceFSIOperator::checkFunction -> FluidInterfaceDisp" << "\n";
#endif

			//operMethod->FluidInterfaceDisp( (LifeV::Vector&) M_FSIOperator->lambdaFluidRepeated() );

	        //M_base = operMethod->bcvFluidInterfaceDisp();

			break;

		case FluidLoadToStructure :

#ifdef DEBUG
			Debug( 5022 ) << "BCInterfaceFSIOperator::checkFunction -> FluidLoadToStructure" << "\n";
#endif

			operMethod->setFluidLoadToStructure( M_FSIOperator->sigmaSolidRepeated() );

	        M_base = operMethod->bcvFluidLoadToStructure();

			break;

		case HarmonicExtensionVelToFluid :

#ifdef DEBUG
			Debug( 5022 ) << "BCInterfaceFSIOperator::checkFunction -> HarmonicExtensionVelToFluid" << "\n";
#endif

			M_FSIOperator->setHarmonicExtensionVelToFluid( M_FSIOperator->veloFluidMesh() );

			M_base = M_FSIOperator->bcvHarmonicExtensionVelToFluid();

			break;

		case SolidLoadToStructure :

#ifdef DEBUG
			Debug( 5022 ) << "BCInterfaceFSIOperator::checkFunction -> SolidLoadToStructure" << "\n";
#endif

			M_FSIOperator->setSolidLoadToStructure( M_FSIOperator->minusSigmaFluidRepeated() );

			M_base = M_FSIOperator->bcvSolidLoadToStructure();

			break;

		case StructureDispToHarmonicExtension :

#ifdef DEBUG
			Debug( 5022 ) << "BCInterfaceFSIOperator::checkFunction -> StructureDispToHarmonicExtension" << "\n";
#endif

			operMethod->setStructureDispToHarmonicExtension( M_FSIOperator->lambdaFluidRepeated() );

	        M_base = operMethod->bcvStructureDispToHarmonicExtension();

			break;

		case StructureDispToSolid :

#ifdef DEBUG
			Debug( 5022 ) << "BCInterfaceFSIOperator::checkFunction -> StructureDispToSolid" << "\n";
#endif

			break;

		case StructureToFluid :

#ifdef DEBUG
			Debug( 5022 ) << "BCInterfaceFSIOperator::checkFunction -> StructureToFluid" << "\n";
#endif

			M_FSIOperator->setStructureToFluid( M_FSIOperator->veloFluidMesh() );
			M_FSIOperator->setStructureToFluidParametres();

			M_base = M_FSIOperator->bcvStructureToFluid();

			break;
	}
}

} // Namespace LifeV

#endif /* __BCInterfaceFSIOperator_H */
