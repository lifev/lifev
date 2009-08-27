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
   \file BCInterfaceFSI.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-04-23
 */

#ifndef __BCInterfaceFSI_H
#define __BCInterfaceFSI_H 1





// ===================================================
//! Include
// ===================================================
#include <life/lifecore/life.hpp>
#include <life/lifefem/bcVector.hpp>

#include <life/lifesolver/exactJacobianBase.hpp>
#include <life/lifesolver/fixedPointBase.hpp>
//#include <life/lifesolver/steklovPoincareBase.hpp>

#include <string>

#include <lifemc/lifefem/BCInterfaceData.hpp>





// ===================================================
//! Namespaces & Enums
// ===================================================
namespace LifeV {




/*!
 * \class BCInterfaceFSI - Fake class for non-FSI problems.
 * \brief Fake class for non-FSI problems
 *
 *  @author Cristiano Malossi
 *  @see
 *
 */
template <class Operator>
class BCInterfaceFSI
//     :
//     public LifeV::Application
{
public:

	/** @name Constructors & Destructor
     */
    //@{

	BCInterfaceFSI( )											{}
	BCInterfaceFSI( const BCInterfaceData<Operator>& /*data*/ ) {}
	BCInterfaceFSI( const BCInterfaceFSI& /*fsi*/ ) 			{}
    ~BCInterfaceFSI() 											{}

    //@}




    /** @name Methods
     */
    //@{

	BCInterfaceFSI& operator=( const BCInterfaceFSI& /*fsi*/ ) 	{}
	void setData( const BCInterfaceData<Operator>& /*data*/ )	{}
	bool compare( const BCInterfaceData<Operator>& /*data*/ )	{ return true; }

    //@}



    /** @name Get functions
     */
    //@{

    BCVectorInterface& getBase()								{ return *M_base; }

    //@}

private:

	boost::shared_ptr<BCVectorInterface>						M_base;
};





/*!
 * \class BCInterfaceFSI - Specialized template implementation for FSI problems.
 * \brief LifeV bcVector wrapper for BCInterface (FSI problems).
 *
 *  @author Cristiano Malossi
 *  @see
 *
 *  This class allows to use impose interface conditions for FSI problems.
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
template <>
class BCInterfaceFSI<FSIOperator>
//     :
//     public LifeV::Application
{
public:

	/** @name Constructors & Destructor
     */
    //@{

    //! Constructor
	BCInterfaceFSI( );

    //! Constructor
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	BCInterfaceFSI( const BCInterfaceData<FSIOperator>& data );

	//! Copy constructor
	/*!
	 * \param fsiOperator		- BCInterfaceFSI
	 */
	BCInterfaceFSI( const BCInterfaceFSI& fsi );

    //! Destructor
    ~BCInterfaceFSI() {}

    //@}




    /** @name Methods
     */
    //@{

	//! Operator =
	/*!
	 * \param fsiOperator		- BCInterfaceFSI
	 */
	BCInterfaceFSI& operator=( const BCInterfaceFSI& fsi );

	//! Set data
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	void setData( const BCInterfaceData<FSIOperator>& data );

	//! Compare function
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	bool compare( const BCInterfaceData<FSIOperator>& data );

    //@}



    /** @name Get functions
     */
    //@{

    BCVectorInterface& getBase() { return *M_base; }

    //@}

private:

	/** @name Private functions
	*/
	//@{

	inline void checkMethod( void );

	template <class method>
	inline void checkFunction( void );

    //@}

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

	boost::shared_ptr<FSIOperator>						M_operator;
	std::string											M_baseString;
	boost::shared_ptr<BCVectorInterface>				M_base;

	std::map<std::string, FSIMethod> 					M_mapMethod;
	std::map<std::string, FSIFunction> 					M_mapFunction;
};





// ===================================================
//! Private functions
// ===================================================
template <class method>
inline void
BCInterfaceFSI<FSIOperator>::checkFunction( void )
{
	method *operMethod = dynamic_cast<method *>(&*M_operator);

	//Set mapFunction
	M_mapFunction["DerFluidLoadToFluid"]				= DerFluidLoadToFluid;
	M_mapFunction["DerFluidLoadToStructure"]			= DerFluidLoadToStructure;
	M_mapFunction["DerHarmonicExtensionVelToFluid"]		= DerHarmonicExtensionVelToFluid;
	M_mapFunction["DerStructureDispToSolid"]			= DerStructureDispToSolid;
	M_mapFunction["FluidInterfaceDisp"]					= FluidInterfaceDisp;
	M_mapFunction["FluidLoadToStructure"]				= FluidLoadToStructure;
	M_mapFunction["HarmonicExtensionVelToFluid"] 		= HarmonicExtensionVelToFluid;
	M_mapFunction["SolidLoadToStructure"]				= SolidLoadToStructure;
	M_mapFunction["StructureDispToHarmonicExtension"]	= StructureDispToHarmonicExtension;
	M_mapFunction["StructureDispToSolid"]				= StructureDispToSolid;
	M_mapFunction["StructureToFluid"]					= StructureToFluid;

	switch ( M_mapFunction[M_baseString] )
	{
		case DerFluidLoadToFluid :

#ifdef DEBUG
			Debug( 5029 ) << "BCInterfaceFSI::checkFunction                          DerFluidLoadToFluid" << "\n";
#endif

			break;

		case DerFluidLoadToStructure :

#ifdef DEBUG
			Debug( 5029 ) << "BCInterfaceFSI::checkFunction                          DerFluidLoadToStructure" << "\n";
#endif
		    if ( !M_operator->isSolid() )
		    	return;

			operMethod->setDerFluidLoadToStructure( M_operator->sigmaSolidRepeated() );

	        M_base = operMethod->bcvDerFluidLoadToStructure();

			break;

		case DerHarmonicExtensionVelToFluid :

#ifdef DEBUG
			Debug( 5029 ) << "BCInterfaceFSI::checkFunction                          DerHarmonicExtensionVelToFluid" << "\n";
#endif

		    if ( !M_operator->isFluid() )
		    	return;

			operMethod->setDerHarmonicExtensionVelToFluid( M_operator->derVeloFluidMesh() );

	        M_base = operMethod->bcvDerHarmonicExtensionVelToFluid();

			break;

		case DerStructureDispToSolid :

#ifdef DEBUG
			Debug( 5029 ) << "BCInterfaceFSI::checkFunction                          DerStructureDispToSolid" << "\n";
#endif

			break;

		case FluidInterfaceDisp :

#ifdef DEBUG
			Debug( 5029 ) << "BCInterfaceFSI::checkFunction                          FluidInterfaceDisp" << "\n";
#endif

			//operMethod->FluidInterfaceDisp( (LifeV::Vector&) M_operator->lambdaFluidRepeated() );

	        //M_base = operMethod->bcvFluidInterfaceDisp();

			break;

		case FluidLoadToStructure :

#ifdef DEBUG
			Debug( 5029 ) << "BCInterfaceFSI::checkFunction                          FluidLoadToStructure" << "\n";
#endif

		    if ( !M_operator->isSolid() )
		    	return;

			operMethod->setFluidLoadToStructure( M_operator->sigmaSolidRepeated() );

	        M_base = operMethod->bcvFluidLoadToStructure();

			break;

		case HarmonicExtensionVelToFluid :

#ifdef DEBUG
			Debug( 5029 ) << "BCInterfaceFSI::checkFunction                          HarmonicExtensionVelToFluid" << "\n";
#endif

		    if ( !M_operator->isFluid() )
		    	return;

			M_operator->setHarmonicExtensionVelToFluid( M_operator->veloFluidMesh() );

			M_base = M_operator->bcvHarmonicExtensionVelToFluid();

			break;

		case SolidLoadToStructure :

#ifdef DEBUG
			Debug( 5029 ) << "BCInterfaceFSI::checkFunction                          SolidLoadToStructure" << "\n";
#endif
		    if ( !M_operator->isFluid() )
		    	return;

			M_operator->setSolidLoadToStructure( M_operator->minusSigmaFluidRepeated() );

			M_base = M_operator->bcvSolidLoadToStructure();

			break;

		case StructureDispToHarmonicExtension :

#ifdef DEBUG
			Debug( 5029 ) << "BCInterfaceFSI::checkFunction                          StructureDispToHarmonicExtension" << "\n";
#endif

		    if ( !M_operator->isFluid() )
		    	return;

			operMethod->setStructureDispToHarmonicExtension( M_operator->lambdaFluidRepeated() );

	        M_base = operMethod->bcvStructureDispToHarmonicExtension();

			break;

		case StructureDispToSolid :

#ifdef DEBUG
			Debug( 5029 ) << "BCInterfaceFSI::checkFunction                          StructureDispToSolid" << "\n";
#endif

			break;

		case StructureToFluid :

#ifdef DEBUG
			Debug( 5029 ) << "BCInterfaceFSI::checkFunction                          StructureToFluid" << "\n";
#endif

		    if ( !M_operator->isFluid() )
		    	return;

			M_operator->setStructureToFluid( M_operator->veloFluidMesh() );
			M_operator->setStructureToFluidParametres();

			M_base = M_operator->bcvStructureToFluid();

			break;
	}
}

} // Namespace LifeV

#endif /* __BCInterfaceFSI_H */
