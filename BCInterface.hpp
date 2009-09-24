/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-04-01

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
   \file BCInterface.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-04-01
 */

#ifndef __BCInterface_H
#define __BCInterface_H 1

#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifefem/bcHandler.hpp>

#include <string>
#include <vector>

#include <lifemc/lifefem/BCInterfaceData.hpp>
#include <lifemc/lifefem/BCInterfaceFunction.hpp>
#include <lifemc/lifefem/BCInterfaceFunctionFile.hpp>
#include <lifemc/lifefem/BCInterfaceOperatorFunction.hpp>
#include <lifemc/lifefem/BCInterfaceOperatorFunctionFile.hpp>
#include <lifemc/lifefem/BCInterfaceFSI.hpp>

namespace LifeV {
//! BCInterface - LifeV Interface to load Boundary Conditions completely from a GetPot file
/*
 *  @author Cristiano Malossi
 *
 *  This class allows to impose boundary conditions completely from a file.
 *
 *  <b>EXAMPLE - DATA FILE</b>
 *
 *  In the GetPot data file, inside each subsection ([fluid], [solid], etc...) BCInterface reads a new section: [boundary_conditions].
 *
 *  Inside the new section there is a list of condition which correspond to other sub-section
 *  with the same name. The list must be inside the apex ' '.
 *
 *  Each condition has a similar structure; here there is an example:
 *
 *  [InFlow]
 *  type       = Essential
 *  flag       = 2
 *  mode       = Full
 *  component  = 3
 *  function   = '(0, 0, 3*0.03*(1/4-(x^2+y^2))'
 *
 *  NOTE: All the parameters are case sensitive.
 *
 *  type - can be: Essential Natural Mixte Flux
 *  flag - contains the flag
 *  mode - can be: Full Component Scalar Tangential Normal.
 *  component - if mode is Scalar, Tangential or Normal it is missing.
 *              if mode is Component it contains the ID of the component (or of the components list inside apex)
 *              if mode is Full it contains the total number of components
 *  function - contains the function. See BCInterfaceFunction for more details about the syntax.
 *
 *  <b>NOTE:</b>
 *
 *  The string "function" represent the base module and can be replaced by other expanded modules.
 *  Up to now we have the following modules for function:
 *
 *  - function
 *  - functionFile
 *  - FSIfunction
 *  - FSIfunctionFile
 *  - OSEENfunction
 *  - OSEENfunctionFile
 *  - FSI
 *
 *  To see some example look at test_fsi.
 *
 *  <b>EXAMPLE - HOW TO USE</b>
 *
 *  Here there is a short example on how to use it.
 *
 *  1) You can define your BCInterface class in a shared pointer:
 *     boost::shared_ptr<BCInterface> 	M_fluidBC;
 *
 *  2) You pass to the BCInterface the GetPot data file and the subsection to read inside [boundary_conditions]:
 *     M_fluidBC.reset( new BCInterface(data_file, "fluid") );
 *
 *  3) If you have operator conditions you have to add an operator
 *     M_fluidBC->setOperator( M_fsi->FSIOper() );
 *
 *  4) Then you can build the handler
 *     M_fluidBC->buildHandler();
 *
 *  5) Finally, to get the handler you can use:
 *     M_fluidBC->Handler_ptr();
 *
 *  NOTE:
 *
 *  You can add manually more conditions by using addBC() after the call to buildHandler() function.
 *  In this case you have to manually set the TOTAL number of boundary conditions
 *  by using setHandlerParameters() function BEFORE building the handler.
 */
template <class Operator>
class BCInterface
{
public:

	typedef singleton< factory< BCInterfaceFunction<Operator>, std::string > >	BCInterfaceFunctionFactory;

	//! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     * \param dataFile - GetPot data file
     * \param dataSection - Subsection inside [boundary_conditions]
     */
	BCInterface( const GetPot& dataFile, const std::string& dataSection );

	//! Copy constructor
	/*!
	 * \param interface	- BCInterface
	 */
	BCInterface( const BCInterface& interface );

    //! Destructor
    ~BCInterface() {}

    //@}



    //! @name Get functions
    //@{

    //! Get the BCHandler
	const BCHandler& 						Handler() 		const	{ return *M_handler; }

	//! Get the shared_ptr to the BCHandler
	const boost::shared_ptr<BCHandler>& 	Handler_ptr() 	const	{ return  M_handler; } //Remove & ??

	//! Get the data container
	      BCInterfaceData<Operator>& 		GetDataContainer()		{ return  M_data; }

    //@}


	//! @name Methods
    //@{

	//! Operator =
	/*!
	 * \param interface - BCInterface
	 */
	BCInterface& operator=( const BCInterface& interface );

	//! Set an operator
    /*!
     * \param Oper - operator
     */
	void SetOperator( const boost::shared_ptr<Operator>& Oper ) { M_data.SetOperator( Oper ); }

	//! Update the variables inside the operator
	void UpdateOperatorVariables( void );

	//! Set manually Handler parameters: you need it only if you are adding manually some parameters by calling addBC
    /*!
     * \param bcNumber - total number of the boundary conditions (files + added manually)
     * \param hint - hint
     */
	void SetHandlerParameters( const ID& bcNumber, const BCHandler::BCHints& hint = BCHandler::HINT_BC_NONE );

	//! Build the bcHandler with the data file provided with the constructor
	void BuildHandler( void );

	//! Read a boundary condition from a different file and add it to the data container
    /*!
     * \param name - name of the boundary condition
     * \param dataSection - section in the data file
     * \param dataFile - external data file
     */
	void ReadExternalBC( const BCName& name, const std::string& dataSection, const GetPot& dataFile );

	//! Insert the external boundary condition in the BChandler
	void InsertExternalBC( void );

    //@}



	//! @name External interface for BCHandler functions
	//@{

	//! Add a Boundary Condition using the standard interface of the BCHandler
    /*!
     * \param name - name of the condition
     * \param flag - list of flags
     * \param type - type of the condition
     * \param mode - mode of the condition
     * \param base - base of the condition
     */
	template <class BCBase>
	void addBC( const BCName& name,
				const BCFlag& flag,
				const BCType& type,
				const BCMode& mode,
					  BCBase& base );

	//! Add a Boundary Condition with component using the standard interface of the BCHandler
    /*!
     * \param name - name of the condition
     * \param flag - list of flags
     * \param type - type of the condition
     * \param mode - mode of the condition
     * \param base - base of the condition
     * \param comp - component of the condition
     */
	template <class BCBase, class BCComp>
	void addBC( const BCName& name,
				const BCFlag& flag,
				const BCType& type,
				const BCMode& mode,
					  BCBase& base,
				const BCComp& comp );

    //@}

private:

	//! @name Private Methods
    //@{

	inline void SetList( const char* conditions );

    inline void BuildBase( void );

    template <class BCInterfaceBase>
    inline bool NewBase( std::map<std::string,size_type>& map, const std::vector< boost::shared_ptr<BCInterfaceBase> >& vector );

	template <class BCInterfaceBase>
	inline void AddBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector );

	template <class BCInterfaceBase>
	inline void AddBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector, const std::string& Oper );

	template <class BCBase>
	inline void AddBCManager( BCBase& base );

    //@}

	// GetPot data file
	GetPot 												M_dataFile;
	std::string	 										M_dataSection;

	std::vector<BCName>									M_list;
	UInt												M_listSize;
	ID													M_bcNumber;

	// Handler and parameters
	BCHandler::BCHints 									M_hint;
	boost::shared_ptr<BCHandler> 						M_handler;

	// Data
	BCInterfaceData<Operator>							M_data;

	// Functions
	static std::map<std::string,size_type>										M_mapFunction;
	static std::vector< boost::shared_ptr<BCInterfaceFunction<Operator> > >		M_vectorFunction;

	// FSI Functions
	static std::map<std::string,size_type>										M_mapFSI;
	static std::vector< boost::shared_ptr<BCInterfaceFSI<Operator> > >			M_vectorFSI;
};





// ===================================================
//! Template implementation
// ===================================================
// Initialize static variables
template <class Operator>
std::map<std::string,size_type>										BCInterface<Operator>::M_mapFunction;
template <class Operator>
std::vector< boost::shared_ptr<BCInterfaceFunction<Operator> > >	BCInterface<Operator>::M_vectorFunction;

template <class Operator>
std::map<std::string,size_type>										BCInterface<Operator>::M_mapFSI;
template <class Operator>
std::vector< boost::shared_ptr<BCInterfaceFSI<Operator> > >			BCInterface<Operator>::M_vectorFSI;



// ===================================================
//! Constructors
// ===================================================
template <class Operator>
BCInterface<Operator>::BCInterface( const GetPot& dataFile, const std::string& dataSection ) :
	M_dataFile						( dataFile ),
	M_dataSection					( dataSection + "/boundary_conditions/" ),
	M_list							( ),
	M_listSize						( 0 ),
	M_bcNumber						( 0 ),
	M_hint							( BCHandler::HINT_BC_NONE ),
	M_handler						( ),
	M_data							( )
{

#ifdef DEBUG
	    Debug( 5020 ) << "BCInterface::BCInterface------------------------------" << "\n";
#endif

	//Set other parameters
	SetList( (M_dataSection + "list").c_str() );

	//Factory registration
	BCInterfaceFunctionFactory::instance().registerProduct( "function", 		&createFunction<Operator> );
	BCInterfaceFunctionFactory::instance().registerProduct( "functionFile", 	&createFunctionFile<Operator> );
	BCInterfaceFunctionFactory::instance().registerProduct( "OPERfunction", 	&createOperatorFunction<Operator> );
	BCInterfaceFunctionFactory::instance().registerProduct( "OPERfunctionFile",	&createOperatorFunctionFile<Operator> );
}

template <class Operator>
BCInterface<Operator>::BCInterface( const BCInterface& interface ) :
	M_dataFile						( interface.M_dataFile ),
	M_dataSection					( interface.M_dataSection ),
	M_list							( interface.M_list ),
	M_listSize						( interface.M_listSize ),
	M_bcNumber						( interface.M_bcNumber ),
	M_hint							( interface.M_hint ),
	M_handler						( interface.M_handler ),
	M_data							( interface.M_data )
{
}



// ===================================================
//! Methods
// ===================================================
template <class Operator>
BCInterface<Operator>&
BCInterface<Operator>::operator=( const BCInterface& interface )
{
    if ( this != &interface )
    {
    	M_dataFile						= interface.M_dataFile;
    	M_dataSection					= interface.M_dataSection;
    	M_list							= interface.M_list;
    	M_listSize						= interface.M_listSize;
    	M_bcNumber						= interface.M_bcNumber;
    	M_hint							= interface.M_hint;
    	M_handler						= interface.M_handler;
    	M_data							= interface.M_data;
    }

	return *this;
}

template <class Operator>
void
BCInterface<Operator>::SetHandlerParameters( const ID& bcNumber, const BCHandler::BCHints& hint )
{
	M_bcNumber 	= bcNumber;
	M_hint 		= hint;

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface::setHandlerParameters          M_bcNumber: " << M_bcNumber << "\n";
    Debug( 5020 ) << "                                               M_hint: " << M_hint << "\n";
#endif

}

template <class Operator>
void
BCInterface<Operator>::BuildHandler( void )
{

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::buildHandler\n";
#endif

	M_handler.reset( new BCHandler ( M_bcNumber, M_hint ) );

	for ( UInt i(0) ; i < M_listSize ; ++i )
	{
		M_data.ReadBC( M_list[i], M_dataSection, M_dataFile );

		BuildBase();
	}
}

template <class Operator>
void
BCInterface<Operator>::ReadExternalBC( const BCName& name, const std::string& dataSection, const GetPot& dataFile )
{

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::readExternalBC\n";
#endif

	M_data.ReadBC( name, dataSection, dataFile );
}

template <class Operator>
void
BCInterface<Operator>::InsertExternalBC( void )
{

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::insertExternalBC\n";
#endif

	BuildBase();
}

template <class Operator>
void
BCInterface<Operator>::UpdateOperatorVariables( void )
{

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::UpdateOperatorVariables\n";
#endif

	for ( UInt i(0) ; i < M_vectorFunction.size() ; ++i )
	{
		BCInterfaceOperatorFunction<Operator> *Oper = dynamic_cast<BCInterfaceOperatorFunction<Operator> *> (&(*M_vectorFunction[i]));

		if ( Oper != 0 )
			Oper->UpdateOperatorVariables();
	}
}

template <class Operator> template <class BCBase>
void
BCInterface<Operator>::addBC( 	const BCName& name,
								const BCFlag& flag,
								const BCType& type,
								const BCMode& mode,
									  BCBase& base )
{

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::addBC (without component)" << "\n\n";
#endif

		M_handler->addBC( name, flag, type, mode, base );
}



template <class Operator> template <class BCBase, class BCComp>
void
BCInterface<Operator>::addBC( 	const BCName& name,
								const BCFlag& flag,
								const BCType& type,
								const BCMode& mode,
									  BCBase& base,
								const BCComp& comp )
{

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::addBC (with component)" << "\n\n";
#endif

		M_handler->addBC( name, flag, type, mode, base, comp );
}



// ===================================================
//! Private Methods
// ===================================================
template <class Operator>
inline void
BCInterface<Operator>::SetList( const char* conditions )
{
    M_listSize = M_dataFile.vector_variable_size( conditions );

    M_list.reserve( M_listSize );
    for ( UInt i(0) ; i < M_listSize ; ++i )
    	M_list.push_back(M_dataFile(conditions, " ", i));

    M_bcNumber = M_listSize;
}

template <class Operator>
inline void
BCInterface<Operator>::BuildBase( void )
{

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::BuildBase\n";
#endif

	switch ( M_data.GetBase().second )
	{
		case function :
		case functionFile :
		case OPERfunction :
		case OPERfunctionFile :

			if ( NewBase( M_mapFunction, M_vectorFunction ) )
				AddBase( M_vectorFunction, M_data.GetBase().first );

			AddBCManager( M_vectorFunction[M_mapFunction[M_data.GetBaseString()]]->GetBase() );

			break;

		case FSI :

			if ( NewBase( M_mapFSI, M_vectorFSI ) )
				AddBase( M_vectorFSI );

			AddBCManager( M_vectorFSI[M_mapFSI[M_data.GetBaseString()]]->GetBase() );

			break;
	}
}



template <class Operator> template <class BCInterfaceBase>
inline bool
BCInterface<Operator>::NewBase( std::map<std::string,size_type>& map, const std::vector< boost::shared_ptr<BCInterfaceBase> >& vector )
{
	//Check if the baseString has been already used
	for ( std::map<std::string, size_type>::iterator j = map.begin() ; j != map.end() ; ++j )
		if( vector[j->second]->Compare( M_data ) )
		{

#ifdef DEBUG
			Debug( 5020 ) << "BCInterface::NewBase                                   NO" << "\n";
#endif

			return false;
		}

	//Add baseString to the map
	size_type size = map.size();
	map[M_data.GetBaseString()] = size;

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::NewBase                                   YES" << "\n";
#endif

	return true;
}


template <class Operator> template <class BCInterfaceBase>
inline void
BCInterface<Operator>::AddBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector )
{
	boost::shared_ptr<BCInterfaceBase> Function( new BCInterfaceBase( M_data ) );
	baseVector.push_back( Function );
}


template <class Operator> template <class BCInterfaceBase>
inline void
BCInterface<Operator>::AddBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector, const std::string& Oper )
{
	boost::shared_ptr<BCInterfaceBase> Function( BCInterfaceFunctionFactory::instance().createObject( Oper ) );

	Function->SetData( M_data );

	baseVector.push_back( Function );
}



template <class Operator> template <class BCBase>
inline void
BCInterface<Operator>::AddBCManager( BCBase& base )
{
	switch ( M_data.GetMode() )
	{
		case Scalar :
		case Normal :
		case Tangential :

#ifdef DEBUG
			Debug( 5020 ) << "BCInterface::AddBCManager                              Scalar, Normal, Tangential" << "\n\n";
#endif

			M_handler->addBC( M_data.GetName(), M_data.GetFlag(), M_data.GetType(), M_data.GetMode(), base );

			break;

		case Full :

#ifdef DEBUG
			Debug( 5020 ) << "BCInterface::AddBCManager                              Full" << "\n\n";
#endif

			M_handler->addBC( M_data.GetName(), M_data.GetFlag(), M_data.GetType(), M_data.GetMode(), base, M_data.GetComN() );

			break;

		case Component :

#ifdef DEBUG
			Debug( 5020 ) << "BCInterface::AddBCManager                              Component" << "\n\n";
#endif

			M_handler->addBC( M_data.GetName(), M_data.GetFlag(), M_data.GetType(), M_data.GetMode(), base, M_data.GetComV() );

			break;
	}
}

} // Namespace LifeV

#endif /* __BCInterface_H */
