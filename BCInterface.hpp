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





// ===================================================
//! Include
// ===================================================
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
#include <lifemc/lifefem/BCInterfaceFSIFunction.hpp>
#include <lifemc/lifefem/BCInterfaceFSIFunctionFile.hpp>
#include <lifemc/lifefem/BCInterfaceFSI.hpp>





// ===================================================
//! Namespaces & Static initializations
// ===================================================
namespace LifeV {





/*!
 * \class BCInterface
 * \brief LifeV Interface to load Boundary Conditions completely
 *        from a GetPot file
 *
 *  @author Cristiano Malossi
 *  @see
 *
 *  This class allows to impose boundary conditions completely from a file.
 *
 *
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
 *  To see some example look at test_cylinder and test_fsi.
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
 *
 */
template <class Operator>
class BCInterface
//     :
//     public LifeV::Application
{
public:

	typedef singleton< factory< BCInterfaceFunction<Operator>, std::string > >	BCInterfaceFunctionFactory;

    /** @name Constructors & Destructor
     */
    //@{

    //! Constructor
    /*!
     * \param dataFile		- GetPot data file
     * \param dataSection 	- Subsection inside [boundary_conditions]
     */
	BCInterface( const GetPot& dataFile, const std::string& dataSection );

	//! Copy constructor
	/*!
	 * \param interface		- BCInterface
	 */
	BCInterface( const BCInterface& interface );

    //! Destructor
    ~BCInterface() {}

    //@}



    /** @name Get functions
     */
    //@{

	const BCHandler& 						Handler() 		const { return *M_handler; }
	const boost::shared_ptr<BCHandler>& 	Handler_ptr() 	const { return  M_handler; } //Remove & ??

    //@}



    /** @name Methods
     */
    //@{

	//! Operator =
	/*!
	 * \param interface		- BCInterface
	 */
	BCInterface& operator=( const BCInterface& interface );

	//! Set an operator
    /*!
     * \param Oper			- operator
     */
	void setOperator( const boost::shared_ptr<Operator>& Oper ) { M_data.set_operator( Oper ); }

	//! Set manually Handler parameters: you need it only if you are adding manually some parameters by calling addBC
    /*!
     * \param bcNumber		- total number of the boundary conditions (files + added manually)
     * \param hint			- hint
     */
	void setHandlerParameters( const ID& bcNumber, const BCHandler::BCHints& hint = BCHandler::HINT_BC_NONE );

	//! Build the bcHandler
	void buildHandler( void );

    //@}



	/** @name External interface for BCHandler functions
	 */
	//@{

	//! Add a Boundary Condition
    /*!
     * \param name			- name of the condition
     * \param flag			- list of flags
     * \param type			- type of the condition
     * \param mode			- mode of the condition
     * \param base			- base of the condition
     */
	template <class BCBase>
	void addBC( const BCName& name,
				const BCFlag& flag,
				const BCType& type,
				const BCMode& mode,
					  BCBase& base );

	//! Add a Boundary Condition
    /*!
     * \param name			- name of the condition
     * \param flag			- list of flags
     * \param type			- type of the condition
     * \param mode			- mode of the condition
     * \param base			- base of the condition
     * \param comp			- component of the condition
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

    /** @name Private functions
     */
    //@{

    template <class BCInterfaceBase>
    inline bool newBase( std::map<std::string,size_type>& map, const std::vector< boost::shared_ptr<BCInterfaceBase> >& vector );

	template <class BCInterfaceBase>
	inline void addBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector );

	template <class BCInterfaceBase>
	inline void addBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector, const std::string& Oper );

	template <class BCBase>
	void addBCManager( BCBase& base );

	inline void setList( const char* conditions );

    inline void autosetHandlerParameters( void );

    inline void readFlag( const char* flag );

    inline void readType( const char* type );

    inline void readMode( const char* mode );

    inline void readComV( const char* component );

    inline void readBase( const std::string& base );

    inline bool isBase( const char* base );

    //@}

	enum BCBaseList{ function, functionFile, OSEENfunction, OSEENfunctionFile, FSIfunction, FSIfunctionFile, FSI };

	// GetPot data file
	GetPot 												M_dataFile;
	std::string	 										M_dataSection;

	std::vector<BCName>									M_list;
	UInt												M_listSize;
	bool												M_autoSetParameters;
	ID													M_bcNumber;

	// Handler and parameters
	BCHandler::BCHints 									M_hint;
	boost::shared_ptr<BCHandler> 						M_handler;

	// Maps
	std::map<std::string, BCType> 						M_mapType;
	std::map<std::string, BCMode> 						M_mapMode;
	std::map<std::string, BCBaseList> 					M_mapBase;

	// Data
	BCInterfaceData<Operator>							M_data;

	// Base
	std::pair<std::string, BCBaseList>					M_base;

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
	M_autoSetParameters				( true ),
	M_bcNumber						( 0 ),
	M_hint							( BCHandler::HINT_BC_ONLY_ESSENTIAL ),
	M_handler						( ),
	M_mapType						( ),
	M_mapMode						( ),
	M_mapBase						( ),
	M_data							( ),
	M_base							( )

{

#ifdef DEBUG
	    Debug( 5020 ) << "BCInterface::BCInterface------------------------------" << "\n";
#endif

	//Set mapType
	M_mapType["Essential"] 			= Essential;
	M_mapType["Natural"] 			= Natural;
	M_mapType["Mixte"] 				= Mixte;
	M_mapType["Flux"] 				= Flux;

	//Set mapMode
	M_mapMode["Scalar"] 			= Scalar;
	M_mapMode["Full"] 				= Full;
	M_mapMode["Component"] 			= Component;
	M_mapMode["Normal"] 			= Normal;
	M_mapMode["Tangential"] 		= Tangential;

	//Set mapBase
	M_mapBase["function"] 			= function;
	M_mapBase["functionFile"] 		= functionFile;
	M_mapBase["OSEENfunction"]		= OSEENfunction;
	M_mapBase["OSEENfunctionFile"]	= OSEENfunctionFile;
	M_mapBase["FSIfunction"]		= FSIfunction;
	M_mapBase["FSIfunctionFile"]	= FSIfunctionFile;
	M_mapBase["FSI"]				= FSI;

	//Set other parameters
	setList( (M_dataSection + "list").c_str() );

	//Factory registration
	BCInterfaceFunctionFactory::instance().registerProduct( "function", 			&createFunction<Operator> );
	BCInterfaceFunctionFactory::instance().registerProduct( "functionFile", 		&createFunctionFile<Operator> );
	BCInterfaceFunctionFactory::instance().registerProduct( "OSEENfunction", 		&createOperatorFunction<Operator> );
	BCInterfaceFunctionFactory::instance().registerProduct( "OSEENfunctionFile",	&createOperatorFunctionFile<Operator> );
	BCInterfaceFunctionFactory::instance().registerProduct( "FSIfunction", 			&createOperatorFunction<Operator> );
	BCInterfaceFunctionFactory::instance().registerProduct( "FSIfunctionFile", 		&createOperatorFunctionFile<Operator> );
}

template <class Operator>
BCInterface<Operator>::BCInterface( const BCInterface& interface ) :
	M_dataFile						( interface.M_dataFile ),
	M_dataSection					( interface.M_dataSection ),
	M_list							( interface.M_list ),
	M_listSize						( interface.M_listSize ),
	M_autoSetParameters				( interface.M_autoSetParameters ),
	M_bcNumber						( interface.M_bcNumber ),
	M_hint							( interface.M_hint ),
	M_handler						( interface.M_handler ),
	M_mapType						( interface.M_mapType ),
	M_mapMode						( interface.M_mapMode ),
	M_mapBase						( interface.M_mapBase ),
	M_data							( interface.M_data ),
	M_base							( interface.M_base )
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
    	M_autoSetParameters				= interface.M_autoSetParameters;
    	M_bcNumber						= interface.M_bcNumber;
    	M_hint							= interface.M_hint;
    	M_handler						= interface.M_handler;
    	M_mapType						= interface.M_mapType;
    	M_mapMode						= interface.M_mapMode;
    	M_mapBase						= interface.M_mapBase;
    	M_data							= interface.M_data;
    	M_base							= interface.M_base;
    }

	return *this;
}



template <class Operator>
void
BCInterface<Operator>::setHandlerParameters( const ID& bcNumber, const BCHandler::BCHints& hint )
{
	M_bcNumber 	= bcNumber;
	M_hint 		= hint;

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface::setHandlerParameters          M_bcNumber: " << M_bcNumber << "\n";
    Debug( 5020 ) << "                                               M_hint: " << M_hint << "\n";
#endif

	M_autoSetParameters = false;
}



template <class Operator>
void
BCInterface<Operator>::buildHandler( void )
{

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::buildHandler         M_autoSetParameters: " << M_autoSetParameters << "\n";
#endif

	if ( M_autoSetParameters )
		autosetHandlerParameters();

	M_handler.reset( new BCHandler ( M_bcNumber, M_hint ) );

	for ( UInt i(0) ; i < M_listSize ; ++i )
	{
		M_data.set_name( M_list[i] );

		readFlag( (M_dataSection + M_data.get_name() + "/flag").c_str() );
		readType( (M_dataSection + M_data.get_name() + "/type").c_str() );
		readMode( (M_dataSection + M_data.get_name() + "/mode").c_str() );
		readComV( (M_dataSection + M_data.get_name() + "/component").c_str() );
		readBase(  M_dataSection + M_data.get_name() + "/" );

		switch ( M_base.second )
		{
			case function :
			case functionFile :
			case OSEENfunction :
			case OSEENfunctionFile :
			case FSIfunction :
			case FSIfunctionFile :

				if ( newBase( M_mapFunction, M_vectorFunction ) )
					addBase( M_vectorFunction, M_base.first );

				addBCManager( M_vectorFunction[M_mapFunction[M_data.get_baseString()]]->getBase() );

				break;

			case FSI :

				if ( newBase( M_mapFSI, M_vectorFSI ) )
					addBase( M_vectorFSI );

				addBCManager( M_vectorFSI[M_mapFSI[M_data.get_baseString()]]->getBase() );

				break;

		}
	}
}



template <class Operator> template <class BCInterfaceBase>
inline bool
BCInterface<Operator>::newBase( std::map<std::string,size_type>& map, const std::vector< boost::shared_ptr<BCInterfaceBase> >& vector )
{
	//Check if the baseString has been already used
	for ( std::map<std::string, size_type>::iterator j = map.begin() ; j != map.end() ; ++j )
		if( vector[j->second]->compare( M_data ) )
		{

#ifdef DEBUG
			Debug( 5020 ) << "BCInterface::newBase                                   NO" << "\n";
#endif

			return false;
		}

	//Add baseString to the map
	size_type size = map.size();
	map[M_data.get_baseString()] = size;

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::newBase                                   YES" << "\n";
#endif

	return true;
}


template <class Operator> template <class BCInterfaceBase>
inline void
BCInterface<Operator>::addBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector )
{
	boost::shared_ptr<BCInterfaceBase> Function( new BCInterfaceBase( M_data ) );
	baseVector.push_back( Function );
}


template <class Operator> template <class BCInterfaceBase>
inline void
BCInterface<Operator>::addBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector, const std::string& Oper )
{
	boost::shared_ptr<BCInterfaceBase> Function( BCInterfaceFunctionFactory::instance().createObject( Oper ) );

	Function->setData( M_data );

	baseVector.push_back( Function );
}



template <class Operator> template <class BCBase>
void
BCInterface<Operator>::addBCManager( BCBase& base )
{
	switch ( M_data.get_mode() )
	{
		case Scalar :
		case Normal :
		case Tangential :

#ifdef DEBUG
			Debug( 5020 ) << "BCInterface::addBCManager                              Scalar, Normal, Tangential" << "\n\n";
#endif

			M_handler->addBC( M_data.get_name(), M_data.get_flag(), M_data.get_type(), M_data.get_mode(), base );

			break;

		case Full :

#ifdef DEBUG
			Debug( 5020 ) << "BCInterface::addBCManager                              Full" << "\n\n";
#endif

			M_handler->addBC( M_data.get_name(), M_data.get_flag(), M_data.get_type(), M_data.get_mode(), base, M_data.get_comN() );

			break;

		case Component :

#ifdef DEBUG
			Debug( 5020 ) << "BCInterface::addBCManager                              Component" << "\n\n";
#endif

			M_handler->addBC( M_data.get_name(), M_data.get_flag(), M_data.get_type(), M_data.get_mode(), base, M_data.get_comV() );

			break;
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
//! Private functions
// ===================================================
template <class Operator>
inline void
BCInterface<Operator>::setList( const char* conditions )
{
    M_listSize = M_dataFile.vector_variable_size( conditions );

    M_list.reserve( M_listSize );
    for ( UInt i(0) ; i < M_listSize ; ++i )
    	M_list.push_back(M_dataFile(conditions, " ", i));
}



template <class Operator>
inline void
BCInterface<Operator>::autosetHandlerParameters( void )
{
	for ( UInt i(0) ; i < M_listSize ; ++i )
	{
		readType( (M_dataSection + M_list[i] + "/type").c_str() );
		if ( M_data.get_type() != Essential )
			M_hint = BCHandler::HINT_BC_NONE;

		M_bcNumber += M_dataFile.vector_variable_size((M_dataSection + M_list[i] + "/flag").c_str());
	}

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface::autosetHandlerParameters      M_bcNumber: " << M_bcNumber << "\n";
    Debug( 5020 ) << "                                               M_hint: " << M_hint << "\n\n";
#endif

}



template <class Operator>
inline void
BCInterface<Operator>::readFlag( const char* flag )
{
	M_data.set_flag( M_dataFile(flag, 0) );

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface::readFlag                            flag: " << static_cast<Real>(M_data.get_flag()) << "\n";
#endif
}



template <class Operator>
inline void
BCInterface<Operator>::readType( const char* type )
{
	M_data.set_type( M_mapType[M_dataFile(type, "Essential")] );

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::readType                            type: " << M_data.get_type() << " (" << M_dataFile(type, "Essential") << ")\n";
#endif
}



template <class Operator>
inline void
BCInterface<Operator>::readMode( const char* mode )
{
	M_data.set_mode( M_mapMode[M_dataFile(mode, "Full")] );

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::readMode                            mode: " << M_data.get_mode() << " (" << M_dataFile(mode, "Full") << ")\n";
#endif
}



template <class Operator>
inline void
BCInterface<Operator>::readComV( const char* component )
{
    UInt componentSize = M_dataFile.vector_variable_size(component);

    M_data.reset_comV( componentSize );

    for (UInt j(0) ; j < componentSize ; ++j)
    	M_data.set_comV( M_dataFile(component, 0, j) );

#ifdef DEBUG
    std::stringstream output;
    output << "BCInterface::readComV                            comV: ";
    for (UInt i(0) ; i < static_cast<UInt>(M_data.get_comV().size()) ; ++i )
    	output << M_data.get_comV()[i] << " ";
    Debug( 5020 ) << output.str() << "\n";
#endif

}



template <class Operator>
inline void
BCInterface<Operator>::readBase( const std::string& base )
{
	for ( typename std::map<std::string, BCBaseList>::iterator j = M_mapBase.begin() ; j != M_mapBase.end() ; ++j )
		if ( isBase( (base + j->first).c_str() ) )
		{
			M_base.first  = j->first;
			M_base.second = M_mapBase[j->first];

#ifdef DEBUG
			Debug( 5020 ) << "BCInterface::readBase                            base: " << M_base.second << " (" << j->first << ")\n";
			Debug( 5020 ) << "                                           baseString: " << M_data.get_baseString() << "\n";
#endif

			break;
		}
}



template <class Operator>
inline bool
BCInterface<Operator>::isBase( const char* base )
{
	M_data.set_baseString( M_dataFile( base, " " ) );

	return M_dataFile.checkVariable( base );
}

} // Namespace LifeV

#endif /* __BCInterface_H */
