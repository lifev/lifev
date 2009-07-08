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

#include <lifemc/lifefem/BCInterfaceFunction.hpp>
#include <lifemc/lifefem/BCInterfaceFSIOperator.hpp>





// ===================================================
//! Namespaces & Enums
// ===================================================
namespace LifeV {

enum BCBaseList{function, fsi};





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
 *  type - can be: Essential Natural Mixte
 *  flag - contains the flag (or the list of flags inside apex '...')
 *  mode - can be: Full Component Scalar Tangential Normal.
 *  component - if mode is Scalar, Tangential or Normal it is missing.
 *              if mode is Component it contains the ID of the component (or of the components list inside apex)
 *              if mode is Full it contains the total number of components
 *  function - contains the function. See BCInterfaceFunction for more details about the syntax.
 *
 *	In the case of FSI problems the syntax is exactly the same but function is replaced by fsi:
 *
 *	[../Interface]
 *	type       = Essential
 *	flag       = 1
 *	mode       = Full
 *	component  = 3
 *	fsi	       = StructureToFluid
 *
 *	where fsi contains the fsi operator. See BCInterfaceFSIOperator for more details about the list of available operators.
 *
 *  NOTE:
 *
 *  To see some example look at test_cylinder and test_fsi.
 *
 *
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
 *  3) If you have fsi conditions you have to add an operator
 *     M_fluidBC->setFSIOperator( M_fsi->operFSI() );
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
 *  \TODO Make static BCInterfaceFunction, to have just one class (and hence one Parser) for the same function strings.
 *  \TODO Find a way to impose component in a more general way (actually to impose 3rd component, we have to provide also the first two: (0,0,'3rd component')
 */
class BCInterface
//     :
//     public LifeV::Application
{
public:

	// ===================================================
	//! Typedef
	// ===================================================

	typedef std::string									BCName;
	typedef std::vector<EntityFlag>						BCFlag;
	//typedef UInt										BCComN;
	typedef std::vector<ID>								BCComV;



	// ===================================================
	//! Public functions
	// ===================================================

    /** @name Constructors & Destructor
     */
    //@{

    //! Constructor
    /*!
     * \param dataFile    - GetPot data file
     * \param dataSection - Subsection inside [conditions]
     */
	BCInterface( GetPot const& dataFile, const std::string dataSection );

	//! Copy constructor
	/*!
	 * \param interface - BCInterface
	 */
	BCInterface( const BCInterface& interface );

	//! Operator =
	/*!
	 * \param interface - BCInterface
	 */
	BCInterface& operator=( const BCInterface& interface );

    //! Destructor
    ~BCInterface() {}

    //@}



    /** @name Get functions
     */
    //@{

	const BCHandler& 						Handler() 		const { return *M_handler; }
	const boost::shared_ptr<BCHandler>& 	Handler_ptr() 	const { return  M_handler; }

    //@}



    /** @name Members functions
     */
    //@{

	//! Set manually Handler parameters: you need it only if you are adding manually some parameters by calling addBC
    /*!
     * \param bcNumber - total number of the boundary conditions (files + added manually)
     * \param hint     - hint
     */
	void setHandlerParameters( const ID bcNumber, const BCHandler::BCHints hint = BCHandler::HINT_BC_NONE );

	//! Build the bcHandler
	void buildHandler( void );

	//! Add a Boundary Condition
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

	//! Add a Boundary Condition
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

	//! Set an FSIOperator (need only for fsi BC)
    /*!
     * \param oper - FSIOperator
     */
	void setFSIOperator( const boost::shared_ptr<FSIOperator>& oper ) { M_FSIOperator = oper; }


    //@}

private:

	// ===================================================
	//! Private variables
	// ===================================================

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

	boost::shared_ptr<FSIOperator>						M_FSIOperator;

	// Operators
	std::vector< boost::shared_ptr<BCInterfaceFunction> > 		M_functionVector;
	std::vector< boost::shared_ptr<BCInterfaceFSIOperator> > 	M_FSIOperatorVector;

	// BC options
	BCName												M_name;
	BCFlag 												M_flag;
	BCType 												M_type;
	BCMode 												M_mode;
	//BCComN 												M_comN;
	BCComV 												M_comV;

	BCBaseList 											M_base;
	std::string											M_baseString;



	// ===================================================
	//! Private functions
	// ===================================================

    /** @name Private functions
     */
    //@{

	/*
	//! addBase
	template <class BCInterfaceBase>
	inline void addBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector );
	*/

	//! addBase
	template <class BCInterfaceBase, class BCparameter>
	inline void addBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector, BCparameter& p );

	//! addBCManager
	template <class BCBase>
	void addBCManager( BCBase& base );

    //! setList
	inline void setList( const char* conditions );

    //! setParameters
    inline void autosetHandlerParameters( void );

    //! readFlag
    inline void readFlag( const char* flag );

    //! readType
    inline void readType( const char* type );

    //! readMode
    inline void readMode( const char* mode );

    //! readComponentNumber
    //inline void readComponentNumber( const char* component );

    //! readComponentVector
    inline void readComponentVector( const char* component );

    //! readBase
    inline void readBase( const std::string base );

    //! isBase
    inline bool isBase( const char* base );

    //@}
};



// ===================================================
//! Template function
// ===================================================
/*
template <class BCInterfaceBase>
inline void
BCInterface::addBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector )
{
	boost::shared_ptr<BCInterfaceBase> Function( new BCInterfaceBase( M_baseString ) );
	baseVector.push_back( Function );
}
*/


template <class BCInterfaceBase, class BCparameter>
inline void
BCInterface::addBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector, BCparameter& p )
{
	boost::shared_ptr<BCInterfaceBase> Function( new BCInterfaceBase( M_baseString, p ) );
	baseVector.push_back( Function );
}



template <class BCBase>
void BCInterface::addBCManager( BCBase& base )
{
	switch ( M_mode )
	{
		case Scalar :
		case Normal :
		case Tangential :

#ifdef DEBUG
			Debug( 5020 ) << "BCInterface::addBCManager (Scalar, Normal, Tangential)" << "\n";
#endif

			addBC( M_name, M_flag, M_type, M_mode, base );

			break;

		case Full :

			//readComponentNumber( (M_dataSection + M_name + "/component").c_str() );
			//readComponentVector( (M_dataSection + M_name + "/component").c_str() );

#ifdef DEBUG
			Debug( 5020 ) << "BCInterface::addBCManager (Full)" << "\n";
#endif

			addBC( M_name, M_flag, M_type, M_mode, base, M_comV.front() );

			break;

		case Component :

			//readComponentVector( (M_dataSection + M_name + "/component").c_str() );

#ifdef DEBUG
			Debug( 5020 ) << "BCInterface::addBCManager (Component)" << "\n";
#endif

			addBC( M_name, M_flag, M_type, M_mode, base, M_comV );

			break;
	}
}



template <class BCBase>
void BCInterface::addBC( 	const BCName& name,
							const BCFlag& flag,
							const BCType& type,
							const BCMode& mode,
								  BCBase& base )
{

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::addBC (without component)" << "\n\n";
#endif

	for ( UInt j(0) ; j < flag.size() ; ++j )
		M_handler->addBC( name, flag[j], type, mode, base );
}



template <class BCBase, class BCComp>
void BCInterface::addBC( 	const BCName& name,
							const BCFlag& flag,
							const BCType& type,
							const BCMode& mode,
								  BCBase& base,
							const BCComp& comp )
{

#ifdef DEBUG
	Debug( 5020 ) << "BCInterface::addBC (with component)" << "\n\n";
#endif

	for ( UInt j(0) ; j < flag.size() ; ++j )
		M_handler->addBC( name, flag[j], type, mode, base, comp );
}

} // Namespace LifeV

#endif /* __BCInterface_H */
