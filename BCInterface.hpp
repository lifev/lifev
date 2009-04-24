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
//! Namespaces
// ===================================================
using namespace LifeV;
enum BCBaseList{function, fsi};





/*!
 * \class BCInterface
 * \brief LifeV Interface to load Boundary Conditions completely
 *        from a GetPot file
 *
 *  @author Cristiano Malossi
 *  @see
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
	typedef UInt										BCComN;
	typedef std::vector<ID>								BCComV;



	// ===================================================
	//! Public functions
	// ===================================================

    /** @name Constructors & Destructor
     */
    //@{

    //! Constructor
    /*!
      \param GetPot data file
      \param GetPot section
    */
	BCInterface( GetPot const& dataFile, std::string dataSection );

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

	// Set Handler parameters
	void setHandlerParameters( const ID bcNumber, const BCHandler::BCHints hint = BCHandler::HINT_BC_NONE );

	// Build BChandler
	void buildHandler( void );

	// Add BC
	template <class BCBase>
	void addBC( const BCName& name,
				const BCFlag& flag,
				const BCType& type,
				const BCMode& mode,
					  BCBase& base );

	template <class BCBase, class BCComp>
	void addBC( const BCName& name,
				const BCFlag& flag,
				const BCType& type,
				const BCMode& mode,
					  BCBase& base,
				const BCComp&  comp );

	// Set FSI operator
	void setFSIOperator(const boost::shared_ptr<FSIOperator>& oper );


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
	ID													M_bcNumber;

	// Handler and parameters
	boost::shared_ptr<BCHandler> 						M_handler;
	bool												M_autoSetParameters;
	BCHandler::BCHints 									M_hint;

	std::map<std::string, BCType> 						M_mapType;
	std::map<std::string, BCMode> 						M_mapMode;
	std::map<std::string, BCBaseList> 					M_mapBase;

	// Operators
	std::vector< boost::shared_ptr<BCInterfaceFunction> > 		M_functionVector;
	std::vector< boost::shared_ptr<BCInterfaceFSIOperator> > 	M_FSIOperatorVector;

	boost::shared_ptr<FSIOperator>						M_FSIOperator;

	// BC options
	BCName												M_name;
	BCFlag 												M_flag;
	BCType 												M_type;
	BCMode 												M_mode;
	BCComN 												M_comN;
	BCComV 												M_comV;

	BCBaseList 											M_base;
	std::string											M_baseString;



	// ===================================================
	//! Private functions
	// ===================================================

    /** @name Private functions
     */
    //@{

	//! addBase
	template <class BCInterfaceBase>
	inline void addBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector );

	//! addBase
	template <class BCInterfaceBase, class BCOperator>
	inline void addBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector, BCOperator& Operator );

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
    inline void readComponentNumber( const char* component );

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
template <class BCInterfaceBase>
inline void
BCInterface::addBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector )
{
	boost::shared_ptr<BCInterfaceBase> Function( new BCInterfaceBase( M_baseString ) );
	baseVector.push_back( Function );
}



template <class BCInterfaceBase, class BCOperator>
inline void
BCInterface::addBase( std::vector< boost::shared_ptr<BCInterfaceBase> >& baseVector, BCOperator& Operator )
{
	boost::shared_ptr<BCInterfaceBase> Function( new BCInterfaceBase( M_baseString, Operator ) );
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

			addBC( M_name, M_flag, M_type, M_mode, base);

			break;

		case Full :

			readComponentNumber( (M_dataSection + M_name + "/component").c_str() );
			addBC( M_name, M_flag, M_type, M_mode, base, M_comN);

			break;

		case Component :

			readComponentVector( (M_dataSection + M_name + "/component").c_str() );
			addBC( M_name, M_flag, M_type, M_mode, base, M_comV);

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
	for (UInt j(0) ; j < flag.size() ; ++j)
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
	for (UInt j(0) ; j < flag.size() ; ++j)
		M_handler->addBC( name, flag[j], type, mode, base, comp );
}

#endif /* __BCInterface_H */
