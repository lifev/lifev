/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-07-17

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
   \file BCInterfaceData.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-07-17
 */

#ifndef __BCInterfaceData_H
#define __BCInterfaceData_H 1





// ===================================================
//! Include
// ===================================================
#include <life/lifecore/life.hpp>
#include <life/lifefem/bcCond.hpp>

#include <boost/algorithm/string.hpp>
#include <string>





// ===================================================
//! Namespaces & Typedef
// ===================================================
namespace LifeV {

typedef std::string									BCName;
typedef EntityFlag									BCFlag;
typedef std::vector<ID>								BCComV;





/*!
 * \class BCInterfaceData
 * \brief A simple container for main BC data read from the GetPot file.
 *
 *  @author Cristiano Malossi
 *  @see
 *
 */
template <class Operator>
class BCInterfaceData
{
public:

    /** @name Constructors & Destructor
     */
    //@{

    //! Constructor
	BCInterfaceData( );

	//! Copy constructor
	/*!
	 * \param data - BCInterfaceData
	 */
	BCInterfaceData( const BCInterfaceData& data );

	//! Operator =
	/*!
	 * \param data - BCInterfaceData
	 */
	BCInterfaceData& operator=( const BCInterfaceData& data );

    //! Destructor
    ~BCInterfaceData() {}

    //@}



    /** @name Methods
     */
    //@{
	/*!
	 * \param components - number of components to reserve
	 */
    void reset_comV( const UInt& components = 0 );

    //@}



    /** @name Set functions
     */
    //@{

	void set_operator( const boost::shared_ptr<Operator>& Oper) { M_operator = Oper; }
	void set_name( const BCName& name ) 						{ M_name = name; }
	void set_flag( const BCFlag& flag ) 						{ M_flag = flag; }
	void set_type( const BCType& type ) 						{ M_type = type; }
	void set_mode( const BCMode& mode ) 						{ M_mode = mode; }
	void set_comV( const BCComV& comV ) 						{ M_comV = comV; }
	void set_comV( const UInt&   comV ) 						{ M_comV.push_back( comV ); }
	void set_comV( const UInt&   comV, const UInt& index ) 		{ M_comV[index] = comV; }
	void set_baseString( const std::string& baseString ) 		{
																	M_baseString = baseString;
																	boost::replace_all( M_baseString, " ",  "" );
																}

    //@}



    /** @name Get functions
     */
    //@{

	const boost::shared_ptr<Operator>& 	get_operator()		const { return M_operator; }
	const BCName& 						get_name() 			const { return M_name; }
	const BCFlag& 						get_flag() 			const { return M_flag; }
	const BCType& 						get_type() 			const { return M_type; }
	const BCMode& 						get_mode() 			const { return M_mode; }
	const BCComV& 						get_comV() 			const { return M_comV; }
	const ID& 							get_comN() 			const { return M_comV.front(); }
	const std::string& 					get_baseString()	const { return M_baseString; }

    //@}

private:

	boost::shared_ptr<Operator>							M_operator;
	BCName												M_name;
	BCFlag 												M_flag;
	BCType 												M_type;
	BCMode 												M_mode;
	BCComV 												M_comV;
	std::string											M_baseString;
};





// ===================================================
//! Constructors
// ===================================================
template <class Operator>
BCInterfaceData<Operator>::BCInterfaceData( ) :
	M_operator		( ),
	M_name			( ),
	M_flag			( ),
	M_type			( ),
	M_mode			( ),
	M_comV			( ),
	M_baseString	( )
{}



template <class Operator>
BCInterfaceData<Operator>::BCInterfaceData( const BCInterfaceData& data ) :
	M_operator		( data.M_operator ),
	M_name			( data.M_name ),
	M_flag			( data.M_flag ),
	M_type			( data.M_type ),
	M_mode			( data.M_mode ),
	M_comV			( data.M_comV ),
	M_baseString	( data.M_baseString )
{}


template <class Operator>
BCInterfaceData<Operator>&
BCInterfaceData<Operator>::operator=( const BCInterfaceData& data )
{
	if ( this != &data )
	{
		M_operator		= data.M_operator;
		M_name			= data.M_name;
		M_flag			= data.M_flag;
		M_type			= data.M_type;
		M_mode			= data.M_mode;
		M_comV			= data.M_comV;
		M_baseString	= data.M_baseString;
	}

	return *this;
}





// ===================================================
//! Methods
// ===================================================
template <class Operator>
void
BCInterfaceData<Operator>::reset_comV( const UInt& components )
{
	M_comV.clear();
	M_comV.reserve( components );
}

} // Namespace LifeV

#endif /* __BCInterfaceData_H */
