//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing the BCInterface3DData class
 *
 *  @date 17-07-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface3DData_H
#define BCInterface3DData_H 1

#include <lifemc/lifesolver/BCInterface3DDefinitions.hpp>

namespace LifeV
{

//! BCInterface3DData - The BCInterface data container
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterface3DData class provides a general container to pass information
 *  to all the BCInterface functions.
 */
class BCInterface3DData
{
public:

    //! @name Type definitions
    //@{

    typedef baseList3D_Type                                                            bcBaseList_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface3DData();

    //! Copy constructor
    /*!
     * @param data BCInterface3DData
     */
    BCInterface3DData( const BCInterface3DData& data );

    //! Destructor
    virtual ~BCInterface3DData() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param data BCInterface3DData
     * @return reference to a copy of the class
     */
    BCInterface3DData& operator=( const BCInterface3DData& data );

    //@}


    //! @name Methods
    //@{
    /*!
     * @param fileName Name of the data file.
     * @param dataSection BC section
     * @param name name of the boundary condition
     */
    void readBC( const std::string& fileName, const std::string& dataSection, const BCName& name );


    //! Display general information about the content of the class
    /*!
     * @param output specify the output format (std::cout by default)
     */
    void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the name of the boundary condition
    /*!
     * @param name Boundary condition name
     */
    void setName( const BCName& name ) { M_name = name; }

    //! Set the flag of the boundary condition
    /*!
     * @param flag Boundary condition flag
     */
    void setFlag( const BCFlag& flag ) { M_flag = flag; }

    //! Set the type of the boundary condition
    /*!
     * @param type Boundary condition type
     */
    void setType( const BCType& type ) { M_type = type; }

    //! Set the mode of the boundary condition
    /*!
     * @param mode Boundary condition mode
     */
    void setMode( const BCMode& mode ) { M_mode = mode; }

    //! Set the components vector of the boundary condition
    /*!
     * @param comV Boundary condition components vector
     */
    void setComV( const BCComV& comV ) { M_comV = comV; }

    //! Set the i-component of the components vector of the boundary condition
    /*!
     * @param comV Boundary condition component
     * @param index Index value
     */
    void setComV( const UInt& comV, const UInt& index ) { M_comV[index] = comV; }

    //! Add a component to the component vector of the boundary condition
    /*!
     * @param comV Boundary condition component
     */
    void addComV( const UInt& comV ) { M_comV.push_back( comV ); }

    //! Set the direction string of the boundary condition
    /*!
     * @param direction Boundary condition direction string
     */
    void setDirection( const std::string& direction ) { M_direction = direction; }

    //! Set the base string of the boundary condition
    /*!
     * @param baseString Boundary condition base string
     */
    void setBaseString( const std::string& baseString );

    //! Set the base type of the boundary condition
    /*!
     * @param base Boundary condition base type
     */
    void setBase( const std::pair< std::string, bcBaseList_Type >& base ) { M_base = base; }

    //@}


    //! @name Get Methods
    //@{

    //! Get the name of the boundary condition
    /*!
     * @return Boundary condition name
     */
    const BCName& name() const { return M_name; }

    //! Get the flag of the boundary condition
    /*!
     * @return Boundary condition flag
     */
    const BCFlag& flag() const { return M_flag; }

    //! Get the type of the boundary condition
    /*!
     * @return Boundary condition type
     */
    const BCType& type() const { return M_type; }

    //! Get the mode of the boundary condition
    /*!
     * @return Boundary condition mode
     */
    const BCMode& mode() const { return M_mode; }

    //! Get the vector of components of the boundary condition
    /*!
     * @return Boundary condition vector of components
     */
    const BCComV& comV() const { return M_comV; }

    //! Get the number of components of the boundary condition
    /*!
     * @return Number of components of the boundary condition
     */
    const ID& comN() const { return M_comV.front(); }

    //! Get the direction string of the boundary condition
    /*!
     * @return Boundary condition direction
     */
    const std::string& direction() const { return M_direction; }

    //! Get the base string of the boundary condition
    /*!
     * @return Boundary condition base string
     */
    const std::string& baseString() const { return M_baseString; }

    //! Get the base type of the boundary condition
    /*!
     * @return Boundary condition base
     */
    const std::pair< std::string, bcBaseList_Type >& base() const { return M_base; }

    //@}

private:

    //! @name Private Methods
    //@{

    void readFlag( const GetPot& dataFile, const char* flag ) { M_flag = dataFile( flag, 0 ); }

    void readType( const GetPot& dataFile, const char* type ) { M_type = M_mapType[dataFile( type, "Essential" )]; }

    void readMode( const GetPot& dataFile, const char* mode ) { M_mode = M_mapMode[dataFile( mode, "Full" )]; }

    void readComV( const GetPot& dataFile, const char* component );

    void readDirection( const GetPot& dataFile, const char* direction ) { M_direction = dataFile( direction, " " ); }

    void readBase( const GetPot& dataFile, const std::string& base );

    bool isBase( const GetPot& dataFile, const char* base );

    //@}

    BCName                                    M_name;
    BCFlag                                    M_flag;
    BCType                                    M_type;
    BCMode                                    M_mode;
    BCComV                                    M_comV;
    std::string                               M_direction;
    std::string                               M_baseString;
    std::pair< std::string, bcBaseList_Type > M_base;

    // Maps
    std::map< std::string, BCType >           M_mapType;
    std::map< std::string, BCMode >           M_mapMode;
    std::map< std::string, bcBaseList_Type >  M_mapBase;
};

} // Namespace LifeV

#endif /* BCInterface3DData_H */
