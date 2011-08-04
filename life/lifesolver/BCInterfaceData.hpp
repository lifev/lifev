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
 *  @brief File containing the BCInterfaceData class
 *
 *  @date 17-07-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceData_H
#define BCInterfaceData_H 1

#include <lifemc/lifesolver/BCInterfaceDefinitions.hpp>

namespace LifeV
{

//! BCInterfaceData - The BCInterface data container
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterfaceData class provides a general container to pass information
 *  to all the BCInterface functions.
 */
class BCInterfaceData
{
public:

    //! @name Type definitions
    //@{

    typedef std::vector< Real >                                                        resistanceContainer_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceData();

    //! Copy constructor
    /*!
     * @param data BCInterfaceData
     */
    BCInterfaceData( const BCInterfaceData& data );

    //! Destructor
    virtual ~BCInterfaceData() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param data BCInterfaceData
     * @return reference to a copy of the class
     */
    BCInterfaceData& operator=( const BCInterfaceData& data );

    //@}


    //! @name Methods
    //@{

    //! Read parameters for all kind of BC
    /*!
     * @param fileName Name of the data file.
     * @param dataSection BC section
     * @param name name of the boundary condition
     */
    void readBC( const std::string& fileName, const std::string& dataSection, const bcName_Type& name );

    //! Set the directional base as the current base
    void setDirectionalBase() { M_base = M_baseDirectional; M_baseString = M_baseStringDirectional; }

    //! Set the Robin Alpha base as the current base
    void setRobinBaseAlpha() { M_base = M_baseRobinAlpha; M_baseString = M_baseStringRobinAlpha; }

    //! Set the Robin Beta base as the current base
    void setRobinBaseBeta() { M_base = M_baseRobinBeta; M_baseString = M_baseStringRobinBeta; }

    //! Display general information about the content of the class
    /*!
     * @param output specify the output format (std::cout by default)
     */
    void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the base string of the boundary condition
    /*!
     * @param baseString Boundary condition base string
     */
    void setBaseString( const std::string& baseString );

    //! Set the base type of the boundary condition
    /*!
     * @param base Boundary condition base type
     */
    void setBase( const std::pair< std::string, baseList_Type >& base ) { M_base = base; }

    //! Set the side of the boundary condition
    /*!
     * @param flag Boundary condition side
     */
    void setSide( const OneDimensional::bcSide_Type& side ) { M_side = side; }

    //! Set the line of the boundary condition
    /*!
     * @param line Boundary condition line
     */
    void setLine( const OneDimensional::bcLine_Type& line ) { M_line = line; }

    //! Set the quantity of the boundary condition
    /*!
     * @param quantity Boundary condition quantity
     */
    void setQuantity( const OneDimensional::bcType_Type& quantity ) { M_quantity = quantity; }

    //! Set the name of the boundary condition
    /*!
     * @param name Boundary condition name
     */
    void setName( const bcName_Type& name ) { M_name = name; }

    //! Set the flag of the boundary condition
    /*!
     * @param flag Boundary condition flag
     */
    void setFlag( const bcFlag_Type& flag ) { M_flag = flag; }

    //! Set the type of the boundary condition
    /*!
     * @param type Boundary condition type
     */
    void setType( const bcType_Type& type ) { M_type = type; }

    //! Set the mode of the boundary condition
    /*!
     * @param mode Boundary condition mode
     */
    void setMode( const bcMode_Type& mode ) { M_mode = mode; }

    //! Set the components vector of the boundary condition
    /*!
     * @param comV Boundary condition components vector
     */
    void setComV( const bcComponentsVec_Type& comV ) { M_comV = comV; }

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

    //@}


    //! @name Get Methods
    //@{

    //! Get the base string of the boundary condition
    /*!
     * @return Boundary condition base string
     */
    const std::string& baseString() const { return M_baseString; }

    //! Get the base type of the boundary condition
    /*!
     * @return Boundary condition base
     */
    const std::pair< std::string, baseList_Type >& base() const { return M_base; }

    //! Get the base map of the boundary condition
    /*!
     * @return Boundary condition base map
     */
    const std::map< std::string, baseList_Type >& mapBase() const { return M_mapBase; }

    //! Get the flag of the boundary condition
    /*!
     * @return Boundary condition side
     */
    const OneDimensional::bcSide_Type& side() const { return M_side; }

    //! Get the mode of the boundary condition
    /*!
     * @return Boundary condition line
     */
    const OneDimensional::bcLine_Type& line() const { return M_line; }

    //! Get the quantity of the boundary condition
    /*!
     * @return Boundary condition quantity
     */
    const OneDimensional::bcType_Type& quantity() const { return M_quantity; }

    //! Get the resistance vector {R1, R2, R3 ...}
    /*!
     * @return Boundary condition resistance vector
     */
    const resistanceContainer_Type& resistance() const { return M_resistance; }

    //! Get the capacitance
    /*!
     * @return Boundary condition capacitance
     */
    const Real& capacitance() const { return M_capacitance; }

    //! Get the name of the boundary condition
    /*!
     * @return Boundary condition name
     */
    const bcName_Type& name() const { return M_name; }

    //! Get the flag of the boundary condition
    /*!
     * @return Boundary condition flag
     */
    const bcFlag_Type& flag() const { return M_flag; }

    //! Get the type of the boundary condition
    /*!
     * @return Boundary condition type
     */
    const bcType_Type& type() const { return M_type; }

    //! Get the mode of the boundary condition
    /*!
     * @return Boundary condition mode
     */
    const bcMode_Type& mode() const { return M_mode; }

    //! Get the vector of components of the boundary condition
    /*!
     * @return Boundary condition vector of components
     */
    const bcComponentsVec_Type& comV() const { return M_comV; }

    //! Get the number of components of the boundary condition
    /*!
     * @return Number of components of the boundary condition
     */
    const ID& comN() const { return M_comV.front(); }

    //@}

private:

    //! @name Private Methods
    //@{

    void readSide( const GetPot& dataFile, const char* side ) {  M_side = M_mapSide[dataFile( side, "left" )]; }

    void readLine( const GetPot& dataFile, const char* line ) { M_line = M_mapLine[dataFile( line, "first" )]; }

    void readQuantity( const GetPot& dataFile, const char* quantity ) { M_quantity = M_mapQuantity[dataFile( quantity, "A" )]; }

    void readResistance( const GetPot& dataFile, const char* resistance );

    void readCapacitance( const GetPot& dataFile, const char* capacitance ) { M_capacitance = dataFile( capacitance, 0 ); }

    void readFlag( const GetPot& dataFile, const char* flag ) { M_flag = dataFile( flag, 0 ); }

    void readType( const GetPot& dataFile, const char* type ) { M_type = M_mapType[dataFile( type, "Essential" )]; }

    void readMode( const GetPot& dataFile, const char* mode ) { M_mode = M_mapMode[dataFile( mode, "Full" )]; }

    void readComV( const GetPot& dataFile, const char* component );

    void readBase( const GetPot& dataFile, const std::string& path, std::pair< std::string, baseList_Type >& base, std::string& baseString );

    bool isBase( const GetPot& dataFile, const char* base, std::string& baseString );

    //@}


    //! @name Common Private Members
    //@{

    std::pair< std::string, baseList_Type >                        M_base;
    std::string                                                    M_baseString;

    std::pair< std::string, baseList_Type >                        M_baseRobinAlpha;
    std::string                                                    M_baseStringRobinAlpha;

    std::pair< std::string, baseList_Type >                        M_baseRobinBeta;
    std::string                                                    M_baseStringRobinBeta;

    std::pair< std::string, baseList_Type >                        M_baseDirectional;
    std::string                                                    M_baseStringDirectional;

    std::map< std::string, baseList_Type >                         M_mapBase;

    //@}


    //! @name 1D Private Members
    //@{

    OneDimensional::bcSide_Type                                    M_side;
    OneDimensional::bcLine_Type                                    M_line;
    OneDimensional::bcType_Type                                    M_quantity;

    resistanceContainer_Type                                       M_resistance;
    Real                                                           M_capacitance;

    // Maps
    std::map< std::string, OneDimensional::bcSide_Type >           M_mapSide;
    std::map< std::string, OneDimensional::bcType_Type >           M_mapQuantity;
    std::map< std::string, OneDimensional::bcLine_Type >           M_mapLine;

    //@}


    //! @name 3D Private Members
    //@{

    bcName_Type                                                    M_name;
    bcFlag_Type                                                    M_flag;
    bcType_Type                                                    M_type;
    bcMode_Type                                                    M_mode;
    bcComponentsVec_Type                                           M_comV;

    // Maps
    std::map< std::string, bcType_Type >                           M_mapType;
    std::map< std::string, bcMode_Type >                           M_mapMode;

    //@}
};

} // Namespace LifeV

#endif /* BCInterfaceData_H */
