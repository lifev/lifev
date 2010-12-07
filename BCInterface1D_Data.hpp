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
 *  @brief File containing the BCInterface1D_Data class
 *
 *  @date 10-05-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface_Data1D_H
#define BCInterface_Data1D_H 1

#include <lifemc/lifesolver/BCInterface1D_Definitions.hpp>

namespace LifeV
{

//! BCInterface1D_Data - The BCInterface data container
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterface1D_Data class provides a general container to pass information
 *  to all the BCInterface functions.
 */
class BCInterface1D_Data
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface1D_BaseList                                                     bcBaseList_Type;
    typedef std::vector< Real >                                                        resistanceContainer_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface1D_Data();

    //! Copy constructor
    /*!
     * @param data BCInterface1D_Data
     */
    BCInterface1D_Data( const BCInterface1D_Data& data );

    //! Destructor
    virtual ~BCInterface1D_Data() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param data BCInterface1D_Data
     * @return reference to a copy of the class
     */
    BCInterface1D_Data& operator=( const BCInterface1D_Data& data );

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

    //! Set the side of the boundary condition
    /*!
     * @param flag Boundary condition side
     */
    void setSide( const OneD_BCSide& side ) { M_side = side; }

    //! Set the line of the boundary condition
    /*!
     * @param line Boundary condition line
     */
    void setLine( const OneD_BCLine& line ) { M_line = line; }

    //! Set the quantity of the boundary condition
    /*!
     * @param quantity Boundary condition quantity
     */
    void setQuantity( const OneD_BC& quantity ) { M_quantity = quantity; }

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

    //! Get the flag of the boundary condition
    const OneD_BCSide& side() const { return M_side; }

    //! Get the mode of the boundary condition
    const OneD_BCLine& line() const { return M_line; }

    //! Get the quantity of the boundary condition
    const OneD_BC& quantity() const { return M_quantity; }

    //! Get the base string of the boundary condition
    const std::string& baseString() const { return M_baseString; }

    //! Get the base type of the boundary condition
    const std::pair< std::string, bcBaseList_Type >& base() const { return M_base; }

    //! Get the resistance vector {R1, R2, R3 ...}
    const resistanceContainer_Type& resistance() const { return M_resistance; }

    //@}

private:

    //! @name Private Methods
    //@{

    void readSide( const GetPot& dataFile, const char* side ) {  M_side = M_mapSide[dataFile( side, "left" )]; }

    void readLine( const GetPot& dataFile, const char* line ) { M_line = M_mapLine[dataFile( line, "first" )]; }

    void readQuantity( const GetPot& dataFile, const char* quantity ) { M_quantity = M_mapQuantity[dataFile( quantity, "A" )]; }

    void readBase( const GetPot& dataFile, const std::string& base );

    bool isBase( const GetPot& dataFile, const char* base );

    void readResistance( const GetPot& dataFile, const char* resistance );

    //@}

    OneD_BCSide                               M_side;
    OneD_BCLine                               M_line;
    OneD_BC                                   M_quantity;
    std::string                               M_baseString;
    std::pair< std::string, bcBaseList_Type > M_base;

    resistanceContainer_Type                  M_resistance;

    // Maps
    std::map< std::string, OneD_BCSide >      M_mapSide;
    std::map< std::string, OneD_BC >          M_mapQuantity;
    std::map< std::string, OneD_BCLine >      M_mapLine;
    std::map< std::string, bcBaseList_Type >  M_mapBase;
};

} // Namespace LifeV

#endif /* BCInterface_Data1D_H */
