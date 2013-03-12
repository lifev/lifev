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
 *  @brief File containing the BCInterfaceData1D class
 *
 *  @date 17-07-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceData1D_H
#define BCInterfaceData1D_H 1

// 1D BCHandler include
#include <lifev/one_d_fsi/fem/OneDFSIBCHandler.hpp>

// BCInterface includes
#include <lifev/bc_interface/fem/BCInterfaceData.hpp>

namespace LifeV
{

//! BCInterfaceData1D - The BCInterface1D data container
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterfaceData1D class provides a general container for all the data
 *  required by the 1D boundary conditions.
 */
class BCInterfaceData1D: public virtual BCInterfaceData
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterfaceData                                                     dataContainer_Type;
    typedef std::vector< Real >                                                 resistanceContainer_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceData1D();

    //! Copy constructor
    /*!
     * @param data BCInterfaceData1D
     */
    BCInterfaceData1D ( const BCInterfaceData1D& data );

    //! Destructor
    virtual ~BCInterfaceData1D() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param data BCInterfaceData1D
     * @return reference to a copy of the class
     */
    BCInterfaceData1D& operator= ( const BCInterfaceData1D& data );

    //@}


    //! @name Methods
    //@{

    //! Read parameters for all kind of BC
    /*!
     * @param fileName Name of the data file.
     * @param dataSection BC section
     * @param name name of the boundary condition
     */
    void readBC ( const std::string& fileName, const std::string& dataSection, const std::string& name );

    //! Display general information about the content of the class
    /*!
     * @param output specify the output format (std::cout by default)
     */
    void showMe ( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the side of the boundary condition
    /*!
     * @param flag Boundary condition side
     */
    void setSide ( const OneDFSI::bcSide_Type& side )
    {
        M_side = side;
    }

    //! Set the line of the boundary condition
    /*!
     * @param line Boundary condition line
     */
    void setLine ( const OneDFSI::bcLine_Type& line )
    {
        M_line = line;
    }

    //! Set the quantity of the boundary condition
    /*!
     * @param quantity Boundary condition quantity
     */
    void setQuantity ( const OneDFSI::bcType_Type& quantity )
    {
        M_quantity = quantity;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Get the flag of the boundary condition
    /*!
     * @return Boundary condition side
     */
    const OneDFSI::bcSide_Type& side() const
    {
        return M_side;
    }

    //! Get the mode of the boundary condition
    /*!
     * @return Boundary condition line
     */
    const OneDFSI::bcLine_Type& line() const
    {
        return M_line;
    }

    //! Get the quantity of the boundary condition
    /*!
     * @return Boundary condition quantity
     */
    const OneDFSI::bcType_Type& quantity() const
    {
        return M_quantity;
    }

    //! Get the resistance vector {R1, R2, R3 ...}
    /*!
     * @return Boundary condition resistance vector
     */
    const resistanceContainer_Type& resistance() const
    {
        return M_resistance;
    }

    //! Get the capacitance
    /*!
     * @return Boundary condition capacitance
     */
    const Real& capacitance() const
    {
        return M_capacitance;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    void readSide ( const GetPot& dataFile, const char* side )
    {
        M_side = M_mapSide[dataFile ( side, "left" )];
    }

    void readLine ( const GetPot& dataFile, const char* line )
    {
        M_line = M_mapLine[dataFile ( line, "first" )];
    }

    void readQuantity ( const GetPot& dataFile, const char* quantity )
    {
        M_quantity = M_mapQuantity[dataFile ( quantity, "A" )];
    }

    void readResistance ( const GetPot& dataFile, const char* resistance );

    void readCapacitance ( const GetPot& dataFile, const char* capacitance )
    {
        M_capacitance = dataFile ( capacitance, 0 );
    }

    //@}


    //! @name Private Members
    //@{

    OneDFSI::bcSide_Type                                           M_side;
    OneDFSI::bcLine_Type                                           M_line;
    OneDFSI::bcType_Type                                           M_quantity;

    resistanceContainer_Type                                       M_resistance;
    Real                                                           M_capacitance;

    // Maps
    std::map< std::string, OneDFSI::bcSide_Type >                  M_mapSide;
    std::map< std::string, OneDFSI::bcType_Type >                  M_mapQuantity;
    std::map< std::string, OneDFSI::bcLine_Type >                  M_mapLine;

    //@}
};

} // Namespace LifeV

#endif /* BCInterfaceData1D_H */
