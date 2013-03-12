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
 *  @brief File containing the BCInterfaceData0D class
 *
 *  @date 17-07-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceData0D_H
#define BCInterfaceData0D_H 1

// 0D BCHandler include
#include <lifev/zero_dimensional/fem/ZeroDimensionalBCHandler.hpp>

// BCInterface includes
#include <lifev/bc_interface/core/bc/BCInterfaceData.hpp>

namespace LifeV
{

//! BCInterfaceData0D - The BCInterface1D data container
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterfaceData0D class provides a general container for all the data
 *  required by the 0D boundary conditions.
 */
class BCInterfaceData0D: public virtual BCInterfaceData
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterfaceData                                                     dataContainer_Type;
    typedef ZeroDimensionalBCHandler::bcType_Type                               bcType_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceData0D();

    //! Copy constructor
    /*!
     * @param data BCInterfaceData0D
     */
    BCInterfaceData0D ( const BCInterfaceData0D& data );

    //! Destructor
    virtual ~BCInterfaceData0D() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param data BCInterfaceData0D
     * @return reference to a copy of the class
     */
    BCInterfaceData0D& operator= ( const BCInterfaceData0D& data );

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

    //! Set the flag of the boundary condition
    /*!
     * @param flag Boundary condition flag
     */
    void setFlag ( const bcFlag_Type& flag )
    {
        M_flag = flag;
    }

    //! Set the type of the boundary condition
    /*!
     * @param type Boundary condition type
     */
    void setType ( const bcType_Type& type )
    {
        M_type = type;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Get the flag of the boundary condition
    /*!
     * @return Boundary condition flag
     */
    const bcFlag_Type& flag() const
    {
        return M_flag;
    }

    //! Get the type of the boundary condition
    /*!
     * @return Boundary condition type
     */
    const bcType_Type& type() const
    {
        return M_type;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    void readFlag ( const GetPot& dataFile, const char* flag )
    {
        M_flag = dataFile ( flag, 0 );
    }

    void readType ( const GetPot& dataFile, const char* type )
    {
        M_type = M_mapType[dataFile ( type, "Current" )];
    }

    //@}


    //! @name Private Members
    //@{


    bcFlag_Type                                                    M_flag;
    bcType_Type                                                    M_type;


    // Maps
    std::map< std::string, bcType_Type >                           M_mapType;

    //@}
};

} // Namespace LifeV

#endif /* BCInterfaceData0D_H */
