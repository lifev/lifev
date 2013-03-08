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

#include <lifev/bc_interface/fem/BCInterfaceDefinitions.hpp>

namespace LifeV
{

//! BCInterfaceData - The BCInterface data container
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterfaceData class provides a general interface for the data container in order to pass information
 *  to all the BCInterface classes.
 */
class BCInterfaceData
{
public:

    //! @name Type definitions
    //@{

    typedef std::vector< Real >                                                 parametersContainer_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceData();

    //! Copy constructor
    /*!
     * @param data BCInterfaceData
     */
    BCInterfaceData ( const BCInterfaceData& data );

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
    BCInterfaceData& operator= ( const BCInterfaceData& data );

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

    //! Set the base string of the boundary condition
    /*!
     * @param baseString Boundary condition base string
     */
    void setBaseString ( const std::string& baseString );

    //! Set the base type of the boundary condition
    /*!
     * @param base Boundary condition base type
     */
    void setBase ( const std::pair< std::string, baseList_Type >& base )
    {
        M_base = base;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Get the base string of the boundary condition
    /*!
     * @return Boundary condition base string
     */
    const std::string& baseString() const
    {
        return M_baseString;
    }

    //! Get the base type of the boundary condition
    /*!
     * @return Boundary condition base
     */
    const std::pair< std::string, baseList_Type >& base() const
    {
        return M_base;
    }

    //! Get the base map of the boundary condition
    /*!
     * @return Boundary condition base map
     */
    const std::map< std::string, baseList_Type >& mapBase() const
    {
        return M_mapBase;
    }

    //! Get the parameters vector {A, B, C, ...}
    /*!
     * @return Boundary condition parameters vector
     */
    const parametersContainer_Type& parameters() const
    {
        return M_parameters;
    }

    //@}

protected:

    //! @name Private Methods
    //@{

    void readBase ( const GetPot& dataFile, const std::string& path, std::pair< std::string, baseList_Type >& base, std::string& baseString );

    bool isBase ( const GetPot& dataFile, const char* base, std::string& baseString );

    void readParameters ( const GetPot& dataFile, const char* parameters );

    //@}


    //! @name Common Private Members
    //@{

    std::pair< std::string, baseList_Type >                        M_base;
    std::string                                                    M_baseString;

    std::map< std::string, baseList_Type >                         M_mapBase;

    parametersContainer_Type                                       M_parameters;

    //@}
};

} // Namespace LifeV

#endif /* BCInterfaceData_H */
