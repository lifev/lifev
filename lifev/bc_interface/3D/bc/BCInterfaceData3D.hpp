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
 *  @brief File containing the BCInterfaceData3D class
 *
 *  @date 17-07-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceData3D_H
#define BCInterfaceData3D_H 1

// 3D BCHandler include
#include <lifev/core/fem/BCHandler.hpp>

// Data interpolator include
#include <lifev/core/fem/BCDataInterpolator.hpp>

// BCInterface includes
#include <lifev/bc_interface/core/bc/BCInterfaceData.hpp>

namespace LifeV
{

//! BCInterfaceData3D - The BCInterface3D data container
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterfaceData3D class provides a general container for all the data
 *  required by the 3D boundary conditions.
 */
class BCInterfaceData3D: public virtual BCInterfaceData
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterfaceData                                                     dataContainer_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceData3D();

    //! Copy constructor
    /*!
     * @param data BCInterfaceData3D
     */
    BCInterfaceData3D ( const BCInterfaceData3D& data );

    //! Destructor
    virtual ~BCInterfaceData3D() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param data BCInterfaceData3D
     * @return reference to a copy of the class
     */
    BCInterfaceData3D& operator= ( const BCInterfaceData3D& data );

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

    //! Set the directional base as the current base
    void setDirectionalBase()
    {
        M_base = M_baseDirectional;
        M_baseString = M_baseStringDirectional;
    }

    //! Set the Robin Alpha base as the current base
    void setRobinBaseAlpha()
    {
        M_base = M_baseRobinAlpha;
        M_baseString = M_baseStringRobinAlpha;
    }

    //! Set the Robin Beta base as the current base
    void setRobinBaseBeta()
    {
        M_base = M_baseRobinBeta;
        M_baseString = M_baseStringRobinBeta;
    }

    //! Display general information about the content of the class
    /*!
     * @param output specify the output format (std::cout by default)
     */
    void showMe ( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the name of the boundary condition
    /*!
     * @param name Boundary condition name
     */
    void setName ( const bcName_Type& name )
    {
        M_name = name;
    }

    //! Set the flag of the boundary condition
    /*!
     * @param flag Boundary condition flag
     */
    void setFlag ( const bcFlag_Type& flag )
    {
        M_boundaryID = flag;
    }

    //! Set the type of the boundary condition
    /*!
     * @param type Boundary condition type
     */
    void setType ( const bcType_Type& type )
    {
        M_type = type;
    }

    //! Set the mode of the boundary condition
    /*!
     * @param mode Boundary condition mode
     */
    void setMode ( const bcMode_Type& mode )
    {
        M_mode = mode;
    }

    //! Set the components vector of the boundary condition
    /*!
     * @param componentsVector Boundary condition components vector
     */
    void setComponentsVector ( const bcComponentsVec_Type& componentsVector )
    {
        M_componentsVector = componentsVector;
    }

    //! Set the i-component of the components vector of the boundary condition
    /*!
     * @param componentsVector Boundary condition component
     * @param index Index value
     */
    void setComponentsVector ( const UInt& componentsVector, const UInt& index )
    {
        M_componentsVector[index] = componentsVector;
    }

    //! Add a component to the component vector of the boundary condition
    /*!
     * @param componentsVector Boundary condition component
     */
    void addComponentsVector ( const UInt& componentsVector )
    {
        M_componentsVector.push_back ( componentsVector );
    }

    //@}


    //! @name Get Methods
    //@{

    //! Get the name of the boundary condition
    /*!
     * @return Boundary condition name
     */
    const bcName_Type& name() const
    {
        return M_name;
    }

    //! Get the flag of the boundary condition
    /*!
     * @return Boundary condition flag
     */
    const bcFlag_Type& flag() const
    {
        return M_boundaryID;
    }

    //! Get the type of the boundary condition
    /*!
     * @return Boundary condition type
     */
    const bcType_Type& type() const
    {
        return M_type;
    }

    //! Get the mode of the boundary condition
    /*!
     * @return Boundary condition mode
     */
    const bcMode_Type& mode() const
    {
        return M_mode;
    }

    //! Get the vector of components of the boundary condition
    /*!
     * @return Boundary condition vector of components
     */
    const bcComponentsVec_Type& componentsVector() const
    {
        return M_componentsVector;
    }

    //! Get the number of components of the boundary condition
    /*!
     * Note that this method should not be called for the case of "Component" boundary conditions,
     * since it does not return the size of the M_componentsVector. It has to be used only for the case
     * of "Full" boundary conditions.
     *
     * @return Number of components of the boundary condition
     */
    const ID& componentsNumber() const
    {
        return M_componentsVector.front();
    }

    //@}

private:

    //! @name Private Methods
    //@{

    void readFlag ( const GetPot& dataFile, const char* flag )
    {
        M_boundaryID = dataFile ( flag, 0 );
    }

    void readType ( const GetPot& dataFile, const char* type )
    {
        M_type = M_mapType[dataFile ( type, "Essential" )];
    }

    void readMode ( const GetPot& dataFile, const char* mode )
    {
        M_mode = M_mapMode[dataFile ( mode, "Full" )];
    }

    void readComponentsVector ( const GetPot& dataFile, const char* component );

    //@}


    //! @name Private Members
    //@{

    std::pair< std::string, baseList_Type >                        M_baseRobinAlpha;
    std::string                                                    M_baseStringRobinAlpha;

    std::pair< std::string, baseList_Type >                        M_baseRobinBeta;
    std::string                                                    M_baseStringRobinBeta;

    std::pair< std::string, baseList_Type >                        M_baseDirectional;
    std::string                                                    M_baseStringDirectional;

    bcName_Type                                                    M_name;
    bcType_Type                                                    M_type;
    bcMode_Type                                                    M_mode;
    bcComponentsVec_Type                                           M_componentsVector;

    // Maps
    std::map< std::string, bcType_Type >                           M_mapType;
    std::map< std::string, bcMode_Type >                           M_mapMode;

    //@}
};

} // Namespace LifeV

#endif /* BCInterfaceData3D_H */
