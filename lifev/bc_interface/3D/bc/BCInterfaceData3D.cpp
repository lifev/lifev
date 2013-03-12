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

#include <lifev/bc_interface/3D/bc/BCInterfaceData3D.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterfaceData3D::BCInterfaceData3D() :
    BCInterfaceData         (),
    M_baseRobinAlpha        (),
    M_baseStringRobinAlpha  (),
    M_baseRobinBeta         (),
    M_baseStringRobinBeta   (),
    M_baseDirectional       (),
    M_baseStringDirectional (),
    M_name                  (),
    M_flag                  (),
    M_type                  (),
    M_mode                  (),
    M_componentsVector      (),
    M_mapType               (),
    M_mapMode               ()
{
    //Set mapType
    M_mapType["Essential"]         = Essential;
    M_mapType["EssentialEdges"]    = EssentialEdges;
    M_mapType["EssentialVertices"] = EssentialVertices;
    M_mapType["Natural"]           = Natural;
    M_mapType["Robin"]             = Robin;
    M_mapType["Flux"]              = Flux;
    M_mapType["Resistance"]        = Resistance;

    //Set mapMode
    M_mapMode["Scalar"]      = Scalar;
    M_mapMode["Full"]        = Full;
    M_mapMode["Component"]   = Component;
    M_mapMode["Normal"]      = Normal;
    M_mapMode["Tangential"]  = Tangential;
    M_mapMode["Directional"] = Directional;
}

BCInterfaceData3D::BCInterfaceData3D ( const BCInterfaceData3D& data ) :
    BCInterfaceData         ( data ),
    M_baseRobinAlpha        ( data.M_baseRobinAlpha ),
    M_baseStringRobinAlpha  ( data.M_baseStringRobinAlpha ),
    M_baseRobinBeta         ( data.M_baseRobinBeta ),
    M_baseStringRobinBeta   ( data.M_baseStringRobinBeta ),
    M_baseDirectional       ( data.M_baseDirectional ),
    M_baseStringDirectional ( data.M_baseStringDirectional ),
    M_name                  ( data.M_name ),
    M_flag                  ( data.M_flag ),
    M_type                  ( data.M_type ),
    M_mode                  ( data.M_mode ),
    M_componentsVector      ( data.M_componentsVector ),
    M_mapType               ( data.M_mapType ),
    M_mapMode               ( data.M_mapMode )
{
}

// ===================================================
// Operators
// ===================================================
BCInterfaceData3D&
BCInterfaceData3D::operator= ( const BCInterfaceData3D& data )
{
    if ( this != &data )
    {
        BCInterfaceData::operator= ( data );
        M_baseRobinAlpha        = data.M_baseRobinAlpha;
        M_baseStringRobinAlpha  = data.M_baseStringRobinAlpha;
        M_baseRobinBeta         = data.M_baseRobinBeta;
        M_baseStringRobinBeta   = data.M_baseStringRobinBeta;
        M_baseDirectional       = data.M_baseDirectional;
        M_baseStringDirectional = data.M_baseStringDirectional;
        M_name                  = data.M_name;
        M_flag                  = data.M_flag;
        M_type                  = data.M_type;
        M_mode                  = data.M_mode;
        M_componentsVector      = data.M_componentsVector;
        M_mapType               = data.M_mapType;
        M_mapMode               = data.M_mapMode;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
BCInterfaceData3D::readBC ( const std::string& fileName, const std::string& dataSection, const std::string& name )
{
    // Call to the base class
    dataContainer_Type::readBC ( fileName, dataSection, name );

    // Read 3D data
    GetPot dataFile ( fileName );

    M_name = name;

    readFlag ( dataFile, ( dataSection + name + "/flag" ).c_str() );
    readType ( dataFile, ( dataSection + name + "/type" ).c_str() );
    readMode ( dataFile, ( dataSection + name + "/mode" ).c_str() );
    readComponentsVector ( dataFile, ( dataSection + name + "/component" ).c_str() );

    // Read base
    if ( M_type == Robin )
    {
        readBase ( dataFile, dataSection + name + "/RobinAlpha/", M_baseRobinAlpha, M_baseStringRobinAlpha );
        readBase ( dataFile, dataSection + name + "/RobinBeta/", M_baseRobinBeta, M_baseStringRobinBeta );
    }
    if ( M_mode == Directional )
    {
        readBase ( dataFile, dataSection + name + "/Directional/", M_baseDirectional, M_baseStringDirectional );
    }
}

void
BCInterfaceData3D::showMe ( std::ostream& output ) const
{
    // Call to the base class
    dataContainer_Type::showMe ( output );

    // Show 3D data
    output << "Flag              = " << static_cast< Real > ( M_flag ) << std::endl;
    output << "Type              = " << M_type << std::endl;
    output << "Mode              = " << M_mode << std::endl;
    output << "Components Vector = ";
    for ( UInt i (0); i < static_cast<UInt> ( M_componentsVector.size() ); ++i )
    {
        output << M_componentsVector[i] << " ";
    }
    output << "\n";
}

// ===================================================
// Private Methods
// ===================================================
void
BCInterfaceData3D::readComponentsVector ( const GetPot& dataFile, const char* component )
{
    UInt componentSize = dataFile.vector_variable_size ( component );

    M_componentsVector.resize ( componentSize );
    for ( UInt j ( 0 ); j < componentSize; ++j )
    {
        M_componentsVector[j] = dataFile ( component, 0, j );
    }
}

} // Namespace LifeV
