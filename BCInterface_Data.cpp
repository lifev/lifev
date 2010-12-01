//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief BCInterface_Data
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 17-07-2009
 */

#include <lifemc/lifesolver/BCInterface_Data.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterface_Data::BCInterface_Data() :
        M_name               (),
        M_flag               (),
        M_type               (),
        M_mode               (),
        M_comV               (),
        M_direction          (),
        M_baseString         (),
        M_base               (),
        M_mapType            (),
        M_mapMode            (),
        M_mapBase            ()
{
    //Set mapType
    M_mapType["Essential"]  = Essential;
    M_mapType["Natural"]    = Natural;
    M_mapType["Mixte"]      = Mixte;
    M_mapType["Flux"]       = Flux;
    M_mapType["Resistance"] = Resistance;

    //Set mapMode
    M_mapMode["Scalar"]      = Scalar;
    M_mapMode["Full"]        = Full;
    M_mapMode["Component"]   = Component;
    M_mapMode["Normal"]      = Normal;
    M_mapMode["Tangential"]  = Tangential;
    M_mapMode["Directional"] = Directional;

    //Set mapBase
    M_mapBase["function"]         = BCInterface_function;
    M_mapBase["functionFile"]     = BCInterface_functionFile;
    M_mapBase["OPERfunction"]     = BCInterface_OPERfunction;
    M_mapBase["OPERfunctionFile"] = BCInterface_OPERfunctionFile;
    M_mapBase["FSI"]              = BCInterface_OPERFSI;
}

BCInterface_Data::BCInterface_Data( const BCInterface_Data& data ) :
        M_name              ( data.M_name ),
        M_flag              ( data.M_flag ),
        M_type              ( data.M_type ),
        M_mode              ( data.M_mode ),
        M_comV              ( data.M_comV ),
        M_direction         ( data.M_direction ),
        M_baseString        ( data.M_baseString ),
        M_base              ( data.M_base ),
        M_mapType           ( data.M_mapType ),
        M_mapMode           ( data.M_mapMode ),
        M_mapBase           ( data.M_mapBase )
{
}

// ===================================================
// Operators
// ===================================================
BCInterface_Data&
BCInterface_Data::operator=( const BCInterface_Data& data )
{
    if ( this != &data )
    {
        M_name              = data.M_name;
        M_flag              = data.M_flag;
        M_type              = data.M_type;
        M_mode              = data.M_mode;
        M_comV              = data.M_comV;
        M_direction         = data.M_direction;
        M_baseString        = data.M_baseString;
        M_base              = data.M_base;
        M_mapType           = data.M_mapType;
        M_mapMode           = data.M_mapMode;
        M_mapBase           = data.M_mapBase;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void BCInterface_Data::ReadBC( const std::string& FileName,
                               const std::string& dataSection,
                               const BCName& name )
{
    M_name = name;

    ReadFlag( FileName, ( dataSection + name + "/flag" ).c_str() );
    ReadType( FileName, ( dataSection + name + "/type" ).c_str() );
    ReadMode( FileName, ( dataSection + name + "/mode" ).c_str() );
    ReadComV( FileName, ( dataSection + name + "/component" ).c_str() );
    ReadDirection( FileName, ( dataSection + name + "/direction" ).c_str() );
    ReadBase( FileName, dataSection + name + "/" );
}

// ===================================================
// Methods
// ===================================================
void BCInterface_Data::SetName( const BCName& name )
{
    M_name = name;
}

void BCInterface_Data::SetFlag( const BCFlag& flag )
{
    M_flag = flag;
}

void BCInterface_Data::SetType( const BCType& type )
{
    M_type = type;
}

void BCInterface_Data::SetMode( const BCMode& mode )
{
    M_mode = mode;
}

void BCInterface_Data::SetComV( const BCComV& comV )
{
    M_comV = comV;
}

void BCInterface_Data::AddComV( const UInt& comV )
{
    M_comV.push_back( comV );
}

void BCInterface_Data::SetComV( const UInt& comV, const UInt& index )
{
    M_comV[index] = comV;
}

void BCInterface_Data::SetDirection( const std::string& direction )
{
    M_direction = direction;
}

void BCInterface_Data::SetBaseString( const std::string& baseString )
{
    M_baseString = baseString;
    boost::replace_all( M_baseString, " ", "" );
}

void BCInterface_Data::SetBase( const std::pair< std::string, BCBaseList_Type >& base )
{
    M_base = base;
}

// ===================================================
// Get Methods
// ===================================================
const BCName&
BCInterface_Data::GetName() const
{
    return M_name;
}

const BCFlag&
BCInterface_Data::GetFlag() const
{
    return M_flag;
}

const BCType&
BCInterface_Data::GetType() const
{
    return M_type;
}

const BCMode&
BCInterface_Data::GetMode() const
{
    return M_mode;
}

const BCComV&
BCInterface_Data::GetComV() const
{
    return M_comV;
}

const ID&
BCInterface_Data::GetComN() const
{
    return M_comV.front();
}

const std::string&
BCInterface_Data::GetDirection() const
{
    return M_direction;
}

const std::string&
BCInterface_Data::GetBaseString() const
{
    return M_baseString;
}

const std::pair< std::string, BCInterface_Data::BCBaseList_Type >&
BCInterface_Data::GetBase() const
{
    return M_base;
}

// ===================================================
// Private Methods
// ===================================================
inline void BCInterface_Data::ReadFlag( const std::string& FileName, const char* flag )
{
    GetPot DataFile( FileName );
    M_flag = DataFile( flag, 0 );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface_Data::ReadFlag                        flag: " << static_cast<Real>( M_flag ) << "\n";
#endif
}

inline void BCInterface_Data::ReadType( const std::string& FileName, const char* type )
{
    GetPot DataFile( FileName );
    M_type = M_mapType[DataFile( type, "Essential" )];

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface_Data::ReadType                        type: " << M_type << " (" << DataFile(type, "Essential") << ")\n";
#endif
}

inline void BCInterface_Data::ReadMode( const std::string& FileName, const char* mode )
{
    GetPot DataFile( FileName );
    M_mode = M_mapMode[DataFile( mode, "Full" )];

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface_Data::ReadMode                        mode: " << M_mode << " (" << DataFile(mode, "Full") << ")\n";
#endif
}

inline void BCInterface_Data::ReadComV( const std::string& FileName, const char* component )
{
    GetPot DataFile( FileName );
    UInt componentSize = DataFile.vector_variable_size( component );

    M_comV.clear();
    M_comV.reserve( componentSize );

    for ( UInt j( 0 ); j < componentSize; ++j )
        M_comV.push_back( DataFile( component, 0, j ) );

#ifdef HAVE_LIFEV_DEBUG
    std::stringstream output;
    output << "BCInterface_Data::ReadComV                        comV: ";
    for ( UInt i(0); i < static_cast<UInt>( M_comV.size() ); ++i )
        output << M_comV[i] << " ";
    Debug( 5020 ) << output.str() << "\n";
#endif

}

inline void BCInterface_Data::ReadDirection( const std::string& FileName, const char* direction )
{
    GetPot DataFile( FileName );
    M_direction = DataFile( direction, " " );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface_Data::ReadDirection              direction: " << M_direction << "\n";
#endif
}

inline void BCInterface_Data::ReadBase( const std::string& FileName, const std::string& base )
{
    for ( std::map< std::string, BCBaseList_Type >::iterator j = M_mapBase.begin(); j
            != M_mapBase.end(); ++j )
        if ( IsBase( FileName, ( base + j->first ).c_str() ) )
        {
            M_base.first = j->first;
            M_base.second = M_mapBase[j->first];

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5020 ) << "BCInterface_Data::ReadBase                        base: " << M_base.second << " (" << j->first << ")\n";
            Debug( 5020 ) << "                                           baseString: " << M_baseString << "\n";
#endif

            break;
        }
}

inline bool BCInterface_Data::IsBase( const std::string& FileName, const char* base )
{
    GetPot DataFile( FileName );
    M_baseString = DataFile( base, " " );

    return DataFile.checkVariable( base );
}

} // Namespace LifeV
