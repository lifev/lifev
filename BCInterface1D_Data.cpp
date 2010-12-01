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
 *  @date 10-05-2010
 */


#include <lifemc/lifesolver/BCInterface1D_Data.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterface1D_Data::BCInterface1D_Data() :
        M_side               (),
        M_line               (),
        M_quantity           (),
        M_baseString         (),
        M_base               (),
        M_resistance         (),
        M_mapSide            (),
        M_mapQuantity        (),
        M_mapLine            (),
        M_mapBase            ()
{
    //Set mapSide
    M_mapSide["left"]   = OneD_left;
    M_mapSide["right"]  = OneD_right;

    //Set mapQuantity
    M_mapQuantity["A"]      = OneD_A;
    M_mapQuantity["Q"]      = OneD_Q;
    M_mapQuantity["W1"]     = OneD_W1;
    M_mapQuantity["W2"]     = OneD_W2;
    M_mapQuantity["P"]      = OneD_P;

    //Set mapLine
    M_mapLine["first"]  = OneD_first;
    M_mapLine["second"] = OneD_second;

    //Set mapBase
    M_mapBase["function"]         = BCInterface1D_function;
    M_mapBase["functionFile"]     = BCInterface1D_functionFile;
    M_mapBase["OPERfunction"]     = BCInterface1D_OPERfunction;
    M_mapBase["OPERfunctionFile"] = BCInterface1D_OPERfunctionFile;
    M_mapBase["Default"]          = BCInterface1D_Default;
}

BCInterface1D_Data::BCInterface1D_Data( const BCInterface1D_Data& data ) :
        M_side              ( data.M_side ),
        M_line              ( data.M_line ),
        M_quantity          ( data.M_quantity ),
        M_baseString        ( data.M_baseString ),
        M_base              ( data.M_base ),
        M_resistance        ( data.M_resistance ),
        M_mapSide           ( data.M_mapSide ),
        M_mapQuantity       ( data.M_mapQuantity ),
        M_mapLine           ( data.M_mapLine ),
        M_mapBase           ( data.M_mapBase )
{
}

// ===================================================
// Operators
// ===================================================
BCInterface1D_Data&
BCInterface1D_Data::operator=( const BCInterface1D_Data& data )
{
    if ( this != &data )
    {
        M_side              = data.M_side;
        M_line              = data.M_line;
        M_quantity          = data.M_quantity;
        M_baseString        = data.M_baseString;
        M_base              = data.M_base;
        M_resistance        = data.M_resistance;
        M_mapSide           = data.M_mapSide;
        M_mapQuantity       = data.M_mapQuantity;
        M_mapLine           = data.M_mapLine;
        M_mapBase           = data.M_mapBase;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void BCInterface1D_Data::ReadBC( const std::string& FileName,
                                 const std::string& dataSection,
                                 const BCName&      name )
{
    ReadSide( FileName, ( dataSection + name + "/side" ).c_str() );
    ReadQuantity( FileName, ( dataSection + name + "/quantity" ).c_str() );
    ReadLine( FileName, ( dataSection + name + "/line" ).c_str() );
    ReadBase( FileName, dataSection + name + "/" );
    ReadResistance( FileName, ( dataSection + name + "/resistance" ).c_str() );
}

// ===================================================
// Methods
// ===================================================
void BCInterface1D_Data::SetSide( const OneD_BCSide& side )
{
    M_side = side;
}

void BCInterface1D_Data::SetQuantity( const OneD_BC& quantity )
{
    M_quantity = quantity;
}

void BCInterface1D_Data::SetLine( const OneD_BCLine& line )
{
    M_line = line;
}

void BCInterface1D_Data::SetBaseString( const std::string& baseString )
{
    M_baseString = baseString;
    boost::replace_all( M_baseString, " ", "" );
}

void BCInterface1D_Data::SetBase( const std::pair< std::string, BCBaseList_Type >& base )
{
    M_base = base;
}

// ===================================================
// Get Methods
// ===================================================
const OneD_BCSide&
BCInterface1D_Data::GetSide() const
{
    return M_side;
}

const OneD_BC&
BCInterface1D_Data::GetQuantity() const
{
    return M_quantity;
}

const OneD_BCLine&
BCInterface1D_Data::GetLine() const
{
    return M_line;
}

const std::string&
BCInterface1D_Data::GetBaseString() const
{
    return M_baseString;
}

const std::pair< std::string, BCInterface1D_Data::BCBaseList_Type >&
BCInterface1D_Data::GetBase() const
{
    return M_base;
}

const BCInterface1D_Data::ResistanceContainer_Type&
BCInterface1D_Data::GetResistance() const
{
    return M_resistance;
}

// ===================================================
// Private Methods
// ===================================================
inline void BCInterface1D_Data::ReadSide( const std::string& FileName, const char* side )
{
    GetPot DataFile( FileName );
    M_side = M_mapSide[DataFile( side, "left" )];

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface1D_Data::ReadSide                        side: " << static_cast<Real>( M_side ) << "\n";
#endif
}

inline void BCInterface1D_Data::ReadQuantity( const std::string& FileName, const char* quantity )
{
    GetPot DataFile( FileName );
    M_quantity = M_mapQuantity[DataFile( quantity, "A" )];

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface1D_Data::ReadQuantity                quantity: " << M_quantity << " (" << DataFile(quantity, "A") << ")\n";
#endif
}

inline void BCInterface1D_Data::ReadLine( const std::string& FileName, const char* line )
{
    GetPot DataFile( FileName );
    M_line = M_mapLine[DataFile( line, "first" )];

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface1D_Data::ReadLine                        line: " << M_line << " (" << DataFile(line, "first") << ")\n";
#endif
}

inline void BCInterface1D_Data::ReadBase( const std::string& FileName, const std::string& base )
{
    for ( std::map< std::string, BCBaseList_Type >::iterator j = M_mapBase.begin(); j
            != M_mapBase.end(); ++j )
        if ( IsBase( FileName, ( base + j->first ).c_str() ) )
        {
            M_base.first = j->first;
            M_base.second = M_mapBase[j->first];

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5020 ) << "BCInterface1D_Data::ReadBase                     base: " << M_base.second << " (" << j->first << ")\n";
            Debug( 5020 ) << "                                           baseString: " << M_baseString << "\n";
#endif

            break;
        }
}

inline bool BCInterface1D_Data::IsBase( const std::string& FileName, const char* base )
{
    GetPot DataFile( FileName );
    M_baseString = DataFile( base, " " );

    return DataFile.checkVariable( base );
}

inline void BCInterface1D_Data::ReadResistance( const std::string& FileName, const char* resistance )
{
    GetPot DataFile( FileName );
    UInt resistanceSize = DataFile.vector_variable_size( resistance );

    M_resistance.clear();
    M_resistance.reserve( resistanceSize );

    for ( UInt j( 0 ); j < resistanceSize; ++j )
        M_resistance.push_back( DataFile( resistance, 0, j ) );

#ifdef HAVE_LIFEV_DEBUG
    std::stringstream output;
    output << "BCInterface1D_Data::ReadResistance                  resistance: ";
    for ( UInt i(0); i < resistanceSize; ++i )
        output << M_resistance[i] << " ";
    Debug( 5020 ) << output.str() << "\n";
#endif
}

} // Namespace LifeV
