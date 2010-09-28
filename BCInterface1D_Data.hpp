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

#ifndef BCInterface_Data1D_H
#define BCInterface_Data1D_H 1

#include <lifemc/lifesolver/BCInterface1D_Definitions.hpp>

namespace LifeV {

//! BCInterface1D_Data - The BCInterface data container
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterface1D_Data class provides a general container to pass information
 *  to all the BCInterface functions.
 */
template< class Operator >
class BCInterface1D_Data
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface1D_BaseList                                                     BCBaseList_Type;
    typedef std::vector< Real >                                                        ResistanceContainer_Type;

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
    ~BCInterface1D_Data() {}

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
     * @param FileName Name of the data file.
     * @param dataSection BC section
     * @param name name of the boundary condition
     */
    void ReadBC( const std::string& FileName, const std::string& dataSection, const BCName& name );

    //@}


    //! @name Set Methods
    //@{

    //! Set the operator
    /*!
     * @param Oper Operator
     */
    void SetOperator( const boost::shared_ptr< Operator >& Oper );

    //! Set the side of the boundary condition
    /*!
     * @param flag Boundary condition side
     */
    void SetSide( const OneD_BCSide& side );

    //! Set the line of the boundary condition
    /*!
     * @param line Boundary condition line
     */
    void SetLine( const OneD_BCLine& line );

    //! Set the quantity of the boundary condition
    /*!
     * @param quantity Boundary condition quantity
     */
    void SetQuantity( const OneD_BC& quantity );

    //! Set the base string of the boundary condition
    /*!
     * @param baseString Boundary condition base string
     */
    void SetBaseString( const std::string& baseString );

    //! Set the base type of the boundary condition
    /*!
     * @param base Boundary condition base type
     */
    void SetBase( const std::pair< std::string, BCBaseList_Type >& base );

    //@}


    //! @name Get Methods
    //@{

    //! Get the shared_ptr to the operator
    const boost::shared_ptr< Operator >& GetOperator() const;

    //! Get the flag of the boundary condition
    const OneD_BCSide& GetSide() const;

    //! Get the mode of the boundary condition
    const OneD_BCLine& GetLine() const;

    //! Get the quantity of the boundary condition
    const OneD_BC& GetQuantity() const;

    //! Get the base string of the boundary condition
    const std::string& GetBaseString() const;

    //! Get the base type of the boundary condition
    const std::pair< std::string, BCBaseList_Type >& GetBase() const;

    //! Get the resistance vector {R1, R2, R3 ...}
    const ResistanceContainer_Type& GetResistance() const;

    //@}

private:

    //! @name Private Methods
    //@{

    inline void ReadSide( const std::string& FileName, const char* side );

    inline void ReadLine( const std::string& FileName, const char* line );

    inline void ReadQuantity( const std::string& FileName, const char* quantity );

    inline void ReadBase( const std::string& FileName, const std::string& base );

    inline bool IsBase( const std::string& FileName, const char* base );

    inline void ReadResistance( const std::string& FileName, const char* resistance );

    //@}

    boost::shared_ptr< Operator >             M_operator;

    OneD_BCSide                               M_side;
    OneD_BCLine                               M_line;
    OneD_BC                                   M_quantity;
    std::string                               M_baseString;
    std::pair< std::string, BCBaseList_Type > M_base;

    ResistanceContainer_Type                  M_resistance;

    // Maps
    std::map< std::string, OneD_BCSide >      M_mapSide;
    std::map< std::string, OneD_BC >          M_mapQuantity;
    std::map< std::string, OneD_BCLine >      M_mapLine;
    std::map< std::string, BCBaseList_Type >  M_mapBase;
};

// ===================================================
// Constructors
// ===================================================
template< class Operator >
BCInterface1D_Data< Operator >::BCInterface1D_Data() :
    M_operator           (),
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

template< class Operator >
BCInterface1D_Data< Operator >::BCInterface1D_Data( const BCInterface1D_Data& data ) :
    M_operator          ( data.M_operator ),
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
template< class Operator >
BCInterface1D_Data< Operator >&
BCInterface1D_Data< Operator >::operator=( const BCInterface1D_Data& data )
{
    if ( this != &data )
    {
        M_operator          = data.M_operator;
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
template< class Operator >
inline void BCInterface1D_Data< Operator >::ReadBC( const std::string& FileName,
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
template< class Operator >
void BCInterface1D_Data< Operator >::SetOperator( const boost::shared_ptr< Operator >& Oper )
{
    M_operator = Oper;
}

template< class Operator >
void BCInterface1D_Data< Operator >::SetSide( const OneD_BCSide& side )
{
    M_side = side;
}

template< class Operator >
void BCInterface1D_Data< Operator >::SetQuantity( const OneD_BC& quantity )
{
    M_quantity = quantity;
}

template< class Operator >
void BCInterface1D_Data< Operator >::SetLine( const OneD_BCLine& line )
{
    M_line = line;
}

template< class Operator >
void BCInterface1D_Data< Operator >::SetBaseString( const std::string& baseString )
{
    M_baseString = baseString;
    boost::replace_all( M_baseString, " ", "" );
}

template< class Operator >
void BCInterface1D_Data< Operator >::SetBase( const std::pair< std::string, BCBaseList_Type >& base )
{
    M_base = base;
}

// ===================================================
// Get Methods
// ===================================================
template< class Operator >
const boost::shared_ptr< Operator >&
BCInterface1D_Data< Operator >::GetOperator() const
{
    return M_operator;
}

template< class Operator >
const OneD_BCSide&
BCInterface1D_Data< Operator >::GetSide() const
{
    return M_side;
}

template< class Operator >
const OneD_BC&
BCInterface1D_Data< Operator >::GetQuantity() const
{
    return M_quantity;
}

template< class Operator >
const OneD_BCLine&
BCInterface1D_Data< Operator >::GetLine() const
{
    return M_line;
}

template< class Operator >
const std::string&
BCInterface1D_Data< Operator >::GetBaseString() const
{
    return M_baseString;
}

template< class Operator >
const std::pair< std::string, typename BCInterface1D_Data< Operator >::BCBaseList_Type >&
BCInterface1D_Data< Operator >::GetBase() const
{
    return M_base;
}

template< class Operator >
const typename BCInterface1D_Data< Operator >::ResistanceContainer_Type&
BCInterface1D_Data< Operator >::GetResistance() const
{
    return M_resistance;
}

// ===================================================
// Private Methods
// ===================================================
template< class Operator >
inline void BCInterface1D_Data< Operator >::ReadSide( const std::string& FileName, const char* side )
{
    GetPot DataFile( FileName );
    M_side = M_mapSide[DataFile( side, "left" )];

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface1D_Data::ReadSide                        side: " << static_cast<Real>( M_side ) << "\n";
#endif
}

template< class Operator >
inline void BCInterface1D_Data< Operator >::ReadQuantity( const std::string& FileName, const char* quantity )
{
    GetPot DataFile( FileName );
    M_quantity = M_mapQuantity[DataFile( quantity, "A" )];

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface1D_Data::ReadQuantity                quantity: " << M_quantity << " (" << DataFile(quantity, "A") << ")\n";
#endif
}

template< class Operator >
inline void BCInterface1D_Data< Operator >::ReadLine( const std::string& FileName, const char* line )
{
    GetPot DataFile( FileName );
    M_line = M_mapLine[DataFile( line, "first" )];

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface1D_Data::ReadLine                        line: " << M_line << " (" << DataFile(line, "first") << ")\n";
#endif
}

template< class Operator >
inline void BCInterface1D_Data< Operator >::ReadBase( const std::string& FileName, const std::string& base )
{
    for ( typename std::map< std::string, BCBaseList_Type >::iterator j = M_mapBase.begin(); j
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

template< class Operator >
inline bool BCInterface1D_Data< Operator >::IsBase( const std::string& FileName, const char* base )
{
    GetPot DataFile( FileName );
    M_baseString = DataFile( base, " " );

    return DataFile.checkVariable( base );
}

template< class Operator >
inline void BCInterface1D_Data< Operator >::ReadResistance( const std::string& FileName, const char* resistance )
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

#endif /* BCInterface_Data1D_H */
