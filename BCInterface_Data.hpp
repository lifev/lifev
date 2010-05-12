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

#ifndef BCInterface_Data_H
#define BCInterface_Data_H 1

#include <lifemc/lifesolver/BCInterface_Definitions.hpp>

namespace LifeV {

//! BCInterface_Data - The BCInterface data container
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterface_Data class provides a general container to pass information
 *  to all the BCInterface functions.
 */
template< class Operator >
class BCInterface_Data
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface_BaseList                                                       BCBaseList_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface_Data();

    //! Copy constructor
    /*!
     * @param data BCInterface_Data
     */
    BCInterface_Data( const BCInterface_Data& data );

    //! Destructor
    ~BCInterface_Data() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param data BCInterface_Data
     * @return reference to a copy of the class
     */
    BCInterface_Data& operator=( const BCInterface_Data& data );

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

    //! Set the name of the boundary condition
    /*!
     * @param name Boundary condition name
     */
    void SetName( const BCName& name );

    //! Set the flag of the boundary condition
    /*!
     * @param flag Boundary condition flag
     */
    void SetFlag( const BCFlag& flag );

    //! Set the type of the boundary condition
    /*!
     * @param type Boundary condition type
     */
    void SetType( const BCType& type );

    //! Set the mode of the boundary condition
    /*!
     * @param mode Boundary condition mode
     */
    void SetMode( const BCMode& mode );

    //! Set the components vector of the boundary condition
    /*!
     * @param comV Boundary condition components vector
     */
    void SetComV( const BCComV& comV );

    //! Add a component to the component vector of the boundary condition
    /*!
     * @param comV Boundary condition component
     */
    void AddComV( const UInt& comV );

    //! Set the i-component of the components vector of the boundary condition
    /*!
     * @param comV Boundary condition component
     * @param index Index value
     */
    void SetComV( const UInt& comV, const UInt& index );

    //! Set the direction string of the boundary condition
    /*!
     * @param direction Boundary condition direction string
     */
    void SetDirection( const std::string& direction );

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

    //! Get the name of the boundary condition
    const BCName& GetName() const;

    //! Get the flag of the boundary condition
    const BCFlag& GetFlag() const;

    //! Get the type of the boundary condition
    const BCType& GetType() const;

    //! Get the mode of the boundary condition
    const BCMode& GetMode() const;

    //! Get the vector of components of the boundary condition
    const BCComV& GetComV() const;

    //! Get the number of components of the boundary condition
    const ID& GetComN() const;

    //! Get the direction string of the boundary condition
    const std::string& GetDirection() const;

    //! Get the base string of the boundary condition
    const std::string& GetBaseString() const;

    //! Get the base type of the boundary condition
    const std::pair< std::string, BCBaseList_Type >& GetBase() const;

    //@}

private:

    //! @name Private Methods
    //@{

    inline void ReadFlag( const std::string& FileName, const char* flag );

    inline void ReadType( const std::string& FileName, const char* type );

    inline void ReadMode( const std::string& FileName, const char* mode );

    inline void ReadComV( const std::string& FileName, const char* component );

    inline void ReadDirection( const std::string& FileName, const char* direction );

    inline void ReadBase( const std::string& FileName, const std::string& base );

    inline bool IsBase( const std::string& FileName, const char* base );

    //@}

    boost::shared_ptr< Operator >             M_operator;

    BCName                                    M_name;
    BCFlag                                    M_flag;
    BCType                                    M_type;
    BCMode                                    M_mode;
    BCComV                                    M_comV;
    std::string                               M_direction;
    std::string                               M_baseString;
    std::pair< std::string, BCBaseList_Type > M_base;

    // Maps
    std::map< std::string, BCType >           M_mapType;
    std::map< std::string, BCMode >           M_mapMode;
    std::map< std::string, BCBaseList_Type >  M_mapBase;
};

// ===================================================
// Constructors
// ===================================================
template< class Operator >
BCInterface_Data< Operator >::BCInterface_Data() :
    M_operator           (),
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

template< class Operator >
BCInterface_Data< Operator >::BCInterface_Data( const BCInterface_Data& data ) :
    M_operator          ( data.M_operator ),
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
template< class Operator >
BCInterface_Data< Operator >&
BCInterface_Data< Operator >::operator=( const BCInterface_Data& data )
{
    if ( this != &data )
    {
        M_operator          = data.M_operator;
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
template< class Operator >
inline void BCInterface_Data< Operator >::ReadBC( const std::string& FileName,
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
template< class Operator >
void BCInterface_Data< Operator >::SetOperator( const boost::shared_ptr< Operator >& Oper )
{
    M_operator = Oper;
}

template< class Operator >
void BCInterface_Data< Operator >::SetName( const BCName& name )
{
    M_name = name;
}

template< class Operator >
void BCInterface_Data< Operator >::SetFlag( const BCFlag& flag )
{
    M_flag = flag;
}

template< class Operator >
void BCInterface_Data< Operator >::SetType( const BCType& type )
{
    M_type = type;
}

template< class Operator >
void BCInterface_Data< Operator >::SetMode( const BCMode& mode )
{
    M_mode = mode;
}

template< class Operator >
void BCInterface_Data< Operator >::SetComV( const BCComV& comV )
{
    M_comV = comV;
}

template< class Operator >
void BCInterface_Data< Operator >::AddComV( const UInt& comV )
{
    M_comV.push_back( comV );
}

template< class Operator >
void BCInterface_Data< Operator >::SetComV( const UInt& comV, const UInt& index )
{
    M_comV[index] = comV;
}

template< class Operator >
void BCInterface_Data< Operator >::SetDirection( const std::string& direction )
{
    M_direction = direction;
}

template< class Operator >
void BCInterface_Data< Operator >::SetBaseString( const std::string& baseString )
{
    M_baseString = baseString;
    boost::replace_all( M_baseString, " ", "" );
}

template< class Operator >
void BCInterface_Data< Operator >::SetBase( const std::pair< std::string, BCBaseList_Type >& base )
{
    M_base = base;
}

// ===================================================
// Get Methods
// ===================================================
template< class Operator >
const boost::shared_ptr< Operator >&
BCInterface_Data< Operator >::GetOperator() const
{
    return M_operator;
}

template< class Operator >
const BCName&
BCInterface_Data< Operator >::GetName() const
{
    return M_name;
}

template< class Operator >
const BCFlag&
BCInterface_Data< Operator >::GetFlag() const
{
    return M_flag;
}

template< class Operator >
const BCType&
BCInterface_Data< Operator >::GetType() const
{
    return M_type;
}

template< class Operator >
const BCMode&
BCInterface_Data< Operator >::GetMode() const
{
    return M_mode;
}

template< class Operator >
const BCComV&
BCInterface_Data< Operator >::GetComV() const
{
    return M_comV;
}

template< class Operator >
const ID&
BCInterface_Data< Operator >::GetComN() const
{
    return M_comV.front();
}

template< class Operator >
const std::string&
BCInterface_Data< Operator >::GetDirection() const
{
    return M_direction;
}

template< class Operator >
const std::string&
BCInterface_Data< Operator >::GetBaseString() const
{
    return M_baseString;
}

template< class Operator >
const std::pair< std::string, typename BCInterface_Data< Operator >::BCBaseList_Type >&
BCInterface_Data< Operator >::GetBase() const
{
    return M_base;
}

// ===================================================
// Private Methods
// ===================================================
template< class Operator >
inline void BCInterface_Data< Operator >::ReadFlag( const std::string& FileName, const char* flag )
{
    GetPot DataFile( FileName );
    M_flag = DataFile( flag, 0 );

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface_Data::ReadFlag                        flag: " << static_cast<Real>( M_flag ) << "\n";
#endif
}

template< class Operator >
inline void BCInterface_Data< Operator >::ReadType( const std::string& FileName, const char* type )
{
    GetPot DataFile( FileName );
    M_type = M_mapType[DataFile( type, "Essential" )];

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface_Data::ReadType                        type: " << M_type << " (" << DataFile(type, "Essential") << ")\n";
#endif
}

template< class Operator >
inline void BCInterface_Data< Operator >::ReadMode( const std::string& FileName, const char* mode )
{
    GetPot DataFile( FileName );
    M_mode = M_mapMode[DataFile( mode, "Full" )];

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface_Data::ReadMode                        mode: " << M_mode << " (" << DataFile(mode, "Full") << ")\n";
#endif
}

template< class Operator >
inline void BCInterface_Data< Operator >::ReadComV( const std::string& FileName, const char* component )
{
    GetPot DataFile( FileName );
    UInt componentSize = DataFile.vector_variable_size( component );

    M_comV.clear();
    M_comV.reserve( componentSize );

    for ( UInt j( 0 ); j < componentSize; ++j )
        M_comV.push_back( DataFile( component, 0, j ) );

#ifdef DEBUG
    std::stringstream output;
    output << "BCInterface_Data::ReadComV                        comV: ";
    for ( UInt i(0); i < static_cast<UInt>( M_comV.size() ); ++i )
        output << M_comV[i] << " ";
    Debug( 5020 ) << output.str() << "\n";
#endif

}

template< class Operator >
inline void BCInterface_Data< Operator >::ReadDirection( const std::string& FileName, const char* direction )
{
    GetPot DataFile( FileName );
    M_direction = DataFile( direction, " " );

#ifdef DEBUG
    Debug( 5020 ) << "BCInterface_Data::ReadDirection              direction: " << M_direction << "\n";
#endif
}

template< class Operator >
inline void BCInterface_Data< Operator >::ReadBase( const std::string& FileName, const std::string& base )
{
    for ( typename std::map< std::string, BCBaseList_Type >::iterator j = M_mapBase.begin(); j
            != M_mapBase.end(); ++j )
        if ( IsBase( FileName, ( base + j->first ).c_str() ) )
        {
            M_base.first = j->first;
            M_base.second = M_mapBase[j->first];

#ifdef DEBUG
            Debug( 5020 ) << "BCInterface_Data::ReadBase                        base: " << M_base.second << " (" << j->first << ")\n";
            Debug( 5020 ) << "                                           baseString: " << M_baseString << "\n";
#endif

            break;
        }
}

template< class Operator >
inline bool BCInterface_Data< Operator >::IsBase( const std::string& FileName, const char* base )
{
    GetPot DataFile( FileName );
    M_baseString = DataFile( base, " " );

    return DataFile.checkVariable( base );
}

} // Namespace LifeV

#endif /* BCInterface_Data_H */
