/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*!
  \file bcCond.cc
  \brief Implementations for bcCond.h
  \version 1.0
  \author M.A. Fernandez
  \date 07/2002

*/


#include "bcCond.hpp"

namespace LifeV
{

// ============ BCBase ================
//! Constructor for BC
BCBase::BCBase( const std::string& name, const EntityFlag& flag,
                  const BCType& type, const BCMode& mode,
                  BCFunctionBase& bcf, const std::vector<ID>& comp )
{
    _name = name;
    _flag = flag;
    _type = type;
    if ( mode != Component )
        ERROR_MSG( "You should use a more specific constructor for this mode" );
    _mode = mode;

    switch ( type )
    {
        case Essential:
        case Natural:
            _bcf = new BCFunctionBase( bcf );
            break;
        case Mixte:
            // We have to hold a pointer to a mixte functor
            _bcf = new BCFunctionMixte( static_cast<BCFunctionMixte&>( bcf ) );
            break;
    }
    _dataVector = 0;

    UInt nComp = comp.size();
    _comp.reserve( nComp );
    _comp.insert( _comp.end(), comp.begin(), comp.end() );

    _finalised = 0;
}

//! Constructor for BC without components for Scalar, Tangential or Normal  mode problems
BCBase::BCBase( const std::string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
                  BCFunctionBase& bcf )
{
    _name = name;
    _flag = flag;
    _type = type;

    UInt nComp;
    switch ( _mode = mode )
    {
    case Scalar:
        nComp = 1;
        _comp.reserve( nComp );
        _comp.push_back( 1 );
        switch ( type )
        {
        case Essential:
        case Natural:
            _bcf = new BCFunctionBase( bcf );
            break;
        case Mixte:
            // We have to hold a pointer to a mixte functor
            _bcf = new BCFunctionMixte( static_cast<BCFunctionMixte&>( bcf ) );
            break;
        }
        break;
    case Tangential:
        nComp = nDimensions - 1;
        _comp.reserve( nComp );
        for ( ID i = 1; i <= nComp; ++i )
            _comp.push_back( i );
        _bcf = new BCFunctionBase( bcf );
        break;
    case Normal:
        nComp = 1;
        _comp.reserve( nComp );
        _comp.push_back( nDimensions );
        _bcf = new BCFunctionBase( bcf );
        break;
    default:
        ERROR_MSG( "You should use a more specific constructor for this mode" );
    }
    _dataVector = 0;
    _finalised = 0;
}


//! Constructor for BC without list of components for Full mode problems
BCBase::BCBase( const std::string& name, const EntityFlag& flag, const BCType& type,
                  const BCMode& mode, BCFunctionBase& bcf, const UInt& nComp )
{
    _name = name;
    _flag = flag;
    _type = type;

    if ( mode != Full )
        ERROR_MSG( "You should use a more specific constructor for this mode" );
    _mode = mode;


    switch ( type )
    {
    case Essential:
    case Natural:
        _bcf = new BCFunctionBase( bcf );
        break;
    case Mixte:
        // We have to hold a pointer to a mixte functor
        _bcf = new BCFunctionMixte( static_cast<BCFunctionMixte&>( bcf ) );
        break;
    }

    _comp.reserve( nComp );
    for ( ID i = 1; i <= nComp; ++i )
        _comp.push_back( i );

    _dataVector = 0;
    _finalised = 0;
}



//! Constructor for BC with data vector
BCBase::BCBase( const std::string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
                  BCVectorBase& bcv, const std::vector<ID>& comp )
{
    _name = name;
    _flag = flag;
    _type = type;
    if ( mode != Component )
        ERROR_MSG( "You should use a more specific constructor for this mode" );
    _mode = mode;

    _dataVector = 1;
    _bcv = &bcv;


    UInt nComp = comp.size();
    _comp.reserve( nComp );
    _comp.insert( _comp.end(), comp.begin(), comp.end() );

    _finalised = 0;
}

//! Constructor for BC with data vector, without components for Scalar, Tangential or Normal  mode problems
BCBase::BCBase( const std::string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
                  BCVectorBase& bcv )
{
    _name = name;
    _flag = flag;
    _type = type;

    UInt nComp;
    switch ( _mode = mode )
    {
    case Scalar:
        nComp = 1;
        _comp.reserve( nComp );
        _comp.push_back( 1 );

        _dataVector = 1;
        _bcv = &bcv;

        break;
    case Tangential:
        nComp = nDimensions - 1;
        _comp.reserve( nComp );
        for ( ID i = 1; i <= nComp; ++i )
            _comp.push_back( i );

        _dataVector = 1;
        _bcv = &bcv;

        break;
    case Normal:
        nComp = 1;
        _comp.reserve( nComp );
        _comp.push_back( nDimensions );

        _dataVector = 1;
        _bcv = &bcv;

        break;
    default:
        ERROR_MSG( "You should use a more specific constructor for this mode" );
    }
    _finalised = 0;
}


//! Constructor for BC with data vector, without list of components for Full mode problems
BCBase::BCBase( const std::string& name, const EntityFlag& flag, const BCType& type,
                  const BCMode& mode, BCVectorBase& bcv, const UInt& nComp )
{
    _name = name;
    _flag = flag;
    _type = type;

    if ( mode != Full )
        ERROR_MSG( "You should use a more specific constructor for this mode" );
    _mode = mode;

    _dataVector = 1;
    _bcv = &bcv;

    _comp.reserve( nComp );
    for ( ID i = 1; i <= nComp; ++i )
        _comp.push_back( i );

    _finalised = 0;
}


//! Destructor (we have a vector of pointers to ID's and Functors)
BCBase::~BCBase()
{

    // \warning Important: The set member list0 is always empty at this
    // point, it is just an auxiliary container used at the moment of
    // the boundary update (see BCHandler::bdUpdate)

    switch ( _type )
    {
    case Essential:
        // Deleting identifiers
        for ( IDIterator i = _idList.begin();i != _idList.end();++i )
            delete static_cast<IdentifierEssential*>( *i );

        break;
    case Natural:
    case Mixte:
        for ( IDIterator i = _idList.begin();i != _idList.end();++i )
            // Deleting identifiers
            delete static_cast<IdentifierNatural*>( *i );
        break;
    }

    if ( _dataVector )
    {
        _dataVector = 0;
    }
    else
    {
        // Deleting user defined functors
        switch ( _type )
        {
        case Natural:
        case Essential:
            delete _bcf;
            break;
        case Mixte:
            delete static_cast<BCFunctionMixte*>( _bcf );
            break;
        }
    }
}



//! Assignment operator for BC (we have a vector of pointers to ID's
//and a pointer to user defined functions)
BCBase & BCBase::operator=( const BCBase& BCb )
{

    _name = BCb._name;
    _flag = BCb._flag;
    _type = BCb._type;
    _mode = BCb._mode;
    _finalised = BCb._finalised;
    _dataVector = BCb._dataVector;

    switch ( _type )
    {
    case Essential:
    case Natural:
        if ( _dataVector )
        {
            _bcv = BCb._bcv;
        }
        else
        {
            _bcf = new BCFunctionBase( *( BCb._bcf ) );
        }
        break;
    case Mixte:
        if ( _dataVector )
        {
            _bcv = BCb._bcv;
        }
        else
        {
            _bcf = new BCFunctionMixte( *( static_cast<BCFunctionMixte*>( BCb._bcf ) ) );
        }
        break;
    }

    _comp = BCb._comp;

    // Important!!: The set member list0 is always empty at this
    // point, it is just an auxiliary container used at the moment of
    // the boundary update (see BCHandler::bdUpdate)

    // The list of ID's must be empty
    if ( !_idList.empty() || !BCb._idList.empty() )
        ERROR_MSG( "The BC assigment operator does not work with lists of identifiers which are not empty" );

    return *this;
}

//! Copy constructor for BC (we have a vector of pointers to ID's and a pointer to user defined functions)
BCBase::BCBase( const BCBase& BCb )
{
    _name = BCb._name;
    _flag = BCb._flag;
    _type = BCb._type;
    _mode = BCb._mode;
    _finalised = BCb._finalised;
    _dataVector = BCb._dataVector;

    switch ( _type )
    {
    case Essential:
    case Natural:
        if ( _dataVector )
        {
            _bcv = BCb._bcv;
        }
        else
        {
            _bcf = new BCFunctionBase( *( BCb._bcf ) );
        }
        break;
    case Mixte:
        if ( _dataVector )
        {
            _bcv = BCb._bcv;
        }
        else
        {
            _bcf = new BCFunctionMixte( *( static_cast<BCFunctionMixte*>( BCb._bcf ) ) );
        }
        break;
    }
    _comp = BCb._comp;

    // Important!!: The set member list0 is always empty at this point, it is just
    // an auxiliary container used at the moment of the boundary update (see BCHandler::bdUpdate)

    // The list of ID's must be empty
    if ( !_idList.empty() || !BCb._idList.empty() )
        ERROR_MSG( "The BC copy constructor does not work whith list of identifiers which are not empty" );
}

//! Returns the BC name
std::string BCBase::name() const
{
    return _name;
}

//! Returns the BC associated flag
EntityFlag BCBase::flag() const
{
    return _flag;
}

//! Returns the BC type
BCType BCBase::type() const
{
    return _type;
}

//! Returns the BC mode
BCMode BCBase::mode() const
{
    return _mode;
}

//! Returns the number of components involved in this boundary condition
UInt BCBase::numberOfComponents() const
{
    return _comp.size();
}

//! eturns the global i-th component involved in the boundary condition
ID BCBase::component( const ID i ) const
{
    ASSERT_BD( i >= 1 && i <= _comp.size() );
    return _comp[ i -1 ];
}

//! Returns wether the list is finalised and the vector of ID's is then accessible
bool BCBase::finalised() const
{
    return _finalised;
}

//! Overloading function operator by calling the (*_bcf)() user specified function
Real BCBase::operator() ( const Real& t, const Real& x, const Real& y,
                           const Real& z, const ID& i ) const
{
    return _bcf->operator() ( t,x, y, z, i );
}



//! Returns a pointer  to the user defined STL functor
const BCFunctionBase* BCBase::pointerToFunctor() const
{
    return _bcf;
}

//! Returns a pointer  to the BCVector
const BCVectorBase* BCBase::pointerToBCVector() const
{
    return _bcv;
}

//! True is a data vector has been provided
bool BCBase::dataVector() const
{
    return _dataVector;
}



Real BCBase::operator() ( const ID& iDof, const ID& iComp ) const
{
    if ( _dataVector )
        return ( *_bcv ) ( iDof, iComp );
    else
    {
        ERROR_MSG( "A data vector must be specified before calling this method" );
        return 0.;
    }
}


//! Returns the value of the mixte coefficient (in BC Vector)
Real BCBase::mixteCoef() const
{
    if ( _dataVector )
        return ( *_bcv ).mixteCoef();
    else
    {
        ERROR_MSG( "A data vector must be specified before calling this method" );
        return 0.;
    }

}

//! Returns a pointer  to the i-th elements in the (finalised) list
//! (counting from 1 ' a la FORTRAN')
const IdentifierBase* BCBase::operator() ( const ID& i ) const
{
    return this->operator[] ( i-1 );
}

//! Returns a pointer to the i-th elements in the (finalised) list
//! (counting from 0 ' a la C')
const IdentifierBase* BCBase::operator[] ( const Index_t& i ) const
{
    ASSERT_PRE( _finalised, "BC List should be finalised before being accessed" );
    ASSERT_BD( i >= 0 && i < _idList.size() );
    return _idList[ i ];
}


//! Add a new identifier to the preliminary list of Identifiers
void BCBase::addIdentifier( IdentifierBase* iden )
{
    list0.insert( iden );
}


//! Transfer between the list and vector containers of ID's
void BCBase::finalise()
{
    if ( ! list0.empty() )
    {
        _idList.clear();
        _idList.reserve( list0.size() );
        std::copy( list0.begin(), list0.end(), std::inserter( _idList, _idList.end() ) );
        list0.clear();
    }
    _finalised = true;
}


//! Returns the liste size
UInt BCBase::list_size() const
{
    return _idList.size();
}

//! Output
std::ostream& BCBase::showMe( bool verbose, std::ostream & out ) const
{
    out << "********************************" << std::endl;
    out << "BC Name: " << _name << std::endl;
    out << "Flag: " << _flag << std::endl;
    out << "Type: " << _type << std::endl;
    out << "Mode: " << _mode << std::endl;
    out << "Number of components: " << _comp.size() << std::endl;
    out << "List of components: ";
    for ( Index_t i = 0; i < _comp.size(); ++i )
        out << _comp[ i ] << " ";
    out << std::endl;
    out << "Number of stored ID's: " << _idList.size() << std::endl;

    if ( verbose && _finalised )
    {
        unsigned int count( 0 ), lines( 10 );
        out << "IDs in list";
        for ( std::vector<IdentifierBase*>::const_iterator i = _idList.begin(); i != _idList.end(); i++ )
        {
            if ( count++ % lines == 0 )
            {
                out << std::endl;
            }
            out << ( *i ) ->id() << "  ";
        }
        if ( count % lines != 0 )
        {
            out << std::endl;
        }
        if ( dataVector() )
        {
            _bcv->showMe( verbose, out );
        }
    }

    out << "********************************" << std::endl;
    return out;
}



}
