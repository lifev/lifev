/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-11

  Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano

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
/**
   \file bcHandler.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-11
*/
#include <sstream>
#include <stdexcept>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/if.hpp>

#include <bcHandler.hpp>

namespace LifeV
{

// ============ BCHandler ================

//! Constructor doing nothing (the user must call setNumber(..))
BCHandler::BCHandler() : _bdUpdateDone( 0 ), _fullEssential( 0 )
{}


//! Constructor taking the number of BC to be stored
BCHandler::BCHandler( const ID& nbc, const bool& fullEssential )
    : _bdUpdateDone( 0 ),
      _fullEssential( fullEssential )
{
    _bcList.reserve( nbc );
}

//! Constructor taking the number of BC to be stored
BCHandler::BCHandler( const ID& nbc )
    : _bdUpdateDone( 0 ),
      _fullEssential( 0 )
{
    _bcList.reserve( nbc );
}

//! How many BC stored?
Index_t BCHandler::size() const
{
    return _bcList.size();
}

//! Is there no BC into the list?
bool BCHandler::empty() const
{
    return _bcList.empty();
}

//! Adding new BC to the list with user defined functions
void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCFunctionBase& bcf,
                       const std::vector<ID>& comp )
{
    // Adding BC
    _bcList.push_back( BCBase( name, flag, type, mode, bcf, comp ) );
    std::sort( _bcList.begin(), _bcList.end() );
}
void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCFunctionBase& bcf )
{
    // Adding BC
    _bcList.push_back( BCBase( name, flag, type, mode, bcf ) );
    std::sort( _bcList.begin(), _bcList.end() );
}

void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCFunctionBase& bcf,
                       const UInt& nComp )
{
    // Adding BC
    _bcList.push_back( BCBase( name, flag, type, mode, bcf, nComp ) );
    std::sort( _bcList.begin(), _bcList.end() );
}


//! Adding new BC to the list with data vectors
void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCVectorBase& bcv,
                       const std::vector<ID>& comp )
{
    // Adding BC
    _bcList.push_back( BCBase( name, flag, type, mode, bcv, comp ) );
    std::sort( _bcList.begin(), _bcList.end() );
}
void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCVectorBase& bcv )
{
    // Adding BC
    _bcList.push_back( BCBase( name, flag, type, mode, bcv ) );
    std::sort( _bcList.begin(), _bcList.end() );
}

void
BCHandler::addBC( const std::string& name,
                  const EntityFlag& flag,
                  const BCType& type,
                  const BCMode& mode,
                  BCVectorBase& bcv,
                  const UInt& nComp )
{
    // Adding BC
    _bcList.push_back( BCBase( name, flag, type, mode, bcv, nComp ) );
    std::sort( _bcList.begin(), _bcList.end() );
}

BCBase*
BCHandler::findBC( std::string const& __name )
{
    BCBase* __bc = 0;
    std::for_each( _bcList.begin(),
                   _bcList.end(),
                   boost::lambda::if_then( boost::lambda::bind( &BCBase::name, boost::lambda::_1 ) == __name,
                                           boost::lambda::var( __bc ) = &boost::lambda::_1 ) );

    //! handle invalid name case: ie we didnot find the name in the _bcList
    if ( !__bc )
    {
        std::ostringstream __ex;
        __ex << "Invalid name for BC to be modified : " << __name << "\n"
             << "The list of available BCs is:\n";
        std::for_each( _bcList.begin(),
                       _bcList.end(),
                       std::cout << boost::lambda::bind( &BCBase::name, boost::lambda::_1 )
                       << boost::lambda::constant( "\n" ) );
        throw std::invalid_argument( __ex.str() );
    }
    return __bc;
}

BCBase*
BCHandler::findBC( int lab)
{
    BCBase* __bc = 0;
    std::for_each( _bcList.begin(),
                   _bcList.end(),
                   boost::lambda::if_then( boost::lambda::bind( &BCBase::flag, boost::lambda::_1 ) == lab,
                                           boost::lambda::var( __bc ) = &boost::lambda::_1 ) );

    return __bc;
}

void
BCHandler::modifyBC( std::string const& __name, BCFunctionBase& __bcf )
{
    BCBase* __bc = findBC( __name );

    __bc->setBCFunction( __bcf );
}
void
BCHandler::modifyBC( std::string const& __name, BCVectorBase& __bcv )
{
    BCBase* __bc = findBC( __name );

    __bc->setBCVector( __bcv );
}

void
BCHandler::modifyBC( int lab, BCFunctionBase& __bcf )
{
    BCBase* __bc = findBC( lab );

    __bc->setBCFunction( __bcf );
}
void
BCHandler::modifyBC( int lab, BCVectorBase& __bcv )
{
    BCBase* __bc = findBC( lab );

    __bc->setBCVector( __bcv );
}

// returns true if the bdUpdate has been done before
bool BCHandler::bdUpdateDone() const
{
    return _bdUpdateDone;
}

//! returns true if all the stored BC are of Essential type
bool BCHandler::hasOnlyEssential() const
{
    if ( empty() )
    {
        return _fullEssential;
    }
    else
    {
        bool listFullEssential = true;
        for( ConstIterator it=_bcList.begin(); it!=_bcList.end(); ++it )
        {
            listFullEssential &=
                ( it->type() == Essential ) &&
                ( it->mode() == Full );
        }
        if ( listFullEssential != _fullEssential )
        {
            std::ostringstream __ex;
            __ex << "BCHandler::hasOnlyEssential(): state is not consistent:\n"
                 << "flag from constructor says    " << _fullEssential << "\n"
                 << "added boundary conditions say " << listFullEssential;
            throw std::logic_error( __ex.str() );
        }
        return listFullEssential;
    }
}

//! Extracting BC from the list
BCBase& BCHandler::operator[] ( const Index_t& i )
{
    return _bcList[ i ];
}


const BCBase& BCHandler::operator[] ( const Index_t& i ) const
{
    return _bcList[ i ];
}

BCBase& BCHandler::GetBCWithFlag(const EntityFlag& aFlag){

  Index_t i;

  for(i = 0; i <= _bcList.size(); i++){
    if(aFlag == _bcList[i].flag()){
      break;
    }
  }

  return _bcList[i];
}

const BCBase& BCHandler::GetBCWithFlag(const EntityFlag& aFlag) const {

  Index_t i;

  for(i = 0; i <= _bcList.size(); i++){
    if(aFlag == _bcList[i].flag()){
      break;
    }
  }
  return _bcList[i];
}


UInt BCHandler::getBCbyName(const std::string _BCName) const
{
    int iBC = -1;

    for (UInt jBC = 0; jBC < _bcList.size(); jBC++)
        if (_bcList[jBC].name() == _BCName) iBC = jBC;

    return iBC;
}


BCType BCHandler::boundaryType(const EntityFlag& aFlag) const
{
  BCType CurrType;

  for(UInt i = 0; i <= _bcList.size(); i++){
    if(aFlag == _bcList[i].flag()){
      CurrType = _bcList[i].type();
      break;
    }
  }

  return CurrType;
}
//! Ouput
std::ostream & BCHandler::showMe( bool verbose, std::ostream & out ) const
{
    out << " Boundary Conditions Handler ====>" << std::endl;
    out << " Number of BC stored " << size() << std::endl;

    out << " List => " << std::endl;
    for ( UInt i = 0; i < _bcList.size(); ++i )
    {
        _bcList[ i ].showMe( verbose, out );
    }
    out << " <===========================>" << std::endl;
    return out;
}
}
