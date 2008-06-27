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

#include <life/lifefem/bcHandler.hpp>

namespace LifeV
{

BCHandler::BCHandler( const ID& nbc, BCHints hint ):
    M_bdUpdateDone( 0 ),
    M_hint( hint )
{
    if ( nbc > 0 )
    {
        M_bcList.reserve( nbc );
    }
}


BCHandler::BCHandler( const BCHandler &BCh):
    M_bdUpdateDone(BCh.M_bdUpdateDone),
    M_hint(BCh.M_hint),
    M_bcList(BCh.M_bcList)
{
}


Index_t BCHandler::size() const
{
    return M_bcList.size();
}


bool BCHandler::empty() const
{
    return M_bcList.empty();
}


void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCFunctionBase& bcf,
                       const std::vector<ID>& comp )
{
    // Adding BC
    M_bcList.push_back( BCBase( name, flag, type, mode, bcf, comp ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}


void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCFunctionBase& bcf )
{
    // Adding BC
    M_bcList.push_back( BCBase( name, flag, type, mode, bcf ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}

void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCFunctionBase& bcf,
                       const UInt& nComp )
{
    // Adding BC
    M_bcList.push_back( BCBase( name, flag, type, mode, bcf, nComp ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}


void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCVectorBase& bcv,
                       const std::vector<ID>& comp )
{
    // Adding BC
    M_bcList.push_back( BCBase( name, flag, type, mode, bcv, comp ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}
void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCVectorBase& bcv )
{
    // Adding BC
    M_bcList.push_back( BCBase( name, flag, type, mode, bcv ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
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
    M_bcList.push_back( BCBase( name, flag, type, mode, bcv, nComp ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}
void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCFunctionUDepBase& bcf )
{
    // Adding BC
    M_bcList.push_back( BCBase( name, flag, type, mode, bcf ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}


BCBase*
BCHandler::findBC( std::string const& __name )
{
    BCBase* __bc = 0;
    std::for_each( M_bcList.begin(),
                   M_bcList.end(),
                   boost::lambda::if_then( boost::lambda::bind( &BCBase::name, boost::lambda::_1 ) == __name,
                                           boost::lambda::var( __bc ) = &boost::lambda::_1 ) );

    //! handle invalid name case: ie we didnot find the name in the M_bcList
    if ( !__bc )
    {
        std::ostringstream __ex;
        __ex << "Invalid name for BC to be modified : " << __name << "\n"
             << "The list of available BCs is:\n";
        std::for_each( M_bcList.begin(),
                       M_bcList.end(),
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
    std::for_each( M_bcList.begin(),
                   M_bcList.end(),
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
BCHandler::modifyBC( std::string const& __name, BCFunctionUDepBase& __bcf )
{
    BCBase* __bc = findBC( __name );

    __bc->setBCFunction( __bcf );
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
void
BCHandler::modifyBC( int lab, BCFunctionUDepBase& __bcf )
{
    BCBase* __bc = findBC( lab );

    __bc->setBCFunction( __bcf );
}


bool BCHandler::bdUpdateDone() const
{
    return M_bdUpdateDone;
}


bool BCHandler::hasOnlyEssential() const
{
    if ( empty() )
    {
        return M_hint == HINT_BC_ONLY_ESSENTIAL;
    }
    else
    {
        bool storedOnlyEssential = ( M_hint == HINT_BC_ONLY_ESSENTIAL );
        bool listOnlyEssential = listHasOnlyEssential();

        if ( listOnlyEssential != storedOnlyEssential )
        {
            std::ostringstream __ex;
            __ex << "BCHandler::hasOnlyEssential(): state is not consistent:"
                 << "\nhint from constructor says    " << storedOnlyEssential
                 << "\nadded boundary conditions say " << listOnlyEssential;
            std::cerr << std::endl << "Throwing exception:\n"
                      << __ex.str() << std::endl;
            throw std::logic_error( __ex.str() );
        }
        return storedOnlyEssential;
    }
}


// operator overloading

BCHandler& BCHandler::operator = (const BCHandler &BCh)
{
    if (this != &BCh)
    {
        M_bdUpdateDone = BCh.M_bdUpdateDone;
        M_hint         = BCh.M_hint;
        M_bcList       = BCh.M_bcList;
    }

    return *this;
}



BCBase& BCHandler::operator[] ( const Index_t& i )
{
    return M_bcList[ i ];
}


const BCBase& BCHandler::operator[] ( const Index_t& i ) const
{
    return M_bcList[ i ];
}


BCBase& BCHandler::GetBCWithFlag(const EntityFlag& aFlag){

  Index_t i;

  for(i = 0; i <= M_bcList.size(); i++){
    if(aFlag == M_bcList[i].flag()){
      break;
    }
  }

  return M_bcList[i];
}


const BCBase& BCHandler::GetBCWithFlag(const EntityFlag& aFlag) const {

    Index_t i;

    for(i = 0; i <= M_bcList.size(); i++){
        if(aFlag == M_bcList[i].flag()){
            break;
        }
    }
    return M_bcList[i];
}


/*const*/ UInt BCHandler::getBCbyName(const std::string __BCName) const
{
    UInt iBC( ( UInt )-1 );

    for (UInt jBC = 0; jBC < M_bcList.size(); jBC++)
        if (M_bcList[jBC].name() == __BCName)
            iBC = jBC;

    if ( iBC == UInt( -1 ) )
    {
        std::ostringstream __ex;
        __ex << __BCName << " was not found in this Boundary conditions set\n"
             << "This set contains \n";
        for ( UInt i = 0; i < M_bcList.size(); ++i )
        {
            M_bcList[ i ].showMe( true, __ex );
        }
        throw std::invalid_argument( __ex.str() );
    }

    return iBC;
}


BCType BCHandler::boundaryType(const EntityFlag& aFlag) const
{
    BCType CurrType(Natural);

    for(UInt i = 0; i <= M_bcList.size(); i++){
        if(aFlag == M_bcList[i].flag()){
            CurrType = M_bcList[i].type();
            break;
        }
    }

  return CurrType;
}


std::ostream & BCHandler::showMe( bool verbose, std::ostream & out ) const
{
    out << " Boundary Conditions Handler ====>" << std::endl;
    out << " Number of BC stored " << size() << std::endl;

    out << " List => " << std::endl;
    for ( UInt i = 0; i < M_bcList.size(); ++i )
    {
        M_bcList[ i ].showMe( verbose, out );
    }
    out << " <===========================>" << std::endl;
    return out;
}


bool BCHandler::listHasOnlyEssential() const
{
    std::map<EntityFlag, EssentialStatus> statusMap;
    for( ConstIterator it = M_bcList.begin(); it != M_bcList.end(); ++it )
    {
        // make sure that this flag is in the map
        EssentialStatus& status = statusMap[it->flag()];
        if ( it->type() == Essential )
        {
            switch ( it->mode() )
            {
                case Scalar:
                    status.setAllComponents();
                    break;
                case Full:
                    status.setAllComponents();
                    break;
                case Component:
                    {
                        UInt nComp = it->numberOfComponents();
                        for (UInt iComp=1; iComp<=nComp; ++iComp)
                        {
                            status.setComponent( it->component(iComp) );
                        }
                    }
                    break;
                case Normal:
                    status.setNormal();
                    break;
                case Tangential:
                    status.setTangential();
                    break;
                default:
                    {
                        std::ostringstream __ex;
                        __ex << "BCHandler::hasOnlyEssential(): BC mode "
                             << "unknown\n";
                        std::cerr << std::endl << "Throwing exception:\n"
                                  << __ex.str() << std::endl;
                        throw std::logic_error( __ex.str() );
                    }
                    break;
            }
        }
    }

    bool listOnlyEssential = true;
    for (std::map<EntityFlag, EssentialStatus>::const_iterator
             it = statusMap.begin(); it != statusMap.end(); ++it)
    {
        listOnlyEssential &= it->second.isEssential();
    }
    return listOnlyEssential;
}

} // namespace LifeV
