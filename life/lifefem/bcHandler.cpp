//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2096 EPFL, Politecnico di Milano, INRIA
               2006-2010 EPFL, Politecnico di Milano

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
/**
   @file
   @brief File containing the implementation of a class for BC handling

   @author Miguel Fernandez
   @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>

   @date 2004-10-11
*/

#include <sstream>
#include <stdexcept>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/if.hpp>

#include <life/lifefem/bcHandler.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
BCHandler::BCHandler( const ID& nbc, const BCHints& hint ):
    M_bdUpdateDone    ( 0 ),
    M_hint            ( hint ),
    M_offset          ( 0 ),
    M_notFoundMarkers ( )
{
    if ( nbc > 0 )
        M_bcList.reserve( nbc );
}

BCHandler::BCHandler( const BCHandler& BCh ):
    M_bdUpdateDone    ( BCh.M_bdUpdateDone ),
    M_hint            ( BCh.M_hint ),
    M_bcList          ( BCh.M_bcList ),
    M_offset          ( BCh.M_offset ),
    M_notFoundMarkers ( BCh.M_notFoundMarkers )
{
}

// ===================================================
// Operators
// ===================================================
BCHandler&
BCHandler::operator = (const BCHandler &BCh)
{
    if (this != &BCh)
    {
        M_bdUpdateDone = BCh.M_bdUpdateDone;
        M_hint         = BCh.M_hint;
        M_bcList       = BCh.M_bcList;
    }

    return *this;
}

BCBase&
BCHandler::operator[] ( const Index_t& i )
{
    return M_bcList[ i ];
}

const BCBase&
BCHandler::operator[] ( const Index_t& i ) const
{
    return M_bcList[ i ];
}

// ===================================================
// Methods
// ===================================================
void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCFunctionBase& bcf,
                       const std::vector<ID>& comp )
{
    M_bcList.push_back( BCBase( name, flag, type, mode, bcf, comp ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}


void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCFunctionBase& bcf )
{
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
    M_bcList.push_back( BCBase( name, flag, type, mode, bcv, comp ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}

void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCVectorBase& bcv )
{
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
    M_bcList.push_back( BCBase( name, flag, type, mode, bcv, nComp ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}

void BCHandler::addBC( const std::string& name,
                       const EntityFlag& flag,
                       const BCType& type,
                       const BCMode& mode,
                       BCFunctionUDepBase& bcf )
{
    M_bcList.push_back( BCBase( name, flag, type, mode, bcf ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
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
BCHandler::modifyBC( Int lab, BCFunctionBase& __bcf )
{
    BCBase* __bc = findBC( lab );

    __bc->setBCFunction( __bcf );
}

void
BCHandler::modifyBC( Int lab, BCVectorBase& __bcv )
{
    BCBase* __bc = findBC( lab );

    __bc->setBCVector( __bcv );
}

void
BCHandler::modifyBC( Int lab, BCFunctionUDepBase& __bcf )
{
    BCBase* __bc = findBC( lab );

    __bc->setBCFunction( __bcf );
}

void
BCHandler::merge( BCHandler& bch )
{
    sumOffsets();
    bch.sumOffsets();
    M_bcList.insert(M_bcList.end(), bch.M_bcList.begin(), bch.M_bcList.end());
    M_bdUpdateDone = M_bdUpdateDone && bch.M_bdUpdateDone;
    M_offset = 0;
}

void
BCHandler::showMe( bool verbose, std::ostream& out ) const
{
    out << " Boundary Conditions Handler ====>" << std::endl;
    out << " Number of BC stored " << size() << std::endl;

    out << " List => " << std::endl;
    for ( UInt i = 0; i < M_bcList.size(); ++i )
        M_bcList[ i ].showMe( verbose, out );
    out << " <===========================>" << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void
BCHandler::setOffset( const UInt& offset )
{
    M_offset = offset;
}

void
BCHandler::setOffset( std::string const& name, Int offset )
{
  BCBase* bc = findBC( name );

  if (bc == 0)
      std::cout << "BCHandler::setOffset : BC " << name << " not found ... ";

  bc->setOffset(offset);
}

// ===================================================
// Get Methods
// ===================================================
BCType
BCHandler::boundaryType(const EntityFlag& aFlag) const
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

BCBase&
BCHandler::GetBCWithFlag(const EntityFlag& aFlag)
{
    Index_t i;

    for(i = 0; i <= M_bcList.size(); i++)
        if(aFlag == M_bcList[i].flag())
            break;

    return M_bcList[i];
}

const BCBase&
BCHandler::GetBCWithFlag(const EntityFlag& aFlag) const
{
    Index_t i;

    for(i = 0; i <= M_bcList.size(); i++)
        if(aFlag == M_bcList[i].flag())
            break;

    return M_bcList[i];
}

std::vector<BCName>
BCHandler::getBCWithType( const BCType& type )
{
    std::vector<BCName> vectorName;

    for( size_type i = 0; i < M_bcList.size(); ++i )
        if( M_bcList[i].type() == type)
            vectorName.push_back( M_bcList[i].name() );

    return vectorName;
}

UInt
BCHandler::getNumberBCWithType( const BCType& type )
{
    UInt typeNumber = 0;

    for( size_type i = 0; i < M_bcList.size(); ++i )
        if( M_bcList[i].type() == type)
            ++typeNumber;

    return typeNumber;
}

UInt
BCHandler::getBCbyName(const std::string __BCName) const
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

UInt
BCHandler::offset() const
{
    return M_offset;
}

BCHandler::BCBase_Iterator
BCHandler::begin()
{
    return M_bcList.begin();
}

BCHandler::BCBase_Iterator
BCHandler::end()
{
    return M_bcList.end();
}

Index_t
BCHandler::size() const
{
    return M_bcList.size();
}

bool
BCHandler::empty() const
{
    return M_bcList.empty();
}

bool
BCHandler::bdUpdateDone() const
{
    return M_bdUpdateDone;
}

bool
BCHandler::hasOnlyEssential() const
{
    if ( empty() )
        return M_hint == HINT_BC_ONLY_ESSENTIAL;
    else if ( M_hint == HINT_BC_ONLY_ESSENTIAL )
    {
        if ( listHasOnlyEssential() != true )
        {
            std::ostringstream __ex;
            __ex << "BCHandler::hasOnlyEssential(): state is not consistent:"
                 << "\nhint from constructor says:    yes"
                 << "\nadded boundary conditions say: no";
            std::cerr << std::endl << "Throwing exception:\n"
                      << __ex.str() << std::endl;
            throw std::logic_error( __ex.str() );
        }
        return true;
    }

    return listHasOnlyEssential();
}

// ===================================================
// Private Methods
// ===================================================
BCBase*
BCHandler::findBC( const std::string& __name )
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
BCHandler::findBC( const Int& lab)
{
    BCBase* __bc = 0;
    std::for_each( M_bcList.begin(),
                   M_bcList.end(),
                   boost::lambda::if_then( boost::lambda::bind( &BCBase::flag, boost::lambda::_1 ) == lab,
                                           boost::lambda::var( __bc ) = &boost::lambda::_1 ) );

    return __bc;
}

bool
BCHandler::listHasOnlyEssential() const
{
    std::map<EntityFlag, EssentialStatus> statusMap;
    for( BCBase_ConstIterator it = M_bcList.begin(); it != M_bcList.end(); ++it )
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
                case Directional:
                    status.setDirectional();
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

void
BCHandler::sumOffsets()
{
    for ( BCBase_Iterator it = M_bcList.begin(); it != M_bcList.end(); ++it )
        it->setOffset(it->offset()+M_offset);
}

} // namespace LifeV
