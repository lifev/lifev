/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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


#include <life/lifefem/bcCond.hpp>

namespace LifeV
{

// ============ BCBase ================
//! Constructor for BC
BCBase::BCBase( const std::string& name, const EntityFlag& flag,
                  const BCType& type, const BCMode& mode,
                  BCFunctionBase& bcf, const std::vector<ID>& comp )
    :
    _M_isUDep(false),
    _M_name( name ),
    _M_flag( flag ),
    _M_type( type ),
    _M_mode( mode ),
    _M_bcf( FactoryCloneBCFunction::instance().createObject( &bcf ) ),
    _M_dataVector( false ),
    _M_comp( comp ),
    _M_offset( -1 ),
    _M_finalised( false )
{
    if ( _M_mode != Component ) {
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
}

//! Constructor for BC without components for Scalar, Tangential or Normal  mode problems
BCBase::BCBase( const std::string& name,
                const EntityFlag&  flag,
                const BCType&      type,
                const BCMode&      mode,
                BCFunctionBase&    bcf ):
    _M_isUDep(false),
    _M_name( name ),
    _M_flag( flag ),
    _M_type( type ),
    _M_mode( mode ),
    _M_bcf( FactoryCloneBCFunction::instance().createObject( &bcf ) ),
    _M_dataVector( false ),
    _M_comp(),
    _M_offset( -1 ),
    _M_finalised( false )
{
    //if(type==3)//flux
    //{
          //	++BCBase::M_fluxes;
    //}
    UInt nComp;
    switch ( _M_mode = mode )
        {
        case Scalar:
            nComp = 1;
            _M_comp.reserve( nComp );
            _M_comp.push_back( 1 );
            break;
        case Tangential:
            nComp = nDimensions;
            _M_comp.reserve( nComp );
            for ( ID i = 1; i <= nComp; ++i )
                _M_comp.push_back( i );
            break;
        case Normal:
            nComp = nDimensions;
            _M_comp.reserve( nComp );
            for ( ID i = 1; i <= nComp; ++i )
                _M_comp.push_back( i );
            break;
        default:
            ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
        }
}


//! Constructor for BC without list of components for Full mode problems
BCBase::BCBase( const std::string& name,
                const EntityFlag&  flag,
                const BCType&      type,
                const BCMode&      mode,
                BCFunctionBase&    bcf,
                const UInt&        nComp )

    :
    _M_isUDep(false),
    _M_name( name ),
    _M_flag( flag ),
    _M_type( type ),
    _M_mode( mode ),
    _M_bcf( FactoryCloneBCFunction::instance().createObject( &bcf ) ),
    _M_dataVector( false ),
    _M_comp(),
    _M_offset( -1 ),
    _M_finalised( false )
{
    if ( _M_mode != Full ) {
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
    _M_comp.reserve( nComp );
    for ( ID i = 1; i <= nComp; ++i )
        _M_comp.push_back( i );

//     if(type==3)//flux
//       {
// 	++BCBase::M_fluxes;
//      }

}



//! Constructor for BC with data vector
BCBase::BCBase( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCVectorBase& bcv,
                const std::vector<ID>& comp )
    :
    _M_isUDep(false),
    _M_name( name ),
    _M_flag( flag ),
    _M_type( type ),
    _M_mode( mode ),
    _M_bcf(),_M_bcfUDep(),
    _M_bcv( FactoryCloneBCVector::instance().createObject( &bcv )  ),
    _M_dataVector( true ),
    _M_comp(comp),
    _M_offset( -1 ),
    _M_finalised( false )
{
    if ( mode != Component ) {
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
//     if(type==3)//flux
//       {
// 	++BCBase::M_fluxes;
//       }
}

//! Constructor for BC with data vector, without components for Scalar, Tangential or Normal  mode problems
BCBase::BCBase( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCVectorBase& bcv )
    :
    _M_isUDep(false),
    _M_name( name ),
    _M_flag( flag ),
    _M_type( type ),
    _M_mode( mode ),
    _M_bcf(),_M_bcfUDep(),
    _M_bcv( FactoryCloneBCVector::instance().createObject( &bcv ) ),
    _M_dataVector( true ),
    _M_comp(),
    _M_finalised( false )
{
    UInt nComp;
    switch ( _M_mode = mode )
    {
    case Scalar:
        nComp = 1;
        _M_comp.reserve( nComp );
        _M_comp.push_back( 1 );

        break;
    case Tangential:
        nComp = nDimensions - 1;
        _M_comp.reserve( nComp );
        for ( ID i = 1; i <= nComp; ++i )
            _M_comp.push_back( i );

        break;
    case Normal:
        nComp = 1;
        _M_comp.reserve( nComp );
        _M_comp.push_back( nDimensions );

        break;
    default:
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
//     if(type==3)//flux
//       {
// 	++BCBase::M_fluxes;
//       }
}


//! Constructor for BC with data vector, without list of components for Full mode problems
BCBase::BCBase( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCVectorBase& bcv,
                const UInt& nComp )
    :
    _M_isUDep(false),
    _M_name( name ),
    _M_flag( flag ),
    _M_type( type ),
    _M_mode( mode ),
    _M_bcf(),_M_bcfUDep(),
    _M_bcv( FactoryCloneBCVector::instance().createObject( &bcv ) ),
    _M_dataVector( true ),
    _M_comp(),
    _M_offset( -1 ),
    _M_finalised( false )
{
    if ( mode != Full ) {
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }

    _M_comp.reserve( nComp );
    for ( ID i = 1; i <= nComp; ++i )
        _M_comp.push_back( i );

//     if(type==3)//flux
//       {
// 	++BCBase::M_fluxes;
//       }
}

BCBase::BCBase( const std::string&     name,
                const EntityFlag&      flag,
                const BCType&          type,
                const BCMode&          mode,
                BCFunctionUDepBase&    bcf,
                const std::vector<ID>& comp ):
    _M_isUDep(true),
    _M_name( name ),
    _M_flag( flag ),
    _M_type( type ),
    _M_mode( mode ),
    _M_bcfUDep( FactoryCloneBCFunctionUDep::instance().createObject( &bcf ) ),
    _M_dataVector( false ),
    _M_comp( comp ),
    _M_finalised( false )
{
    if ( _M_mode != Component ) {
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
//     if(type==3)//flux
//       {
// 	++BCBase::M_fluxes;
//       }
}
BCBase::BCBase( const std::string&  name,
                const EntityFlag&   flag,
                const BCType&       type,
                const BCMode&       mode,
                BCFunctionUDepBase& bcf):
    _M_isUDep(true),
    _M_name( name ),
    _M_flag( flag ),
    _M_type( type ),
    _M_mode( mode ),
    _M_bcfUDep( FactoryCloneBCFunctionUDep::instance().createObject( &bcf ) ),
    _M_dataVector( false ),
    _M_comp(),
    _M_offset( -1 ),
    _M_finalised( false )
{

    UInt nComp;
    switch ( _M_mode = mode )
    {
        case Scalar:
            nComp = 1;
            _M_comp.reserve( nComp );
            _M_comp.push_back( 1 );
            break;
        case Tangential:
            nComp = nDimensions - 1;
            _M_comp.reserve( nComp );
            for ( ID i = 1; i <= nComp; ++i )
                _M_comp.push_back( i );
            break;
        case Normal:
            nComp = 1;
            _M_comp.reserve( nComp );
            _M_comp.push_back( nDimensions );
            break;
        default:
            ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
//     if(type==3)//flux
//       {
// 	++BCBase::M_fluxes;
//       }
}
BCBase::BCBase( const std::string&  name,
                const EntityFlag&   flag,
                const BCType&       type,
                const BCMode&       mode,
                BCFunctionUDepBase& bcf,
                const UInt&         nComp )
    :
    _M_isUDep(true),
    _M_name( name ),
    _M_flag( flag ),
    _M_type( type ),
    _M_mode( mode ),
    _M_bcfUDep( FactoryCloneBCFunctionUDep::instance().createObject( &bcf ) ),
    _M_dataVector( false ),
    _M_comp(),
    _M_offset( -1 ),
    _M_finalised( false )
{
    if ( _M_mode != Full ) {
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }

    _M_comp.reserve( nComp );
    for ( ID i = 1; i <= nComp; ++i )
        _M_comp.push_back( i );

//     if(type==3)//flux
//       {
// 	++BCBase::M_fluxes;
//       }
}





//! Destructor (we have a vector of pointers to ID's and Functors)
BCBase::~BCBase()
{
}



//! Assignment operator for BC (we have a vector of pointers to ID's
//and a pointer to user defined functions)
BCBase & BCBase::operator=( const BCBase& BCb )
{

    _M_name = BCb._M_name;
    _M_flag = BCb._M_flag;
    _M_type = BCb._M_type;
    _M_mode = BCb._M_mode;
    _M_finalised = BCb._M_finalised;
    _M_dataVector = BCb._M_dataVector;
    _M_bcfUDep=BCb._M_bcfUDep;
    _M_bcv = BCb._M_bcv;
    _M_bcf = BCb._M_bcf;
    _M_offset  = BCb._M_offset;
    _M_comp = BCb._M_comp;
    _M_isUDep=BCb._M_isUDep;

    // Important!!: The set member list0 is always empty at this
    // point, it is just an auxiliary container used at the moment of
    // the boundary update (see BCHandler::bdUpdate)

    // The list of ID's must be empty
    if ( !_M_idList.empty() || !BCb._M_idList.empty() ) {
        ERROR_MSG( "BCBase::operator= : The BC assigment operator does not work with lists of identifiers which are not empty" );
    }

    return *this;
}

//! Copy constructor for BC (we have a vector of pointers to ID's and a pointer to user defined functions)
BCBase::BCBase( const BCBase& BCb )
    :
    _M_isUDep(BCb._M_isUDep),
    _M_name( BCb._M_name ),
    _M_flag( BCb._M_flag ),
    _M_type( BCb._M_type ),
    _M_mode( BCb._M_mode ),
    _M_bcf( BCb._M_bcf ),
    _M_bcfUDep(BCb._M_bcfUDep),
    _M_bcv( BCb._M_bcv ),
    _M_dataVector( BCb._M_dataVector ),
    _M_comp( BCb._M_comp ),
    _M_offset   ( BCb._M_offset ),
    _M_finalised( BCb._M_finalised )
{
    // Important!!: The set member list0 is always empty at this point, it is just
    // an auxiliary container used at the moment of the boundary update (see BCHandler::bdUpdate)

    // The list of ID's must be empty
    if ( !_M_idList.empty() || !BCb._M_idList.empty() ) {
        ERROR_MSG( "BCBase::BCBase : The BC copy constructor does not work whith list of identifiers which are not empty" );
    }
}

//! Returns the BC name
std::string BCBase::name() const
{
    return _M_name;
}

//! Returns the BC associated flag
EntityFlag BCBase::flag() const
{
    return _M_flag;
}

//! Returns the BC type
BCType BCBase::type() const
{
    return _M_type;
}

//! Returns the BC mode
BCMode BCBase::mode() const
{
    return _M_mode;
}
bool BCBase::isUDep() const
{
  return _M_isUDep;
}

//! Returns the number of components involved in this boundary condition
UInt BCBase::numberOfComponents() const
{
    return _M_comp.size();
}

//! eturns the global i-th component involved in the boundary condition
ID BCBase::component( const ID i ) const
{
    ASSERT_BD( i >= 1 && i <= _M_comp.size() );
    return _M_comp[ i -1 ];
}

//! Returns wether the list is finalised and the vector of ID's is then accessible
bool BCBase::finalised() const
{
    return _M_finalised;
}

//! Overloading function operator by calling the (*_M_bcf)() user specified function
Real BCBase::operator() ( const Real& t, const Real& x, const Real& y,
                           const Real& z, const ID& i ) const
{
    return _M_bcf->operator() ( t,x, y, z, i );
}

//! Overloading function operator by calling the (*_M_bcf)() user specified function
/* new overloading for BCFunctionUDepending */
Real BCBase::operator() ( const Real& t, const Real& x, const Real& y,
                           const Real& z, const ID& i, const Real& u ) const
{
    /* is there a better way ? */
Debug(800)<<"debug800 in BCBase::operator(6x)\n";
   return _M_bcfUDep->operator()(t,x, y, z, i, u);
Debug(800)<<"debug800 out BCBase::operator(6x)\n";
}



//! Returns a pointer  to the user defined STL functor
const BCFunctionBase* BCBase::pointerToFunctor() const
{
    return _M_bcf.get();
}

//! Returns a pointer  to the BCVector
const BCVectorBase* BCBase::pointerToBCVector() const
{
    return _M_bcv.get();
}

//! Returns a pointer  to the user defined STL functor
const BCFunctionUDepBase* BCBase::pointerToFunctorUDep() const
{
    return _M_bcfUDep.get();
}

//! True is a data vector has been provided
bool BCBase::dataVector() const
{
    return _M_dataVector;
}

void
BCBase::setBCVector( BCVectorBase& __v )
{
    _M_bcv = boost::shared_ptr<BCVectorBase >( FactoryCloneBCVector::instance().createObject( &__v ) );
    _M_dataVector = true;
    _M_isUDep=false;
}

void
BCBase::setBCFunction( BCFunctionBase& __f )
{
    _M_bcf = boost::shared_ptr<BCFunctionBase>( FactoryCloneBCFunction::instance().createObject( &__f ) );
    _M_dataVector = false;
    _M_isUDep=false;
}

void
BCBase::setBCFunction( BCFunctionUDepBase& __f )
{
    _M_bcfUDep = boost::shared_ptr<BCFunctionUDepBase>( FactoryCloneBCFunctionUDep::instance().createObject( &__f ) );
    _M_dataVector = false;
    _M_isUDep=true;
}

Real BCBase::operator() ( const ID& iDof, const ID& iComp ) const
{
    if ( _M_dataVector )
        return ( *_M_bcv ) ( iDof, iComp );
    else
    {
        ERROR_MSG( "BCBase::operator() : A data vector must be specified before calling this method" );
        return 0.;
    }
}

//! Return true if mixte coefficient is bcVector, false otherside
bool  BCBase::ismixteVec()  const
{
    if ( _M_dataVector )
      {
     	return  (*_M_bcv).ismixteVec();
      }
    else
    {
        ERROR_MSG( "BCBase::mixte : A data vector must be specified before calling this method" );
        return 0.;
    }
}

//! Return true if beta coefficient  is bcVector, false otherside
bool BCBase::isbetaVec()   const
{
    if ( _M_dataVector )
      {

	return   (*_M_bcv).isbetaVec();
      }
    else
    {
        ERROR_MSG( "BCBase::beta: A data vector must be specified before calling this method" );
        return 0.;
    }
}

//! Return true if gamma coefficient is bcVector, false otherside
bool BCBase::isgammaVec()  const
{
    if ( _M_dataVector )
      {

	return (*_M_bcv).isgammaVec();
      }
    else
    {
        ERROR_MSG( "BCBase::gamma : A data vector must be specified before calling this method" );
        return 0.;
    }
}

//! Returns the value of the mixte coefficient (in BC Vector)
Real BCBase::mixteCoef() const
{
    if ( _M_dataVector )
        return ( *_M_bcv ).mixteCoef();
    else
    {
        ERROR_MSG( "BCBase::mixteCoef : A data vector must be specified before calling this method" );
        return 0.;
    }

}

//! Returns the value of the beta coefficient (in BC Vector)
Real BCBase::betaCoef() const
{
    if ( _M_dataVector )
        return ( *_M_bcv ).betaCoef();
    else
    {
        ERROR_MSG( "BCBase::mixteCoef : A data vector must be specified before calling this method" );
        return 0.;
    }

}

//! Returns the value of the gamma coefficient (in BC Vector)
Real BCBase::gammaCoef() const
{
    if ( _M_dataVector )
        return ( *_M_bcv ).gammaCoef();
    else
    {
        ERROR_MSG( "BCBase::mixteCoef : A data vector must be specified before calling this method" );
        return 0.;
    }

}

//! Returns the value of the mixte coefficient vector (in BC Vector) M.Prosi
Real BCBase::MixteVec( const ID& iDof, const ID& iComp ) const
{
    if ( _M_dataVector )
        return ( *_M_bcv).MixteVec( iDof, iComp );
    else
    {
        ERROR_MSG( "BCBase::MixteVec : A data vector must be specified before calling this method" );
        return 0.;
    }

}


//! Returns the value of the beta coefficient vector (in BC Vector)
Real BCBase::BetaVec( const ID& iDof, const ID& iComp ) const
{
    if ( _M_dataVector )
        return ( *_M_bcv).BetaVec( iDof, iComp );
    else
    {
        ERROR_MSG( "BCBase::MixteVec : A data vector must be specified before calling this method" );
        return 0.;
    }

}

//! Returns the value of the gamma coefficient vector (in BC Vector)
Real BCBase::GammaVec( const ID& iDof, const ID& iComp ) const
{
    if ( _M_dataVector )
        return ( *_M_bcv).GammaVec( iDof, iComp );
    else
    {
        ERROR_MSG( "BCBase::MixteVec : A data vector must be specified before calling this method" );
        return 0.;
    }

}


//! Returns a pointer  to the i-th elements in the (finalised) list
//! (counting from 1 ' a la FORTRAN')
const IdentifierBase*
BCBase::operator() ( const ID& i ) const
{
    return this->operator[] ( i-1 );
}

//! Returns a pointer to the i-th elements in the (finalised) list
//! (counting from 0 ' a la C')
const IdentifierBase*
BCBase::operator[] ( const Index_t& i ) const
{
    ASSERT_PRE( _M_finalised, "BC List should be finalised before being accessed" );
    ASSERT_BD( i < _M_idList.size() );
    return _M_idList[ i ].get();
}


//! Add a new identifier to the preliminary list of Identifiers
void
BCBase::addIdentifier( IdentifierBase* iden )
{
    list0.insert( boost::shared_ptr<IdentifierBase>( iden ) );
}


//! Transfer between the list and vector containers of ID's
void
BCBase::finalise()
{
    if ( ! list0.empty() )
    {
        _M_idList.clear();
        _M_idList.reserve( list0.size() );
        std::copy( list0.begin(), list0.end(), std::inserter( _M_idList, _M_idList.end() ) );
        list0.clear();
    }
    _M_finalised = true;
}


//! Returns the liste size
UInt
BCBase::list_size() const
{
    return _M_idList.size();
}

//! Output
std::ostream&
BCBase::showMe( bool verbose, std::ostream & out ) const
{
    out << "********************************" << std::endl;
    out << "BC Name              : " << _M_name << std::endl;
    out << "Flag                 : " << _M_flag << std::endl;
    out << "Type                 : " << _M_type << std::endl;
    out << "Mode                 : " << _M_mode << std::endl;
    out << "Number of components : " << _M_comp.size() << std::endl;
    out << "List of components   : ";
    for ( Index_t i = 0; i < _M_comp.size(); ++i )
        out << _M_comp[ i ] << " ";
    out << std::endl;
    out << "Offset               : " << _M_offset << std::endl;
    out << "Number of stored ID's: " << _M_idList.size() << std::endl;

    if ( verbose && _M_finalised )
    {
        unsigned int count( 0 ), lines( 10 );
        out << "IDs in list";
        for ( std::vector<boost::shared_ptr<IdentifierBase> >::const_iterator i = _M_idList.begin();
              i != _M_idList.end(); i++ )
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
            _M_bcv->showMe( verbose, out );
        }
    }

    out << "********************************" << std::endl;
    return out;
}

//  UInt BCBase::M_fluxes(0);

}
