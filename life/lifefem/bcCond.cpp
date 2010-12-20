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
    @file
    @brief classes to handle boundary conditions

    @author M.A. Fernandez
    @contributor Lucia Mirabella <lucia.mirabell@gmail.com>
    @contributor Mauro Perego <perego.mauro@gmail.com>
    @maintainer Mauro Perego <perego.mauro@gmail.com>

    @date 06-2002

    @date 11-2002

  This file contains the classes which may be used to store boundary
  conditions. A boundary condition object will have the following
  elements:
<ol>
  <li> a name identifying a specific BC,

  <li> a flag identifying a specific part of the mesh boundary,

  <li> a type (Essential, Natural, Mixte, Flux, Resistance),

  <li> a mode of implementation (Scalar, Full, Component, Normal,
     Tangential, Resistance, Directional),

  <li> a functor holding the data function,

  <li> a bool vector  describing the components involved in this boundary condition

  <li> a list of pointers to identifiers allowing the user to know to
     which DOF the boundary condition applies.
</ol>

 */

#include <life/lifecore/life.hpp>
#include <life/lifefem/bcCond.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

BCBase::BCBase()
{
}

BCBase::BCBase( const std::string& name, const bcFlag_Type& flag,
                const bcType_Type& type, const bcMode_Type& mode,
                BCFunctionBase& bcFunction, const bcComponentsVec_Type& components )
        :
        M_name( name ),
        M_flag( flag ),
        M_type( type ),
        M_mode( mode ),
        M_components( components ),
        M_bcFunction( FactoryCloneBCFunction::instance().createObject( &bcFunction ) ),
        M_bcFunctionFEVectorDependent(),
        M_bcVector(),
        M_isStored_BcVector( false ),
        M_isStored_BcFunctionVectorDependent(false),
        M_offset( -1 ),
        M_finalized( false )
{
    if ( M_mode != Component )
    {
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
}

BCBase::BCBase( const bcName_Type& name,
                const bcFlag_Type&  flag,
                const bcType_Type&      type,
                const bcMode_Type&      mode,
                BCFunctionBase&    bcFunction ):
        M_name( name ),
        M_flag( flag ),
        M_type( type ),
        M_mode( mode ),
        M_components(),
        M_bcFunction( FactoryCloneBCFunction::instance().createObject( &bcFunction ) ),
        M_bcFunctionFEVectorDependent(),
        M_bcVector(),
        M_isStored_BcVector( false ),
        M_isStored_BcFunctionVectorDependent(false),
        M_offset( -1 ),
        M_finalized( false )
{
    UInt numberOfComponents;
    switch ( M_mode = mode )
    {
    case Scalar:
        numberOfComponents = 1;
        M_components.reserve( numberOfComponents );
        M_components.push_back( 1 );
        break;
    case Tangential:
        numberOfComponents = nDimensions;
        M_components.reserve( numberOfComponents );
        for ( ID i = 1; i <= numberOfComponents; ++i )
            M_components.push_back( i );
        break;
    case Normal:
        // Normal Essential boundary condition (cf Gwenol Grandperrin Master Thesis)
        if (type == Essential)
        {
            numberOfComponents = 1;
            M_components.reserve( numberOfComponents );
            M_components.push_back( nDimensions );
        }
        else
        {
            numberOfComponents = nDimensions;
            M_components.reserve( numberOfComponents );
            for ( ID i = 1; i <= numberOfComponents; ++i )
                M_components.push_back( i );
        }
        break;
    case Directional:
        numberOfComponents = 1;
        M_components.reserve( numberOfComponents );
        M_components.push_back( nDimensions );
        break;
    default:
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
}

BCBase::BCBase( const bcName_Type& name,
                const bcFlag_Type&  flag,
                const bcType_Type&      type,
                const bcMode_Type&      mode,
                BCFunctionBase&    bcFunction,
                const UInt&        numberOfComponents )
        :
        M_name( name ),
        M_flag( flag ),
        M_type( type ),
        M_mode( mode ),
        M_components(),
        M_bcFunction( FactoryCloneBCFunction::instance().createObject( &bcFunction ) ),
        M_bcFunctionFEVectorDependent(),
        M_bcVector(),
        M_isStored_BcVector( false ),
        M_isStored_BcFunctionVectorDependent(false),
        M_offset( -1 ),
        M_finalized( false )
{
    if ( M_mode != Full )
    {
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
    M_components.reserve( numberOfComponents );
    for ( ID i = 1; i <= numberOfComponents; ++i )
        M_components.push_back( i );

}


BCBase::BCBase( const bcName_Type& name,
                const bcFlag_Type& flag,
                const bcType_Type& type,
                const bcMode_Type& mode,
                BCVectorBase& bcVector,
                const bcComponentsVec_Type& components )
        :
        M_name( name ),
        M_flag( flag ),
        M_type( type ),
        M_mode( mode ),
        M_components(components),
        M_bcFunction(),
        M_bcFunctionFEVectorDependent(),
        M_bcVector( FactoryCloneBCVector::instance().createObject( &bcVector )  ),
        M_isStored_BcVector( true ),
        M_isStored_BcFunctionVectorDependent(false),
        M_offset( -1 ),
        M_finalized( false )
{
    if ( mode != Component )
    {
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
}

BCBase::BCBase( const bcName_Type& name,
                const bcFlag_Type& flag,
                const bcType_Type& type,
                const bcMode_Type& mode,
                BCVectorBase& bcVector )
        :
        M_name( name ),
        M_flag( flag ),
        M_type( type ),
        M_mode( mode ),
        M_components(),
        M_bcFunction(),
        M_bcFunctionFEVectorDependent(),
        M_bcVector( FactoryCloneBCVector::instance().createObject( &bcVector ) ),
        M_isStored_BcVector( true ),
        M_isStored_BcFunctionVectorDependent(false),
        M_finalized( false )
{
    UInt numberOfComponents;
    switch ( M_mode = mode )
    {
    case Scalar:
        numberOfComponents = 1;
        M_components.reserve( numberOfComponents );
        M_components.push_back( 1 );

        break;
    case Tangential:
        numberOfComponents = nDimensions - 1;
        M_components.reserve( numberOfComponents );
        for ( ID i = 1; i <= numberOfComponents; ++i )
            M_components.push_back( i );

        break;
    case Normal:
        numberOfComponents = 1;
        M_components.reserve( numberOfComponents );
        M_components.push_back( nDimensions );

        break;
    case Directional:
        numberOfComponents = 1;
        M_components.reserve( numberOfComponents );
        M_components.push_back( nDimensions );
        break;
    default:
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
}


BCBase::BCBase( const bcName_Type& name,
                const bcFlag_Type& flag,
                const bcType_Type& type,
                const bcMode_Type& mode,
                BCVectorBase& bcVector,
                const UInt& numberOfComponents )
        :
        M_name( name ),
        M_flag( flag ),
        M_type( type ),
        M_mode( mode ),
        M_components(),
        M_bcFunction(),
        M_bcFunctionFEVectorDependent(),
        M_bcVector( FactoryCloneBCVector::instance().createObject( &bcVector ) ),
        M_isStored_BcVector( true ),
        M_isStored_BcFunctionVectorDependent(false),
        M_offset( -1 ),
        M_finalized( false )
{
    if ( mode != Full )
    {
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }

    M_components.reserve( numberOfComponents );
    for ( ID i = 1; i <= numberOfComponents; ++i )
        M_components.push_back( i );

}

BCBase::BCBase( const bcName_Type&     name,
                const bcFlag_Type&      flag,
                const bcType_Type&          type,
                const bcMode_Type&          mode,
                BCFunctionUDepBase&    bcFunctionFEVectorDependent,
                const bcComponentsVec_Type& components ):
        M_name( name ),
        M_flag( flag ),
        M_type( type ),
        M_mode( mode ),
        M_components( components ),
        M_bcFunction(),
        M_bcFunctionFEVectorDependent( FactoryCloneBCFunctionUDep::instance().createObject( &bcFunctionFEVectorDependent ) ),
        M_bcVector(),
        M_isStored_BcVector( false ),
        M_isStored_BcFunctionVectorDependent(true),
        M_finalized( false )
{
    if ( M_mode != Component )
    {
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }
}

BCBase::BCBase( const bcName_Type&  name,
                const bcFlag_Type&   flag,
                const bcType_Type&       type,
                const bcMode_Type&       mode,
                BCFunctionUDepBase& bcFunctionFEVectorDependent):
        M_name( name ),
        M_flag( flag ),
        M_type( type ),
        M_mode( mode ),
        M_components(),
        M_bcFunction(),
        M_bcFunctionFEVectorDependent( FactoryCloneBCFunctionUDep::instance().createObject( &bcFunctionFEVectorDependent ) ),
        M_bcVector(),
        M_isStored_BcVector( false ),
        M_isStored_BcFunctionVectorDependent(true),
        M_offset( -1 ),
        M_finalized( false )
{

    UInt numberOfComponents;
    switch ( M_mode = mode )
    {
    case Scalar:
        numberOfComponents = 1;
        M_components.reserve( numberOfComponents );
        M_components.push_back( 1 );
        break;
    case Tangential:
        numberOfComponents = nDimensions - 1;
        M_components.reserve( numberOfComponents );
        for ( ID i = 1; i <= numberOfComponents; ++i )
            M_components.push_back( i );
        break;
    case Normal:
        numberOfComponents = 1;
        M_components.reserve( numberOfComponents );
        M_components.push_back( nDimensions );
        break;
    case Directional:
        numberOfComponents = 1;
        M_components.reserve( numberOfComponents );
        M_components.push_back( nDimensions );
        break;
    default:
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }

}

BCBase::BCBase( const bcName_Type&  name,
                const bcFlag_Type&   flag,
                const bcType_Type&       type,
                const bcMode_Type&       mode,
                BCFunctionUDepBase& bcFunctionFEVectorDependent,
                const UInt&         numberOfComponents )
        :
        M_name( name ),
        M_flag( flag ),
        M_type( type ),
        M_mode( mode ),
        M_components(),
        M_bcFunction(),
        M_bcFunctionFEVectorDependent( FactoryCloneBCFunctionUDep::instance().createObject( &bcFunctionFEVectorDependent ) ),
        M_bcVector(),
        M_isStored_BcVector( false ),
        M_isStored_BcFunctionVectorDependent(true),
        M_offset( -1 ),
        M_finalized( false )
{
    if ( M_mode != Full )
    {
        ERROR_MSG( "BCBase::BCBase: You should use a more specific constructor for this mode" );
    }

    M_components.reserve( numberOfComponents );
    for ( ID i = 1; i <= numberOfComponents; ++i )
        M_components.push_back( i );

}


BCBase::BCBase( const BCBase& bcBase )
        :
        M_name( bcBase.M_name ),
        M_flag( bcBase.M_flag ),
        M_type( bcBase.M_type ),
        M_mode( bcBase.M_mode ),
        M_components( bcBase.M_components ),
        M_bcFunction( bcBase.M_bcFunction ),
        M_bcFunctionFEVectorDependent(bcBase.M_bcFunctionFEVectorDependent),
        M_bcVector( bcBase.M_bcVector ),
        M_isStored_BcVector( bcBase.M_isStored_BcVector ),
        M_isStored_BcFunctionVectorDependent(bcBase.M_isStored_BcFunctionVectorDependent),
        M_idSet( ),
        M_idList( ),
        M_offset   ( bcBase.M_offset ),
        M_finalized( bcBase.M_finalized )
{
    // Important!!: The set member M_idSet is always empty at this point, it is just
    // an auxiliary container used at the moment of the boundary update (see BCHandler::bcUpdate)

    // The list of ID's must be empty
    if ( !M_idList.empty() || !bcBase.M_idList.empty() )
    {
        ERROR_MSG( "BCBase::BCBase : The BC copy constructor does not work with list of identifiers which are not empty" );
    }
}


BCBase::~BCBase()
{
}

// ===================================================
// Methods
// ===================================================

ID BCBase::component( const ID i ) const
{
    ASSERT_BD( i >= 1 && i <= M_components.size() );
    return M_components[ i -1 ];
}

bool  BCBase::ismixteVec()  const
{
    if ( M_isStored_BcVector )
    {
        return  (*M_bcVector).ismixteVec();
    }
    else
    {
        ERROR_MSG( "BCBase::mixte : A data vector must be specified before calling this method" );
        return 0.;
    }
}

bool BCBase::isbetaVec()   const
{
    if ( M_isStored_BcVector )
    {

        return   (*M_bcVector).isbetaVec();
    }
    else
    {
        ERROR_MSG( "BCBase::beta: A data vector must be specified before calling this method" );
        return 0.;
    }
}


Real BCBase::MixteVec( const ID& iDof, const ID& iComponent ) const
{
    if ( M_isStored_BcVector )
        return ( *M_bcVector).MixteVec( iDof, iComponent );
    else
    {
        ERROR_MSG( "BCBase::MixteVec : A data vector must be specified before calling this method" );
        return 0.;
    }
}

Real BCBase::BetaVec( const ID& iDof, const ID& iComponent ) const
{
    if ( M_isStored_BcVector )
        return ( *M_bcVector).BetaVec( iDof, iComponent );
    else
    {
        ERROR_MSG( "BCBase::MixteVec : A data vector must be specified before calling this method" );
        return 0.;
    }

}


const BCFunctionBase* BCBase::pointerToFunctor() const
{
    return M_bcFunction.get();
}

const BCFunctionUDepBase* BCBase::pointerToFunctorUDep() const
{
    return M_bcFunctionFEVectorDependent.get();
}

const BCVectorBase* BCBase::pointerToBCVector() const
{
    return M_bcVector.get();
}

void
BCBase::addIdentifier( IdentifierBase* identifierToAddPtr )
{
    M_idSet.insert( boost::shared_ptr<IdentifierBase>( identifierToAddPtr ) );
}


UInt
BCBase::list_size() const
{
    return M_idList.size();
}

std::ostream&
BCBase::showMe( bool verbose, std::ostream & out ) const
{
    out << "********************************" << std::endl;
    out << "BC Name              : " << M_name << std::endl;
    out << "Flag                 : " << M_flag << std::endl;
    out << "Type                 : " << M_type << std::endl;
    out << "Mode                 : " << M_mode << std::endl;
    out << "Number of components : " << M_components.size() << std::endl;
    out << "List of components   : ";
    for ( ID i = 0; i < M_components.size(); ++i )
        out << M_components[ i ] << " ";
    out << std::endl;
    out << "Offset               : " << M_offset << std::endl;
    out << "Number of stored ID's: " << M_idList.size() << std::endl;

    if ( verbose && M_finalized )
    {
        unsigned int count( 0 ), lines( 10 );
        out << "IDs in list";
        for ( std::vector<boost::shared_ptr<IdentifierBase> >::const_iterator i = M_idList.begin();
                i != M_idList.end(); i++ )
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
            M_bcVector->showMe( verbose, out );
        }
    }

    out << "********************************" << std::endl;
    return out;
}


// ===================================================
// Operators
// ===================================================

BCBase & BCBase::operator=( const BCBase& BCb )
{

    M_name = BCb.M_name;
    M_flag = BCb.M_flag;
    M_type = BCb.M_type;
    M_mode = BCb.M_mode;
    M_finalized = BCb.M_finalized;
    M_isStored_BcVector = BCb.M_isStored_BcVector;
    M_bcFunctionFEVectorDependent=BCb.M_bcFunctionFEVectorDependent;
    M_bcVector = BCb.M_bcVector;
    M_bcFunction = BCb.M_bcFunction;
    M_offset  = BCb.M_offset;
    M_components = BCb.M_components;
    M_isStored_BcFunctionVectorDependent=BCb.M_isStored_BcFunctionVectorDependent;

    // Important!!: The set member M_idSet is always empty at this
    // point, it is just an auxiliary container used at the moment of
    // the boundary update (see BCHandler::bcUpdate)

    // The list of ID's must be empty
    if ( !M_idList.empty() || !BCb.M_idList.empty() )
    {
        ERROR_MSG( "BCBase::operator= : The BC assigment operator does not work with lists of identifiers which are not empty" );
    }

    return *this;
}

const IdentifierBase*
BCBase::operator[] ( const ID& i ) const
{
    ASSERT_PRE( M_finalized, "BC List should be finalized before being accessed" );
    ASSERT_BD( i < M_idList.size() );
    return M_idList[ i ].get();
}

const IdentifierBase*
BCBase::operator() ( const ID& i ) const
{
    return this->operator[] ( i-1 );
}

Real BCBase::operator() ( const Real& t, const Real& x, const Real& y,
                          const Real& z, const ID& iComponent ) const
{
    return M_bcFunction->operator() ( t,x, y, z, iComponent );
}

Real BCBase::operator() ( const Real& t, const Real& x, const Real& y,
                          const Real& z, const ID& iComponent, const Real& u ) const
{
    /* is there a better way ? */
    Debug(800)<<"debug800 in BCBase::operator(6x)\n";
    return M_bcFunctionFEVectorDependent->operator()(t,x, y, z, iComponent, u);
    Debug(800)<<"debug800 out BCBase::operator(6x)\n";
}

Real BCBase::operator() ( const ID& iDof, const ID& iComponent ) const
{
    if ( M_isStored_BcVector )
        return ( *M_bcVector ) ( iDof, iComponent );
    else
    {
        ERROR_MSG( "BCBase::operator() : A data vector must be specified before calling this method" );
        return 0.;
    }
}


// ===================================================
// Set Methods
// ===================================================

void
BCBase::setBCVector( const BCVectorBase& bcVector )
{
    M_bcVector = boost::shared_ptr<BCVectorBase >( FactoryCloneBCVector::instance().createObject( &bcVector ) );
    M_isStored_BcVector = true;
    M_isStored_BcFunctionVectorDependent=false;
}

void
BCBase::setBCFunction( const BCFunctionBase& bcFunction )
{
    M_bcFunction = boost::shared_ptr<BCFunctionBase>( FactoryCloneBCFunction::instance().createObject( &bcFunction ) );
    M_isStored_BcVector = false;
    M_isStored_BcFunctionVectorDependent=false;
}

void
BCBase::setBCFunction( const BCFunctionUDepBase& bcFunctionFEVectorDependent )
{
    M_bcFunctionFEVectorDependent = boost::shared_ptr<BCFunctionUDepBase>( FactoryCloneBCFunctionUDep::instance().createObject( &bcFunctionFEVectorDependent ) );
    M_isStored_BcVector = false;
    M_isStored_BcFunctionVectorDependent=true;
}


// ===================================================
// Get Methods
// ===================================================

std::string BCBase::name() const
{
    return M_name;
}

bcFlag_Type BCBase::flag() const
{
    return M_flag;
}

bcType_Type BCBase::type() const
{
    return M_type;
}

bcMode_Type BCBase::mode() const
{
    return M_mode;
}

UInt BCBase::numberOfComponents() const
{
    return M_components.size();
}

Real BCBase::mixteCoef() const
{
    if ( M_isStored_BcVector )
        return ( *M_bcVector ).mixteCoef();
    else
    {
        ERROR_MSG( "BCBase::mixteCoef : A data vector must be specified before calling this method" );
        return 0.;
    }

}

Real BCBase::resistanceCoef() const
{
    if ( M_isStored_BcVector )
        return ( *M_bcVector ).resistanceCoef();
    else
    {
        ERROR_MSG( "BCBase::resistanceCoef : A data vector must be specified before calling this method" );
        return 0.;
    }

}

Real BCBase::betaCoef() const
{
    if ( M_isStored_BcVector )
        return ( *M_bcVector ).betaCoef();
    else
    {
        ERROR_MSG( "BCBase::mixteCoef : A data vector must be specified before calling this method" );
        return 0.;
    }

}

bool BCBase::dataVector() const
{
    return M_isStored_BcVector;
}

bool BCBase::finalized() const
{
    return M_finalized;
}

bool BCBase::isUDep() const
{
    return M_isStored_BcFunctionVectorDependent;
}


// ===================================================
// Private Methods
// ===================================================

void
BCBase::finalize()
{
    if ( ! M_idSet.empty() )
    {
        M_idList.clear();
        M_idList.reserve( M_idSet.size() );
        std::copy( M_idSet.begin(), M_idSet.end(), std::inserter( M_idList, M_idList.end() ) );
        M_idSet.clear();
    }
    M_finalized = true;
}


}
