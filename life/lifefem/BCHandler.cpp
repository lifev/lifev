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
    @brief File containing BCHandler class for handling boundary conditions

    @author Miguel Fernandez <miguel.fernandez@inria.fr>
    @contributor Christophe Prud'homme <christophe.prudhomme@epfl.ch>
    @contributor Mauro Perego <perego.mauro@gmail.com>
    @maintainer Mauro Perego <perego.mauro@gmail.com>

    @date 10-11-2004
 *///@HEADER

#include <sstream>
#include <stdexcept>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/if.hpp>

#include <life/lifefem/bcHandler.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
BCHandler::BCHandler():
        M_bcUpdateDone    ( 0 ),
        M_offset          ( 0 )
{
}

BCHandler::BCHandler( const BCHandler& BCh ):
        M_bcUpdateDone    ( false ), //TODO change this! (related with BCBase copy constructor)
        M_bcList          ( BCh.M_bcList ),
        M_offset          ( BCh.M_offset ),
        M_notFoundMarkers ( BCh.M_notFoundMarkers )
{
}


BCHandler::~BCHandler()
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
        M_bcUpdateDone    = BCh.M_bcUpdateDone;
        M_bcList          = BCh.M_bcList;
        M_notFoundMarkers = BCh.M_notFoundMarkers;
    }

    return *this;
}

BCBase&
BCHandler::operator[] ( const ID& i )
{
    return M_bcList[ i ];
}

const BCBase&
BCHandler::operator[] ( const ID& i ) const
{
    return M_bcList[ i ];
}

// ===================================================
// Methods
// ===================================================

void
BCHandler::addBC( const bcName_Type& name,
                  const bcFlag_Type& flag,
                  const bcType_Type& type,
                  const bcMode_Type& mode,
                  BCFunctionBase& bcFunction,
                  const bcComponentsVec_Type& components )
{
    M_bcList.push_back( BCBase( name, flag, type, mode, bcFunction, components ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}


void
BCHandler::addBC( const bcName_Type& name,
                  const bcFlag_Type& flag,
                  const bcType_Type& type,
                  const bcMode_Type& mode,
                  BCFunctionBase& bcFunction )
{
    M_bcList.push_back( BCBase( name, flag, type, mode, bcFunction ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}

void
BCHandler::addBC( const bcName_Type& name,
                  const bcFlag_Type& flag,
                  const bcType_Type& type,
                  const bcMode_Type& mode,
                  BCFunctionBase& bcFunction,
                  const UInt& numComponents )
{
    M_bcList.push_back( BCBase( name, flag, type, mode, bcFunction, numComponents ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}

void
BCHandler::addBC( const bcName_Type& name,
                  const bcFlag_Type& flag,
                  const bcType_Type& type,
                  const bcMode_Type& mode,
                  BCVectorBase& bcVector,
                  const bcComponentsVec_Type& numComponents )
{
    M_bcList.push_back( BCBase( name, flag, type, mode, bcVector, numComponents ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}

void
BCHandler::addBC( const bcName_Type& name,
                  const bcFlag_Type& flag,
                  const bcType_Type& type,
                  const bcMode_Type& mode,
                  BCVectorBase& bcVector )
{
    M_bcList.push_back( BCBase( name, flag, type, mode, bcVector ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}

void
BCHandler::addBC( const bcName_Type& name,
                  const bcFlag_Type& flag,
                  const bcType_Type& type,
                  const bcMode_Type& mode,
                  BCVectorBase& bcVector,
                  const UInt& numComponents )
{
    M_bcList.push_back( BCBase( name, flag, type, mode, bcVector, numComponents ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}

void
BCHandler::addBC( const bcName_Type& name,
                  const bcFlag_Type& flag,
                  const bcType_Type& type,
                  const bcMode_Type& mode,
                  BCFunctionUDepBase& bcUDepFunction )
{
    M_bcList.push_back( BCBase( name, flag, type, mode, bcUDepFunction ) );
    std::sort( M_bcList.begin(), M_bcList.end() );
}

void
BCHandler::modifyBC( bcName_Type const& name, BCFunctionBase const& bcFunction )
{
    BCBase* bcBasePtr = M_findBC( name );

    bcBasePtr->setBCFunction( bcFunction );
}

void
BCHandler::modifyBC( bcName_Type const& name, BCVectorBase const& bcVector )
{
    BCBase* bcBasePtr = M_findBC( name );

    bcBasePtr->setBCVector( bcVector );
}

void
BCHandler::modifyBC( std::string const& name, BCFunctionUDepBase const& bcUDepFunction )
{
    BCBase* bcBasePtr = M_findBC( name );

    bcBasePtr->setBCFunction( bcUDepFunction );
}

void
BCHandler::modifyBC( bcFlag_Type const& aFlag, BCFunctionBase const& bcFunction )
{
    BCBase* bcBasePtr = M_findBC( aFlag );

    bcBasePtr->setBCFunction( bcFunction );
}

void
BCHandler::modifyBC( bcFlag_Type const& aFlag, BCVectorBase const& bcVector )
{
    BCBase* bcBasePtr = M_findBC( aFlag );

    bcBasePtr->setBCVector( bcVector );
}

void
BCHandler::modifyBC( bcFlag_Type const& aFlag, BCFunctionUDepBase const& bcFunction )
{
    BCBase* bcBasePtr = M_findBC( aFlag );

    bcBasePtr->setBCFunction( bcFunction );
}

void
BCHandler::merge( BCHandler& bcHandler )
{
    M_sumOffsets();
    bcHandler.M_sumOffsets();
    M_bcList.insert(M_bcList.end(), bcHandler.M_bcList.begin(), bcHandler.M_bcList.end());
    M_bcUpdateDone = M_bcUpdateDone && bcHandler.M_bcUpdateDone;
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
BCHandler::setOffset( const bcName_Type& name, Int offset )
{
    BCBase* bc = M_findBC( name );

    if (bc == 0)
        std::cout << "BCHandler::setOffset : BC " << name << " not found ... ";

    bc->setOffset(offset);
}


BCBase&
BCHandler::findBCWithFlag(const bcFlag_Type& aFlag)
{
    ID i;

    for (i = 0; i <= M_bcList.size(); i++)
        if (aFlag == M_bcList[i].flag())
            break;

    return M_bcList[i];
}

const BCBase&
BCHandler::findBCWithFlag(const bcFlag_Type& aFlag) const
{
    ID i;

    for (i = 0; i <= M_bcList.size(); i++)
        if (aFlag == M_bcList[i].flag())
            break;

    return M_bcList[i];
}

std::vector<bcName_Type>
BCHandler::findAllBCWithType( const bcType_Type& type ) const
{
    std::vector<bcName_Type> vectorName;

    for ( std::size_t i = 0; i < M_bcList.size(); ++i )
        if ( M_bcList[i].type() == type)
            vectorName.push_back( M_bcList[i].name() );

    return vectorName;
}

UInt
BCHandler::numberOfBCWithType( const bcType_Type& type ) const
{
    UInt typeNumber = 0;

    for ( std::size_t i = 0; i < M_bcList.size(); ++i )
        if ( M_bcList[i].type() == type)
            ++typeNumber;

    return typeNumber;
}

UInt
BCHandler::findBCWithName(bcName_Type const & name) const
{
    UInt iBC( 0 );

    for ( ; iBC < M_bcList.size(); iBC++ )
        if (M_bcList[iBC].name() == name)
            break;

    if ( iBC == M_bcList.size() )
    {
        std::ostringstream __ex;
        __ex << name << " was not found in this Boundary conditions set\n"
        << "This set contains \n";
        for ( UInt i = 0; i < M_bcList.size(); ++i )
        {
            M_bcList[ i ].showMe( true, __ex );
        }
        throw std::invalid_argument( __ex.str() );
    }

    return iBC;
}



bool
BCHandler::hasOnlyEssential() const
{
    std::map<bcFlag_Type, std::set<ID> > nonEssentialConditions;
    std::set<ID> nonEssentialComponents;
    for (UInt i=1; i<=nDimensions; i++)
        nonEssentialComponents.insert(i);

    for ( bcBaseConstIterator_Type it = M_bcList.begin(); it != M_bcList.end(); ++it )
        nonEssentialConditions.insert(std::make_pair(it->flag(), nonEssentialComponents) );

    for ( bcBaseConstIterator_Type it = begin(); it != end(); ++it )
    {
        if ( it->type() == Essential )
        {
            switch (it->mode())
            {
            case Full:
            case Normal:
                nonEssentialConditions.erase(it->flag());
                break;
            case Scalar:
                nonEssentialConditions.find(it->flag())->second.erase(1);
                break;
            case Component:
                for ( UInt iComp = 1; iComp <= it->numberOfComponents(); ++iComp )
                    nonEssentialConditions.find(it->flag())->second.erase( it->component(iComp) );
                if ( nonEssentialConditions.find(it->flag())->second.empty() )
                    nonEssentialConditions.erase(it->flag());
                break;
            default:
                break;
            }
        }
    }
    return ( nonEssentialConditions.empty() );
}

// ===================================================
// Private Methods
// ===================================================
BCBase*
BCHandler::M_findBC( bcName_Type const& name )
{
    BCBase* bcBasePtr = 0;
    std::for_each( M_bcList.begin(),
                   M_bcList.end(),
                   boost::lambda::if_then( boost::lambda::bind( &BCBase::name, boost::lambda::_1 ) == name,
                                           boost::lambda::var( bcBasePtr ) = &boost::lambda::_1 ) );

    //! handle invalid name case: ie we didnot find the name in the M_bcList
    if ( !bcBasePtr )
    {
        std::ostringstream __ex;
        __ex << "Invalid name for BC to be modified : " << name << "\n"
        << "The list of available BCs is:\n";
        std::for_each( M_bcList.begin(),
                       M_bcList.end(),
                       std::cout << boost::lambda::bind( &BCBase::name, boost::lambda::_1 )
                                 << boost::lambda::constant( "\n" ) );
        throw std::invalid_argument( __ex.str() );
    }
    return bcBasePtr;
}

BCBase*
BCHandler::M_findBC( bcFlag_Type const& aFlag)
{
    BCBase* bcBasePtr = 0;
    std::for_each( M_bcList.begin(),
                   M_bcList.end(),
                   boost::lambda::if_then( boost::lambda::bind( &BCBase::flag, boost::lambda::_1 ) == aFlag,
                                           boost::lambda::var( bcBasePtr ) = &boost::lambda::_1 ) );

    return bcBasePtr;
}



void
BCHandler::M_sumOffsets()
{
    for ( bcBaseIterator_Type it = M_bcList.begin(); it != M_bcList.end(); ++it )
        it->setOffset(it->offset()+M_offset);
}

} // namespace LifeV
