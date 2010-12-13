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
    @brief File contains classes to holds the FE vectors used for prescribing boundary conditions

	@author Miguel Fernandez <miguel.fernandez@inria.fr>
    @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
    @author Vincent Martin <vincent.martin@inria.fr>
    @contributor Mauro Perego <perego.mauro@gmail.com>
    @maintainer Mauro Perego <perego.mauro@gmail.com>

    */


#include <life/lifefem/bcVector.hpp>
namespace LifeV
{
// ===================================
// Implementation for BCVectorBase
// ===================================



// ===================================
// Constructors
//====================================


BCVectorBase::BCVectorBase()
        :
        M_robinBoundaryMassCoeff( 0.0 ),
        M_resistanceCoeff( 0.0 ),
        M_betaCoeff(1.0),
        M_gammaCoeff(1.0),
        M_isRobinBdMassCoeffAVector( false ),
        M_isBetaCoeffAVector( false  ),
        M_isGammaCoeffAVector( false ),
        M_type( 0 ),
        M_finalized( false )
{}


BCVectorBase::BCVectorBase( const EpetraVector& rightHandSideVector , const UInt numberOfTotalDof, UInt type )
        :
        M_rightHandSideVectorPtr       ( &rightHandSideVector  ),
        M_numberOfTotalDof( numberOfTotalDof ),
        M_robinBoundaryMassCoeff ( 0.0 ),
        M_resistanceCoeff( 0.0 ),
        M_betaCoeff(1.0),
        M_gammaCoeff(1.0),
        M_isRobinBdMassCoeffAVector( false ),
        M_isBetaCoeffAVector( false  ),
        M_isGammaCoeffAVector( false ),
        M_type      ( type ),
        M_finalized ( false )
{}

BCVectorBase::BCVectorBase( BCVectorBase const& bcVectorBase ):
        M_rightHandSideVectorPtr ( bcVectorBase.M_rightHandSideVectorPtr ),
        M_numberOfTotalDof ( bcVectorBase.M_numberOfTotalDof ),
        M_robinBoundaryMassCoeff ( bcVectorBase.M_robinBoundaryMassCoeff ),
        M_resistanceCoeff ( bcVectorBase.M_resistanceCoeff ),
        M_betaCoeff ( bcVectorBase.M_betaCoeff ),
        M_gammaCoeff ( bcVectorBase.M_gammaCoeff ),
        M_isRobinBdMassCoeffAVector ( bcVectorBase.M_isRobinBdMassCoeffAVector ),
        M_isBetaCoeffAVector ( bcVectorBase.M_isBetaCoeffAVector ),
        M_isGammaCoeffAVector ( bcVectorBase.M_isGammaCoeffAVector ),
        M_type ( bcVectorBase.M_type ),
        M_finalized ( bcVectorBase.M_finalized )
{

}

// ===================================
// Operators
//====================================


BCVectorBase&
BCVectorBase::operator=( BCVectorBase const& bcVectorBase )
{
    if ( this != &bcVectorBase )
    {
        M_rightHandSideVectorPtr        = bcVectorBase.M_rightHandSideVectorPtr;
        M_numberOfTotalDof = bcVectorBase.M_numberOfTotalDof;
        M_robinBoundaryMassCoeff  = bcVectorBase.M_robinBoundaryMassCoeff;
        M_resistanceCoeff  = bcVectorBase.M_resistanceCoeff;
        M_betaCoeff   = bcVectorBase.M_betaCoeff;
        M_gammaCoeff  = bcVectorBase.M_gammaCoeff;
        M_type       = bcVectorBase.M_type;
        M_isRobinBdMassCoeffAVector = bcVectorBase.M_isRobinBdMassCoeffAVector;
        M_isBetaCoeffAVector  = bcVectorBase.M_isBetaCoeffAVector;
        M_isGammaCoeffAVector = bcVectorBase.M_isGammaCoeffAVector;
        M_finalized  = bcVectorBase.M_finalized;
    }
    return *this;
}


Real
BCVectorBase::operator() ( const ID& globalDofId, const ID& component ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    return ( *M_rightHandSideVectorPtr ) ( ( component - 1 ) * M_numberOfTotalDof + globalDofId );
}


// ===================================
// Methods
//====================================

Real
BCVectorBase::MixteVec ( const ID& globalDofId, const ID& component ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    return ( *M_rightHandSideVectorPtr ) ( ( component - 1 ) * M_numberOfTotalDof + globalDofId );
}


Real
BCVectorBase::BetaVec ( const ID& globalDofId, const ID& component ) const
{


    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );

    return ( *M_betaCoeffVectorPtr ) ( ( component - 1 ) * M_numberOfTotalDof + globalDofId );
}

Real
BCVectorBase::GammaVec ( const ID& globalDofId, const ID& component ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    return ( *M_gammaCoeffVectorPtr ) ( ( component - 1 ) * M_numberOfTotalDof + globalDofId );
}


// ===================================
// Set Methods
//====================================

void
BCVectorBase::setVector( const EpetraVector& rightHandSideVector , UInt numberOfTotalDof, UInt type )
{
    M_rightHandSideVectorPtr = &rightHandSideVector  ;
    M_numberOfTotalDof = numberOfTotalDof;
    M_type = type;
    M_finalized = true;
}

void
BCVectorBase::setMixteVec( const EpetraVector& robinBoundaryMassCoeffVector )
{
    M_isRobinBdMassCoeffAVector = true;
    M_robinBoundaryMassCoeffVectorPtr= &robinBoundaryMassCoeffVector;
}

void
BCVectorBase::setBetaVec( const EpetraVector& betaCoeffVector )
{
    M_isBetaCoeffAVector = true;
    M_betaCoeffVectorPtr= &betaCoeffVector;
}

void
BCVectorBase::setGammaVec( const EpetraVector& gammaCoeffVector )
{
    M_isGammaCoeffAVector = true;
    M_gammaCoeffVectorPtr= &gammaCoeffVector;
}




// ===================================
// Implementation for BCVector
// ===================================



// ===================================
// Constructors
//====================================


BCVector::BCVector( EpetraVector& rightHandSideVector, UInt const numberOfTotalDof, UInt type )
        :
        BCVectorBase( rightHandSideVector, numberOfTotalDof, type )
{
    this->setFinalized( true );
}


BCVector::BCVector( const BCVector& bcVector ):
        BCVectorBase( bcVector ) {}


// ===================================
// Operators
//====================================


//! Assignment operator for BCVector
BCVector&
BCVector::operator=( const BCVector& bcVector )
{
    if ( this != &bcVector )
        bcVectorBase_Type::operator=( bcVector );

    return *this;
}

// ===================================
// Methods
//====================================


std::ostream&
BCVector::showMe( bool /* verbose */, std::ostream & out ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    out << "+++++++++++++++++++++++++++++++" << std::endl;
    out << "BC Vector Interface: " << std::endl;
    out << "number of interface vector Dof : " << this->nbTotalDOF() << std::endl;
    out << "==>Interface Dof :\n";
    out << "+++++++++++++++++++++++++++++++" << std::endl;
    return out;
}





BCVectorBase*
createBCVector( BCVectorBase const* __bc )
{
    return new BCVector( ( BCVector const& )*__bc );
}
// register BCFunctionBase in factory for cloning
const bool __bcvec = FactoryCloneBCVector::instance().registerProduct( typeid(BCVector), &createBCVector );


// ====================================
// Implementation for BCVectorInterface
// ====================================


// ===================================
// Constructors
//====================================


BCVectorInterface::BCVectorInterface( const EpetraVector& rightHandSideVector, UInt numberOfTotalDof,
                                      const dofInterfacePtr_Type& interfaceDofPtr, UInt type )
        :
        BCVectorBase( rightHandSideVector, numberOfTotalDof, type ),
        M_interfaceDofPtr( interfaceDofPtr )
{
    this->setFinalized( true );
}


BCVectorInterface::BCVectorInterface( const BCVectorInterface & bcVectorInterface ):
        BCVectorBase( bcVectorInterface ),
        M_interfaceDofPtr( bcVectorInterface.M_interfaceDofPtr ) {}


// ===================================
// Operators
//====================================

BCVectorInterface&
BCVectorInterface::operator=( const BCVectorInterface & bcVectorInterface )
{
    if ( this != &bcVectorInterface )
    {
        bcVectorBase_Type::operator=( bcVectorInterface );
        M_interfaceDofPtr = bcVectorInterface.M_interfaceDofPtr;
    }

    return *this;
}


Real
BCVectorInterface::operator() ( const ID& globalDofId, const ID& component ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    return ( *M_rightHandSideVectorPtr ) (( component - 1 ) * M_numberOfTotalDof + M_interfaceDofPtr->getInterfaceDof( globalDofId ));
}

// ===================================
// methods
//====================================

void BCVectorInterface::setup( const EpetraVector& rightHandSideVector, UInt numberOfTotalDof, const dofInterfacePtr_Type& interfaceDofPtr, UInt type )
{
    M_rightHandSideVectorPtr        = &rightHandSideVector;
    M_numberOfTotalDof = numberOfTotalDof;
    M_robinBoundaryMassCoeff  = 0.0;
    M_resistanceCoeff  = 0.0;
    M_betaCoeff   = 1.0;
    M_gammaCoeff  = 1.0;
    M_type       = type;
    M_interfaceDofPtr      = interfaceDofPtr;
    M_isRobinBdMassCoeffAVector = false;
    M_isBetaCoeffAVector  = false;
    M_isGammaCoeffAVector = false;
    this->setFinalized( true );
}


void
BCVectorInterface::setVector( const EpetraVector& rightHandSideVector, UInt numberOfTotalDof, const dofInterfacePtr_Type& interfaceDofPtr, UInt type )
{
    ASSERT_PRE( !this->isFinalized(), "BC Vector cannot be set twice." );

    bcVectorBase_Type::setVector( rightHandSideVector, numberOfTotalDof, type );

    M_interfaceDofPtr = interfaceDofPtr;

}


Real
BCVectorInterface::MixteVec( const ID& globalDofId, const ID& component ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    return ( *M_robinBoundaryMassCoeffVectorPtr ) (( component - 1 ) * M_numberOfTotalDof + M_interfaceDofPtr->getInterfaceDof( globalDofId ));
}


Real
BCVectorInterface::BetaVec( const ID& globalDofId, const ID& component ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    return ( *M_betaCoeffVectorPtr ) (( component - 1 ) * M_numberOfTotalDof + M_interfaceDofPtr->getInterfaceDof( globalDofId ));
}

Real
BCVectorInterface::GammaVec( const ID& globalDofId, const ID& component ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    return ( *M_gammaCoeffVectorPtr ) (( component - 1 ) * M_numberOfTotalDof + M_interfaceDofPtr->getInterfaceDof( globalDofId ));
}


std::ostream&
BCVectorInterface::showMe( bool verbose, std::ostream & out ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    out << "+++++++++++++++++++++++++++++++" << std::endl;
    out << "BC Vector Interface: " << std::endl;
    out << "number of interface vector Dof : " << M_numberOfTotalDof << std::endl;
    out << "==>Interface Dof :\n";
    M_interfaceDofPtr->showMe( verbose, out );  // no showMe(..) in Miguel's DofInterface
    out << "+++++++++++++++++++++++++++++++" << std::endl;
    return out;
}

BCVectorBase*
createBCVectorInterface( BCVectorBase const* bcVectorBase )
{
    return new BCVectorInterface( ( BCVectorInterface const& )* bcVectorBase );
}

// register BCFunctionBase in factory for cloning
const bool bcVectorBase = FactoryCloneBCVector::instance().registerProduct( typeid(BCVectorInterface), &createBCVectorInterface );

}
