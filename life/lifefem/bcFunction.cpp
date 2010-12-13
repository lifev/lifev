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
    @brief File contains BCNormalManager class for handling normal essential boundary conditions

	@author Miguel Fernandez <miguel.fernandez@inria.fr>
    @contributor Christophe Prud'homme <christophe.prudhomme@epfl.ch>
    @contributor Mauro Perego <perego.mauro@gmail.com>
    @maintainer Mauro Perego <perego.mauro@gmail.com>

    @date 10-12-2004
 *///@HEADER


#include <life/lifecore/life.hpp>
#include <life/lifefem/bcFunction.hpp>

namespace LifeV
{

//==================================================
// BCFunctionBase
//==================================================


//==================================================
// Constructors
//==================================================


BCFunctionBase::BCFunctionBase( function_Type userDefinedFunction )
        :
        M_userDefinedFunction( userDefinedFunction )
{
}

BCFunctionBase::BCFunctionBase( const BCFunctionBase& bcFunctionBase )
        :
        M_userDefinedFunction( bcFunctionBase.M_userDefinedFunction )
{
}


//==================================================
// Operators
//==================================================


BCFunctionBase&
BCFunctionBase::operator= ( const BCFunctionBase& bcFunctionBase )
{
    if (this != &bcFunctionBase)
    {
        M_userDefinedFunction = bcFunctionBase.M_userDefinedFunction;
    }
    return *this;
}


BCFunctionBase*
createBCFunctionBase( BCFunctionBase const* bcFunctionBase )
{
    return new BCFunctionBase( ( BCFunctionBase const& )* bcFunctionBase );
}
// register BCFunctionBase in factory for cloning
const bool bcBaseFactory = FactoryCloneBCFunction::instance().registerProduct( typeid(BCFunctionBase), &createBCFunctionBase );




//==================================================
// BCFunctionMixte
//==================================================

//==================================================
// Constructors
//==================================================

BCFunctionMixte::BCFunctionMixte( const BCFunctionMixte& bcFunctionMixte )
        :
        BCFunctionBase( bcFunctionMixte ),
        M_robinBoundaryMassCoeffFunction( bcFunctionMixte.M_robinBoundaryMassCoeffFunction )
{
}

BCFunctionMixte::BCFunctionMixte( const function_Type& rightHandSideFunction, const function_Type& massCoeffFunction )
        :
        BCFunctionBase( rightHandSideFunction ),
        M_robinBoundaryMassCoeffFunction( massCoeffFunction )

{
}


//==================================================
// Assignment Operator
//==================================================


BCFunctionMixte&
BCFunctionMixte::operator=( const BCFunctionMixte& bcFunctionMixte )
{
    if ( this != &bcFunctionMixte )
    {
        this->BCFunctionBase::operator=( bcFunctionMixte );
        M_robinBoundaryMassCoeffFunction = bcFunctionMixte.M_robinBoundaryMassCoeffFunction;
    }
    return *this;
}


//==================================================
// Set Methods
//==================================================

void
BCFunctionMixte::setFunctions_Mixte( const function_Type& rightHandSideFunction, const function_Type& massCoeffFunction )
{
    setFunction( rightHandSideFunction );
    M_robinBoundaryMassCoeffFunction = massCoeffFunction;
}



BCFunctionBase*
createBCFunctionMixte( BCFunctionBase const* __bc )
{
    return new BCFunctionMixte( ( BCFunctionMixte const& )*__bc );
}
// register BCFunctionMixte in factory for cloning
const bool __bcmixte = FactoryCloneBCFunction::instance().registerProduct( typeid(BCFunctionMixte), &createBCFunctionMixte );


//==================================================
// BCFunctionUDepBase
//==================================================

//==================================================
// Constructors
//==================================================


BCFunctionUDepBase::BCFunctionUDepBase(const function_Type& userDefinedFunction ): M_userDefinedFunction(userDefinedFunction)
{
}


BCFunctionUDepBase::BCFunctionUDepBase(const BCFunctionUDepBase& bcFunctionUDepBase ):
        M_userDefinedFunction(bcFunctionUDepBase.M_userDefinedFunction)
{
}


//==================================================
// Operators
//==================================================

BCFunctionUDepBase&
BCFunctionUDepBase::operator= ( const BCFunctionUDepBase& bcFunctionUDepBase )
{
    if (this != &bcFunctionUDepBase)
    {
        M_userDefinedFunction = bcFunctionUDepBase.M_userDefinedFunction;
    }
    return *this;
}


BCFunctionUDepBase*
createBCFunctionUDep( BCFunctionUDepBase const* bcFunctionUDepBase )
{
    return new BCFunctionUDepBase( ( BCFunctionUDepBase const& )* bcFunctionUDepBase );
}
const bool bcFunctionUDepBase = FactoryCloneBCFunctionUDep::instance().registerProduct(
                                    typeid(BCFunctionUDepBase), &createBCFunctionUDep );


//==================================================
// BCFunctionUDepMixte
//==================================================


//==================================================
// Constructor
//==================================================

BCFunctionUDepMixte::BCFunctionUDepMixte( const BCFunctionUDepMixte& bcFunctionUDepMixte )
        :
        BCFunctionUDepBase( bcFunctionUDepMixte ),
        M_robinBoundaryMassCoeffFunction( bcFunctionUDepMixte.M_robinBoundaryMassCoeffFunction )
{
}

BCFunctionUDepMixte::BCFunctionUDepMixte( const function_Type& rightHandSideFunction, const function_Type& massCoeffFunction )
        :
        BCFunctionUDepBase( rightHandSideFunction ),
        M_robinBoundaryMassCoeffFunction( massCoeffFunction )

{
}


//==================================================
// Operators
//==================================================

BCFunctionUDepMixte&
BCFunctionUDepMixte::operator=( const BCFunctionUDepMixte& bcFunctionUDepMixte )
{
    if ( this != &bcFunctionUDepMixte )
    {
        this->BCFunctionUDepBase::operator=( bcFunctionUDepMixte );
        M_robinBoundaryMassCoeffFunction = bcFunctionUDepMixte.M_robinBoundaryMassCoeffFunction;
    }
    return *this;
}


//==================================================
// Set Methods
//==================================================

void
BCFunctionUDepMixte::setFunctions_Mixte( const function_Type& rightHandSideFunction, const function_Type& massCoeffFunction )
{
    setFunction( rightHandSideFunction );
    M_robinBoundaryMassCoeffFunction = massCoeffFunction;
}


BCFunctionUDepBase*
createBCFunctionUDepMixte( BCFunctionUDepBase const* bcFunctionUDepMixte )
{
    return new BCFunctionUDepMixte( ( BCFunctionUDepMixte const& )* bcFunctionUDepMixte );
}
const bool bcFunctionUDepMixte = FactoryCloneBCFunctionUDep::instance().registerProduct(
                                     typeid(BCFunctionUDepMixte), &createBCFunctionUDepMixte );


//==================================================
// BCFunctionDirectional
//==================================================

//==================================================
// Constructors
//==================================================
BCFunctionDirectional::BCFunctionDirectional( const BCFunctionDirectional& bcFunctionDirectional )
        :
        BCFunctionBase( bcFunctionDirectional ),
        M_userDefinedVersorsFunction( bcFunctionDirectional.M_userDefinedVersorsFunction )
{
}

BCFunctionDirectional::BCFunctionDirectional( const function_Type& userDefinedFunction, const function_Type& userDefinedVersorsFunction )
        :
        BCFunctionBase( userDefinedFunction ),
        M_userDefinedVersorsFunction( userDefinedVersorsFunction )

{
}

//==================================================
// Operators
//==================================================

BCFunctionDirectional&
BCFunctionDirectional::operator=( const BCFunctionDirectional& bcFunctionDirectional )
{
    if ( this != &bcFunctionDirectional )
    {
        this->BCFunctionBase::operator=( bcFunctionDirectional );
        M_userDefinedVersorsFunction = bcFunctionDirectional.M_userDefinedVersorsFunction;
    }
    return *this;
}


//==================================================
// Set Method
//==================================================
void
BCFunctionDirectional::setFunctions_Directional( const function_Type& userDefinedFunction, const function_Type& userDefinedVersorsFunction )
{
    setFunction( userDefinedFunction );
    M_userDefinedVersorsFunction = userDefinedVersorsFunction;
}


BCFunctionBase*
createBCFunctionDirectional( BCFunctionBase const* bcFunctionDirectional )
{
    return new BCFunctionDirectional( ( BCFunctionDirectional const& )* bcFunctionDirectional);
}
// register BCFunctionMixte in factory for cloning
const bool bcFunctionDirectional = FactoryCloneBCFunction::instance().registerProduct( typeid(BCFunctionDirectional), &createBCFunctionDirectional );

} //End of namespace LifeV
