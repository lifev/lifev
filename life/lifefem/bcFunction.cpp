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
    @brief File contains BCManageNormal class for handling normal essential boundary conditions

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
// BCFunctionRobin
//==================================================

//==================================================
// Constructors
//==================================================

BCFunctionRobin::BCFunctionRobin( const BCFunctionRobin& bcFunctionRobin )
        :
        BCFunctionBase( bcFunctionRobin ),
        M_robinBoundaryMassCoeffFunction( bcFunctionRobin.M_robinBoundaryMassCoeffFunction )
{
}

BCFunctionRobin::BCFunctionRobin( const function_Type& rightHandSideFunction, const function_Type& massCoeffFunction )
        :
        BCFunctionBase( rightHandSideFunction ),
        M_robinBoundaryMassCoeffFunction( massCoeffFunction )

{
}


//==================================================
// Assignment Operator
//==================================================


BCFunctionRobin&
BCFunctionRobin::operator=( const BCFunctionRobin& bcFunctionRobin )
{
    if ( this != &bcFunctionRobin )
    {
        this->BCFunctionBase::operator=( bcFunctionRobin );
        M_robinBoundaryMassCoeffFunction = bcFunctionRobin.M_robinBoundaryMassCoeffFunction;
    }
    return *this;
}


//==================================================
// Set Methods
//==================================================

void
BCFunctionRobin::setFunctions_Robin( const function_Type& rightHandSideFunction, const function_Type& massCoeffFunction )
{
    setFunction( rightHandSideFunction );
    M_robinBoundaryMassCoeffFunction = massCoeffFunction;
}



BCFunctionBase*
createBCFunctionRobin( BCFunctionBase const* __bc )
{
    return new BCFunctionRobin( ( BCFunctionRobin const& )*__bc );
}
// register BCFunctionRobin in factory for cloning
const bool __bcRobin = FactoryCloneBCFunction::instance().registerProduct( typeid(BCFunctionRobin), &createBCFunctionRobin );


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
// BCFunctionUDepRobin
//==================================================


//==================================================
// Constructor
//==================================================

BCFunctionUDepRobin::BCFunctionUDepRobin( const BCFunctionUDepRobin& bcFunctionUDepRobin )
        :
        BCFunctionUDepBase( bcFunctionUDepRobin ),
        M_robinBoundaryMassCoeffFunction( bcFunctionUDepRobin.M_robinBoundaryMassCoeffFunction )
{
}

BCFunctionUDepRobin::BCFunctionUDepRobin( const function_Type& rightHandSideFunction, const function_Type& massCoeffFunction )
        :
        BCFunctionUDepBase( rightHandSideFunction ),
        M_robinBoundaryMassCoeffFunction( massCoeffFunction )

{
}


//==================================================
// Operators
//==================================================

BCFunctionUDepRobin&
BCFunctionUDepRobin::operator=( const BCFunctionUDepRobin& bcFunctionUDepRobin )
{
    if ( this != &bcFunctionUDepRobin )
    {
        this->BCFunctionUDepBase::operator=( bcFunctionUDepRobin );
        M_robinBoundaryMassCoeffFunction = bcFunctionUDepRobin.M_robinBoundaryMassCoeffFunction;
    }
    return *this;
}


//==================================================
// Set Methods
//==================================================

void
BCFunctionUDepRobin::setFunctions_Robin( const function_Type& rightHandSideFunction, const function_Type& massCoeffFunction )
{
    setFunction( rightHandSideFunction );
    M_robinBoundaryMassCoeffFunction = massCoeffFunction;
}


BCFunctionUDepBase*
createBCFunctionUDepRobin( BCFunctionUDepBase const* bcFunctionUDepRobin )
{
    return new BCFunctionUDepRobin( ( BCFunctionUDepRobin const& )* bcFunctionUDepRobin );
}
const bool bcFunctionUDepRobin = FactoryCloneBCFunctionUDep::instance().registerProduct(
                                     typeid(BCFunctionUDepRobin), &createBCFunctionUDepRobin );


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
// register BCFunctionRobin in factory for cloning
const bool bcFunctionDirectional = FactoryCloneBCFunction::instance().registerProduct( typeid(BCFunctionDirectional), &createBCFunctionDirectional );

} //End of namespace LifeV
