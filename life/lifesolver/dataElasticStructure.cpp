//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano

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

/*!
 *  @file
 *  @brief DataElasticStructure - File containing a data container for solid problems with elastic structure
 *
 *  @version 1.0
 *  @author M.A. Fernandez
 *  @date 01-10-2003
 *
 *  @version 1.18
 *  @author Cristiano Malossi
 *  @date 10-06-2010
 *
 *  @version 1.19
 *  @author Gilles Fourestey
 */

#include <life/lifesolver/dataElasticStructure.hpp>

namespace LifeV {

DataElasticStructure::DataElasticStructure() :
        M_time                             ( ),
        M_density                          ( ),
        M_thickness                        ( ),
        M_poisson                          ( ),
        M_young                            ( ),
        M_order                            ( ),
        M_factor                           ( ),
        M_verbose                          ( )
{
}

DataElasticStructure::DataElasticStructure( const DataElasticStructure& dataElasticStructure ):
    DataTime                           ( dataElasticStructure ),
    M_time                             ( dataElasticStructure.M_time ),
    M_density                          ( dataElasticStructure.M_density ),
    M_thickness                        ( dataElasticStructure.M_thickness ),
    M_poisson                          ( dataElasticStructure.M_poisson ),
    M_young                            ( dataElasticStructure.M_young ),
    M_order                            ( dataElasticStructure.M_order ),
    M_factor                           ( dataElasticStructure.M_factor ),
    M_verbose                          ( dataElasticStructure.M_verbose )
{
}

// ===================================================
// Operators
// ===================================================
DataElasticStructure&
DataElasticStructure::operator=( const DataElasticStructure& dataElasticStructure )
{
    if ( this != &dataElasticStructure )
    {
        M_time                             = dataElasticStructure.M_time;
        M_density                          = dataElasticStructure.M_density;
        M_thickness                        = dataElasticStructure.M_thickness;
        M_poisson                          = dataElasticStructure.M_poisson;
        M_young                            = dataElasticStructure.M_young;
        M_order                            = dataElasticStructure.M_order;
        M_factor                           = dataElasticStructure.M_factor;
        M_verbose                          = dataElasticStructure.M_verbose;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
DataElasticStructure::setup( const GetPot& dataFile, const std::string& section )
{
    // If data time has not been set
    if ( !M_time.get() )
        M_time.reset( new Time_Type( dataFile, section + "/time_discretization" ) );

    // physics
    M_density   = dataFile( ( section + "/physics/density" ).data(), 1. );
    M_thickness = dataFile( ( section + "/physics/thickness" ).data(), 0.1 );
    M_useExactJacobian = dataFile( ( section + "/useExactJacobian" ).data(), false );

    UInt materialsNumber = dataFile.vector_variable_size( ( section + "/physics/material_flag" ).data() );
    if ( materialsNumber == 0 )
    {

        //WARNING("The material flag was not set from data file. Its value will be dedced from the first volume marker.");
//         M_young[1]   = dataFile( ( section + "/physics/young" ).data(), 0. );
//         M_poisson[1] = dataFile( ( section + "/physics/poisson" ).data(), 0. );
    }
    else
    {
        ASSERT( materialsNumber == dataFile.vector_variable_size( ( section + "/physics/young" ).data()),   "!!! ERROR: Inconsistent size for Young Modulus !!!");
        ASSERT( materialsNumber == dataFile.vector_variable_size( ( section + "/physics/poisson" ).data() ), "!!! ERROR: Inconsistent size for Poisson Coeff. !!!");

        UInt material(0);
        for ( UInt i(0) ; i < materialsNumber ; ++i )
        {
            material            = dataFile( ( section + "/physics/material_flag" ).data(), 0., i );
            M_young[material]   = dataFile( ( section + "/physics/young" ).data(), 0., i );
            M_poisson[material] = dataFile( ( section + "/physics/poisson" ).data(), 0., i );
        }
    }

    // space_discretization
    M_order     = dataFile( "solid/space_discretization/order", "P1" );

    // miscellaneous
    M_factor  = dataFile( "solid/miscellaneous/factor", 1.0 );
    M_verbose = dataFile( "solid/miscellaneous/verbose", 1 );
    M_solidType = dataFile( "solid/physics/solidType", "linearVenantKirchhof" );
}

void
DataElasticStructure::showMe( std::ostream& output ) const
{
    // physics
    output << "\n*** Values for data [solid/physics]\n\n";
    output << "density                          = " << M_density << std::endl;
    output << "thickness                        = " << M_thickness << std::endl;
    for ( MaterialContainer_ConstIterator i = M_young.begin() ; i != M_young.end() ; ++i )
        output << "young[" << i->first << "]                         = " << i->second << std::endl;
    for ( MaterialContainer_ConstIterator i = M_poisson.begin() ; i != M_poisson.end() ; ++i )
        output << "poisson[" << i->first << "]                       = " << i->second << std::endl;

    for ( MaterialContainer_ConstIterator i = M_poisson.begin() ; i != M_poisson.end() ; ++i )
    {
        output << "Lame - lambda[" << i->first << "]                 = " << lambda( i->first ) << std::endl;
        output << "Lame - mu[" << i->first << "]                     = " << mu( i->first ) << std::endl;
    }

    output << "\n*** Values for data [solid/miscellaneous]\n\n";
    output << "deformation factor               = " << M_factor << std::endl;
    output << "verbose                          = " << M_verbose << std::endl;

    output << "\n*** Values for data [solid/space_discretization]\n\n";
    output << "FE order                         = " << M_order << std::endl;

    output << "\n*** Values for data [solid/time_discretization]\n\n";
    M_time->showMe( output );
}

// ===================================================
// Set Method
// ===================================================
void
DataElasticStructure::setDataTime( const Time_ptrType DataTime )
{
    M_time = DataTime;
}

void
DataElasticStructure::setDensity( const Real& density )
{
    M_density = density;
}

void
DataElasticStructure::setThickness( const Real& thickness )
{
    M_thickness = thickness;
}

void
DataElasticStructure::setPoisson( const Real& poisson, const UInt& material )
{
    M_poisson[material] = poisson;
}

void
DataElasticStructure::setYoung( const Real& young, const UInt& material )
{
    M_young[material] = young;
}

// ===================================================
// Get Method
// ===================================================
DataElasticStructure::Time_ptrType
DataElasticStructure::dataTime() const
{
    return M_time;
}

const Real&
DataElasticStructure::rho() const
{
    return M_density;
}

const Real&
DataElasticStructure::thickness() const
{
    return M_thickness;
}

const Real&
DataElasticStructure::poisson( const UInt& material ) const
{
    MaterialContainer_Type::const_iterator IT = M_poisson.find( material );
    if(IT != M_poisson.end())
        return M_poisson.find( material )->second;
    else
    {
        //WARNING("the Poisson modulus has not been set");
        return 0;
    }
}

const Real&
DataElasticStructure::young( const UInt& material ) const
{
    MaterialContainer_Type::const_iterator IT = M_young.find( material );
    if(IT != M_young.end())
        return IT->second;
    else
    {
        //WARNING("the Young modulus has not been set");
        return 0;
    }
}

Real
DataElasticStructure::lambda( const UInt& material ) const
{
    return M_young.find( material )->second * M_poisson.find( material )->second /
           ( ( 1.0 + M_poisson.find( material )->second ) * ( 1.0 - 2.0 * M_poisson.find( material )->second ) );
}

Real
DataElasticStructure::mu( const UInt& material ) const
{
    return M_young.find( material )->second/( 2.0 * ( 1.0 + M_poisson.find( material )->second ) );
}

const std::string&
DataElasticStructure::order() const
{
    return M_order;
}

const Real&
DataElasticStructure::factor() const
{
    return M_factor;
}

const UInt&
DataElasticStructure::verbose() const
{
    return M_verbose;
}

const std::string&
DataElasticStructure::solidType()
{
    return M_solidType;
}

const bool&
DataElasticStructure::useExactJacobian() const
{
    return M_useExactJacobian;
}

}
