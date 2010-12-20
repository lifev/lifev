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
 *  @file
 *  @brief DataElasticStructure - File containing a data container for solid problems with elastic structure
 *
 *  @version 1.0
 *  @date 01-10-2003
 *  @author M.A. Fernandez
 *
 *  @version 1.18
 *  @date 10-06-2010
 *  @author Cristiano Malossi
 *
 *  @version 1.19
 *  @author Gilles Fourestey
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#include <life/lifecore/life.hpp>
#include <life/lifesolver/dataElasticStructure.hpp>

//=====================================================
// Constructors
//=====================================================

namespace LifeV
{
VenantKirchhoffElasticData::VenantKirchhoffElasticData():
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

VenantKirchhoffElasticData::VenantKirchhoffElasticData( const VenantKirchhoffElasticData& VenantKirchhoffElasticData ):
        DataTime                           ( VenantKirchhoffElasticData ),
        M_time                             ( VenantKirchhoffElasticData.M_time ),
        M_density                          ( VenantKirchhoffElasticData.M_density ),
        M_thickness                        ( VenantKirchhoffElasticData.M_thickness ),
        M_poisson                          ( VenantKirchhoffElasticData.M_poisson ),
        M_young                            ( VenantKirchhoffElasticData.M_young ),
        M_order                            ( VenantKirchhoffElasticData.M_order ),
        M_factor                           ( VenantKirchhoffElasticData.M_factor ),
        M_verbose                          ( VenantKirchhoffElasticData.M_verbose )
{
}

// ===================================================
// Operators
// ===================================================
VenantKirchhoffElasticData&
VenantKirchhoffElasticData::operator=( const VenantKirchhoffElasticData& VenantKirchhoffElasticData )
{
    if ( this != &VenantKirchhoffElasticData )
    {
        M_time                             = VenantKirchhoffElasticData.M_time;
        M_density                          = VenantKirchhoffElasticData.M_density;
        M_thickness                        = VenantKirchhoffElasticData.M_thickness;
        M_poisson                          = VenantKirchhoffElasticData.M_poisson;
        M_young                            = VenantKirchhoffElasticData.M_young;
        M_order                            = VenantKirchhoffElasticData.M_order;
        M_factor                           = VenantKirchhoffElasticData.M_factor;
        M_verbose                          = VenantKirchhoffElasticData.M_verbose;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
VenantKirchhoffElasticData::setup( const GetPot& dataFile, const std::string& section )
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
VenantKirchhoffElasticData::showMe( std::ostream& output ) const
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
        output << "Lame - lambda[" << i->first << "]                 = " << getLambda( i->first ) << std::endl;
        output << "Lame - mu[" << i->first << "]                     = " << getMu( i->first ) << std::endl;
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
VenantKirchhoffElasticData::setDataTime( const TimePtr_Type DataTime )
{
    M_time = DataTime;
}

void
VenantKirchhoffElasticData::setDensity( const Real& density )
{
    M_density = density;
}

void
VenantKirchhoffElasticData::setThickness( const Real& thickness )
{
    M_thickness = thickness;
}

void
VenantKirchhoffElasticData::setPoisson( const Real& poisson, const UInt& material )
{
    M_poisson[material] = poisson;
}

void
VenantKirchhoffElasticData::setYoung( const Real& young, const UInt& material )
{
    M_young[material] = young;
}

// ===================================================
// Get Method
// ===================================================
const VenantKirchhoffElasticData::TimePtr_Type
VenantKirchhoffElasticData::getDataTime() const
{
    return M_time;
}

const Real&
VenantKirchhoffElasticData::getRho() const
{
    return M_density;
}

const Real&
VenantKirchhoffElasticData::getThickness() const
{
    return M_thickness;
}

const Real&
VenantKirchhoffElasticData::getPoisson( const UInt& material ) const
{
    MaterialContainer_Type::const_iterator IT = M_poisson.find( material );
    if (IT != M_poisson.end())
        return M_poisson.find( material )->second;
    else
    {
        //WARNING("the Poisson modulus has not been set");
        return 0;
    }
}

const Real&
VenantKirchhoffElasticData::getYoung( const UInt& material ) const
{
    MaterialContainer_Type::const_iterator IT = M_young.find( material );
    if (IT != M_young.end())
        return IT->second;
    else
    {
        //WARNING("the Young modulus has not been set");
        return 0;
    }
}

Real
VenantKirchhoffElasticData::getLambda( const UInt& material ) const
{
    return M_young.find( material )->second * M_poisson.find( material )->second /
           ( ( 1.0 + M_poisson.find( material )->second ) * ( 1.0 - 2.0 * M_poisson.find( material )->second ) );
}

Real
VenantKirchhoffElasticData::getMu( const UInt& material ) const
{
    return M_young.find( material )->second/( 2.0 * ( 1.0 + M_poisson.find( material )->second ) );
}

const std::string&
VenantKirchhoffElasticData::getOrder() const
{
    return M_order;
}

const Real&
VenantKirchhoffElasticData::getFactor() const
{
    return M_factor;
}

const UInt&
VenantKirchhoffElasticData::getVerbose() const
{
    return M_verbose;
}

const std::string&
VenantKirchhoffElasticData::getSolidType()
{
    return M_solidType;
}

const bool&
VenantKirchhoffElasticData::getUseExactJacobian() const
{
    return M_useExactJacobian;
}

}
