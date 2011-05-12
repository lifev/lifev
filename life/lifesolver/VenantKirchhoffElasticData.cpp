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

#include <life/lifesolver/VenantKirchhoffElasticData.hpp>

namespace LifeV
{

//=====================================================
// Constructors
//=====================================================
VenantKirchhoffElasticData::VenantKirchhoffElasticData():
        M_time                             ( ),
        M_density                          ( ),
        M_thickness                        ( ),
        M_materialsFlagSet                 ( false ),
        M_poisson                          ( ),
        M_young                            ( ),
        M_order                            ( ),
        M_factor                           ( ),
        M_verbose                          ( )
{
}

VenantKirchhoffElasticData::VenantKirchhoffElasticData( const VenantKirchhoffElasticData& venantKirchhoffElasticData ):
        TimeData                           ( venantKirchhoffElasticData ),
        M_time                             ( venantKirchhoffElasticData.M_time ),
        M_density                          ( venantKirchhoffElasticData.M_density ),
        M_thickness                        ( venantKirchhoffElasticData.M_thickness ),
        M_materialsFlagSet                 ( venantKirchhoffElasticData.M_materialsFlagSet ),
        M_poisson                          ( venantKirchhoffElasticData.M_poisson ),
        M_young                            ( venantKirchhoffElasticData.M_young ),
        M_order                            ( venantKirchhoffElasticData.M_order ),
        M_factor                           ( venantKirchhoffElasticData.M_factor ),
        M_verbose                          ( venantKirchhoffElasticData.M_verbose )
{
}

// ===================================================
// Operators
// ===================================================
VenantKirchhoffElasticData&
VenantKirchhoffElasticData::operator=( const VenantKirchhoffElasticData& venantKirchhoffElasticData )
{
    if ( this != &venantKirchhoffElasticData )
    {
        M_time                             = venantKirchhoffElasticData.M_time;
        M_density                          = venantKirchhoffElasticData.M_density;
        M_thickness                        = venantKirchhoffElasticData.M_thickness;
        M_materialsFlagSet                 = venantKirchhoffElasticData.M_materialsFlagSet;
        M_poisson                          = venantKirchhoffElasticData.M_poisson;
        M_young                            = venantKirchhoffElasticData.M_young;
        M_order                            = venantKirchhoffElasticData.M_order;
        M_factor                           = venantKirchhoffElasticData.M_factor;
        M_verbose                          = venantKirchhoffElasticData.M_verbose;
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
        M_time.reset( new time_Type( dataFile, section + "/time_discretization" ) );

    // physics
    M_solidType = dataFile( ( section + "/physics/solidType" ).data(), "linearVenantKirchhoff" );
    M_density   = dataFile( ( section + "/physics/density"   ).data(), 1. );
    M_thickness = dataFile( ( section + "/physics/thickness" ).data(), 0.1 );

    UInt materialsNumber = dataFile.vector_variable_size( ( section + "/physics/material_flag" ).data() );
    if ( materialsNumber == 0 )
    {
        // If no material is specified in the data file the code assume that there is just one material
        // and by default it is memorized with ID 1. Getters and Setters have been designed to deal with thic choice.
        M_materialsFlagSet = false;

        M_young[1]   = dataFile( ( section + "/physics/young"   ).data(), 0. );
        M_poisson[1] = dataFile( ( section + "/physics/poisson" ).data(), 0. );
    }
    else
    {
        M_materialsFlagSet = true;

        // These asserts are commented because in some cases we need to initialize the materials with default Young and Poisson, setting the correct values a posteriori.
        //ASSERT( M_materialsFlagSet == dataFile.vector_variable_size( ( section + "/physics/young"   ).data()),  "!!! ERROR: Inconsistent size for Young Modulus !!!");
        //ASSERT( M_materialsFlagSet == dataFile.vector_variable_size( ( section + "/physics/poisson" ).data() ), "!!! ERROR: Inconsistent size for Poisson Coef. !!!");

        UInt material(0);
        for ( UInt i(0) ; i < materialsNumber ; ++i )
        {
            material            = dataFile( ( section + "/physics/material_flag" ).data(), 0, i );
            M_young[material]   = dataFile( ( section + "/physics/young"         ).data(), 0., i );
            M_poisson[material] = dataFile( ( section + "/physics/poisson"       ).data(), 0., i );
        }
    }

    // space_discretization
    M_order            = dataFile( ( section + "/space_discretization/order" ).data(), "P1" );

    // miscellaneous
    M_factor           = dataFile( ( section + "/miscellaneous/factor"  ).data(), 1.0 );
    M_verbose          = dataFile( ( section + "/miscellaneous/verbose" ).data(), 1 );
    M_useExactJacobian = dataFile( ( section + "/useExactJacobian"      ).data(), false );
}

void
VenantKirchhoffElasticData::showMe( std::ostream& output ) const
{
    // physics
    output << "\n*** Values for data [solid/physics]\n\n";
    output << "density                          = " << M_density << std::endl;
    output << "thickness                        = " << M_thickness << std::endl;
    for ( materialContainerIterator_Type i = M_young.begin() ; i != M_young.end() ; ++i )
        output << "young[" << i->first << "]                         = " << i->second << std::endl;
    for ( materialContainerIterator_Type i = M_poisson.begin() ; i != M_poisson.end() ; ++i )
        output << "poisson[" << i->first << "]                       = " << i->second << std::endl;

    for ( materialContainerIterator_Type i = M_poisson.begin() ; i != M_poisson.end() ; ++i )
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
// Get Method
// ===================================================
Real
VenantKirchhoffElasticData::poisson( const UInt& material ) const
{
    materialContainer_Type::const_iterator IT;

    if ( M_materialsFlagSet )
        IT = M_poisson.find( material );
    else
        IT = M_poisson.find( 1 );

    if ( IT != M_poisson.end() )
        return M_poisson.find( material )->second;
    else
    {
        std::cout << " !!! Warning: the Poisson modulus has not been set !!!" << std::endl;
        return 0;
    }
}

Real
VenantKirchhoffElasticData::young( const UInt& material ) const
{
    materialContainer_Type::const_iterator IT;
    if ( M_materialsFlagSet )
        IT = M_young.find( material );
    else
        IT = M_young.find( 1 );

    if ( IT != M_young.end() )
        return IT->second;
    else
    {
        std::cout << " !!! Warning: the Young modulus has not been set !!!" << std::endl;
        return 0;
    }
}

Real
VenantKirchhoffElasticData::lambda( const UInt& material ) const
{
    Real youngC, poissonC;

    youngC   = this->young( material );
    poissonC = this->poisson( material );

    return youngC * poissonC / ( ( 1.0 + poissonC ) * ( 1.0 - 2.0 * poissonC ) );
}

Real
VenantKirchhoffElasticData::mu( const UInt& material ) const
{
    Real youngC, poissonC;

    youngC   = this->young( material );
    poissonC = this->poisson( material );

    return youngC / ( 2.0 * ( 1.0 + poissonC ) );
}

}
