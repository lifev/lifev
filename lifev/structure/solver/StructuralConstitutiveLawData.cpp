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

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

namespace LifeV
{

//=====================================================
// Constructors
//=====================================================
StructuralConstitutiveLawData::StructuralConstitutiveLawData() :
    M_time                             ( ),
    M_timeAdvance                      ( ),
    M_density                          ( ),
    M_thickness                        ( ),
    M_externalPressure                 ( ),
    M_materialsFlagSet                 ( false ),
    M_poisson                          ( ),
    M_young                            ( ),
    M_bulk                             ( ),
    M_alpha                            ( ),
    M_gamma                            ( ),
    M_order                            ( ),
    M_verbose                          ( ),
    M_vectorMaterialFlags              ( )
{
}

StructuralConstitutiveLawData::StructuralConstitutiveLawData ( const StructuralConstitutiveLawData& structuralConstitutiveLawData ) :
    M_time                             ( structuralConstitutiveLawData.M_time ),
    M_timeAdvance                      ( structuralConstitutiveLawData.M_timeAdvance ),
    M_density                          ( structuralConstitutiveLawData.M_density ),
    M_thickness                        ( structuralConstitutiveLawData.M_thickness ),
    M_externalPressure                 ( structuralConstitutiveLawData.M_externalPressure ),
    M_materialsFlagSet                 ( structuralConstitutiveLawData.M_materialsFlagSet ),
    M_poisson                          ( structuralConstitutiveLawData.M_poisson ),
    M_young                            ( structuralConstitutiveLawData.M_young ),
    M_bulk                             ( structuralConstitutiveLawData.M_bulk ),
    M_alpha                            ( structuralConstitutiveLawData.M_alpha ),
    M_gamma                            ( structuralConstitutiveLawData.M_gamma ),
    M_order                            ( structuralConstitutiveLawData.M_order ),
    M_verbose                          ( structuralConstitutiveLawData.M_verbose ),
    M_vectorMaterialFlags              ( structuralConstitutiveLawData.M_vectorMaterialFlags )
{
}

// ===================================================
// Operators
// ===================================================
StructuralConstitutiveLawData&
StructuralConstitutiveLawData::operator= ( const StructuralConstitutiveLawData& structuralConstitutiveLawData )
{
    if ( this != &structuralConstitutiveLawData )
    {
        M_time                             = structuralConstitutiveLawData.M_time;
        M_timeAdvance                      = structuralConstitutiveLawData.M_timeAdvance;
        M_density                          = structuralConstitutiveLawData.M_density;
        M_thickness                        = structuralConstitutiveLawData.M_thickness;
        M_externalPressure                 = structuralConstitutiveLawData.M_externalPressure;
        M_materialsFlagSet                 = structuralConstitutiveLawData.M_materialsFlagSet;
        M_poisson                          = structuralConstitutiveLawData.M_poisson;
        M_young                            = structuralConstitutiveLawData.M_young;
        M_bulk                             = structuralConstitutiveLawData.M_bulk;
        M_alpha                            = structuralConstitutiveLawData.M_alpha;
        M_gamma                            = structuralConstitutiveLawData.M_gamma;
        M_order                            = structuralConstitutiveLawData.M_order;
        M_verbose                          = structuralConstitutiveLawData.M_verbose;
        M_vectorMaterialFlags              = structuralConstitutiveLawData.M_vectorMaterialFlags;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
StructuralConstitutiveLawData::setup ( const GetPot& dataFile, const std::string& section )
{
    // If data time has not been set
    if ( !M_time.get() )
    {
        M_time.reset ( new time_Type ( dataFile, section + "/time_discretization" ) );
    }

    if ( !M_timeAdvance.get() )
    {
        M_timeAdvance.reset ( new timeAdvance_Type ( dataFile, section + "/time_discretization" ) );
    }

    // physics
    M_solidType = dataFile ( ( section + "/physics/solidType" ).data(), "NO_DEFAULT_SOLID_TYPE" );
    M_externalPressure = dataFile ( ( section + "/physics/externalPressure" ).data(), 0. );
    M_density   = dataFile ( ( section + "/physics/density"   ).data(), 1. );
    M_thickness = dataFile ( ( section + "/physics/thickness" ).data(), 0.1 );

    UInt materialsNumber = dataFile.vector_variable_size ( ( section + "/physics/material_flag" ).data() );

    ASSERT ( materialsNumber, "Set the materrial_flag variable in [solid]/physics");

    if ( materialsNumber == 0 )
    {
        // If no material is specified in the data file the code assume that there is just one material
        // and by default it is memorized with ID 1. Getters and Setters have been designed to deal with thic choice.
        M_materialsFlagSet = false;
	M_vectorMaterialFlags.resize(1);

	M_vectorMaterialFlags[0] = 1;
        M_young[1]   = dataFile( ( section + "/physics/young"   ).data(), 0. );
        M_poisson[1] = dataFile( ( section + "/physics/poisson" ).data(), 0. );

        M_bulk[1] = dataFile ( ( section + "/physics/bulk"   ).data(), 1e9 );
        M_alpha[1] = dataFile ( ( section + "/physics/alpha" ).data(), 3e6 );
        M_gamma[1] = dataFile ( ( section + "/physics/gamma" ).data(), 0.8 );
    }
    else
    {
        M_materialsFlagSet = true;

        // These asserts are commented because in some cases we need to initialize the materials with default Young and Poisson, setting the correct values a posteriori.
        //ASSERT( M_materialsFlagSet == dataFile.vector_variable_size( ( section + "/physics/young"   ).data()),  "!!! ERROR: Inconsistent size for Young Modulus !!!");
        //ASSERT( M_materialsFlagSet == dataFile.vector_variable_size( ( section + "/physics/poisson" ).data() ), "!!! ERROR: Inconsistent size for Poisson Coef. !!!");

        UInt material (0);
        for ( UInt i (0) ; i < materialsNumber ; ++i )
        {

  	    M_vectorMaterialFlags.resize( materialsNumber );
            material            = dataFile( ( section + "/physics/material_flag" ).data(), 0., i );

	    M_vectorMaterialFlags[i] = material;
            M_young[material]   = dataFile( ( section + "/physics/young"         ).data(), 0., i );
            M_poisson[material] = dataFile( ( section + "/physics/poisson"       ).data(), 0., i );

            M_bulk[material] = dataFile( ( section + "/physics/bulk"         ).data(), 1e9, i );
            M_alpha[material] = dataFile( ( section + "/physics/alpha"       ).data(), 3e6, i );
            M_gamma[material] = dataFile( ( section + "/physics/gamma"       ).data(), 0.8, i );
        }
    }

    // space_discretization
    M_order            = dataFile ( ( section + "/space_discretization/order" ).data(), "P1" );

    // miscellaneous
    M_verbose          = dataFile ( ( section + "/miscellaneous/verbose" ).data(), 0 );
    M_useExactJacobian = dataFile ( ( section + "/useExactJacobian"      ).data(), false );
}

void
StructuralConstitutiveLawData::showMe ( std::ostream& output ) const
{
    // physics
    output << "\n*** Values for data [solid/physics]\n\n";
    output << "external pressure                = " << M_externalPressure << std::endl;
    output << "density                          = " << M_density << std::endl;
    output << "thickness                        = " << M_thickness << std::endl;
    for ( materialContainerIterator_Type i = M_young.begin() ; i != M_young.end() ; ++i )
    {
        output << "young[" << i->first << "]                         = " << i->second << std::endl;
    }
    for ( materialContainerIterator_Type i = M_poisson.begin() ; i != M_poisson.end() ; ++i )
    {
        output << "poisson[" << i->first << "]                       = " << i->second << std::endl;
    }

    for ( materialContainerIterator_Type i = M_bulk.begin() ; i != M_bulk.end() ; ++i )
    {
        output << "bulk[" << i->first << "]                       = " << i->second << std::endl;
    }

    for ( materialContainerIterator_Type i = M_alpha.begin() ; i != M_alpha.end() ; ++i )
    {
        output << "alpha[" << i->first << "]                       = " << i->second << std::endl;
    }

    for ( materialContainerIterator_Type i = M_gamma.begin() ; i != M_gamma.end() ; ++i )
    {
        output << "gamma[" << i->first << "]                       = " << i->second << std::endl;
    }

    for ( materialContainerIterator_Type i = M_poisson.begin() ; i != M_poisson.end() ; ++i )
    {
        output << "Lame - lambda[" << i->first << "]                 = " << lambda ( i->first ) << std::endl;
        output << "Lame - mu[" << i->first << "]                     = " << mu ( i->first ) << std::endl;
    }

    for ( UInt i (0); i < M_vectorMaterialFlags.size(); i++ )
    {
        output << "Position:" << i << " -> Material Flag:            = " << M_vectorMaterialFlags[i] << std::endl;
    }


    output << "\n*** Values for data [solid/miscellaneous]\n\n";
    output << "verbose                          = " << M_verbose << std::endl;

    output << "\n*** Values for data [solid/space_discretization]\n\n";
    output << "FE order                         = " << M_order << std::endl;

    output << "\n*** Values for data [solid/time_discretization]\n\n";
    M_time->showMe ( output );
    M_timeAdvance->showMe ( output );
}

// ===================================================
// Get Method
// ===================================================
Real
StructuralConstitutiveLawData::poisson ( const UInt& material ) const
{
    materialContainer_Type::const_iterator IT;

    if ( M_materialsFlagSet )
    {
        IT = M_poisson.find ( material );
    }
    else
    {
        IT = M_poisson.find ( 1 );
    }

    if ( IT != M_poisson.end() )
    {
        return IT->second;
    }
    else
    {
        std::cout << " !!! Warning: the Poisson modulus has not been set !!!" << std::endl;
        return 0;
    }
}

Real
StructuralConstitutiveLawData::young ( const UInt& material ) const
{
    materialContainer_Type::const_iterator IT;
    if ( M_materialsFlagSet )
    {
        IT = M_young.find ( material );
    }
    else
    {
        IT = M_young.find ( 1 );
    }

    if ( IT != M_young.end() )
    {
        return IT->second;
    }
    else
    {
        std::cout << " !!! Warning: the Young modulus has not been set !!!" << std::endl;
        return 0;
    }
}

Real
StructuralConstitutiveLawData::bulk ( const UInt& material ) const
{
    materialContainer_Type::const_iterator IT;
    if ( M_materialsFlagSet )
    {
        IT = M_bulk.find ( material );
    }
    else
    {
        IT = M_bulk.find ( 1 );
    }

    if ( IT != M_bulk.end() )
    {
        return IT->second;
    }
    else
    {
        std::cout << " !!! Warning: the bulk modulus has not been set !!!" << std::endl;
        return 0;
    }
}

Real
StructuralConstitutiveLawData::alpha ( const UInt& material ) const
{
    materialContainer_Type::const_iterator IT;
    if ( M_materialsFlagSet )
    {
        IT = M_alpha.find ( material );
    }
    else
    {
        IT = M_alpha.find ( 1 );
    }

    if ( IT != M_alpha.end() )
    {
        return IT->second;
    }
    else
    {
        std::cout << " !!! Warning: the alpha modulus has not been set !!!" << std::endl;
        return 0;
    }
}

Real
StructuralConstitutiveLawData::gamma ( const UInt& material ) const
{
    materialContainer_Type::const_iterator IT;
    if ( M_materialsFlagSet )
    {
        IT = M_gamma.find ( material );
    }
    else
    {
        IT = M_gamma.find ( 1 );
    }

    if ( IT != M_gamma.end() )
    {
        return IT->second;
    }
    else
    {
        std::cout << " !!! Warning: the gamma modulus has not been set !!!" << std::endl;
        return 0;
    }
}


Real
StructuralConstitutiveLawData::lambda ( const UInt& material ) const
{
    Real youngC, poissonC;

    youngC   = this->young ( material );
    poissonC = this->poisson ( material );

    return youngC * poissonC / ( ( 1.0 + poissonC ) * ( 1.0 - 2.0 * poissonC ) );
}

Real
StructuralConstitutiveLawData::mu ( const UInt& material ) const
{
    Real youngC, poissonC;

    youngC   = this->young ( material );
    poissonC = this->poisson ( material );

    return youngC / ( 2.0 * ( 1.0 + poissonC ) );
}

}
