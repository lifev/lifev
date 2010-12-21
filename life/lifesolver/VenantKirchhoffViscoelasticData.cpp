//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

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
  @brief dataSecondOrder - Class to secondorder problem (S. Venant Kirchhoff Viscoelastic)

  @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>

  @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
  @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 */

#include <life/lifesolver/dataSecondOrder.hpp>

namespace LifeV
{

//===================================================
// Constructor
//===================================================

DataSecondOrder::DataSecondOrder() :
        M_time                           ( ),
        M_density                        ( ),
        M_thickness                      ( ),
        M_poisson                        ( ),
        M_young                          ( ),
        M_gamma                          ( ),
        M_beta                           ( ),
        M_factor                         ( ),
	M_verbose                        ( ),
	M_order                          ( ),                       
      	M_damping                        ( )
{
}

DataSecondOrder::DataSecondOrder( const DataSecondOrder& DataSecondOrder):
        DataTime               ( DataSecondOrder),
        M_density              ( DataSecondOrder.M_density ),
	M_thickness            ( DataSecondOrder.M_thickness ),
        M_poisson              ( DataSecondOrder.M_poisson ),
        M_young                ( DataSecondOrder.M_young ),
        M_gamma                ( DataSecondOrder.M_gamma ),
        M_beta                 ( DataSecondOrder.M_beta ),
        M_factor               ( DataSecondOrder.M_factor ),
        M_verbose              ( DataSecondOrder.M_verbose ),
        M_order                ( DataSecondOrder.M_order ),
	M_damping              ( DataSecondOrder.M_damping )
{
}

//==================================================
// Operators
//==================================================

DataSecondOrder&
DataSecondOrder::operator=( const DataSecondOrder& DataSecondOrder )
{
    if ( this != &DataSecondOrder )
    {
        M_time                    = DataSecondOrder.M_time;
        M_density                 = DataSecondOrder.M_density;
	M_thickness               = DataSecondOrder.M_thickness;
	M_poisson                 = DataSecondOrder.M_poisson;
        M_young                   = DataSecondOrder.M_young;
        M_gamma                   = DataSecondOrder.M_gamma;
        M_beta                    = DataSecondOrder.M_beta;
        M_factor                  = DataSecondOrder.M_factor;
        M_verbose                 = DataSecondOrder.M_verbose;
	M_order                   = DataSecondOrder.M_order;
	M_damping                 = DataSecondOrder.M_damping;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================

void
DataSecondOrder::setup( const GetPot& dataFile, const std::string& section )
{
    // If data time has not been set
    if ( !M_time.get() )
        M_time.reset( new Time_Type( dataFile, section + "/time_discretization" ) );

    // physics
    M_density   = dataFile( ( section + "/physics/density" ).data(), 1. );
    M_thickness = dataFile( ( section + "/physics/thickness" ).data(), 0.1 );
    
    UInt materialsNumber = dataFile.vector_variable_size( ( section + "/physics/material_flag" ).data() );
 std::cout<<"materialNumber "<<materialsNumber<<"\n";
    if ( materialsNumber == 0 )
      {
	std::cout<<"The material flag was not set from data file. Its value will be dedced from the first volume marker."<<"\n";
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
	    M_gamma[material]   = dataFile( ( section + "/physics/young" ).data(), 0., i );
            M_beta[material] = dataFile( ( section + "/physics/poisson" ).data(), 0., i );
	  }
    }

    M_damping     = dataFile( (section+"/damping").data(), false);

    // space_discretization
    M_order     = dataFile( (section+"/space_discretization/order").data(), "P1" );

    // miscellaneous
    M_factor  = dataFile( (section + "/miscellaneous/factor").data(), 1.0 );
    M_verbose = dataFile( (section + "/miscellaneous/verbose").data(), 1 );

}

void
DataSecondOrder::showMe( std::ostream& output ) const
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

  for ( MaterialContainer_ConstIterator i = M_gamma.begin() ; i != M_gamma.end() ; ++i )
        output << "gamma[" << i->first << "]                         = " << i->second << std::endl;

  for ( MaterialContainer_ConstIterator i = M_beta.begin() ; i != M_beta.end() ; ++i )
        output << "beta[" << i->first << "]                         = " << i->second << std::endl;

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
DataSecondOrder::setDataTime( const Time_ptrType DataTime )
{
    M_time = DataTime;
}

void
DataSecondOrder::setDensity( const Real& density )
{
    M_density = density;
}

void
DataSecondOrder::setThickness( const Real& thickness )
{
    M_thickness = thickness;
}

void
DataSecondOrder::setPoisson( const Real& poisson, const UInt& material )
{
    M_poisson[material] = poisson;
}

void
DataSecondOrder::setYoung( const Real& young, const UInt& material )
{
    M_young[material] = young;
}

void
DataSecondOrder::setGamma( const Real& gamma, const UInt& material )
{
    M_gamma[material] = gamma;
}


void
DataSecondOrder::setBeta( const Real& beta, const UInt& material )
{
    M_beta[material] = beta;
}

// ===================================================
// Get Method
// ===================================================

DataSecondOrder::Time_ptrType
DataSecondOrder::dataTime() const
{
    return M_time;
}

const Real&
DataSecondOrder::rho() const
{
    return M_density;
}

const Real&
DataSecondOrder::thickness() const
{
    return M_thickness;
}

const Real&
DataSecondOrder::poisson( const UInt& material ) const
{
    MaterialContainer_Type::const_iterator IT = M_poisson.find( material );
    if (IT != M_poisson.end())
        return M_poisson.find( material )->second;
    else
      return 0;
}

const Real&
DataSecondOrder::young( const UInt& material ) const
{
    MaterialContainer_Type::const_iterator IT = M_young.find( material );
    if (IT != M_young.end())
        return IT->second;
    else
       return 0;
}

Real
 DataSecondOrder::lambda( const UInt& material ) const
{
    return M_young.find( material )->second * M_poisson.find( material )->second /
           ( ( 1.0 + M_poisson.find( material )->second ) * ( 1.0 - 2.0 * M_poisson.find( material )->second ) );
}

Real 
DataSecondOrder::mu( const UInt& material ) const
{
    return M_young.find( material )->second/( 2.0 * ( 1.0 + M_poisson.find( material )->second ) );
}


const Real&
DataSecondOrder::gamma( const UInt& material ) const
{
    return M_gamma.find( material )->second;
}

const Real&
DataSecondOrder::beta( const UInt& material ) const
{
    return M_beta.find( material )->second;
}

const Real&
DataSecondOrder::factor() const
{
    return M_factor;
}

const std::string&
DataSecondOrder::order() const
{
    return M_order;
}

}
