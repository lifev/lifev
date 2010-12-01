/*
  This file is part of the LifeV library
  Copyright (C) 2010 EPFL, INRIA, Politecnico di Milano and Emory University

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*!
  \file dataADR.cpp
  \date 08/2010
  \version 2.0

  \brief Implementation of DataADR

*/

#include <life/lifesolver/dataADR.hpp>

namespace LifeV
{


// Empty constructor
DataADR::DataADR():
        M_diffusionCoefficient      ( 0. ),
        M_reactionCoefficient       ( 0. ),
        M_steady                    ( 1 ),
        M_solutionFieldDimension            ( 1 ),
        M_dataTimePtr               ( new dataTime_type() ),
        M_verbose                   ( false ),
        M_solFEType                 ( "P1" ),
//                  M_advectionFieldFEType      ( "P1" ),
        M_stabilization_list        ( "Stab. list" ),
        M_stabilizationMethod       ( ADR_NO_STABILIZATION ),
        M_stabilizationCoefficient  ( 0. )
{
    M_stabilization_list.add( "ip",  ADR_IP_STABILIZATION, "interior penalty " );
    M_stabilization_list.add( "sd",  ADR_SD_STABILIZATION, "stream-line difussion" );
    M_stabilization_list.add( "none", ADR_NO_STABILIZATION,  "none (default)" );
}


// Copy constructor
DataADR::DataADR( const DataADR& dataADR ) :
        M_diffusionCoefficient      ( dataADR.M_diffusionCoefficient ),
        M_reactionCoefficient       ( dataADR.M_reactionCoefficient ),
        M_steady                    ( dataADR.M_steady ),
        M_solutionFieldDimension            ( dataADR.M_solutionFieldDimension ),
        M_dataTimePtr               ( dataADR.M_dataTimePtr ),
        M_verbose                   ( dataADR.M_verbose ),
        M_solFEType                 ( dataADR.M_solFEType ),
//                  M_advectionFieldFEType      ( dataADR.M_advectionFieldFEType ),
        M_stabilization_list        ( dataADR.M_stabilization_list ),
        M_stabilizationMethod       ( dataADR.M_stabilizationMethod ),
        M_stabilizationCoefficient  ( dataADR.M_stabilizationCoefficient )
{
}


// Set up the class reading from data file
void
DataADR::setup( const GetPot& dfile, const std::string& section )
{
    // We want a slash dividing the data file section from the variable name but only
    // when not provided by the user or when not looking in the root of the data file
    std::string corrected_section( section );
    if ( ( ! section.empty() ) && ( section[section.length()-1] != '/' ) )
        corrected_section = section + '/';
    // physics
    M_diffusionCoefficient = dfile( (corrected_section+"physics/diffusionCoefficient").data(), 1. );
    M_reactionCoefficient  = dfile( (corrected_section+"physics/reactionCoefficient").data(), 1. );
    M_steady               = dfile( (corrected_section+"physics/steady").data(), 1  );
    M_solutionFieldDimension       = dfile( (corrected_section+"physics/solutionFieldDimension").data(), 1  );

    // miscellaneous
    M_verbose = dfile( (corrected_section+"miscellaneous/verbose").data(), 1 );

    // type of finite element (P1, P2, ...)
    M_solFEType = dfile( (corrected_section+"space_discretization/sol_FEtype").data(), "P1");
//    M_advectionFieldFEType = dfile( (corrected_section+"space_discretization/advectionField_FEtype").data(), "P1");

    // stabilization
    M_stabilizationMethod = ADRStabilization ( M_stabilization_list.value( dfile( (corrected_section+"space_discretization/stabilization").data(), "none") ) );
    M_stabilizationCoefficient = dfile( (corrected_section+"space_discretization/stabilizationCoefficient").data(), 0. );
}


// Output
void
DataADR::showMe( std::ostream& c ) const
{
    c << "\n*** DataADR: values for user-defined data\n";

    c << "\n[/physics]" << std::endl;
    c << "diffusionCoefficient\t\t= " << M_diffusionCoefficient << std::endl;
    c << "reactionCoefficient\t\t= " << M_reactionCoefficient << std::endl;
    c << "steady\t\t= " << M_steady << std::endl;
    c << "solutionFieldDimension\t\t= " << M_solutionFieldDimension << std::endl;

    c << "\n[/miscellaneous]" << std::endl;
    c << "verbose\t\t= " << M_verbose << std::endl;

    c << "\n[/space_discretization]" << std::endl;
    c << "sol_FEtype\t\t= " << M_solFEType << std::endl;
//    c << "advectionField_FEtype\t\t= " << M_advectionFieldFEType << std::endl;

    c << "stabilization\t\t= ";
    switch ( M_stabilizationMethod )
    {
    case ADR_NO_STABILIZATION:
        c << "none" ;
        break;
    case ADR_IP_STABILIZATION:
        c << "ip" ;
        break;
    case ADR_SD_STABILIZATION:
        c << "sd" ;
        break;
    }
    c << std::endl;
    c << "stabilizationCoefficient\t\t= " << M_stabilizationCoefficient << std::endl;
    c << std::endl;
}

}
