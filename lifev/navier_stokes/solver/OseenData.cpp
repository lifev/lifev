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
    @brief File containing the implementation of the file OseenData.hpp

    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @contributor Alexis Aposporidis <aapospo@emory.edu>
    @maintainer

    @date 01-09-2009

 */


#include <lifev/navier_stokes/solver/OseenData.hpp>
#include <lifev/core/LifeV.hpp>


namespace LifeV
{


// ===================================================
// Constructors
// ===================================================

OseenData::OseenData( ) :
    M_time                             ( ),
    M_timeAdvance                      ( ),
    M_density                          ( ),
    M_viscosity                        ( ),
    M_uOrder                           ( ),
    M_pOrder                           ( ),
    M_verbose                          ( ),
    M_dumpInit                         ( ),
    M_dumpPeriod                       ( ),
    M_factor                           ( ),
    M_stokes                           ( false ),
    M_stabMethod                       ( ),
    M_semiImplicit                     ( false ),
    M_shapeDerivatives                 ( false ),
    M_domainVelImplicit                ( false ),
    M_convectiveImplicit               ( false ),
    M_computeMeanValuesPerSection      ( ),
    M_NbZSections                      ( ),
    M_ToleranceSection                 ( ),
    M_XSectionFrontier                 ( ),
    M_ZSectionInit                     ( ),
    M_ZSectionFinal                    ( ),
    M_NbPolygonEdges                   ( ),
    M_stabilizationList                ( "fluid/space_discretization/stabilization" ),
    M_conservativeFormulation          (true)
{
}


OseenData::OseenData ( const OseenData& oseenData ) :
    M_time                             ( oseenData.M_time ),
    M_timeAdvance                      ( oseenData.M_timeAdvance ),
    M_fluidNumber                      ( oseenData.M_fluidNumber ),
    M_density                          ( oseenData.M_density ),
    M_viscosity                        ( oseenData.M_viscosity ),
    M_uOrder                           ( oseenData.M_uOrder ),
    M_pOrder                           ( oseenData.M_pOrder ),
    M_verbose                          ( oseenData.M_verbose ),
    M_dumpInit                         ( oseenData.M_dumpInit ),
    M_dumpPeriod                       ( oseenData.M_dumpPeriod ),
    M_factor                           ( oseenData.M_factor ),
    M_stokes                           ( oseenData.M_stokes ),
    M_stabMethod                       ( oseenData.M_stabMethod ),
    M_semiImplicit                     ( false ),
    M_shapeDerivatives                 ( false ),
    M_domainVelImplicit                ( false ),
    M_convectiveImplicit               ( false ),
    M_computeMeanValuesPerSection      ( oseenData.M_computeMeanValuesPerSection ),
    M_NbZSections                      ( oseenData.M_NbZSections ),
    M_ToleranceSection                 ( oseenData.M_ToleranceSection ),
    M_XSectionFrontier                 ( oseenData.M_XSectionFrontier ),
    M_ZSectionInit                     ( oseenData.M_ZSectionInit ),
    M_ZSectionFinal                    ( oseenData.M_ZSectionFinal ),
    M_NbPolygonEdges                   ( oseenData.M_NbPolygonEdges ),
    M_stabilizationList                ( oseenData.M_stabilizationList ),
    M_conservativeFormulation          ( false )
{
}






// ===================================================
// Methods
// ===================================================

OseenData&
OseenData::operator= ( const OseenData& oseenData )
{
    if ( this != &oseenData )
    {
        M_time                             = oseenData.M_time;
        M_timeAdvance                      = oseenData.M_timeAdvance;
        M_fluidNumber                      = oseenData.M_fluidNumber;
        M_density                          = oseenData.M_density;
        M_viscosity                        = oseenData.M_viscosity;
        M_uOrder                           = oseenData.M_uOrder;
        M_pOrder                           = oseenData.M_pOrder;
        M_verbose                          = oseenData.M_verbose;
        M_dumpInit                         = oseenData.M_dumpInit;
        M_dumpPeriod                       = oseenData.M_dumpPeriod;
        M_factor                           = oseenData.M_factor;
        M_stokes                           = oseenData.M_stokes;
        M_stabMethod                       = oseenData.M_stabMethod;
        M_semiImplicit                     = oseenData.M_semiImplicit;
        M_shapeDerivatives                 = oseenData.M_shapeDerivatives;
        M_domainVelImplicit                = oseenData.M_domainVelImplicit;
        M_convectiveImplicit               = oseenData.M_convectiveImplicit;
        M_computeMeanValuesPerSection      = oseenData.M_computeMeanValuesPerSection;
        M_NbZSections                      = oseenData.M_NbZSections;
        M_ToleranceSection                 = oseenData.M_ToleranceSection;
        M_XSectionFrontier                 = oseenData.M_XSectionFrontier;
        M_ZSectionInit                     = oseenData.M_ZSectionInit;
        M_ZSectionFinal                    = oseenData.M_ZSectionFinal;
        M_NbPolygonEdges                   = oseenData.M_NbPolygonEdges;
        M_stabilizationList                = oseenData.M_stabilizationList;
        M_conservativeFormulation          = oseenData.M_conservativeFormulation;
    }

    return *this;
}


void
OseenData::setup ( const GetPot& dataFile, const std::string& section )
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

    M_stabilizationList.add ( "ip", IP_STABILIZATION,   "interior penalty " );
    M_stabilizationList.add ( "sd", SD_STABILIZATION,   "stream-line diffusion" );
    M_stabilizationList.add ( "none", NO_STABILIZATION, "none (default)" );

    // Physics
    UInt temp = dataFile ( (section + "/physics/fluid_number" ).data(), 0 );

    if (temp == 0) // Old fashion of declaring fluids
    {
        M_fluidNumber = 1;
        M_density.push_back ( dataFile ( ( section + "/physics/density" ).data(), 1. ) );
        M_viscosity.push_back ( dataFile ( ( section + "/physics/viscosity" ).data(), 1. ) );
    }
    else   // New fashion of declaring fluids
    {
        M_fluidNumber = temp;
        M_density = std::vector<Real> (temp, 0);
        M_viscosity = std::vector<Real> (temp, 0);

        for (UInt iter_fluid (0); iter_fluid < temp; ++iter_fluid)
        {
            // build the section name
            std::string iter_fluid_section ( section + "/physics/fluid_");
            iter_fluid_section += number2string (iter_fluid);

            // Read the quantities
            M_density[iter_fluid] = dataFile ( (iter_fluid_section + "/density").c_str() , 1.0);
            M_viscosity[iter_fluid] = dataFile ( (iter_fluid_section + "/viscosity").c_str() , 1.0);
        }
    }

    // FE Order
    M_uOrder       = dataFile ( ( section + "/space_discretization/vel_order" ).data(), "P1");
    M_pOrder       = dataFile ( ( section + "/space_discretization/press_order" ).data(), "P1");

    // Miscellaneous
    M_verbose      = dataFile ( ( section + "/miscellaneous/verbose" ).data(), 1 );
    M_dumpInit    = dataFile ( ( section + "/miscellaneous/dump_init" ).data(), M_time->initialTime() );
    M_dumpPeriod  = dataFile ( ( section + "/miscellaneous/dump_period" ).data(), 1 );
    M_factor       = dataFile ( ( section + "/miscellaneous/factor" ).data(), 0. );
    M_stokes       = dataFile ( ( section + "/miscellaneous/Stokes" ).data(), false );

    M_stabMethod  = NSStabilization ( M_stabilizationList.value (
                                          dataFile ( ( section + "/space_discretization/stabilization" ).data(), "none") ) );

    // Semi-implicit and shape derivatives
    M_shapeDerivatives = dataFile ( ( section + "/useShapeDerivatives" ).data(), false ) ;
    setSemiImplicit ( dataFile ( ( section + "/semiImplicit" ).data(), false ) );
    M_domainVelImplicit = dataFile ( (section + "/domainVelImplicit").data(), false );
    M_convectiveImplicit = dataFile ( (section + "/convectiveImplicit").data(), false );
    M_conservativeFormulation   = dataFile ( ( section + "/conservativeFormulation" ).data(), true );

    // Mean values per section
    M_computeMeanValuesPerSection = dataFile ( ( section + "/valuespersection/computeMeanValuesPerSection" ).data(), 0 );
    M_NbZSections      = dataFile ( ( section + "/valuespersection/nb_z_section" ).data(), 2 );
    M_ToleranceSection = dataFile ( ( section + "/valuespersection/tol_section" ).data(), 2e-2 );
    M_XSectionFrontier = dataFile ( ( section + "/valuespersection/x_section_frontier" ).data(), 0. );
    M_ZSectionInit     = dataFile ( ( section + "/valuespersection/z_section_init" ).data(), -1. );
    M_ZSectionFinal    = dataFile ( ( section + "/valuespersection/z_section_final" ).data(), 0. );
    M_NbPolygonEdges   = dataFile ( ( section + "/valuespersection/nb_polygon_edges" ).data(), 10 );
}


void
OseenData::showMe ( std::ostream& output ) const
{
    if (M_fluidNumber == 1)
    {
        output << "\n*** Values for data [fluid/physics]\n\n";
        output << "density     = " << M_density[0] << std::endl;
        output << "viscosity   = " << M_viscosity[0] << std::endl;
    }
    else
    {
        output << "\n*** Values for data [fluid/physics]\n\n";
        for (UInt iter_fluid (0); iter_fluid < M_fluidNumber; ++iter_fluid)
        {
            output << "fluid " << iter_fluid << std::endl;
            output << "density     = " << M_density[iter_fluid] << std::endl;
            output << "viscosity   = " << M_viscosity[iter_fluid] << std::endl;
        }
    }
    output << "\n*** Values for data [fluid/miscellaneous]\n\n";
    output << "verbose     = " << M_verbose << std::endl;
    output << "initial time for writing solution  = " << M_dumpInit << std::endl;
    output << "number of time steps between two consecutive dumps of the solution = " << M_dumpPeriod << std::endl;
    output << "amplification factor = " << M_factor << std::endl;
    output << "Stokes simulation    = " << M_stokes << std::endl;


    output << "\n*** Values for data [fluid/time_discretization]\n\n";
    M_time->showMe ( output );
    M_timeAdvance->showMe ( output );

    output << "stabilization = ";
    switch ( M_stabMethod )
    {
        case NO_STABILIZATION:
            output << "none" ;
            break;
        case IP_STABILIZATION:
            output << "ip" ;
            break;
        case SD_STABILIZATION:
            output << "sd" ;
            break;
    }
    output << std::endl;

    output << "\n*** Values for data [fluid/valuespersection]\n\n";
    output << "computeMeanValuesPerSection (switch 0: don't compute, 1: compute)  = "
           << M_computeMeanValuesPerSection << std::endl;
    output << "nb_z_section       = " << M_NbZSections << std::endl;
    output << "tol_section        = " << M_ToleranceSection << std::endl;
    output << "x_section_frontier = " << M_XSectionFrontier << std::endl;
    output << "z_section_init     = " << M_ZSectionInit << std::endl;
    output << "z_section_final    = " << M_ZSectionFinal << std::endl;
    output << "nb_polygon_edges   = " << M_NbPolygonEdges << std::endl;
}


} //end namespace LifeV
