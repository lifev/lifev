/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file dataADR.cpp

*/

#include <life/lifesolver/dataADR.hpp>
#include <life/lifecore/dataString.hpp>


//
// IMPLEMENTATION
//


// Constructor
namespace LifeV{


DataADR::DataADR():
    DataTime                     ( ),
    M_diffusivity                ( 0. ),
    M_react                      ( 0. ),
    M_verbose                    ( false ),
    M_factor                     ( 0. ),
    M_order                      ( 0 ),
    M_stab_method                ( ),
    M_stabilization_list         ( "Stab. list" )
{
}




DataADR::DataADR( const DataADR& dataADR ) :
    DataTime                     ( dataADR ),
    M_diffusivity                ( dataADR.M_diffusivity),
    M_react                      ( dataADR.M_react),
    M_verbose                    ( dataADR.M_verbose),
    M_factor                     ( dataADR.M_factor),
    M_order                      ( dataADR.M_order),
    M_stab_method                ( dataADR.M_stab_method),
    M_stabilization_list         ( dataADR.M_stabilization_list)
{
}




void
DataADR::setup(  const GetPot& dfile,  const std::string& section )
{
    M_stabilization_list.add( "ip",  ADR_IP_STABILIZATION, "interior penalty " );
    M_stabilization_list.add( "sd",  ADR_SD_STABILIZATION, "stream-line difussion" );
    M_stabilization_list.add( "none", ADR_NO_STABILIZATION,  "none (default)" );

    // physics
    M_diffusivity = dfile( "adr/physics/diffusivity", 1. );
    M_react = dfile( "adr/physics/react", 1. );

    //miscellaneous
    M_verbose = dfile( "adr/miscellaneous/verbose", 1 );
    // _dump_init = dfile( "adr/miscellaneous/dump_init", getInitialTime() );
    // _dump_period = dfile( "adr/miscellaneous/dump_period", 1 );
    M_factor = dfile( "adr/miscellaneous/factor", 0. );

    M_order = dfile( "adr/space_discretization/order", "P1");

    M_stab_method = ADRStabilization ( M_stabilization_list.value( dfile( "adr/space_discretization/stabilization", "none") ) );

    //mean values per section
    // M_computeMeanValuesPerSection =
    //     dfile( "adr/valuespersection/computeMeanValuesPerSection", 0 );
    // M_NbZSections =
    //     dfile( "adr/valuespersection/nb_z_section", 2 );
    // M_ToleranceSection =
    //     dfile( "adr/valuespersection/tol_section", 2e-2 );
    // M_XSectionFrontier =
    //     dfile( "adr/valuespersection/x_section_frontier", 0. );
    // M_ZSectionInit =
    //     dfile( "adr/valuespersection/z_section_init", -1. );
    // M_ZSectionFinal =
    //     dfile( "adr/valuespersection/z_section_final", 0. );
    // M_NbPolygonEdges =
    //     dfile( "adr/valuespersection/nb_polygon_edges", 10 );
}

// Output
void
DataADR::showMe( std::ostream& c )
{
    // physics
    c << "\n*** Values for data [adr/miscellaneous]\n\n";
    c << "verbose     = " << M_verbose << std::endl;
    //    c << "initial time for writing solution  = " << _dump_init << std::endl;
    //c << "number of time steps between two consecutive dumps of the solution = " << _dump_period << std::endl;
    c << "amplification factor = " << M_factor << std::endl;


    //c << "\n*** Values for data [adr/space_discretization]\n\n";
    //DataMesh::showMe( c );
    c << "\n*** Values for data [adr/time_discretization]\n\n";
    DataTime::showMe( c );
    c << "stabilization = ";
    switch( M_stab_method )
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


    // c << "\n*** Values for data [adr/valuespersection]\n\n";
    // c << "computeMeanValuesPerSection (switch 0: don't compute, 1: compute)  = "
    //   << M_computeMeanValuesPerSection << std::endl;
    // c << "nb_z_section        = " << M_NbZSections << std::endl;
    // c << "tol_section         = " << M_ToleranceSection << std::endl;
    // c << "x_section_frontier  = " << M_XSectionFrontier << std::endl;
    // c << "z_section_init      = " << M_ZSectionInit << std::endl;
    // c << "z_section_final     = " << M_ZSectionFinal << std::endl;
    // c << "nb_polygon_edges    = " << M_NbPolygonEdges << std::endl;

}

////////////////////
// The viscosity

Real
DataADR::diffusivity() const
{
    return M_diffusivity;
}

Real
DataADR::reaction() const
{
    return M_react;
}
////////////////

// verbose variable
UInt
DataADR::verbose() const
{
    return M_verbose;
}

// Amplification factor
Real
DataADR::factor() const
{
    return M_factor;
}

std::string
DataADR::order() const
{
    return M_order;
}


// Stabilization method
ADRStabilization DataADR::stabilization() const
{
    return M_stab_method;
}


//! Mean Values per Sections
//! compute (0) or not (1) the mean values per sections

// UInt
// DataADR::computeMeanValuesPerSection() const
// {
//     return M_computeMeanValuesPerSection;
// }

// //! number of sections
// UInt
// DataADR::NbZSections() const
// {
//     return M_NbZSections;
// }

// //! tolerance for point proximity
// Real
// DataADR::ToleranceSection() const
// {
//     return M_ToleranceSection;
// }
// //! x (see NSHandler): point at the frontier for computation of the area
// //! with a polygonal formula (x -> displacement of the boundary)
// Real
// DataADR::XSectionFrontier() const
// {
//     return M_XSectionFrontier;
// }
// //! lower section
// Real
// DataADR::ZSectionInit() const
// {
//     return M_ZSectionInit;
// }
// //! upper section
// Real
// DataADR::ZSectionFinal() const
// {
//     return M_ZSectionFinal;
// }


// //! number of edges of the polygon (in mesh) describing the circle

// UInt
// DataADR::NbPolygonEdges() const
// {
//     return M_NbPolygonEdges;
// }

// void
// DataADR::setSemiImplicit(const int& SI)
// { M_semiImplicit = SI; }

// int
// DataADR::semiImplicit() const
//  { return M_semiImplicit; }


}
