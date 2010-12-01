/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file dataNavierStokes.cpp
*/


#include <life/lifesolver/dataNavierStokes.hpp>


namespace LifeV
{


// ===================================================
// Constructors
// ===================================================
//template <typename Mesh>
DataNavierStokes::DataNavierStokes( ) :
        M_time                             ( ),
        //M_mesh                             ( ),
        M_density                          ( ),
        M_viscosity                        ( ),
        M_uOrder                           ( ),
        M_pOrder                           ( ),
        M_verbose                          ( ),
        M_dump_init                        ( ),
        M_dump_period                      ( ),
        M_factor                           ( ),
        M_Stokes                           ( false ),
        M_stab_method                      ( ),
        M_semiImplicit                     ( false ),
        M_shapeDerivatives                 ( false ),
        M_computeMeanValuesPerSection      ( ),
        M_NbZSections                      ( ),
        M_ToleranceSection                 ( ),
        M_XSectionFrontier                 ( ),
        M_ZSectionInit                     ( ),
        M_ZSectionFinal                    ( ),
        M_NbPolygonEdges                   ( ),
        M_stabilization_list               ( "fluid/space_discretization/stabilization" )
{
}

// template <typename Mesh>
// DataNavierStokes<Mesh>::DataNavierStokes( const GetPot& dataFile,
//                                           const Time_ptrType DataTime,
//                                           const Mesh_ptrType DataMesh ) :
//         M_time                             ( DataTime ),
//         M_mesh                             ( DataMesh ),
//         M_density                          ( ),
//         M_viscosity                        ( ),
//         M_uOrder                           ( ),
//         M_pOrder                           ( ),
//         M_verbose                          ( ),
//         M_dump_init                        ( ),
//         M_dump_period                      ( ),
//         M_factor                           ( ),
//         M_Stokes                           ( false ),
//         M_stab_method                      ( ),
//         M_semiImplicit                     ( false ),
//         M_shapeDerivatives                 ( false ),
//         M_computeMeanValuesPerSection      ( ),
//         M_NbZSections                      ( ),
//         M_ToleranceSection                 ( ),
//         M_XSectionFrontier                 ( ),
//         M_ZSectionInit                     ( ),
//         M_ZSectionFinal                    ( ),
//         M_NbPolygonEdges                   ( ),
//         M_stabilization_list               ( "fluid/space_discretization/stabilization" )
// {
//     setup( dataFile );
// }

//template <typename Mesh>
DataNavierStokes::DataNavierStokes( const DataNavierStokes& dataNavierStokes ) :
        M_time                             ( dataNavierStokes.M_time ),
        //M_mesh                             ( dataNavierStokes.M_mesh ),
        M_fluid_number                     ( dataNavierStokes.M_fluid_number ),
        M_density                          ( dataNavierStokes.M_density ),
        M_viscosity                        ( dataNavierStokes.M_viscosity ),
        M_uOrder                           ( dataNavierStokes.M_uOrder ),
        M_pOrder                           ( dataNavierStokes.M_pOrder ),
        M_verbose                          ( dataNavierStokes.M_verbose ),
        M_dump_init                        ( dataNavierStokes.M_dump_init ),
        M_dump_period                      ( dataNavierStokes.M_dump_period ),
        M_factor                           ( dataNavierStokes.M_factor ),
        M_Stokes                           ( dataNavierStokes.M_Stokes ),
        M_stab_method                      ( dataNavierStokes.M_stab_method ),
        M_semiImplicit                     ( false ),
        M_shapeDerivatives                 ( false ),
        M_computeMeanValuesPerSection      ( dataNavierStokes.M_computeMeanValuesPerSection ),
        M_NbZSections                      ( dataNavierStokes.M_NbZSections ),
        M_ToleranceSection                 ( dataNavierStokes.M_ToleranceSection ),
        M_XSectionFrontier                 ( dataNavierStokes.M_XSectionFrontier ),
        M_ZSectionInit                     ( dataNavierStokes.M_ZSectionInit ),
        M_ZSectionFinal                    ( dataNavierStokes.M_ZSectionFinal ),
        M_NbPolygonEdges                   ( dataNavierStokes.M_NbPolygonEdges ),
        M_stabilization_list               ( dataNavierStokes.M_stabilization_list )
{
}






// ===================================================
// Methods
// ===================================================
//template <typename Mesh>

DataNavierStokes&
DataNavierStokes::operator=( const DataNavierStokes& dataNavierStokes )
{
    if ( this != &dataNavierStokes )
    {
        M_time                             = dataNavierStokes.M_time;
        //M_mesh                             = dataNavierStokes.M_mesh;
        M_fluid_number                     = dataNavierStokes.M_fluid_number;
        M_density                          = dataNavierStokes.M_density;
        M_viscosity                        = dataNavierStokes.M_viscosity;
        M_uOrder                           = dataNavierStokes.M_uOrder;
        M_pOrder                           = dataNavierStokes.M_pOrder;
        M_verbose                          = dataNavierStokes.M_verbose;
        M_dump_init                        = dataNavierStokes.M_dump_init;
        M_dump_period                      = dataNavierStokes.M_dump_period;
        M_factor                           = dataNavierStokes.M_factor;
        M_Stokes                           = dataNavierStokes.M_Stokes;
        M_stab_method                      = dataNavierStokes.M_stab_method;
        M_semiImplicit                     = dataNavierStokes.M_semiImplicit;
        M_shapeDerivatives                 = dataNavierStokes.M_shapeDerivatives;
        M_computeMeanValuesPerSection      = dataNavierStokes.M_computeMeanValuesPerSection;
        M_NbZSections                      = dataNavierStokes.M_NbZSections;
        M_ToleranceSection                 = dataNavierStokes.M_ToleranceSection;
        M_XSectionFrontier                 = dataNavierStokes.M_XSectionFrontier;
        M_ZSectionInit                     = dataNavierStokes.M_ZSectionInit;
        M_ZSectionFinal                    = dataNavierStokes.M_ZSectionFinal;
        M_NbPolygonEdges                   = dataNavierStokes.M_NbPolygonEdges;
        M_stabilization_list               = dataNavierStokes.M_stabilization_list;
    }

    return *this;
}

// template <typename Mesh>
void
DataNavierStokes::setup( const GetPot& dataFile, const std::string& section )
{
    // If data time has not been set
    if ( !M_time.get() )
        M_time.reset( new Time_Type( dataFile, section + "/time_discretization" ) );

    // If data mesh has not been set
    // if ( !M_mesh.get() )
    //     M_mesh.reset( new Mesh_Type( dataFile, section + "/space_discretization" ) );

    M_stabilization_list.add( "ip", IP_STABILIZATION,   "interior penalty " );
    M_stabilization_list.add( "sd", SD_STABILIZATION,   "stream-line diffusion" );
    M_stabilization_list.add( "none", NO_STABILIZATION, "none (default)" );

    // Physics
    UInt temp = dataFile( (section + "/physics/fluid_number" ).data(), 0 );

    if (temp == 0) // Old fashion of declaring fluids
    {
        M_fluid_number = 1;
        M_density.push_back( dataFile( ( section + "/physics/density" ).data(), 1. ) );
        M_viscosity.push_back ( dataFile( ( section + "/physics/viscosity" ).data(), 1. ) );
    }
    else   // New fashion of declaring fluids
    {
        M_fluid_number = temp;
        M_density = std::vector<Real>(temp,0);
        M_viscosity = std::vector<Real>(temp,0);

        for (UInt iter_fluid(0); iter_fluid<temp; ++iter_fluid)
        {
            // build the section name
            std::string iter_fluid_section( section + "/physics/fluid_");
            iter_fluid_section += number2string(iter_fluid);

            // Read the quantities
            M_density[iter_fluid]= dataFile((iter_fluid_section+"/density").c_str() , 1.0);
            M_viscosity[iter_fluid]= dataFile((iter_fluid_section+"/viscosity").c_str() , 1.0);
        }
    }

    // FE Order
    M_uOrder       = dataFile( ( section + "/space_discretization/vel_order" ).data(), "P1");
    M_pOrder       = dataFile( ( section + "/space_discretization/press_order" ).data(), "P1");

    // Miscellaneous
    M_verbose      = dataFile( ( section + "/miscellaneous/verbose" ).data(), 1 );
    M_dump_init    = dataFile( ( section + "/miscellaneous/dump_init" ).data(), M_time->getInitialTime() );
    M_dump_period  = dataFile( ( section + "/miscellaneous/dump_period" ).data(), 1 );
    M_factor       = dataFile( ( section + "/miscellaneous/factor" ).data(), 0. );
    M_Stokes       = dataFile( ( section + "/miscellaneous/Stokes" ).data(), false );

    M_stab_method  = NSStabilization ( M_stabilization_list.value(
                                           dataFile( ( section + "/space_discretization/stabilization" ).data(), "none") ) );

    // Semi-implicit and shape derivatives
    M_semiImplicit     = dataFile( ( section + "/semiImplicit" ).data(), false ) ;
    M_shapeDerivatives = dataFile( ( section + "/useShapeDerivatives" ).data(), false ) ;
    setSemiImplicit( M_semiImplicit );

    // Mean values per section
    M_computeMeanValuesPerSection = dataFile( ( section + "/valuespersection/computeMeanValuesPerSection" ).data(), 0 );
    M_NbZSections      = dataFile( ( section + "/valuespersection/nb_z_section" ).data(), 2 );
    M_ToleranceSection = dataFile( ( section + "/valuespersection/tol_section" ).data(), 2e-2 );
    M_XSectionFrontier = dataFile( ( section + "/valuespersection/x_section_frontier" ).data(), 0. );
    M_ZSectionInit     = dataFile( ( section + "/valuespersection/z_section_init" ).data(), -1. );
    M_ZSectionFinal    = dataFile( ( section + "/valuespersection/z_section_final" ).data(), 0. );
    M_NbPolygonEdges   = dataFile( ( section + "/valuespersection/nb_polygon_edges" ).data(), 10 );
}

//template <typename Mesh>
void
DataNavierStokes::showMe( std::ostream& output ) const
{
    if (M_fluid_number == 1)
    {
        output << "\n*** Values for data [fluid/physics]\n\n";
        output << "density     = " << M_density[0] << std::endl;
        output << "viscosity   = " << M_viscosity[0] << std::endl;
    }
    else
    {
        output << "\n*** Values for data [fluid/physics]\n\n";
        for (UInt iter_fluid(0); iter_fluid<M_fluid_number; ++iter_fluid)
        {
            output << "fluid " << iter_fluid << std::endl;
            output << "density     = " << M_density[iter_fluid] << std::endl;
            output << "viscosity   = " << M_viscosity[iter_fluid] << std::endl;
        }
    }
    output << "\n*** Values for data [fluid/miscellaneous]\n\n";
    output << "verbose     = " << M_verbose << std::endl;
    output << "initial time for writing solution  = " << M_dump_init << std::endl;
    output << "number of time steps between two consecutive dumps of the solution = " << M_dump_period << std::endl;
    output << "amplification factor = " << M_factor << std::endl;
    output << "Stokes simulation    = " << M_Stokes << std::endl;

    // output << "\n*** Values for data [fluid/space_discretization]\n\n";
    // M_mesh->showMe( output );
    output << "\n*** Values for data [fluid/time_discretization]\n\n";
    M_time->showMe( output );

    output << "stabilization = ";
    switch ( M_stab_method )
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
}
