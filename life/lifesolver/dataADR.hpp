/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
  \file dataNavierStokes.h
  \author M.A. Fernandez
  \date 01/2003
  \version 1.0

  \brief File containing a class for handling NavierStokes data with GetPot

*/
#ifndef _DATANADR_H_
#define _DATANADR_H_
#include <string>
#include <iostream>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifecore/dataString.hpp>

namespace LifeV
{


/*!
  \typedef enum
*/
enum ADRStabilization
{
    ADR_NO_STABILIZATION,       //!< No stabilization
    ADR_IP_STABILIZATION,       //!< Interior penalty
    ADR_SD_STABILIZATION        //!< Stream-line diffusion
};

/*!
  \class DataADR

  Base class which holds usual data for the NavierStokes equations solvers

*/


template <typename Mesh>
class DataADR:
        public DataMesh<Mesh>,
        public DataTime
{
public:

    //! Constructor

    DataADR( const GetPot& dfile, 	const std::string& mesh_section= "adr/discretization",
									const std::string& time_section= "adr/time" );

    DataADR( const DataADR& dataNavierStokes );

    //! Ouptut
    void showMe( std::ostream& c = std::cout );


    //! external setup

    void setup( const GetPot& dfile );

    //! End time
    Real density() const;
    Real diffusivity() const;
    Real reaction() const;

    UInt verbose() const;
    Real dump_init() const;
    UInt dump_period() const;
    Real factor() const;

    std::string order() const;

    ADRStabilization stabilization() const;

    //! a way to obtain the Mean Flux, Mean Pressure and
    //! Mean Area at a section defined by z=z_data.
    UInt computeMeanValuesPerSection() const;
    UInt NbZSections() const;
    Real ToleranceSection() const;
    Real XSectionFrontier() const;
    Real ZSectionInit() const;
    Real ZSectionFinal() const;
    UInt NbPolygonEdges() const;
    int semiImplicit() const;
    void setSemiImplicit(const int& SI);

protected:
    //! Physics

    Real M_diffusivity; // Diffusivity
    Real M_react; // Reaction coefficient

    //! Miscellaneous
    UInt M_verbose; // temporal output verbose
    Real _dump_init; // time for starting the dumping of the results (Alex December 2003)
    UInt _dump_period; // frequency of the dumping (one dump after _dump_period time steps) (Alex December 2003)
    Real M_factor; // amplification factor for moving domains

    std::string  M_order;

    //! Discretization
    ADRStabilization M_stab_method;
    int M_semiImplicit;

private:

    //! To extract Mean Values at a given section z
    UInt M_computeMeanValuesPerSection; //! switch: 0 don't compute it, 1 compute
    UInt M_NbZSections;
    Real M_ToleranceSection;
    Real M_XSectionFrontier;
    Real M_ZSectionInit;
    Real M_ZSectionFinal;
    UInt M_NbPolygonEdges; //! number of edges of the polygon (in mesh) describing the circle

    DataStringList M_stabilization_list;
};


//
// IMPLEMENTATION
//


// Constructor


template <typename Mesh>
DataADR<Mesh>::
DataADR( const GetPot& dfile, const std::string& mesh_section, const std::string& time_section ) :
       DataMesh<Mesh>      ( dfile, mesh_section ),
       DataTime            ( dfile, time_section ),
       M_semiImplicit      ( 0 ),
       M_stabilization_list( (mesh_section + "/stabilization").data() )
       {
            setup(dfile);
       }



template <typename Mesh>
DataADR<Mesh>::
DataADR( const DataADR& dataADR ) :
    DataMesh<Mesh>               ( dataADR ),
    DataTime                     ( dataADR ),
    M_diffusivity                ( dataADR.M_diffusivity),
    M_react                      ( dataADR.M_react),
    M_verbose                    ( dataADR.M_verbose),
    M_factor                     ( dataADR.M_factor),
    M_order                      ( dataADR.M_order),
    M_stab_method                ( dataADR.M_stab_method)
{
}




template <typename Mesh>
void
DataADR<Mesh>::
setup(  const GetPot& dfile )
{
    M_stabilization_list.add( "ip",  ADR_IP_STABILIZATION, "interior penalty " );
    M_stabilization_list.add( "sd",  ADR_SD_STABILIZATION, "stream-line difussion" );
    M_stabilization_list.add( "none", ADR_NO_STABILIZATION,  "none (default)" );

    // physics
    M_diffusivity = dfile( "adr/physics/diffusivity", 1. );
    M_react = dfile( "adr/physics/react", 1. );

    //miscellaneous
    M_verbose = dfile( "adr/miscellaneous/verbose", 1 );
    _dump_init = dfile( "adr/miscellaneous/dump_init", getInitialTime() );
    _dump_period = dfile( "adr/miscellaneous/dump_period", 1 );
    M_factor = dfile( "adr/miscellaneous/factor", 0. );

    M_order = dfile( "adr/discretization/order", "P1");

    M_stab_method = ADRStabilization ( M_stabilization_list.value( dfile( "adr/discretization/stabilization", "none") ) );

    // IP needs boundary faces
    bool ipfaces =  ( M_stab_method == ADR_IP_STABILIZATION ) && (this->meshFaces() != "all" ) ;
    if ( ipfaces ) {
        ERROR_MSG("ERROR: IP requires boundary faces. Put mesh_faces = all in data file." ); }

    //mean values per section
    M_computeMeanValuesPerSection =
        dfile( "adr/valuespersection/computeMeanValuesPerSection", 0 );
    M_NbZSections =
        dfile( "adr/valuespersection/nb_z_section", 2 );
    M_ToleranceSection =
        dfile( "adr/valuespersection/tol_section", 2e-2 );
    M_XSectionFrontier =
        dfile( "adr/valuespersection/x_section_frontier", 0. );
    M_ZSectionInit =
        dfile( "adr/valuespersection/z_section_init", -1. );
    M_ZSectionFinal =
        dfile( "adr/valuespersection/z_section_final", 0. );
    M_NbPolygonEdges =
        dfile( "adr/valuespersection/nb_polygon_edges", 10 );
}

// Output
template <typename Mesh>
void DataADR<Mesh>::
showMe( std::ostream& c )
{
    // physics
    c << "\n*** Values for data [adr/miscellaneous]\n\n";
    c << "verbose     = " << M_verbose << std::endl;
    c << "initial time for writing solution  = " << _dump_init << std::endl;
    c << "number of time steps between two consecutive dumps of the solution = " << _dump_period << std::endl;
    c << "amplification factor = " << M_factor << std::endl;


    c << "\n*** Values for data [adr/discretization]\n\n";
    DataMesh<Mesh>::showMe( c );
    c << "\n*** Values for data [adr/time]\n\n";
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


    c << "\n*** Values for data [adr/valuespersection]\n\n";
    c << "computeMeanValuesPerSection (switch 0: don't compute, 1: compute)  = "
      << M_computeMeanValuesPerSection << std::endl;
    c << "nb_z_section        = " << M_NbZSections << std::endl;
    c << "tol_section         = " << M_ToleranceSection << std::endl;
    c << "x_section_frontier  = " << M_XSectionFrontier << std::endl;
    c << "z_section_init      = " << M_ZSectionInit << std::endl;
    c << "z_section_final     = " << M_ZSectionFinal << std::endl;
    c << "nb_polygon_edges    = " << M_NbPolygonEdges << std::endl;

}

////////////////////
// The viscosity

template <typename Mesh>
Real DataADR<Mesh>::
diffusivity() const
{
    return M_diffusivity;
}

template <typename Mesh>
Real DataADR<Mesh>::
reaction() const
{
    return M_react;
}
////////////////

// verbose variable
template <typename Mesh>
UInt DataADR<Mesh>::
verbose() const
{
    return M_verbose;
}

// Amplification factor
template <typename Mesh>
Real DataADR<Mesh>::
factor() const
{
    return M_factor;
}

template <typename Mesh>
std::string
DataADR<Mesh>::order() const
{
    return M_order;
}


// Stabilization method
template <typename Mesh>
ADRStabilization DataADR<Mesh>::stabilization() const
{
    return M_stab_method;
}


//! Mean Values per Sections
//! compute (0) or not (1) the mean values per sections
template <typename Mesh>
UInt DataADR<Mesh>::
computeMeanValuesPerSection() const
{
    return M_computeMeanValuesPerSection;
}
//! number of sections
template <typename Mesh>
UInt DataADR<Mesh>::
NbZSections() const
{
    return M_NbZSections;
}
//! tolerance for point proximity
template <typename Mesh>
Real DataADR<Mesh>::
ToleranceSection() const
{
    return M_ToleranceSection;
}
//! x (see NSHandler): point at the frontier for computation of the area
//! with a polygonal formula (x -> displacement of the boundary)
template <typename Mesh>
Real DataADR<Mesh>::
XSectionFrontier() const
{
    return M_XSectionFrontier;
}
//! lower section
template <typename Mesh>
Real DataADR<Mesh>::
ZSectionInit() const
{
    return M_ZSectionInit;
}
//! upper section
template <typename Mesh>
Real DataADR<Mesh>::
ZSectionFinal() const
{
    return M_ZSectionFinal;
}
//! number of edges of the polygon (in mesh) describing the circle
template <typename Mesh>
UInt DataADR<Mesh>::
NbPolygonEdges() const
{
    return M_NbPolygonEdges;
}

template <typename Mesh>
void DataADR<Mesh>::
setSemiImplicit(const int& SI)
{ M_semiImplicit = SI; }

template <typename Mesh>
int DataADR<Mesh>::
semiImplicit() const
 { return M_semiImplicit; }
}
#endif





