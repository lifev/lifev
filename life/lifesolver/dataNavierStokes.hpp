/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
#ifndef _DATANAVIERSTOKES_H_
#define _DATANAVIERSTOKES_H_
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
typedef enum NSStabilization
{
    NO_STABILIZATION,     //!< No stabilization
    IP_STABILIZATION,       //!< Interior penalty
    SD_STABILIZATION        //!< Stream-line diffusion
};

/*!
  \class DataNavierStokes

  Base class which holds usual data for the NavierStokes equations solvers

*/
template <typename Mesh>
class DataNavierStokes:
        public DataMesh<Mesh>,
        public DataTime
{
public:

    //! Constructor

    DataNavierStokes( const GetPot& dfile );

    DataNavierStokes( const DataNavierStokes& dataNavierStokes );

    //! Ouptut
    void showMe( std::ostream& c = std::cout );


    //! external setup

    void setup( const GetPot& dfile );

    //! End time
    Real density() const;
    Real viscosity() const;
    Real inittime() const;
    Real endtime() const;

    UInt verbose() const;
    Real dump_init() const;
    UInt dump_period() const;
    Real factor() const;

    std::string uOrder() const;
    std::string pOrder() const;

    NSStabilization stabilization() const;

    //! a way to obtain the Mean Flux, Mean Pressure and
    //! Mean Area at a section defined by z=z_data.
    UInt computeMeanValuesPerSection() const;
    UInt NbZSections() const;
    Real ToleranceSection() const;
    Real XSectionFrontier() const;
    Real ZSectionInit() const;
    Real ZSectionFinal() const;
    UInt NbPolygonEdges() const;

protected:
    //! Physics
    Real _rho; // density
    Real _mu; // viscosity
    Real _inittime; // initial time (Alex December 2003)
    Real _endtime; // end time

    //! Miscellaneous
    UInt _verbose; // temporal output verbose
    Real _dump_init; // time for starting the dumping of the results (Alex December 2003)
    UInt _dump_period; // frequency of the dumping (one dump after _dump_period time steps) (Alex December 2003)
    Real M_factor; // amplification factor for moving domains

    std::string  M_uOrder;
    std::string  M_pOrder;

    //! Discretization
    NSStabilization M_stab_method;

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
DataNavierStokes<Mesh>::
DataNavierStokes( const GetPot& dfile ) :
    DataMesh<Mesh>( dfile, "fluid/discretization" ),
    DataTime( dfile, "fluid/discretization" ),
    M_stabilization_list( "fluid/discretization/stabilization" )
{
    setup(dfile);
}

template <typename Mesh>
DataNavierStokes<Mesh>::
DataNavierStokes( const DataNavierStokes& dataNavierStokes ) :
    DataMesh<Mesh>               ( dataNavierStokes ),
    DataTime                     ( dataNavierStokes ),
    _rho                         (dataNavierStokes._rho),
    _mu                          (dataNavierStokes._mu),
    _inittime                    (dataNavierStokes._inittime),
    _endtime                     (dataNavierStokes._endtime),
    _verbose                     (dataNavierStokes._verbose),
    _dump_init                   (dataNavierStokes._dump_init),
    _dump_period                 (dataNavierStokes._dump_period),
    M_factor                     (dataNavierStokes.M_factor),
    M_uOrder                     (dataNavierStokes.M_uOrder),
    M_pOrder                     (dataNavierStokes.M_pOrder),
    M_stab_method                (dataNavierStokes.M_stab_method),
    M_computeMeanValuesPerSection(dataNavierStokes.M_computeMeanValuesPerSection),
    M_NbZSections                (dataNavierStokes.M_NbZSections),
    M_ToleranceSection           (dataNavierStokes.M_ToleranceSection),
    M_XSectionFrontier           (dataNavierStokes.M_XSectionFrontier),
    M_ZSectionInit               (dataNavierStokes.M_ZSectionInit),
    M_ZSectionFinal              (dataNavierStokes.M_ZSectionFinal),
    M_NbPolygonEdges             (dataNavierStokes.M_NbPolygonEdges),
    M_stabilization_list         (dataNavierStokes.M_stabilization_list)
{
}




template <typename Mesh>
void
DataNavierStokes<Mesh>::
setup(  const GetPot& dfile )
{
    M_stabilization_list.add( "ip", IP_STABILIZATION, "interior penalty " );
    M_stabilization_list.add( "sd", SD_STABILIZATION, "stream-line difussion" );
    M_stabilization_list.add( "none", NO_STABILIZATION,  "none (default)" );


    // physics
    _rho = dfile( "fluid/physics/density", 1. );
    _mu = dfile( "fluid/physics/viscosity", 1. );
    _inittime = dfile( "fluid/physics/inittime", 0. );
    _endtime = dfile( "fluid/physics/endtime", 1. );

    //miscellaneous
    _verbose = dfile( "fluid/miscellaneous/verbose", 1 );
    _dump_init = dfile( "fluid/miscellaneous/dump_init", _inittime );
    _dump_period = dfile( "fluid/miscellaneous/dump_period", 1 );
    M_factor = dfile( "fluid/miscellaneous/factor", 0. );

    M_uOrder = dfile( "fluid/discretization/vel_order", "P1");
    M_pOrder = dfile( "fluid/discretization/press_order", "P1");

    M_stab_method = NSStabilization ( M_stabilization_list.value( dfile( "fluid/discretization/stabilization", "none") ) );

    // IP needs boundary faces
    bool ipfaces =  ( M_stab_method == IP_STABILIZATION ) && (this->meshFaces() != "all" ) ;
    if ( ipfaces )
        ERROR_MSG("ERROR: IP requires boundary faces. Put mesh_faces = all in data file." );

    //mean values per section
    M_computeMeanValuesPerSection =
        dfile( "fluid/valuespersection/computeMeanValuesPerSection", 0 );
    M_NbZSections =
        dfile( "fluid/valuespersection/nb_z_section", 2 );
    M_ToleranceSection =
        dfile( "fluid/valuespersection/tol_section", 2e-2 );
    M_XSectionFrontier =
        dfile( "fluid/valuespersection/x_section_frontier", 0. );
    M_ZSectionInit =
        dfile( "fluid/valuespersection/z_section_init", -1. );
    M_ZSectionFinal =
        dfile( "fluid/valuespersection/z_section_final", 0. );
    M_NbPolygonEdges =
        dfile( "fluid/valuespersection/nb_polygon_edges", 10 );
}

// Output
template <typename Mesh>
void DataNavierStokes<Mesh>::
showMe( std::ostream& c )
{
    // physics
    c << "\n*** Values for data [fluid/physics]\n\n";
    c << "density   = " << _rho << std::endl;
    c << "viscosity = " << _mu << std::endl;
    c << "initial time = " << _inittime << std::endl;
    c << "endtime   = " << _endtime << std::endl;

    c << "\n*** Values for data [fluid/miscellaneous]\n\n";
    c << "verbose   = " << _verbose << std::endl;
    c << "initial time for writing solution  = " << _dump_init << std::endl;
    c << "number of time steps between two consecutive dumps of the solution = " << _dump_period << std::endl;
    c << "amplification factor = " << M_factor << std::endl;


    c << "\n*** Values for data [fluid/discretization]\n\n";
    DataMesh<Mesh>::showMe( c );
    DataTime::showMe( c );
    c << "stabilization = ";
    switch( M_stab_method )
    {
        case NO_STABILIZATION:
            c << "none" ;
            break;
        case IP_STABILIZATION:
            c << "ip" ;
            break;
        case SD_STABILIZATION:
            c << "sd" ;
            break;
    }
    c << std::endl;


    c << "\n*** Values for data [fluid/valuespersection]\n\n";
    c << "computeMeanValuesPerSection (switch 0: don't compute, 1: compute)  = "
      << M_computeMeanValuesPerSection << std::endl;
    c << "nb_z_section  = " << M_NbZSections << std::endl;
    c << "tol_section  = " << M_ToleranceSection << std::endl;
    c << "x_section_frontier  = " << M_XSectionFrontier << std::endl;
    c << "z_section_init  = " << M_ZSectionInit << std::endl;
    c << "z_section_final  = " << M_ZSectionFinal << std::endl;
    c << "nb_polygon_edges  = " << M_NbPolygonEdges << std::endl;

}
////////////////////
// The density
template <typename Mesh>
Real DataNavierStokes<Mesh>::
density() const
{
    return _rho;
}

// The viscosity
template <typename Mesh>
Real DataNavierStokes<Mesh>::
viscosity() const
{
    return _mu;
}


// The initial time
template <typename Mesh>
Real DataNavierStokes<Mesh>::
inittime() const
{
    return _inittime;
}

// The end time
template <typename Mesh>
Real DataNavierStokes<Mesh>::
endtime() const
{
    return _endtime;
}
////////////////

// verbose variable
template <typename Mesh>
UInt DataNavierStokes<Mesh>::
verbose() const
{
    return _verbose;
}
// Dumping start
template <typename Mesh>
Real DataNavierStokes<Mesh>::
dump_init() const
{
    return _dump_init;
}
// Period of dumping
template <typename Mesh>
UInt DataNavierStokes<Mesh>::
dump_period() const
{
    return _dump_period;
}

// Amplification factor
template <typename Mesh>
Real DataNavierStokes<Mesh>::
factor() const
{
    return M_factor;
}

template <typename Mesh>
std::string
DataNavierStokes<Mesh>::uOrder() const
{
    return M_uOrder;
}

template <typename Mesh>
std::string DataNavierStokes<Mesh>::pOrder() const
{
    return M_pOrder;
}


// Stabilization method
template <typename Mesh>
NSStabilization DataNavierStokes<Mesh>::stabilization() const
{
    return M_stab_method;
}


//! Mean Values per Sections
//! compute (0) or not (1) the mean values per sections
template <typename Mesh>
UInt DataNavierStokes<Mesh>::
computeMeanValuesPerSection() const
{
    return M_computeMeanValuesPerSection;
}
//! number of sections
template <typename Mesh>
UInt DataNavierStokes<Mesh>::
NbZSections() const
{
    return M_NbZSections;
}
//! tolerance for point proximity
template <typename Mesh>
Real DataNavierStokes<Mesh>::
ToleranceSection() const
{
    return M_ToleranceSection;
}
//! x (see NSHandler): point at the frontier for computation of the area
//! with a polygonal formula (x -> displacement of the boundary)
template <typename Mesh>
Real DataNavierStokes<Mesh>::
XSectionFrontier() const
{
    return M_XSectionFrontier;
}
//! lower section
template <typename Mesh>
Real DataNavierStokes<Mesh>::
ZSectionInit() const
{
    return M_ZSectionInit;
}
//! upper section
template <typename Mesh>
Real DataNavierStokes<Mesh>::
ZSectionFinal() const
{
    return M_ZSectionFinal;
}
//! number of edges of the polygon (in mesh) describing the circle
template <typename Mesh>
UInt DataNavierStokes<Mesh>::
NbPolygonEdges() const
{
    return M_NbPolygonEdges;
}


}
#endif





