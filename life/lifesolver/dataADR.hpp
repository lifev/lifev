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
//#include <life/lifemesh/dataMesh.hpp>
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


//template <typename Mesh>
class DataADR:
        //        public DataMesh<Mesh>,
        public DataTime
{
public:

    //! Constructor

    DataADR();

    DataADR( const DataADR& dataNavierStokes );

    //! Ouptut
    void showMe( std::ostream& c = std::cout );


    //! external setup

    void setup( const GetPot& dfile, const std::string& section = "adr" );

    //! End time
    Real density()     const;
    Real diffusivity() const;
    Real reaction()    const;

    UInt verbose()     const;
    Real dump_init()   const;
    UInt dump_period() const;
    Real factor()      const;

    std::string order() const;

    ADRStabilization stabilization() const;

    //! a way to obtain the Mean Flux, Mean Pressure and
    //! Mean Area at a section defined by z=z_data.
    // UInt computeMeanValuesPerSection() const;
    // UInt NbZSections() const;
    // Real ToleranceSection() const;
    // Real XSectionFrontier() const;
    // Real ZSectionInit() const;
    // Real ZSectionFinal() const;
    // UInt NbPolygonEdges() const;
    // int  semiImplicit() const;
    // void setSemiImplicit(const int& SI);

private:
    //! Physics

    Real M_diffusivity; // Diffusivity
    Real M_react; // Reaction coefficient

    //! Miscellaneous
    UInt M_verbose; // temporal output verbose

    // Real _dump_init; // time for starting the dumping of the results (Alex December 2003)
    // UInt _dump_period; // frequency of the dumping (one dump after _dump_period time steps) (Alex December 2003
    //    )
    Real M_factor; // amplification factor for moving domains

    std::string  M_order;

    //! Discretization
    ADRStabilization M_stab_method;
    int M_semiImplicit;

//private:

    //! To extract Mean Values at a given section z
    // UInt M_computeMeanValuesPerSection; //! switch: 0 don't compute it, 1 compute
    // UInt M_NbZSections;
    // Real M_ToleranceSection;
    // Real M_XSectionFrontier;
    // Real M_ZSectionInit;
    // Real M_ZSectionFinal;
    // UInt M_NbPolygonEdges; //! number of edges of the polygon (in mesh) describing the circle

    DataStringList M_stabilization_list;
};


}
#endif





