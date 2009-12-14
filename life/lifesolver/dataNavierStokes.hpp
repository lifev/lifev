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
  \file dataNavierStokes.hpp

  \version 1.0
  \date 01/2003
  \author M.A. Fernandez

  \version 1.25
  \date 09/2009
  \author Cristiano Malossi<cristiano.malossi@epfl.ch>

  \version 1.26
  \date 09/2009
  \author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  \brief File containing a class for handling NavierStokes data with GetPot

  The data is now able to store multiple fluids, but the use of the old interface (for only one fluid) is still available and fully compatible.

  The old way is to declare one fluid with its density and viscosity in [fluid/physics/density] and [fluid/physics/viscosity]. The set and get functions can then be used without precising which fluid is concerned. Internally, the informations are stored inside vectors, but with only one value (the 0 value: density is stored in M_density[0]).

  In the new way, one has to declare first of all in [fluid/physics/fluid_number] the number of fluids that will be used. Then every fluid is specified in a section [fluid/physics/fluid_k] where k is the number of the fluid (k from 0 to fluid_number-1, the same numeration is used internally). In this section, the density and the viscosity are declared.

  Example: To declare two fluids, one should use in the data file:
  [fluid]
     [./physics]
       fluid_number = 2

       [./fluid_0]
         density = 1;
	 viscosity = 1;

       [../fluid_1]
         density = 10;
	 viscosity = 100;

  Remark: in case both ways of declaring the fluids are used, the new one has the priority.
  Remark: do not use "fluid_number = 1" with the old way, it will not work properly.

*/
#ifndef _DATANAVIERSTOKES_H_
#define _DATANAVIERSTOKES_H_

#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifecore/dataString.hpp>
#include <life/lifecore/util_string.hpp>

#include <string>
#include <iostream>

namespace LifeV {

enum NSStabilization
{
    NO_STABILIZATION, //!< No stabilization
    IP_STABILIZATION, //!< Interior penalty
    SD_STABILIZATION  //!< Stream-line diffusion
};

//! DataNavierStokes - LifeV Base class which holds usual data for the NavierStokes equations solvers
/*
 *  @author M.A. Fernandez, Cristiano Malossi, Samuel Quinodoz
 *
 */
template <typename Mesh>
class DataNavierStokes : public DataMesh<Mesh>,
                         public DataTime
{
public:

	//! @name Constructors & Destructor
	//@{

    //! Constructor
	DataNavierStokes( 	const GetPot& dataFile,
                        const bool         verbose      = false,
						const std::string& mesh_section = "fluid/space_discretization",
						const std::string& time_section = "fluid/time_discretization");

	//! Copy constructor
	/*!
	 * \param dataNavierStokes - DataNavierStokes
	 */
	DataNavierStokes( const DataNavierStokes& dataNavierStokes );

	//@}



    //! @name Methods
    //@{

	//! Operator=
	/*!
	 * \param dataNavierStokes - DataNavierStokes
	 */
	DataNavierStokes& operator=( const DataNavierStokes& dataNavierStokes );

    //! Ouptut
    void showMe( std::ostream& output = std::cout );


    //! external setup

    void setup( const GetPot& dataFile );

    //@}



    //! @name Set functions
    //@{

  INLINE void density ( const Real& density, const UInt nfluid=0 )
  {
    ASSERT(nfluid< M_fluid_number,"Undeclared fluid");
    M_density[nfluid] = density;
  };

  INLINE void viscosity ( const Real& viscosity, const UInt nfluid=0 )
  {
    ASSERT(nfluid< M_fluid_number,"Undeclared fluid");
    M_viscosity[nfluid] = viscosity;
  };

    void setSemiImplicit       ( const bool SI )         {
                                                           M_semiImplicit = SI;
                                                           if ( M_semiImplicit )
                                                        	   setUseShapeDerivatives(false);
                                                         }

    void setUseShapeDerivatives( const bool SD )         { M_shapeDerivatives = SD; }

    //@}



    //! @name Get functions
    //@{

    INLINE UInt fluid_number() const { return M_fluid_number; };

    INLINE Real density(const UInt& n=0) const
    {
      ASSERT(n<M_fluid_number,"Undelared fluid");
      return M_density[n];
    };
    INLINE Real viscosity(const UInt& n=0) const
    {
      ASSERT(n<M_fluid_number,"Undelared fluid");
      return M_viscosity[n];
    };

    std::string     uOrder()                      const { return M_uOrder; }
    std::string     pOrder()                      const { return M_pOrder; }

    UInt            verbose()                     const { return M_verbose; }
    Real            dump_init()                   const { return M_dump_init; }
    UInt            dump_period()                 const { return M_dump_period; }
    Real            factor()                      const { return M_factor; }

    NSStabilization stabilization()               const { return M_stab_method; }

    bool            isSemiImplicit()              const { return M_semiImplicit; }
    bool            useShapeDerivatives()         const { return M_shapeDerivatives; }

    UInt            computeMeanValuesPerSection() const { return M_computeMeanValuesPerSection; }
    UInt            NbZSections()                 const { return M_NbZSections; }
    Real            ToleranceSection()            const { return M_ToleranceSection; }
    Real            XSectionFrontier()            const { return M_XSectionFrontier; }
    Real            ZSectionInit()                const { return M_ZSectionInit; }
    Real            ZSectionFinal()               const { return M_ZSectionFinal; }
    UInt            NbPolygonEdges()              const { return M_NbPolygonEdges; }

    //@}

protected:

    //! Physics
    UInt             M_fluid_number;
    std::vector<Real> M_density;
    std::vector<Real> M_viscosity;

    //! FE order
    std::string      M_uOrder;
    std::string      M_pOrder;

    //! Miscellaneous
    UInt             M_verbose;     // temporal output verbose
    Real             M_dump_init;   // time for starting the dumping of the results (Alex December 2003)
    UInt             M_dump_period; // frequency of the dumping (one dump after _dump_period time steps) (Alex December 2003)
    Real             M_factor;      // amplification factor for moving domains

    //! Discretization
    NSStabilization  M_stab_method;

private:

    //! To extract Mean Values at a given section z
    bool             M_semiImplicit;
    bool             M_shapeDerivatives;
    UInt             M_computeMeanValuesPerSection; //! switch: 0 don't compute it, 1 compute
    UInt             M_NbZSections;
    Real             M_ToleranceSection;
    Real             M_XSectionFrontier;
    Real             M_ZSectionInit;
    Real             M_ZSectionFinal;
    UInt             M_NbPolygonEdges; //! number of edges of the polygon (in mesh) describing the circle

    DataStringList   M_stabilization_list;
};



// ===================================================
//! Constructors
// ===================================================
template <typename Mesh>
DataNavierStokes<Mesh>::DataNavierStokes( const GetPot&      dataFile,
                                          const bool         verbose,
                                          const std::string& mesh_section,
                                          const std::string& time_section) :
        DataMesh<Mesh>                     ( dataFile, mesh_section, verbose ),
        DataTime                           ( dataFile, time_section ),
        M_density                          ( ),
        M_viscosity                        ( ),
        M_uOrder                           ( ),
        M_pOrder                           ( ),
        M_verbose                          ( ),
        M_dump_init                        ( ),
        M_dump_period                      ( ),
        M_factor                           ( ),
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
    setup( dataFile );
}

template <typename Mesh>
DataNavierStokes<Mesh>::DataNavierStokes( const DataNavierStokes& dataNavierStokes ) :
    DataMesh<Mesh>                     ( dataNavierStokes ),
    DataTime                           ( dataNavierStokes ),
    M_fluid_number                     ( dataNavierStokes.M_fluid_number ),
    M_density                          ( dataNavierStokes.M_density ),
    M_viscosity                        ( dataNavierStokes.M_viscosity ),
    M_uOrder                           ( dataNavierStokes.M_uOrder ),
    M_pOrder                           ( dataNavierStokes.M_pOrder ),
    M_verbose                          ( dataNavierStokes.M_verbose ),
    M_dump_init                        ( dataNavierStokes.M_dump_init ),
    M_dump_period                      ( dataNavierStokes.M_dump_period ),
    M_factor                           ( dataNavierStokes.M_factor ),
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
//! Methods
// ===================================================
template <typename Mesh>
DataNavierStokes<Mesh>&
DataNavierStokes<Mesh>::operator=( const DataNavierStokes& dataNavierStokes )
{
    if ( this != &dataNavierStokes )
    {
        DataMesh<Mesh>::operator=( dataNavierStokes );
        DataTime::operator=( dataNavierStokes );

	M_fluid_number                     = dataNavierStokes.M_fluid_number;
        M_density                          = dataNavierStokes.M_density;
        M_viscosity                        = dataNavierStokes.M_viscosity;
        M_uOrder                           = dataNavierStokes.M_uOrder;
        M_pOrder                           = dataNavierStokes.M_pOrder;
        M_verbose                          = dataNavierStokes.M_verbose;
        M_dump_init                        = dataNavierStokes.M_dump_init;
        M_dump_period                      = dataNavierStokes.M_dump_period;
        M_factor                           = dataNavierStokes.M_factor;
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

template <typename Mesh>
void
DataNavierStokes<Mesh>::setup( const GetPot& dataFile )
{
    M_stabilization_list.add( "ip", IP_STABILIZATION,   "interior penalty " );
    M_stabilization_list.add( "sd", SD_STABILIZATION,   "stream-line diffusion" );
    M_stabilization_list.add( "none", NO_STABILIZATION, "none (default)" );

    // Physics
    UInt temp = dataFile("fluid/physics/fluid_number", 0 );

    if (temp == 0) // Old fashion of declaring fluids
      {
	M_fluid_number = 1;
	M_density.push_back( dataFile( "fluid/physics/density", 1. ) );
	M_viscosity.push_back ( dataFile( "fluid/physics/viscosity", 1. ) );
      }
    else   // New fashion of declaring fluids
      {
	M_fluid_number = temp;
	M_density = std::vector<Real>(temp,0);
	M_viscosity = std::vector<Real>(temp,0);

	for (UInt iter_fluid(0); iter_fluid<temp; ++iter_fluid)
	  {
	    // build the section name
	    std::string iter_fluid_section("fluid/physics/fluid_");
	    iter_fluid_section += number2string(iter_fluid);

	    // Read the quantities
	    M_density[iter_fluid]= dataFile((iter_fluid_section+"/density").c_str() , 1.0);
	    M_viscosity[iter_fluid]= dataFile((iter_fluid_section+"/viscosity").c_str() , 1.0);
	  };
      };

    // FE Order
    M_uOrder      = dataFile( "fluid/space_discretization/vel_order", "P1");
    M_pOrder      = dataFile( "fluid/space_discretization/press_order", "P1");

    // Miscellaneous
    M_verbose      = dataFile( "fluid/miscellaneous/verbose", 1 );
    M_dump_init    = dataFile( "fluid/miscellaneous/dump_init", getInitialTime() );
    M_dump_period  = dataFile( "fluid/miscellaneous/dump_period", 1 );
    M_factor       = dataFile( "fluid/miscellaneous/factor", 0. );

    M_stab_method  = NSStabilization ( M_stabilization_list.value( dataFile( "fluid/space_discretization/stabilization", "none") ) );

    // Semi-implicit and shape derivatives
    M_semiImplicit     = dataFile( "fluid/semiImplicit", false ) ;
    M_shapeDerivatives = dataFile( "fluid/useShapeDerivatives", false ) ;
    setSemiImplicit( M_semiImplicit );

    // Mean values per section
    M_computeMeanValuesPerSection = dataFile( "fluid/valuespersection/computeMeanValuesPerSection", 0 );
    M_NbZSections      = dataFile( "fluid/valuespersection/nb_z_section", 2 );
    M_ToleranceSection = dataFile( "fluid/valuespersection/tol_section", 2e-2 );
    M_XSectionFrontier = dataFile( "fluid/valuespersection/x_section_frontier", 0. );
    M_ZSectionInit     = dataFile( "fluid/valuespersection/z_section_init", -1. );
    M_ZSectionFinal    = dataFile( "fluid/valuespersection/z_section_final", 0. );
    M_NbPolygonEdges   = dataFile( "fluid/valuespersection/nb_polygon_edges", 10 );
}

template <typename Mesh>
void
DataNavierStokes<Mesh>::showMe( std::ostream& output )
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
      };
    };
    output << "\n*** Values for data [fluid/miscellaneous]\n\n";
    output << "verbose     = " << M_verbose << std::endl;
    output << "initial time for writing solution  = " << M_dump_init << std::endl;
    output << "number of time steps between two consecutive dumps of the solution = " << M_dump_period << std::endl;
    output << "amplification factor = " << M_factor << std::endl;


    output << "\n*** Values for data [fluid/space_discretization]\n\n";
    DataMesh<Mesh>::showMe( output );
    output << "\n*** Values for data [fluid/time_discretization]\n\n";
    DataTime::showMe( output );

    output << "stabilization = ";
    switch( M_stab_method )
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

} // end namespace LifeV

#endif
