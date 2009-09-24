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

  \brief File containing a class for handling NavierStokes data with GetPot

*/
#ifndef _DATANAVIERSTOKES_H_
#define _DATANAVIERSTOKES_H_

#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifecore/dataString.hpp>

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
 *  @author M.A. Fernandez, Cristiano Malossi
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
						const std::string& mesh_section= "fluid/space_discretization",
						const std::string& time_section= "fluid/time_discretization" );

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

    void density               ( const Real& density )   { M_density = density; }

    void viscosity             ( const Real& viscosity ) { M_viscosity = viscosity; }

    void setSemiImplicit       ( const bool SI )         {
                                                           M_semiImplicit = SI;
                                                           if ( M_semiImplicit )
                                                        	   setUseShapeDerivatives(false);
                                                         }

    void setUseShapeDerivatives( const bool SD )         { M_shapeDerivatives = SD; }

    //@}



    //! @name Get functions
    //@{

    Real            density()                     const { return M_density; }
    Real            viscosity()                   const { return M_viscosity; }

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
    Real             M_density;   // density
    Real             M_viscosity; // viscosity

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
                                          const std::string& mesh_section,
                                          const std::string& time_section ) :
    DataMesh<Mesh>                     ( dataFile, mesh_section ),
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
    M_density     = dataFile( "fluid/physics/density", 1. );
    M_viscosity   = dataFile( "fluid/physics/viscosity", 1. );

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
    M_semiImplicit     = dataFile( "problem/semiImplicit", false ) ;
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
    output << "\n*** Values for data [fluid/physics]\n\n";
    output << "density     = " << M_density << std::endl;
    output << "viscosity   = " << M_viscosity << std::endl;

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
