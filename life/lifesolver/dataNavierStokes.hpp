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
//#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifecore/dataString.hpp>
#include <life/lifecore/util_string.hpp>

#include <boost/shared_ptr.hpp>

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
/*!
 *  @author M.A. Fernandez, Cristiano Malossi, Samuel Quinodoz
 */
//

class DataNavierStokes
{
public:

    //! @name Type definitions
    //@{

    typedef DataTime                                                  Time_Type;
    typedef boost::shared_ptr< Time_Type >                            Time_ptrType;

    //    typedef DataMesh                                            Mesh_Type;
    //typedef boost::shared_ptr< Mesh_Type >                            Mesh_ptrType;

    //@}


	//! @name Constructors & Destructor
	//@{

    //! Empty Constructor
    DataNavierStokes();

	//! Copy constructor
	/*!
	 * @param dataNavierStokes DataNavierStokes
	 */
	DataNavierStokes( const DataNavierStokes& dataNavierStokes );

	//@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param dataNavierStokes DataNavierStokes
     */
    DataNavierStokes& operator=( const DataNavierStokes& dataNavierStokes );

    //@}


	//! @name Methods
    //@{

    //! Read the dataFile and set all the quantities
    /*!
     * @param dataFile data file
     * @param section section of the file
     */
    void setup( const GetPot& dataFile, const std::string& section = "fluid" );

    //! Display the values
    void showMe( std::ostream& output = std::cout ) const;

    //@}



    //! @name Set methods
    //@{

    //! Set data time container
    /*!
     * @param DataTime shared_ptr to dataTime container
     */
    inline void setDataTime( const Time_ptrType DataTime ) { M_time = DataTime; }

    //! Set mesh container
    /*!
     * @param DataMesh shared_ptr to dataMesh container
     */
    //inline void setDataMesh( const Mesh_ptrType DataMesh ) { M_mesh = DataMesh; }

    inline void density ( const Real& density, const UInt nfluid=0 )
    {
        ASSERT(nfluid< M_fluid_number,"Undeclared fluid");
        M_density[nfluid] = density;
    }

    inline void viscosity ( const Real& viscosity, const UInt nfluid=0 )
    {
        ASSERT(nfluid< M_fluid_number,"Undeclared fluid");
        M_viscosity[nfluid] = viscosity;
    }

    void setStokes             ( const bool Stokes )     { M_Stokes = Stokes; }

    void setSemiImplicit       ( const bool SI )         {
                                                           M_semiImplicit = SI;
                                                           if ( M_semiImplicit )
                                                        	   setUseShapeDerivatives(false);
                                                         }

    void setUseShapeDerivatives( const bool SD )         { M_shapeDerivatives = SD; }

    //@}



    //! @name Get methods
    //@{

    //! Get data time container
    /*!
     * @return shared_ptr to dataTime container
     */
    inline Time_ptrType dataTime( void ) const { return M_time; }

    //! Get mesh container
    /*!
     * @return shared_ptr to dataMesh container
     */
    //inline Mesh_ptrType dataMesh( void ) const { return M_mesh; }

    inline UInt fluid_number() const { return M_fluid_number; };

    inline Real density(const UInt& n=0) const
    {
        ASSERT(n<M_fluid_number,"Undelared fluid");
        return M_density[n];
    }

    inline Real viscosity(const UInt& n=0) const
    {
        ASSERT(n<M_fluid_number,"Undelared fluid");
        return M_viscosity[n];
    }

    std::string     uOrder()                      const { return M_uOrder; }
    std::string     pOrder()                      const { return M_pOrder; }

    UInt            verbose()                     const { return M_verbose; }
    Real            dump_init()                   const { return M_dump_init; }
    UInt            dump_period()                 const { return M_dump_period; }
    Real            factor()                      const { return M_factor; }

    NSStabilization stabilization()               const { return M_stab_method; }

    inline bool     Stokes()                      const { return M_Stokes; }
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

    //! Data containers for time and mesh
    Time_ptrType      M_time;
    //Mesh_ptrType      M_mesh;

    //! Physics
    UInt              M_fluid_number;
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
    bool             M_Stokes;      // true: Stokes problem; false: Navier-Stokes problem

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



} // end namespace LifeV

#endif
