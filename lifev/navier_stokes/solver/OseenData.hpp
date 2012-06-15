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
    @brief File containing a class for handling Navier-Stokes data with GetPot

    @author M.A. Fernandez
            Christiano Malossi <cristiano.malossi@epfl.ch>
            Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @contributor Alexis Aposporidis <aapospo@emory.edu>
    @maintainer

    @date 01-09-2009

 */


#ifndef OSEENDATA_H
#define OSEENDATA_H

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/StringData.hpp>
#include <lifev/core/util/StringUtility.hpp>
#include <lifev/core/fem/TimeData.hpp>
#include <lifev/core/fem/TimeAdvanceData.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <iostream>

namespace LifeV
{

enum NSStabilization
{
    NO_STABILIZATION, //!< No stabilization
    IP_STABILIZATION, //!< Interior penalty
    SD_STABILIZATION  //!< Stream-line diffusion
};


//! @class OseenData
//! OseenData - LifeV Base class which holds usual data for the NavierStokes equations solvers
/*!
 *  @author M.A. Fernandez, Cristiano Malossi, Samuel Quinodoz
 *
 *  The data is now able to store multiple fluids, but the use of the old interface (for only one fluid) is still available and
 *  fully compatible.
 *  The old way is to declare one fluid with its density and viscosity in [fluid/physics/density] and [fluid/physics/viscosity].
 *  The set and get functions can then be used without specifying which fluid is concerned.
 *  Internally, the information is stored in vectors, but with only one value (the 0 value: density is stored in M_density[0]).
 *  In the new way, one has to declare first of all in [fluid/physics/fluid_number] the number of fluids that will be used.
 *  Then every fluid is specified in a section [fluid/physics/fluid_k] where k is the number of the fluid (k from 0 to fluid_number-1,
 *  the same enumeration is used internally). In this section, the density and the viscosity are declared.
 *
 *  Example: To declare two fluids, one should use in the data file:
 *  [fluid]
 *     [./physics]
 *       fluid_number = 2
 *
 *     [./fluid_0]
 *       density = 1;
 *         viscosity = 1;
 *
 *     [../fluid_1]
 *       density = 10;
 *         viscosity = 100;
 *
 *  Remark: in case both ways of declaring the fluids are used, the new one has priority.
 *  Remark: do not use "fluid_number = 1" with the old way, it will not work properly.
 *
 */

class OseenData
{

public:

    //! @name Public Types
    //@{

    typedef TimeData                               time_Type;
    typedef boost::shared_ptr<time_Type>           timePtr_Type;

    typedef TimeAdvanceData                        timeAdvance_Type;
    typedef boost::shared_ptr<timeAdvance_Type>    timeAdvancePtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    OseenData();

    //! Copy constructor
    /*!
     * @param oseenData OseenData
     */
    OseenData( const OseenData& oseenData );

    //! Virtual destructor
    virtual ~OseenData() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param oseenData OseenData
     */
    OseenData& operator=( const OseenData& oseenData );

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
     * @param TimeData shared_ptr to TimeData container
     */
    void setTimeData( const timePtr_Type timeData ) { M_time = timeData; }

    //! Set data time advance container
    /*!
     * @param timeAdvanceData shared_ptr to TimeAdvanceData container
     */
    void setTimeAdvanceData( const timeAdvancePtr_Type timeAdvanceData ) { M_timeAdvance = timeAdvanceData; }

    //! Set the density for the specified fluid
    /*!
     * @param density
     * @param nfluid the fluid number
     */
    void setDensity ( const Real& density, const UInt nfluid=0 )
    {
        ASSERT(nfluid< M_fluidNumber,"Undeclared fluid");
        M_density[nfluid] = density;
    }

    //! Set the viscosity of the fluid
    /*!
     * @param viscosity
     * @param nfluid the fluid number
     */
    void setViscosity ( const Real& viscosity, const UInt nfluid=0 )
    {
        ASSERT(nfluid< M_fluidNumber,"Undeclared fluid");
        M_viscosity[nfluid] = viscosity;
    }

    //! Set this instance of OseenData to either a Stokes or a Navier-Stokes problem
    /*!
     * @param Stokes a boolean that is "true" for a Stokes problem and "false" for a Navier-Stokes
     *                problem
     */
    void setStokes( const bool stokes ) { M_stokes = stokes; }

    //! Set the flag for the semi-implicit treatment of the shape derivatives in FSI simulations
    /*!
     * @param SI the flag for semi implicit treatment
     */
    void setSemiImplicit( const bool SI )
    {
        M_semiImplicit = SI;
        if ( M_semiImplicit )
            setUseShapeDerivatives(false);
    }

    //! Set the flag for using shape derivatives
    /*!
     * @param SD the flag for using shape derivatives
     */
    void setUseShapeDerivatives( const bool SD ) { M_shapeDerivatives = SD; }

    //@}



    //! @name Get methods
    //@{

    //! Get data time container
    /*!
     * @return shared_ptr to TimeData container
     */
    timePtr_Type dataTime() const { return M_time; }

    //! Get data time advance container
    /*!
     * @return shared_ptr to TimeAdvanceData container
     */
    timeAdvancePtr_Type dataTimeAdvance() const { return M_timeAdvance; }

    //! Get the number of the fluid
    /*!
     * @return M_fluidNumber the number of the current fluid
     */
    const UInt& fluidNumber() const { return M_fluidNumber; };

    //! Get the density of the fluid
    /*!
     * @param n the fluid number
     * @return M_density the density of fluid n
     */
    const Real& density(const UInt& n=0) const
    {
        ASSERT(n<M_fluidNumber,"Undeclared fluid");
        return M_density[n];
    }

    //! Get the viscosity of the fluid
    /*!
     * @param n the fluid number
     * @return M_viscosity the viscosity of the fluid n
     */
    const Real& viscosity(const UInt& n=0) const
    {
        ASSERT(n<M_fluidNumber,"Undeclared fluid");
        return M_viscosity[n];
    }

    //! Get the order of the finite elements used for velocity
    /*!
     * @return M_uOrder a string specifying the order of finite elements for velocity
     */
    std::string uOrder() const { return M_uOrder; }

    //! Get the order of the finite elements used for pressure
    /*!
     * @return M_pOrder a string specifying the order of finite elements for pressure
     */
    std::string pOrder() const { return M_pOrder; }

    //! Temporal output verbose
    /*!
     * @return M_verbose
     */
    UInt verbose() const { return M_verbose; }

    //! Dumping of the results
    /*!
     * @return M_dumpInit the time for dumping of the results
     */
    Real dumpInit() const { return M_dumpInit; }

    //! Get the frequency of the dumping
    /*!
     * @return M_dumpPeriod number of time steps after which one dump is performed
     */
    UInt dumpPeriod() const { return M_dumpPeriod; }

    //! Get the amplification factor
    /*!
     * @return M_factor The amplification factor
     */
    Real factor() const { return M_factor; }

    //! Get the stabilization method
    /*!
     * @return M_stabMethod The method used for stabilizing
     */
    NSStabilization stabilization() const { return M_stabMethod; }

    //!
    /*! find out if this is a Stokes or Navier-Stokes problem
     * @return M_stokes Boolean that is "true" for a Stokes and "false" for a Navier-Stokes
     *                  problem
     */
    bool isStokes() const { return M_stokes; }

    //! Find out if a semi-implicit scheme is used
    /*!
     * @return M_semiImplicit "true" if a semi-implicit scheme is used, "false" otherwise
     */
    bool isSemiImplicit() const { return M_semiImplicit; }

    //!Get the flag for using shape derivatives
    /*!
     *@return M_shapeDerivatives Flag for shape derivatives
     *
     */
    bool useShapeDerivatives() const { return M_shapeDerivatives; }


    //!Get the flag for considering implicitly the fluid domain (when it is moving, e.g. ALE)
    /*!
     *@return M_domainVelImplicit Flag for shape derivatives
     *
     */
    bool domainVelImplicit() const { return M_domainVelImplicit; }


    //!Get the flag for considering implicitly the fluid convective term
    /*!
     *@return M_convectiveImplicit Flag for shape derivatives
     *
     */
    bool convectiveImplicit() const { return M_convectiveImplicit; }

    //! Get the number of mean valuNes per section
    /*!
     * @return M_computeMeanValuesPerSection number of mean values
     */
    UInt computeMeanValuesPerSection() const { return M_computeMeanValuesPerSection; }

    //! Get the number of NBZ-sections
    /*!
     * @return M_NbZSections Number of NBZ-sections
     */
    UInt nbZSections() const { return M_NbZSections; }

    //! Tolerance section
    /*!
     * @return M_ToleranceSection The tolerance section
     */
    Real toleranceSection() const { return M_ToleranceSection; }

    //! X-Section frontier
    /*!
     * @return M_XSectionFrontier The x-section frontier
     */
    Real xSectionFrontier() const { return M_XSectionFrontier; }

    //! Z section init
    /*!
     * @return M_ZSectionInit The initial z-section
     */
    Real zSectionInit() const { return M_ZSectionInit; }

    //! Z section final
    /*!
     * @return M_ZSectionFinal The final z-section
     */
    Real zSectionFinal() const { return M_ZSectionFinal; }

    //! Number of edges of the polygon (in the mesh) describing the circle
    /*!
     * @return M_NbPolygonEdges The number of polygon edges
     */
    UInt nbPolygonEdges() const { return M_NbPolygonEdges; }

    //! Returns wether the formulation of the momentum conservation equation is written in conservative form or not.
    /*!
     * @return M_conservativeFormulation the output flag
     */
    bool conservativeFormulation() const { return M_conservativeFormulation; }

    //@}

protected:

    //! Data containers for time and mesh
    timePtr_Type        M_time;
    timeAdvancePtr_Type M_timeAdvance;

    //! @name Physics
    //@{

    //! number of this fluid
    UInt                M_fluidNumber;

    //! density of each fluid
    std::vector<Real>   M_density;

    //! viscosity of each fluid
    std::vector<Real>   M_viscosity;

    //@}


    //! @name FE order
    //@{

    //! order of finite elements for velocity
    std::string         M_uOrder;

    //! order of finite elements for pressure
    std::string         M_pOrder;

    //@}


    //! @name Miscellaneous
    //@{

    //! temporal output verbose
    UInt                M_verbose;

    //! time for starting the dumping of the results (Alex December 2003)
    Real                M_dumpInit;

    //! frequency of the dumping (one dump after _dump_period time steps) (Alex December 2003)
    UInt                M_dumpPeriod;

    //! amplification factor for moving domains
    Real                M_factor;

    //! true: Stokes problem; false: Navier-Stokes problem
    bool                M_stokes;

    //@}

    //! @name Discretization
    //@{

    //! stabilization method
    NSStabilization     M_stabMethod;

    //@}

private:

    //! To extract Mean Values at a given section z
    bool                M_semiImplicit;
    bool                M_shapeDerivatives;
    bool                M_domainVelImplicit;
    bool                M_convectiveImplicit;
    UInt                M_computeMeanValuesPerSection; // switch: 0 don't compute it, 1 compute
    UInt                M_NbZSections;
    Real                M_ToleranceSection;
    Real                M_XSectionFrontier;
    Real                M_ZSectionInit;
    Real                M_ZSectionFinal;
    UInt                M_NbPolygonEdges; // number of edges of the polygon (in mesh) describing the circle

    StringDataList      M_stabilizationList;
    bool                M_conservativeFormulation;
};



} // end namespace LifeV

#endif /* OSEENDATA_H */
