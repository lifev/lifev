/*
  This file is part of the LifeV library
  Copyright (C) 2010 EPFL, INRIA, Politecnico di Milano and Emory University

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
  \file dataADR.hpp
  \date 08/2010
  \version 2.0

  \brief File containing a class for handling data and parameters
         for advection-diffusion-reaction problems

*/
#ifndef _DATAADR_H_
#define _DATAADR_H_

#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifecore/dataString.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{

/*!
  \typedef The kind of stabilization (if any) applied
           to the numerical scheme
*/
enum ADRStabilization
{
    ADR_NO_STABILIZATION,       //!< No stabilization
    ADR_IP_STABILIZATION,       //!< Interior penalty
    ADR_SD_STABILIZATION        //!< Stream-line diffusion
};

/*!
  \class DataADR

  Base class which holds usual data for the ADR equations solvers

*/
class DataADR
{
public:

    //! @name Type definitions
    //@{
    //! The time data handler
    /*!
       The class DataADR "has a" time handler object, to manage time dependent
       problems.
     */
    typedef DataTime                              dataTime_type;
    typedef dataTime_type*                        dataTime_rawptr_type;
    typedef boost::shared_ptr<dataTime_type>      dataTime_ptr_type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    DataADR();

    //! Copy constructor
    /*!
     * @param dataADR an object of type DataADR
     */
    DataADR( const DataADR& dataADR );

    //@}

    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param dataADR an object of type DataADR
     */
    DataADR& operator=( const DataADR& dataADR );

    //@}

    //! @name Methods
    //@{

    //! Read the dataFile and set all the internal quantities
    /*!
     * @param dataFile data file
     * @param section section of the file
     */
    void setup( const GetPot& dataFile, const std::string& section = "adr" );

    //! Display the internal values
    void showMe( std::ostream& output = std::cout ) const;

    //@}

    //! @name Set methods
    //@{

    inline void setDiffusionCoefficient( const Real& diffusionCoefficient )
    {
        M_diffusionCoefficient = diffusionCoefficient;
    }

    inline void setReactionCoefficient( const Real& reactionCoefficient )
    {
        M_reactionCoefficient = reactionCoefficient;
    }

    inline void setSteady( const Real& steady )
    {
        M_steady = steady;
    }

    //! Set data time container
    /*!
     * @param DataTime shared_ptr to dataTime container
     */
    inline void setDataTimePtr( const dataTime_ptr_type dataTimePtr ) { M_dataTimePtr = dataTimePtr; }

    inline void setFieldDimension( const UInt& fieldDim ) { M_solutionFieldDimension = fieldDim; }

    inline void setStabilizationMethod( const ADRStabilization& stabMethod ) { M_stabilizationMethod = stabMethod; }

    inline void setStabilizationCoefficient( const Real& stabCoeff ) { M_stabilizationCoefficient = stabCoeff; }
    //@}

    //! @name Get methods
    //@{

    inline const Real& diffusionCoefficient() const { return M_diffusionCoefficient; }

    inline const Real& reactionCoefficient() const { return M_reactionCoefficient; }

    inline const Real& steady() const { return M_steady; }

    //! Get data time container
    inline dataTime_ptr_type dataTimePtr( void ) const { return M_dataTimePtr; }

    inline const UInt& solutionFieldDimension() const { return M_solutionFieldDimension; }

    inline const UInt& verbose() const { return M_verbose; }

    inline const std::string& solFEType() const { return M_solFEType; }

    // inline const std::string& advectionFieldFEType() const { return M_advectionFieldFEType; }

    inline const ADRStabilization& stabilizationMethod() const { return M_stabilizationMethod; }

    inline const Real& stabilizationCoefficient() const { return M_stabilizationCoefficient; }

    //@}

private:

    //! Physics
    Real M_diffusionCoefficient;     // Diffusion coefficient
    Real M_reactionCoefficient;      // Reaction coefficient
    Real M_steady;                   // is the problem constant in time?
    UInt M_solutionFieldDimension;   // number of components of the unknown field

    //! Data container for time parameters
    dataTime_ptr_type               M_dataTimePtr;

    //! Miscellaneous parameters
    UInt M_verbose;                  // Class output verbose

    //! Type of finite element (P1, P2, ...) for the solution
    std::string  M_solFEType;
    //! Type of finite element (P1, P2, ...) for the advection field
    // std::string  M_advectionFieldFEType;

    //! Discretization
    DataStringList M_stabilization_list;
    ADRStabilization M_stabilizationMethod;
    // stabilization coefficient for the advection field
    Real M_stabilizationCoefficient;

};


}
#endif

