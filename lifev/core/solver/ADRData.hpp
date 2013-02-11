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
    @brief File containing a class for handling data and parameters
         for advection-diffusion-reaction problems

    @author
    @date 00-08-2010

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

 */
#ifndef _DATAADR_H_
#define _DATAADR_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/fem/TimeData.hpp>
#include <lifev/core/util/StringData.hpp>

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

    //! @name Public Types
    //@{

    //! The time data handler
    /*!
       The class DataADR "has a" time handler object, to manage time dependent
       problems.
     */
    typedef TimeData                              TimeData_type;
    typedef boost::shared_ptr<TimeData_type>      TimeData_ptr_type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    DataADR();

    //! Copy constructor
    /*!
     * @param dataADR an object of type DataADR
     */
    DataADR ( const DataADR& dataADR );

    //@}

    //! @name Methods
    //@{

    //! Read the dataFile and set all the internal quantities
    /*!
     * @param dataFile data file
     * @param section section of the file
     */
    void setup ( const GetPot& dataFile, const std::string& section = "adr" );

    //! Display the internal values
    void showMe ( std::ostream& output = std::cout ) const;

    //@}

    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param dataADR an object of type DataADR
     */
    DataADR& operator= ( const DataADR& dataADR );

    //@}

    //! @name Set methods
    //@{

    inline void setDiffusionCoefficient ( const Real& diffusionCoefficient )
    {
        M_diffusionCoefficient = diffusionCoefficient;
    }

    inline void setReactionCoefficient ( const Real& reactionCoefficient )
    {
        M_reactionCoefficient = reactionCoefficient;
    }

    inline void setSteady ( const Real& steady )
    {
        M_steady = steady;
    }

    //! Set data time container
    /*!
     * @param TimeData shared_ptr to TimeData container
     */
    inline void setTimeDataPtr ( const TimeData_ptr_type TimeDataPtr )
    {
        M_TimeDataPtr = TimeDataPtr;
    }

    inline void setFieldDimension ( const UInt& fieldDim )
    {
        M_solutionFieldDimension = fieldDim;
    }

    inline void setStabilizationMethod ( const ADRStabilization& stabMethod )
    {
        M_stabilizationMethod = stabMethod;
    }

    inline void setStabilizationCoefficient ( const Real& stabCoeff )
    {
        M_stabilizationCoefficient = stabCoeff;
    }
    //@}

    //! @name Get methods
    //@{

    inline const Real& diffusionCoefficient() const
    {
        return M_diffusionCoefficient;
    }

    inline const Real& reactionCoefficient() const
    {
        return M_reactionCoefficient;
    }

    inline const Real& steady() const
    {
        return M_steady;
    }

    //! Get data time container
    inline TimeData_ptr_type TimeDataPtr ( void ) const
    {
        return M_TimeDataPtr;
    }

    inline const UInt& solutionFieldDimension() const
    {
        return M_solutionFieldDimension;
    }

    inline const UInt& verbose() const
    {
        return M_verbose;
    }

    inline const std::string& solFEType() const
    {
        return M_solFEType;
    }

    // inline const std::string& advectionFieldFEType() const { return M_advectionFieldFEType; }

    inline const ADRStabilization& stabilizationMethod() const
    {
        return M_stabilizationMethod;
    }

    inline const Real& stabilizationCoefficient() const
    {
        return M_stabilizationCoefficient;
    }

    //@}

private:

    //! Physics
    Real M_diffusionCoefficient;     // Diffusion coefficient
    Real M_reactionCoefficient;      // Reaction coefficient
    Real M_steady;                   // is the problem constant in time?
    UInt M_solutionFieldDimension;   // number of components of the unknown field

    //! Data container for time parameters
    TimeData_ptr_type               M_TimeDataPtr;

    //! Miscellaneous parameters
    UInt M_verbose;                  // Class output verbose

    //! Type of finite element (P1, P2, ...) for the solution
    std::string  M_solFEType;
    //! Type of finite element (P1, P2, ...) for the advection field
    // std::string  M_advectionFieldFEType;

    //! Discretization
    StringDataList M_stabilization_list;
    ADRStabilization M_stabilizationMethod;
    // stabilization coefficient for the advection field
    Real M_stabilizationCoefficient;

};


}
#endif

