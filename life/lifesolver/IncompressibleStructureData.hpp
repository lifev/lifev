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
    @brief File containing the implementation of the file IncompressibleStructureData.hpp

    @author Simone Rossi <simone.rossi@epfl.ch>
    @contributor
    @maintainer Simone Rossi <simone.rossi@epfl.ch>

    @date 04-2011

 */


#ifndef INCOMPRESSIBLESTRUCTUREDATA_H
#define INCOMPRESSIBLESTRUCTUREDATA_H

#include <life/lifefilters/GetPot.hpp>
#include <life/lifecore/LifeV.hpp>
#include <life/lifecore/StringData.hpp>
#include <life/lifecore/StringUtility.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <iostream>

namespace LifeV
{


class IncompressibleStructureData
{

public:

    //! @name Public Types
    //@{


    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    IncompressibleStructureData();


    //! Copy constructor
    /*!
     * @param incompressibleStructureData IncompressibleStructureData
     */
    IncompressibleStructureData( const IncompressibleStructureData& incompressibleStructureData );

    //! Virtual destructor
    virtual ~IncompressibleStructureData() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param incompressibleStructureData IncompressibleStructureData
     */
    IncompressibleStructureData& operator=( const IncompressibleStructureData& incompressibleStructureData );

    //@}


    //! @name Methods
    //@{

    //! Read the dataFile and set all the quantities
    /*!
     * @param dataFile data file
     * @param section section of the file
     */
    void setup( const GetPot& dataFile, const std::string& section = "structure" );

    //! Display the values
    void showMe( std::ostream& output = std::cout ) const;

    //@}



    //! @name Set methods
    //@{



    //! @name Get methods
    //@{

    //! Get the viscosity of the fluid
    /*!
     * @param n the fluid number
     * @return M_viscosity the viscosity of the fluid n
     */
    const Real& viscosity(const UInt& n=0) const
    {
  //      ASSERT(n<M_fluidNumber,"Undeclared fluid");
        return M_viscosity[n];
    }

    //! Get the order of the finite elements used for velocity
    /*!
     * @return M_uOrder a string specifying the order of finite elements for velocity
     */
    std::string     uOrder()                      const { return M_uOrder; }

    //! Get the order of the finite elements used for pressure
    /*!
     * @return M_pOrder a string specifying the order of finite elements for pressure
     */
    std::string     pOrder()                      const { return M_pOrder; }

    //! Temporal output verbose
    /*!
     * @return M_verbose
     */
    UInt            verbose()                     const { return M_verbose; }


    //@}

protected:

    //! @name Physics
    //@{


    //! viscosity of each fluid
    std::vector<Real> M_viscosity;

    //@}


    //! @name FE order
    //@{

    //! order of finite elements for velocity
    std::string      M_uOrder;

    //! order of finite elements for pressure
    std::string      M_pOrder;

    //@}


    //! @name Miscellaneous
    //@{

    //! temporal output verbose
    UInt             M_verbose;



private:


};



} // end namespace LifeV



#endif /* INCOMPRESSIBLESTRUCTUREDATA_H */
