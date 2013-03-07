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
    @brief File containing a class for handling temporal discretization.

    @author Cristiano Malossi <cristiano.malossi@epfl.ch>
    @date 11-06-2012

    @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 */

#ifndef TimeAdvanceData_H
#define TimeAdvanceData_H 1

#include <ostream>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/LifeV.hpp>

namespace LifeV
{

//! TimeAdvanceData - Class for handling temporal discretization.
/*!
 *  @author Cristiano Malossi
 *
 *  The class is a data container for the time discretization.
 */
class TimeAdvanceData
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    TimeAdvanceData();

    //! Constructor
    /*!
     * @param dataFile GetPot data file
     * @param section the section on the data file that contains the information on the time discretization
     */
    TimeAdvanceData ( const GetPot& dataFile, const std::string& section = "time_discretization" );

    //! Copy constructor
    /*!
     * @param TimeAdvanceData - TimeAdvanceData class
     */
    TimeAdvanceData ( const TimeAdvanceData& TimeAdvanceData);

    //! Virtual destructor
    virtual ~TimeAdvanceData() {}

    //@}


    //! @name Methods
    //@{

    //! Read the dataFile and set all the internal quantities
    /*!
     * @param dataFile data file
     * @param section section of the file
     */
    void setup ( const GetPot& dfile, const std::string& section = "time_discretization" );

    //! Display general information about the content of the class
    /*!
        @param output - specify the output format (std::cout by default)
     */
    void showMe ( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the BDF odert to use
    /*!
     * @param order BDF order
     */
    void setOrderBDF ( const UInt& orderBDF )
    {
        M_orderBDF = orderBDF;
    }

    //! Set the theta of TimeAdvanceNewmark scheme
    /*!
     * @param theta - coefficient of TimeAdvanceNewmark scheme
     */
    void setTheta ( const Real& theta )
    {
        M_theta = theta;
    }

    //! Set the theta of TimeAdvanceNewmark scheme
    /*!
     * @param gamma- coefficient of TimeAdvanceNewmark scheme
     */
    void setGamma ( const Real& gamma )
    {
        M_gamma = gamma;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Get the BDF order used
    /*!
     * @return BDF order value
     */
    const UInt& orderBDF() const
    {
        return M_orderBDF;
    }

    //! Return theta parameter of TimeAdvanceNewmark scheme
    /*!
     * @return theta value
     */
    const Real& theta() const
    {
        return M_theta;
    }

    //! Return gamma of TimeAdvanceNewmark scheme
    /*!
     * @return gamma value
     */
    const Real& gamma() const
    {
        return M_gamma;
    }

    //! Return TimeAdvanceNewmark parameters (\f$theta\f$, \f$\gamma\f$)
    /*!
     * @return TimeAdvanceNewmark parameters (\f$theta\f$, \f$\gamma\f$)
     */
    std::vector<Real> coefficientsNewmark();

    //@}

private:

    UInt                    M_orderBDF;       // order of the time discretization formula
    Real                    M_theta;          // TimeAdvanceNewmark parameter
    Real                    M_gamma;          // TimeAdvanceNewmark parameter
};

} // namespace LifeV

#endif // TimeAdvanceData_H
