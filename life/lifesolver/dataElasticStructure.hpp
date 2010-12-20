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
 *  @file
 *  @brief DataElasticStructure - File containing a data container for solid problems with elastic structure
 *
 *  @version 1.0
 *  @date 01-10-2003
 *  @author M.A. Fernandez
 *
 *  @version 1.18
 *  @date 10-06-2010
 *  @author Cristiano Malossi
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef _DATAELASTICSTRUCTURE_H_
#define _DATAELASTICSTRUCTURE_H_

#include <string>
#include <iostream>
#include <map>
#include <life/lifecore/GetPot.hpp>
#include <boost/shared_ptr.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifecore/util_string.hpp>
#include <life/lifefem/dataTime.hpp>

namespace LifeV
{

//! DataElasticStructure - Data container for solid problems with elastic structure
class DataElasticStructure: public DataTime
{
public:

    //! @name Type definitions
    //@{

    typedef DataTime                                                  Time_Type;
    typedef boost::shared_ptr< Time_Type >                            TimePtr_Type;

    typedef std::map<UInt, Real>                                      MaterialContainer_Type;
    typedef MaterialContainer_Type::const_iterator                    MaterialContainer_ConstIterator;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    DataElasticStructure();

    //! Copy constructor
    /*!
     * @param dataElasticStructure - DataElasticStructure
     */
    DataElasticStructure( const DataElasticStructure& dataElasticStructure );

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param dataElasticStructure - DataElasticStructure
     */
    DataElasticStructure& operator=( const DataElasticStructure& dataElasticStructure );

    //@}


    //! @name Methods
    //@{

    //! Read the dataFile and set all the quantities
    /*!
     * @param dataFile data file
     * @param section section of the file
     */
    void setup( const GetPot& dataFile, const std::string& section = "solid" );

    //! Display the values
    void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set methods
    //@{

    //! Set data time container
    /*!
     * @param DataTime shared_ptr to dataTime container
     */
    void setDataTime( const TimePtr_Type DataTime );

    //! Set density
    /*!
     * @param density solid density value
     */
    void setDensity( const Real& density );

    //! Set thickness
    /*!
     * @param thickness solid thickness value
     */
    void setThickness( const Real& thickness );

    //! Set poisson
    /*!
     * @param poisson solid poisson value
     * @param material material ID (1 by default)
     */
    void setPoisson( const Real& poisson, const UInt& material = 1 );

    //! Set Young modulus
    /*!
     * @param Young solid young modulus value
     * @param material material ID (1 by default)
     */
    void setYoung( const Real& young, const UInt& material = 1 );

    //@}


    //! @name Get methods
    //@{

    //! Get data time container
    /*!
     * @return shared_ptr to dataTime container
     */
    const TimePtr_Type getDataTime() const;

    //! Get solid density
    /*!
     * @return Solid density
     */
    const Real& getRho() const;

    //! Get solid thickness
    /*!
     * @return Solid thickness
     */
    const Real& getThickness() const;

    //! Get solid poisson coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid poisson coefficient
     */
    const Real& getPoisson( const UInt& material = 1 ) const;

    //! Get solid young modulus
    /*!
     * @param material material ID (1 by default)
     * @return Solid young modulus
     */
    const Real& getYoung( const UInt& material = 1 ) const;

    //! Get solid first lame coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid first Lame coefficient
     */
    Real getLambda( const UInt& material = 1 ) const;

    //! Get solid second Lame coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid second Lame coefficient
     */
    Real getMu( const UInt& material = 1 ) const;

    //! Get FE order
    /*!
     * @return FE order
     */
    const std::string& getOrder() const;

    //! Get solid amplification factor
    /*!
     * @return Solid amplification factor
     */
    const Real& getFactor() const;

    //! Get verbose level
    /*!
     * @return verbose level
     */
    const UInt& getVerbose() const;

    //! Get solid type
    /*!
     * @return solid type
     */
    const std::string& getSolidType();

    //! Get whether to use or not exact Jacobian
    /*!
     * @return true: if using exact Jacobian, false: otherwise
     */
    const bool& getUseExactJacobian() const;

    //@}

private:

    //! Data containers for time and mesh
    TimePtr_Type           M_time;

    //! Physics
    Real                   M_density;
    Real                   M_thickness;

    MaterialContainer_Type M_poisson;
    MaterialContainer_Type M_young;

    //! Space discretization
    std::string            M_order;

    //! Miscellaneous
    Real                   M_factor;  // amplification factor for deformed mesh
    UInt                   M_verbose; // temporal output verbose

    std::string            M_solidType;
    bool                   M_useExactJacobian;

};

} // end namespace LifeV

#endif // end _DATAELASTICSTRUCTURE_H_
