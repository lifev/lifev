//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
  @file
  @brief dataSecondOrder - Class to secondorder problem (S. Venant Kirchhoff Viscoelastic)
 
  @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
  
  @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
  @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 */

#ifndef VENANTKIRCHHOFFVISCOELASTICDATA_H
#define VENANTKIRCHHOFFVISCOELASTICDATA_H 1

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <string>
#include <iostream>
#include <map>
#include <list>
#include <boost/shared_ptr.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifefem/dataTime.hpp>


namespace LifeV
{
/*!
  \class DataSecondOrder
*/
class DataSecondOrder:
      public DataTime
{
public:
    //! @name Type definitions
    //@{

    typedef DataTime                                 Time_Type;
    typedef boost::shared_ptr< Time_Type >           Time_ptrType;

    typedef std::map<UInt, Real>                    MaterialContainer_Type;
    typedef MaterialContainer_Type::const_iterator  MaterialContainer_ConstIterator;

    //@}

    //! Constructor
    DataSecondOrder();
    DataSecondOrder( const DataSecondOrder& DataSecondOrder);

    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param DataSecondOrder - DataSecondOrder
     */
    DataSecondOrder& operator=( const DataSecondOrder& DataSecondOrder );

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
    void setDataTime( const Time_ptrType DataTime );

    //! Set density
    /*!
     * @param density solid density value
     */
    void setDensity( const Real& density );

    //! Set gamma
    /*!
     * @param gamma damping coefficient.
     */
    void setGamma( const Real& gamma );

    void setGamma( const Real& gamma, const UInt& material );

    //! Set beta
    /*!
     * @param beta damping coefficient.
     */
    void setBeta( const Real& beta);

    void setBeta( const Real& beta, const UInt& material);


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
    Time_ptrType dataTime() const;

    //! Get solid density
    /*!
     * @return Solid density
     */
    const Real& rho() const;

    //! Get solid thickness
    /*!
     * @return Solid thickness
     */
    const Real& thickness() const;

    //! Get solid poisson coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid poisson coefficient
     */
    const Real& poisson( const UInt& material = 1 ) const;

    //! Get solid young modulus
    /*!
     * @param material material ID (1 by default)
     * @return Solid young modulus
     */
    const Real& young( const UInt& material = 1 ) const;

    //! Get solid first lame coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid first Lame coefficient
     */
    Real lambda( const UInt& material = 1 ) const;

    //! Get solid second Lame coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid second Lame coefficient
     */
    Real mu( const UInt& material = 1 ) const;

    //! Get damping coefficients
    /*!
     * @return gamma damping coefficient
     */
    const Real& gamma( const UInt& material = 1 ) const;

    /*!
    * @return beta damping coefficient
    */
    const Real& beta( const UInt& material = 1 ) const;

    //! Get FE order
    /*!
     * @return FE order
     */
    const std::string& order()     const;

    //! Get solid amplification factor
    /*!
     * @return Solid amplification factor
     */
    const Real& factor()    const;

    //! Get verbose level
    /*!
     * @return verbose level
     */
    const UInt& verbose()   const;

    //!
    /*
     * @ return true if Damping coefficient is not zero.
    */
    bool  damping() const
    {
        return M_damping;
    }
    //@}

private:

    //! Data containers for time and mesh
    Time_ptrType           M_time;

    //@name Physics 
    //@{
    //! density
    Real                   M_density;

    //!thickness 
    Real                   M_thickness;

    //!poisson 
    MaterialContainer_Type M_poisson;

    //!young
    MaterialContainer_Type M_young;

    //!damping coefficient (Stiffness)
    MaterialContainer_Type M_gamma;
  
    //! damping coefficint (mass)
    MaterialContainer_Type  M_beta;
    //@}
    
    //! Miscellaneous
    Real                   M_factor; // amplification factor for deformed mesh
    //!verbose
    UInt                   M_verbose; // temporal output verbose
    //! order
    std::string            M_order;

    //!damping true when there is damping term
    bool                   M_damping;  
};

}  //namespace 
#endif /* VENANTKIRCHHOFFVISCOELASTICDATA_H*/
