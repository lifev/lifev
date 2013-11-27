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
  @brief VenantKirchhoffViscoelasticData - Class to secondorder problem (S. Venant Kirchhoff Viscoelastic)

  @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>

  @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
  @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 */

#ifndef VENANTKIRCHHOFFVISCOELASTICDATA_H
#define VENANTKIRCHHOFFVISCOELASTICDATA_H 1


#include <string>
#include <iostream>
#include <map>
#include <list>
#include <boost/shared_ptr.hpp>


#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/TimeData.hpp>
#include <lifev/core/fem/TimeAdvanceData.hpp>


namespace LifeV
{
/*!
  \class VenantKirchhoffViscoelasticData
*/
class VenantKirchhoffViscoelasticData
{
public:
    //! @name Type definitions
    //@{

    typedef TimeData                               time_Type;
    typedef boost::shared_ptr<time_Type>           timePtr_Type;

    typedef TimeAdvanceData                        timeAdvance_Type;
    typedef boost::shared_ptr<timeAdvance_Type>    timeAdvancePtr_Type;

    typedef std::map<UInt, Real>                   MaterialContainer_Type;
    typedef MaterialContainer_Type::const_iterator MaterialContainer_ConstIterator;

    //@}

    //! Constructor
    VenantKirchhoffViscoelasticData();
    VenantKirchhoffViscoelasticData ( const VenantKirchhoffViscoelasticData& VenantKirchhoffViscoelasticData);

    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param VenantKirchhoffViscoelasticData - VenantKirchhoffViscoelasticData
     */
    VenantKirchhoffViscoelasticData& operator= ( const VenantKirchhoffViscoelasticData& VenantKirchhoffViscoelasticData );

    //@}
    //! @name Methods
    //@{

    //! Read the dataFile and set all the quantities
    /*!
     * @param dataFile data file
     * @param section section of the file
     */
    void setup ( const GetPot& dataFile, const std::string& section = "solid" );

    //! Display the values
    void showMe ( std::ostream& output = std::cout ) const;

    //@}

    //! @name Set methods
    //@{

    //! Set data time container
    /*!
     * @param TimeData shared_ptr to TimeData container
     */
    void setTimeData ( const timePtr_Type timeData )
    {
        M_time = timeData;
    }

    //! Set data time advance container
    /*!
     * @param timeAdvanceData shared_ptr to TimeAdvanceData container
     */
    void setTimeAdvanceData ( const timeAdvancePtr_Type timeAdvanceData )
    {
        M_timeAdvance = timeAdvanceData;
    }

    //! Set density
    /*!
     * @param density solid density value
     */
    void setDensity ( const Real& density );

    //! Set gamma
    /*!
     * @param gamma damping coefficient.
     */
    void setGamma ( const Real& gamma );

    void setGamma ( const Real& gamma, const UInt& material );

    //! Set beta
    /*!
     * @param beta damping coefficient.
     */
    void setBeta ( const Real& beta);

    void setBeta ( const Real& beta, const UInt& material);


    //! Set thickness
    /*!
     * @param thickness solid thickness value
     */
    void setThickness ( const Real& thickness );

    //! Set poisson
    /*!
     * @param poisson solid poisson value
     * @param material material ID (1 by default)
     */
    void setPoisson ( const Real& poisson, const UInt& material = 1 );

    //! Set Young modulus
    /*!
     * @param Young solid young modulus value
     * @param material material ID (1 by default)
     */
    void setYoung ( const Real& young, const UInt& material = 1 );

    //@}

    //! @name Get methods
    //@{

    //! Get data time container
    /*!
     * @return shared_ptr to TimeData container
     */
    timePtr_Type dataTime() const
    {
        return M_time;
    }

    //! Get data time container
    /*!
     * @return shared_ptr to TimeAdvanceData container
     */
    timeAdvancePtr_Type dataTimeAdvance() const
    {
        return M_timeAdvance;
    }

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
    Real poisson ( const UInt& material = 1 ) const;

    //! Get solid young modulus
    /*!
     * @param material material ID (1 by default)
     * @return Solid young modulus
     */
    Real young ( const UInt& material = 1 ) const;

    //! Get solid first lame coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid first Lame coefficient
     */
    Real lambda ( const UInt& material = 1 ) const;

    //! Get solid second Lame coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid second Lame coefficient
     */
    Real mu ( const UInt& material = 1 ) const;

    //! Get damping coefficients
    /*!
     * @return gamma damping coefficient
     */
    const Real& gamma ( const UInt& material = 1 ) const;

    /*!
    * @return beta damping coefficient
    */
    const Real& beta ( const UInt& material = 1 ) const;

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
    timePtr_Type        M_time;
    timeAdvancePtr_Type M_timeAdvance;

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

    //!verbose
    UInt                   M_verbose; // temporal output verbose
    //! order
    std::string            M_order;

    //!damping true when there is damping term
    bool                   M_damping;
};

}  //namespace
#endif /* VENANTKIRCHHOFFVISCOELASTICDATA_H*/
