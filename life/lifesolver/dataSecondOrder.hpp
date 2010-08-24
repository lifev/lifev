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
  \file dataSecondOrder.h
  \author F. Nobile, M. Pozzoli, C. Vergara
  \date 02/2010
  \version 1.0
  \brief
*/

#ifndef _DATASECONDORDER_H_
#define _DATASECONDORDER_H_
#include <string>
#include <iostream>
#include <map>
#include <life/lifecore/GetPot.hpp>
#include <boost/shared_ptr.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifefem/dataTime.hpp>
#include <list>

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

    typedef DataTime                                                  Time_Type;
    typedef boost::shared_ptr< Time_Type >                            Time_ptrType;

  //@}

    //! Constructor
    DataSecondOrder();
    DataSecondOrder( const GetPot& dfile );
    DataSecondOrder( const DataSecondOrder& dataSecondOrder);

 //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param dataSecondOrder - DataSecondOrder
     */
    DataSecondOrder& operator=( const DataSecondOrder& dataSecondOrder );

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

    //! Set alpha
    /*!
     * @param alpha damping coefficient.
     */
     void setAlpha( const Real& alpha );

     void setAlpha( const Real& alpha, const UInt& material );

    //! Set beta
    /*!
     * @param beta damping coefficient. 
     */
     void setBeta( const Real& beta);

  void setBeta( const Real& beta, const UInt& material);

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

    //! Get damping coefficients
    /*!
     * @return alpha damping coefficient
     */
     const Real& alpha( const UInt& material = 1 ) const;

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

   bool  isDamping() const
  {
    if (M_beta.size() == 1 || M_alpha.size()==1)
      return  _beta +_alpha;
    return false;
  }
   //@}

private:

    //! Data containers for time and mesh
    Time_ptrType           M_time;

    //! Physics
    Real                   M_density; // densisty
    Real                   _alpha;  // Damping coefficient (Mass)
    Real                   _beta; // Damping coefficient (Stiffness)
    bool                   M_isDamping;  //true if damping else false

    std::map<int, double>  M_alpha;
    std::map<int, double>  M_beta;

    //! Miscellaneous
    Real                  M_factor; // amplification factor for deformed mesh
    UInt                   M_verbose; // temporal output verbose
    std::string        M_order;

};



//
// IMPLEMENTATION
//
DataSecondOrder::DataSecondOrder() :
        M_time                             ( ),
        M_density                        ( ),
        M_alpha                           ( ),
        M_beta                             ( ),
        _alpha                              ( ),
        _beta                                ( ),
        M_factor                           ( ),
        M_verbose                       ( )
{
}

// Constructor
DataSecondOrder::DataSecondOrder( const GetPot& dfile ) :
         DataTime( dfile, "problem/time_discretization" )
{
    // physics
       M_density        = dfile( "problem/physics/density", 1. );
       _alpha    = dfile( "problem/physics/alpha" , 0.0 );
       _beta     = dfile( "problem/physics/beta" , 0.0 );

    // miscellaneous
    M_factor  = dfile( "problem/miscellaneous/factor", 1.0 );
    M_verbose = dfile( "problem/miscellaneous/verbose", 1 );

    M_order  = dfile( "problem/space_discretization/order", "P1");

    std::string flagList;
    flagList = dfile("problem/physics/material_flag", "0");
    std::list<int> fList;
    parseList(flagList, fList);

}

DataSecondOrder::DataSecondOrder( const DataSecondOrder& dataSecondOrder):
  DataTime               ( dataSecondOrder),
  M_density               ( dataSecondOrder.M_density ),
  _alpha                     ( dataSecondOrder._alpha ),
  _beta                       ( dataSecondOrder._beta ),
  M_isDamping         ( dataSecondOrder.M_isDamping),
  M_alpha                  ( dataSecondOrder.M_alpha ),
  M_beta                    ( dataSecondOrder.M_beta ),
  M_factor                  ( dataSecondOrder.M_factor ),
  M_verbose              ( dataSecondOrder.M_verbose ),
  M_order                   ( dataSecondOrder.M_order )
{
}

DataSecondOrder&
DataSecondOrder::operator=( const DataSecondOrder& dataSecondOrder )
{
    if ( this != &dataSecondOrder )
    {
        M_time                       = dataSecondOrder.M_time;
        M_density                  = dataSecondOrder.M_density;
	M_isDamping            = dataSecondOrder.M_isDamping;
        M_alpha                     = dataSecondOrder.M_alpha;
        M_beta                       = dataSecondOrder.M_beta;
        _alpha                        = dataSecondOrder._alpha;
	_beta                          = dataSecondOrder._beta;
        M_order                      = dataSecondOrder.M_order;
        M_factor                     = dataSecondOrder.M_factor;
        M_verbose                 = dataSecondOrder.M_verbose;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================

void
DataSecondOrder::setup( const GetPot& dataFile, const std::string& section )
{
    // If data time has not been set
    if ( !M_time.get() )
        M_time.reset( new Time_Type( dataFile, section + "/time_discretization" ) );

    // physics
    M_density   = dataFile( ( section + "/physics/density" ).data(), 1. );

   M_alpha[1]   = dataFile( ( section + "/physics/young" ).data(), 0. );
   M_beta[1] = dataFile( ( section + "/physics/poisson" ).data(), 0. );

     
    // space_discretization
    M_order     = dataFile( "solid/space_discretization/order", "P1" );

    // miscellaneous
    M_factor  = dataFile( "solid/miscellaneous/factor", 1.0 );
    M_verbose = dataFile( "solid/miscellaneous/verbose", 1 );
}


 
// ===================================================
// Set Method
// ===================================================

 void
DataSecondOrder::setDataTime( const Time_ptrType DataTime )
{
    M_time = DataTime;
}

 void
DataSecondOrder::setDensity( const Real& density )
{
    M_density = density;
}

void
 DataSecondOrder::setAlpha( const Real& alpha )
{
    _alpha = alpha;
}

 void
 DataSecondOrder::setAlpha( const Real& alpha, const UInt& material )
{
    M_alpha[material] = alpha;
}

void
 DataSecondOrder::setBeta( const Real& beta )
{
    _beta = beta;
}


 void
DataSecondOrder::setBeta( const Real& beta, const UInt& material )
{
    M_beta[material] = beta;
}



// ===================================================
// Get Method
// ===================================================

DataSecondOrder::Time_ptrType
DataSecondOrder::dataTime() const
{
    return M_time;
}

const Real&
DataSecondOrder::rho() const
{
    return M_density;
}

 const Real&
DataSecondOrder::alpha( const UInt& material ) const
{
    return M_alpha.find( material )->second;
}


 const Real&
DataSecondOrder::beta( const UInt& material ) const
{
    return M_beta.find( material )->second;
}


 const Real&
DataSecondOrder::factor() const
{
    return M_factor;
}

const std::string&
DataSecondOrder::order() const
{
    return M_order;
}


 
  // Output
void DataSecondOrder::showMe( std::ostream& c ) const
{

    // physics
    c << "\n*** Values for data [problem/physics]\n\n";
    c << "density                          = " << M_density << std::endl;
    c << "alpha                            = " << _alpha << std::endl;
    c << "beta                             = " << _beta << std::endl;
    
    c << "\n*** Values for data [problem/miscellaneous]\n\n";
    c << "deformation factor               = " << M_factor << std::endl;
    c << "verbose                          = " << M_verbose << std::endl;

    c << "\n*** Values for data [problem/space_discretization]\n\n";
   
    c << "\n*** Values for data [problem/time_discretization]\n\n";
    DataTime::showMe( c );
}

}

#endif
