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
  \file dataMonodomain.h
  \author L. Mirabella M. Perego
  \date 11/2007
  \version 1.0

  \brief File containing a class for handling Monodomain data with GetPot

*/
#ifndef _DATAIONIC_H_
#define _DATAIONIC_H_
#include <string>
#include <iostream>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifecore/dataString.hpp>
#include <life/lifearray/tab.hpp>

namespace LifeV
{
using namespace std;

/*!
  \class DataMonodomain

  Base class which holds usual data for the ionic model solvers

*/
template <typename Mesh>
class DataIonic:
        public DataMesh<Mesh>,
        public DataTime
{
public:

    //! Constructors
    DataIonic( const GetPot& dfile );

    DataIonic( const DataIonic& dataIonic );

    //! Ouptut
    void showMe( std::ostream& c = std::cout );

    //! external setup
    void setup( const GetPot& dfile );

    //! End time
    Real endtime() const;

    //! FE space order
    std::string wOrder() const;

    UInt verbose; 
    string mesh_file;
    Real a;
    Real b;
    Real c1;
    Real c2;
    Real d;
    Real T;
    Real A;
    Real u0;
    Real winit;
    string mesh_dir;
        
    
private:


};


//
// IMPLEMENTATION
//


//! Constructors
template <typename Mesh>
DataIonic<Mesh>::
DataIonic( const GetPot& dfile ) :
    DataMesh<Mesh>( dfile, "electric/discretization" ),
    DataTime( dfile, "electric/discretization" )
{
    setup(dfile);
}

template <typename Mesh>
DataIonic<Mesh>::
DataIonic( const DataIonic& dataIonic ) :
    DataMesh<Mesh>               ( dataIonic ),
    DataTime                     ( dataIonic ),
    a(dataIonic.a),
    b(dataIonic.b),
    c1(dataIonic.c1),
    c2(dataIonic.c2),
    d(dataIonic.d),
    T(dataIonic.T),
    A(dataIonic.A),
    u0(dataIonic.u0),   
    winit(dataIonic.winit)
{
}


template <typename Mesh>
void
DataIonic<Mesh>::
setup(  const GetPot& dfile )
{
    a   		= dfile("electric/physics/a",0.13);   // 0.13  adim  //RogersMcCulloch1994
    b   		= dfile("electric/physics/b",0.013);  // 0.013 adim //RogersMcCulloch1994
    c1   		= dfile("electric/physics/c1",0.26);  // 0.26  adim //RogersMcCulloch1994
    c2   		= dfile("electric/physics/c2",0.1);   //0.1    adim //RogersMcCulloch1994
    d   		= dfile("electric/physics/d",1);      //1      adim //RogersMcCulloch1994
    T   		= dfile("electric/physics/T",0.63);    //0.63ms    //RogersMcCulloch1994
    A   		= dfile("electric/physics/A",130);    //130mV    //RogersMcCulloch1994
    u0   		= dfile("electric/physics/u0",-84.0);	  //-84mV    //RogersMcCulloch1994
    winit		= dfile("electric/physics/winit", 0);  
    
}

// Output
template <typename Mesh>
void DataIonic<Mesh>::
showMe( std::ostream& c )
{

}


}
#endif
