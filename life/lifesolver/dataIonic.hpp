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

	\modif J.Castelneau (INRIA)
	\date 06/09
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
	//Mitchell & Schaeffer
	Real		tau_in;   // = 0.8
	Real		tau_out;  // = 18.0
	Real		tau_open; // = 300.0
	Real		tau_close;// = 100.0
	Real		vcrit;    // =  -67.0
	Real 		v_min;
	Real 		v_max;
	Real 		reac_amp;
	Real 		tinit;
	Real 		tend;
	Real 		order_bdf;       //= 1
    bool		has_HeteroTauClose;
private:


};


//
// IMPLEMENTATION
//


//! Constructors
template <typename Mesh>
DataIonic<Mesh>::
DataIonic( const GetPot& dfile ) :
    DataMesh<Mesh>( dfile, "electric/space_discretization" ),
    DataTime( dfile, "electric/time_discretization" )
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
    winit(dataIonic.winit),
	// Mitchell & Schaeffer
    tau_in(dataIonic.tau_in),
    v_min(dataIonic.v_min),
    v_max(dataIonic.v_max),
    reac_amp(dataIonic.reac_amp),
    tau_out(dataIonic.tau_out),
    tau_open(dataIonic.tau_open),
    tau_close(dataIonic.M_tau_close),
    vcrit(dataIonic.vcrit),
    tinit(dataIonic.tinit),
    tend(dataIonic.tend),
    order_bdf(dataIonic.order_bdf),
    has_HeteroTauClose(dataIonic.has_HeteroTauClose)
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
    A   		= dfile("electric/physics/A",110);    //130mV    //RogersMcCulloch1994
    u0   		= dfile("electric/physics/u0",-84.0);	  //-84mV    //RogersMcCulloch1994
    winit		= dfile("electric/physics/winit", 0);
	// Mitchell & Schaeffer
    tau_in    = dfile("electric/physics/tau_in",0.8);
    v_min    = dfile("electric/physics/v_min",-80.0);
    v_max    = dfile("electric/physics/v_max", 20.0);
    reac_amp    = dfile("electric/physics/reac_amp", 0.2);
    tau_out   = dfile("electric/physics/tau_out",18.0);
    tau_open  = dfile("electric/physics/tau_open",100.0);
    tau_close = dfile("electric/physics/tau_close",100.0);
    vcrit     = dfile("electric/physics/vcrit",-67.0);
    tinit     = dfile("electric/physics/init_time",0.0);
    tend      = dfile("electric/physics/end_time",1000.0);
    order_bdf       = dfile("electric/time_discretization/BDF_order",1);
    has_HeteroTauClose = dfile("electric/physics/hasHeteroTauClose",1);
}

// Output
template <typename Mesh>
void DataIonic<Mesh>::
showMe( std::ostream& c )
{

}

}
#endif
