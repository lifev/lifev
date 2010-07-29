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
  \file dataIonic.cpp
  \author L. Mirabella M. Perego
  \date 11/2007
  \version 1.0

  \brief File containing a class for handling Monodomain data with GetPot

	\modif J.Castelneau (INRIA)
	\date 06/09
*/
#include <life/lifesolver/dataIonic.hpp>

namespace LifeV
{
//
// IMPLEMENTATION
//


//! Constructors
DataIonic::
DataIonic( const GetPot& dfile ) :
    DataMesh( dfile, "electric/space_discretization" ),
    DataTime( dfile, "electric/time_discretization" )
{
    setup(dfile);
}

DataIonic::
DataIonic( const DataIonic& dataIonic ) :
    DataMesh		              ( dataIonic ),
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
    tau_out(dataIonic.tau_out),
    tau_open(dataIonic.tau_open),
    tau_close(dataIonic.tau_close),
    vcrit(dataIonic.vcrit),
    v_min(dataIonic.v_min),
    v_max(dataIonic.v_max),
    reac_amp(dataIonic.reac_amp),
    tinit(dataIonic.tinit),
    tend(dataIonic.tend),
    order_bdf(dataIonic.order_bdf),
    has_HeteroTauClose(dataIonic.has_HeteroTauClose)
{
}


void
DataIonic::
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
void DataIonic::
showMe( std::ostream& /*c*/ )
{

}

}

