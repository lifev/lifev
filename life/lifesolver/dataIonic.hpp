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
class DataIonic:
        public DataMesh,
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

}
#endif
