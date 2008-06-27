/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
#ifndef _DATAMONODOMAIN_H_
#define _DATAMONODOMAIN_H_
#include <string>
#include <iostream>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifecore/dataString.hpp>
#include <life/lifearray/tab.hpp>
#include "heart_functors.hpp"



namespace LifeV
{
using namespace std;

/*!
  \class DataMonodomain

  Base class which holds usual data for the Monodomain model solvers

*/
template <typename Mesh>
class DataMonodomain:
        public DataMesh<Mesh>,
        public DataTime
{
public:

    //! Constructors
    DataMonodomain( boost::shared_ptr<HeartFunctors> heart_fct);

    DataMonodomain( const DataMonodomain& dataMonodomain );
    
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const&, Real const&)> fct_type;

    //! Ouptut
    void showMe( std::ostream& c = std::cout );


    //! external setup
    void setup( const GetPot& dfile );

    //! End time
    Real endtime() const {return M_endtime;};
    
    //! verbose
    Real verbose() const {return M_verbose;};        
        
    //! FE space order
    std::string uOrder() const {return M_uOrder;};
    //! Chi
	inline   Real Chi() const {return M_Chi;}
	 //! Chi
	inline   string fibers_file() const {return M_fibers_file;}
	
	inline int heart_diff_fct() const {return M_heart_diff_fct;}
	
	inline bool has_fibers() const {return M_has_fibers;}
	
	//! sigma_l
    inline   Real sigmal() const 	{return M_sigmal;}
    
    //! sigma_t
    inline   Real sigmat() const 	{return M_sigmat;}
    
    //! lambda
    inline   Real lambda() const 	{return M_lambda;}
        
	//! Cm
	inline   Real Cm() const 	{return M_Cm;}
	//! D
	inline   Real D() const 	{return M_D;}
	//! Post_dir
	inline   string Post_dir() const {return M_post_dir;}
	
	 fct_type red_sigma_sphere;
	 fct_type red_sigma_cyl;
	 fct_type red_sigma_box;
    
protected:
    std::string  M_uOrder;
    string M_fibers_file;        
    Real M_Chi;
    Real M_Cm;
    Real M_D;
    Real M_sigmal;
    Real M_sigmat;
    Real M_lambda;
    int M_heart_diff_fct; 
    UInt M_verbose;
    string M_post_dir;   
    Real M_endtime;
    string M_fibers_dir;
    bool M_has_fibers;
    
private:


};
//
// IMPLEMENTATION
//


// Constructors
template <typename Mesh>
DataMonodomain<Mesh>::
DataMonodomain(   boost::shared_ptr<HeartFunctors> heart_fct ) :
    DataMesh<Mesh>( heart_fct->_dataFile, "electric/discretization" ),
    DataTime( heart_fct->_dataFile, "electric/discretization" ),
    red_sigma_sphere(heart_fct->get_reduced_sigma_sphere() ),
    red_sigma_cyl(heart_fct->get_reduced_sigma_cylinder() ),
    red_sigma_box(heart_fct->get_reduced_sigma_box() )
{
    setup(heart_fct->_dataFile);
}

template <typename Mesh>
DataMonodomain<Mesh>::
DataMonodomain( const DataMonodomain& dataMonodomain ) :
    DataMesh<Mesh>               ( dataMonodomain ),
    DataTime                     ( dataMonodomain ),
    M_fibers_file(dataMonodomain.M_fibers_file),
    M_Chi(dataMonodomain.M_Chi),
    M_Cm(dataMonodomain.M_Cm),
    M_D(dataMonodomain.M_D),
    M_sigmal(dataMonodomain.M_sigmal),
    M_sigmat(dataMonodomain.M_sigmat),
    M_lambda(dataMonodomain.M_lambda),
    M_heart_diff_fct(dataMonodomain.M_heart_diff_fct),
    M_endtime(dataMonodomain.M_endtime),
    M_post_dir(dataMonodomain.M_post_dir),
    M_uOrder(dataMonodomain.M_uOrder),
    M_has_fibers(dataMonodomain.M_has_fibers)
{
}




template <typename Mesh>
void
DataMonodomain<Mesh>::
setup(  const GetPot& dfile )
{	
	M_Chi = dfile("electric/physics/Chi", 1e3);   // [1e-3 1/cm]    ColliPavarinoTaccardi2005
	M_Cm = dfile("electric/physics/Cm", 1e-3);	  // [1e-3 mF/cm2]  ColliPavarinoTaccardi2005

	if(dfile("electric/physics/ion_model",1) == 1)
		{
			M_D = dfile("electric/physics/D" , 0.0156);   // 0.0156   [1/Ohm/cm]    L^2/T*D,  L=0.099 cm, T=0.63 ms D=1,  //RogersMcCulloch1994
			M_sigmal =  dfile("electric/physics/sigmal", 0.0328);    // 0.0328   [1/Ohm/cm]   sigmal_LR * D_RM/D_LR
			M_sigmat   = dfile("electric/physics/sigmat", 0.00699);    // 0.00699  [1/Ohm/cm]   sigmat_LR * D_RM/D_LR
		}
	else if (dfile("electric/physics/ion_model",1) == 2)
	{
		M_D = dfile("electric/physics/D" , 5.7e-4)/M_Chi/M_Cm;  // 5.7e-4 [1/Ohm/cm]              sigmal/3 + sigmat*2/3
		M_sigmal =  dfile("electric/physics/sigmal", 1.2e-3)/M_Chi/M_Cm;  // 1.2e-3  [1/Ohm/cm]   sigmal_i*sigmal_e/(sigmal_i+sigmal_e)    ColliPavarinoTaccardi2005
		M_sigmat   = dfile("electric/physics/sigmat", 2.56e-4)/M_Chi/M_Cm; // 2.56e-4 [1/Ohm/cm]   sigmat_i*sigmat_e/(sigmat_i+sigmat_e)    ColliPavarinoTaccardi2005
	}
	M_lambda  =  dfile("electric/physics/lambda", 0.66667); // 0.66667 [adim]       sigmal_e/sigmal_i 
    M_heart_diff_fct = dfile("electric/physics/heart_diff_fct",0);
    M_endtime = dfile( "electric/physics/endtime", 300. );
    M_post_dir  = dfile("electric/miscellaneous/post_dir","./");
    M_uOrder = dfile( "electric/discretization/u_order", "P1");
    M_has_fibers = dfile( "electric/discretization/has_fibers", 0);
    if (M_has_fibers)
    {	
    	std::string fibers_dir = dfile( "electric/discretization/fibers_dir", this->meshDir().c_str() );
    	std::string fibers_file = this->meshFile();
    	fibers_file.replace(fibers_file.find(".mesh"), 5, "fibers");
    	M_fibers_file = fibers_dir + dfile( "electric/discretization/fibers_file", fibers_file.c_str() );
      	std::cout<<"Fibers File: "<<M_fibers_file<<std::endl;    	
    }
    else
    {
    	M_fibers_file="";
    	std::cout<<"Fibers not included!"<<std::endl;
    }
}

// Output
template <typename Mesh>
void DataMonodomain<Mesh>::
showMe( std::ostream& c )
{
    c << "\n*** Values for data [fluid/physics]\n\n";
    c << "endtime   = " << M_endtime << std::endl;
    c << "\n*** Values for data [fluid/miscellaneous]\n\n";
    c << "verbose   = " << M_verbose << std::endl;
}



}
#endif





