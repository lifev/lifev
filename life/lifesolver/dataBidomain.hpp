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
  \mantainer R. Ruiz
   \date 2010-04
*/
#ifndef _DATABIDOMAIN_H_
#define _DATABIDOMAIN_H_
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifecore/dataString.hpp>
#include <life/lifearray/tab.hpp>
#include <life/lifesolver/heartFunctors.hpp>


namespace LifeV
{
using namespace std;

/*!
  \class DataBidomain

  Base class which holds usual data for the Bidomain model solvers

*/
template <typename Mesh>
class DataBidomain:
        public DataMesh<Mesh>,
        public DataTime
{
public:

    //! Constructors
    DataBidomain( boost::shared_ptr<HeartFunctors> heart_fct );
    DataBidomain( const DataBidomain& dataBidomain );

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const&, Real const&)> fct_type;


    //! Ouptut
    void showMe( std::ostream& c = std::cout );


    //! external setup
    void setup( const GetPot& dfile );

    //! verbose
    Real verbose() const {return M_verbose;};


    //! FE space order
    std::string uOrder() const{return M_uOrder;};
    //! Chi
	inline   Real Chi() const {return M_Chi;}
	 //! Chi
	inline   string fibers_file() const {return M_fibers_file;}

	inline int heart_diff_fct() const {return M_heart_diff_fct;}

	inline bool has_fibers() const {return M_has_fibers;}

    //! format vct
	inline bool fibers_format() const{return M_fibers_format;}

	//! sigma_l
	inline   Real sigmal_i() const 	{return M_sigmal_i;}
	inline   Real sigmal_e() const 	{return M_sigmal_e;}

	//! sigma_t
	inline   Real sigmat_i() const 	{return M_sigmat_i;}
	inline   Real sigmat_e() const 	{return M_sigmat_e;}

	//! Cm
	inline   Real Cm() const 	{return M_Cm;}
	//! D
	inline   Real D_i() const 	{return M_D_i;}
	//! Post_dir
	inline   Real D_e() const 	{return M_D_e;}
		//! Post_dir
	inline   string Post_dir() const {return M_post_dir;}

	inline	 UInt order_bdf() const {return M_order_bdf;}

	fct_type red_sigma_sphere;
	fct_type red_sigma_cyl;
	fct_type red_sigma_box;


protected:

    std::string  M_uOrder;
    string M_fibers_file;
    Real M_Chi;
    Real M_Cm;
    Real M_D_i;
    Real M_D_e;
    Real M_sigmal_i;
    Real M_sigmat_i;
    Real M_sigmal_e;
    Real M_sigmat_e;
    int M_heart_diff_fct;
    UInt M_verbose;
    string M_post_dir; //! full name
    string M_fibers_dir;
    bool M_has_fibers;
    // format of fibers file
    bool M_fibers_format;
    // order of time discretization BDF order
    UInt M_order_bdf;

private:


};

//
// IMPLEMENTATION
//


// Constructors
template <typename Mesh>
DataBidomain<Mesh>::
DataBidomain( boost::shared_ptr<HeartFunctors> heart_fct ) :
    DataMesh<Mesh>( heart_fct->_dataFile, "electric/space_discretization" ),
    DataTime( heart_fct->_dataFile, "electric/time_discretization" ),
    red_sigma_sphere(heart_fct->get_reduced_sigma_sphere() ),
    red_sigma_cyl(heart_fct->get_reduced_sigma_cylinder() ),
    red_sigma_box(heart_fct->get_reduced_sigma_box() )
{
    setup(heart_fct->_dataFile);
}


template <typename Mesh>
DataBidomain<Mesh>::
DataBidomain( const DataBidomain& dataBidomain ) :
    DataMesh<Mesh>               ( dataBidomain ),
    DataTime                     ( dataBidomain ),
    M_fibers_file(dataBidomain.M_fibers_file),
    M_Chi(dataBidomain.M_Chi),
    M_Cm(dataBidomain.M_Cm),
    M_D_i(dataBidomain.M_D_i),
    M_D_e(dataBidomain.M_D_e),
    M_sigmal_i(dataBidomain.M_sigmal_i),
    M_sigmat_i(dataBidomain.M_sigmat_i),
    M_sigmal_e(dataBidomain.M_sigmal_e),
    M_sigmat_e(dataBidomain.M_sigmat_e),
    M_heart_diff_fct(dataBidomain.M_heart_diff_fct),
    M_verbose(dataBidomain.M_verbose),
    M_post_dir(dataBidomain.M_post_dir),
    M_uOrder(dataBidomain.M_uOrder),
    M_has_fibers(dataBidomain.M_has_fibers),
    M_fibers_format(dataBidomain.M_fibers_format),
    M_order_bdf(dataBidomain.M_order_bdf)
{
}


template <typename Mesh>
void
DataBidomain<Mesh>::
setup(  const GetPot& dfile )
{
	M_Chi = dfile("electric/physics/Chi",1e3); 	// [1e-3 1/cm]   ColliPavarinoTaccardi2005
	M_Cm = dfile("electric/physics/Cm",1e-3);  	// [1e-3 mF/cm2]   ColliPavarinoTaccardi2005
    	if(dfile("electric/physics/ion_model",1) == 1)
		{
    		M_D_i = dfile("electric/physics/D_i" ,3.3e-2);       // 3.3e-2  [1/Ohm/cm]   D_i_LR * D_RM/D_LR  see dataMonodomain
    		M_D_e = dfile("electric/physics/D_e" ,4.29e-2);      // 4.29e-2 [1/Ohm/cm]	D_e_LR * D_RM/D_LR  see dataMonodomain
    		M_sigmal_i =  dfile("electric/physics/sigmal_i", 8.19e-2); // 8.19e-2 [1/Ohm/cm]   sigmal_i_LR * D_RM/D_LR
    		M_sigmat_i   = dfile("electric/physics/sigmat_i", 8.6e-3); 	// 8.6e-3  [1/Ohm/cm]   sigmat_i_LR * D_RM/D_LR
    		M_sigmal_e =  dfile("electric/physics/sigmal_e", 5.46e-2); // 5.46e-2 [1/Ohm/cm]	sigmal_e_LR * D_RM/D_LR
    		M_sigmat_e   = dfile("electric/physics/sigmat_e",3.69e-2); // 3.69e-2 [1/Ohm/cm]	sigmat_e_LR * D_RM/D_LR
		}
	else if(dfile("electric/physics/ion_model",1) == 2)
		{
		M_D_i = dfile("electric/physics/D_i" , 1.21e-3);  		// sigmal_i/3 + sigmat_i*2/3
		M_D_e = dfile("electric/physics/D_e" , 1.57e-3);  		// sigmal_e/3 + sigmat_e*2/3
		M_sigmal_i =  dfile("electric/physics/sigmal_i", 3e-3);  		// 3e-3      [1/Ohm/cm]   ColliPavarinoTaccardi2005
		M_sigmat_i   = dfile("electric/physics/sigmat_i", 3.1525e-4); 	// 3.1525e-4 [1/Ohm/cm]   ColliPavarinoTaccardi2005
		M_sigmal_e =  dfile("electric/physics/sigmal_e", 2e-3); 		// 2e-3      [1/Ohm/cm]   ColliPavarinoTaccardi2005
		M_sigmat_e   = dfile("electric/physics/sigmat_e",1.3514e-3); 	// 1.3514e-3 [1/Ohm/cm]   ColliPavarinoTaccardi2005
		}
	else if(dfile("electric/physics/ion_model",1) == 3)
                {
                M_D_i = dfile("electric/physics/D_i" ,3.3e-2);       // 3.3e-2  [1/Ohm/cm]   D_i_LR * D_RM/D_LR  see dataMonodomain
                M_D_e = dfile("electric/physics/D_e" ,4.29e-2);      // 4.29e-2 [1/Ohm/cm]      D_e_LR * D_RM/D_LR  see dataMonodomain
                M_sigmal_i =  dfile("electric/physics/sigmal_i", 8.19e-2); // 8.19e-2 [1/Ohm/cm]   sigmal_i_LR * D_RM/D_LR
                M_sigmat_i   = dfile("electric/physics/sigmat_i", 8.6e-3);      // 8.6e-3  [1/Ohm/cm]   sigmat_i_LR * D_RM/D_LR
                M_sigmal_e =  dfile("electric/physics/sigmal_e", 5.46e-2); // 5.46e-2 [1/Ohm/cm]        sigmal_e_LR * D_RM/D_LR
                M_sigmat_e   = dfile("electric/physics/sigmat_e",3.69e-2); // 3.69e-2 [1/Ohm/cm]        sigmat_e_LR * D_RM/D_LR
                }
    M_heart_diff_fct 	= dfile("electric/physics/heart_diff_fct",0);
    M_verbose 		= dfile( "electric/miscellaneous/verbose", 1 );
    M_post_dir  	= dfile("electric/miscellaneous/post_dir","./");
    M_uOrder 		= dfile( "electric/space_discretization/u_order", "P1");
    M_fibers_format 	= dfile("electric/space_discretization/fibers_format",0);
    M_has_fibers 	= dfile( "electric/space_discretization/has_fibers", 0);
    M_order_bdf       	= dfile("electric/time_discretization/BDF_order",1);
    if (M_has_fibers)
    {
        std::string fibers_dir = dfile( "electric/space_discretization/fibers_dir", this->meshDir().c_str() );
        std::string fibers_file = this->meshFile();
        fibers_file.replace(fibers_file.find(".mesh"), 5, "fibers");
        M_fibers_file = fibers_dir + dfile( "electric/space_discretization/fibers_file", fibers_file.c_str() );
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
void DataBidomain<Mesh>::
showMe( std::ostream& c )
{
    c << "\n*** Values for data [fluid/physics]\n\n";
    c << "endtime   = " << getEndTime() << std::endl;
    c << "\n*** Values for data [fluid/miscellaneous]\n\n";
    c << "verbose   = " << M_verbose << std::endl;
}



}
#endif
