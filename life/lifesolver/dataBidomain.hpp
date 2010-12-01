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

class DataBidomain:
        public DataMesh,
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
    std::string uOrder() const {return M_uOrder;};
    //! Chi
    inline   Real Chi() const {return M_Chi;}
    //! Chi
    inline   string fibers_file() const {return M_fibers_file;}

    inline int heart_diff_fct() const {return M_heart_diff_fct;}

    inline bool has_fibers() const {return M_has_fibers;}

    //! format vct
    inline bool fibers_format() const {return M_fibers_format;}

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

}
#endif
