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
#ifndef _DATAMONODOMAIN_H_
#define _DATAMONODOMAIN_H_
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
  \class DataMonodomain

  Base class which holds usual data for the Monodomain model solvers

*/
class DataMonodomain:
        public DataMesh,
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
    string M_fibers_dir;
    bool M_has_fibers;

private:


};

}
#endif
