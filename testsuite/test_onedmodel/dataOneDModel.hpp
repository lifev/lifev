/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/*!
  \file dataOneDModel.hpp
  \author Vincent Martin
  \date 07/2004
  \version 1.0

  \brief File containing a class for data in oned model

*/
#ifndef _DATAONEDMODEL_H_
#define _DATAONEDMODEL_H_
#include <string>
#include <iostream>
#include "life.hpp"
#include "GetPot.hpp"

namespace LifeV
{
/*!
  \class DataOneDModel

  Base class which holds usual data for the Blood flow One Dimensional Model solvers

*/
class DataOneDModel
{
public:

    //! Constructor
    DataOneDModel(const GetPot& dfile);

    // return the different time data
    Real timestep() const;
    Real inittime() const;
    Real endtime()  const;

    // return the two extremal points
    Real xLeft() const;
    Real xRight() const;

    const std::string PostDirectory() const {return _M_post_dir;};

    //! Output
    void showMeData(std::ostream& c=std::cout) const;

protected:

    //! boundcond (boundary conditions)
    int _M_test_case; //! test case

    //
    //! Time
    //
    Real _M_time_beg;  //! starting time
    Real _M_time_end; //! finishing time
    Real _M_time_step; //! time step
    //
    //! Discretization
    //
    std::string _M_mesh_file;
    std::string _M_mesh_dir;
    Real _M_x_left;  //! left coordinate
    Real _M_x_right; //! right coordinate
    UInt _M_nb_elem;    //! number of elements
    //
    //! Miscellaneous
    //
    std::string _M_post_dir; //! full directory name (including path)
    std::string _M_post_file; //! output file name
    int _M_verbose;
};
}
#endif
