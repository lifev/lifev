/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

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
  \file dataOneDModel.cpp
  \author Vincent Martin
  \date 07/2004
  \version 1.0

  \brief Implementation of the class for the data in oned model

*/
#include <lifemc/lifesolver/dataOneDModel.hpp>


namespace LifeV
{

//! constructor using a data file.
DataOneDModel::DataOneDModel(const GetPot& dfile):
        _M_time_beg           (dfile("time/timebeg",0.0)),
        _M_time_end           (dfile("time/timeend",1.0)),
        _M_time_step          (dfile("time/timestep",0.1)),
    //! Discretization
        _M_mesh_file          (dfile("discretization/mesh_file","mesh.mesh")),
        _M_mesh_dir           (dfile("discretization/mesh_dir","./")),
        _M_x_left             (dfile("discretization/x_left",0.)),
        _M_x_right            (dfile("discretization/x_right",1.)),
        _M_nb_elem            (dfile("discretization/nb_elem",10)),
//! Miscellaneous
        _M_post_dir           (dfile("miscellaneous/post_dir","./")),
        _M_post_file          (dfile("miscellaneous/post_file","sol")),
        _M_verbose            (dfile("miscellaneous/verbose",1)),
        _M_postProcessTimeStep(dfile("miscellaneous/postprocess_timestep",0.01)),
        _M_adaptiveMesh       (dfile("problem/adaptiveMesh", 0)),
        _M_mesh               ( )
{


    std::cout << "    Mesh setup ... " << std::endl;

    _M_mesh.reset(new DataOneDModel::mesh_raw_type);

    _M_mesh->setUp( _M_x_left, _M_x_right, _M_nb_elem );
//                     dfile("problem/alpha",50),
//                     dfile("problem/delta",10),
//                     dfile("problem/order",5),
//                     dfile("problem/minDeltaX",1.));
    std::cout << "    ok." << std::endl;
    //! Time
//     _M_time_beg  = dfile("time/timebeg",0.0);
//     _M_time_end  = dfile("time/timeend",1.0);
//     _M_time_step = dfile("time/timestep",0.1);
//     //! Discretization
//     _M_mesh_file = dfile("discretization/mesh_file","mesh.mesh");
//     _M_mesh_dir  = dfile("discretization/mesh_dir","./");
//     _M_x_left    = dfile("discretization/x_left",0.);
//     _M_x_right   = dfile("discretization/x_right",1.);
//     _M_nb_elem   = dfile("discretization/nb_elem",10);
//     //! Miscellaneous
//     _M_post_dir  = dfile("miscellaneous/post_dir","./");
//     _M_post_file = dfile("miscellaneous/post_file","sol");
//     _M_verbose   = dfile("miscellaneous/verbose",1);
//     _M_postProcessTimeStep = dfile("miscellaneous/postprocess_timestep",0.01);
}


// DataOneDModel::DataOneDModel(const DataOneDModel& dataOneDModel):
//     //! Time
//         _M_time_beg            = dataOneDModel._M_time_beg,
//         _M_time_end            = dataOneDModel._time_end,
//         _M_time_step           = dataOneDModel._M_time_step,
//         //! Discretization
//         _M_mesh_file           = dataOneDModel._M_mesh_file,
//         _M_mesh_dir            = dataOneDModel._M_mesh_dir,
//         _M_x_left              = dataOneDModel._M_x_left,
//         _M_x_right             = dataOneDModel._M_x_right,
//         _M_nb_elem             = dataOneDModel._M_nb_elem,
//     //! Miscellaneous
//         _M_post_dir            = dataOneDModel._M_post_dir,
//         _M_post_file           = dataOneDModel._M_post_file,
//         _M_verbose             = dataOneDModel._M_verbose,
//         _M_postProcessTimeStep = dataOneDModel._M_postProcessTimeStep
// }
// {



Real const& DataOneDModel::timestep() const
{
    return _M_time_step;
}
Real const& DataOneDModel::inittime() const
{
    return _M_time_beg;
}
Real const& DataOneDModel::endtime() const
{
    return _M_time_end;
}

void DataOneDModel::settimestep( const Real& dt )
{
    _M_time_step = dt;
}
void DataOneDModel::setinittime( const Real& t0 )
{
    _M_time_beg = t0;
}
void DataOneDModel::setendtime( const Real& T )
{
    _M_time_end = T;
}

Real DataOneDModel::xLeft() const
{
    return _M_x_left;
}
Real DataOneDModel::xRight() const
{
    return _M_x_right;
}
Real DataOneDModel::nbElem() const
{
    return _M_nb_elem;
}


void DataOneDModel::showMe(std::ostream& c) const
{
    //! Time
    c << "\n*** Values for data [time]\n";
    c << "time_beg = " << _M_time_beg << "\n";
    c << "time_end = " << _M_time_end << "\n";
    c << "time_step = " << _M_time_step << "\n";

    //! Discretization
    c << "\n*** Values for data [discretization]\n";
    c << "mesh_file = " << _M_mesh_file << "\n";
    c << "mesh_dir  = " << _M_mesh_dir << "\n";
    c << "x_left = " << _M_x_left << "\n";
    c << "x_right = " << _M_x_right << "\n";
    c << "nb_elem = " << _M_nb_elem << "\n";

    //! Miscellaneous
    c << "\n*** Values for data [miscellaneous]\n\n";
    c << "post_dir         = " << _M_post_dir << "\n";
    c << "post_file        = " << _M_post_file << "\n";
    c << "verbose          = " << _M_verbose << std::endl;
}



Real dot(const DataOneDModel::Vec2D& vec1,
         const DataOneDModel::Vec2D& vec2)
{
    ASSERT_PRE( vec1.size() == 2 && vec2.size() == 2,
                "dot works only for 2D vectors" );

    //boost::numeric::ublas::vector<double>   vec;

    return vec1[0]*vec2[0] + vec1[1] * vec2[1];

}


}
