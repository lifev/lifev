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
  DataOneDModel::DataOneDModel(const GetPot& dfile, const std::string& section):
    M_time_beg               (dfile((section + "/time/timebeg").data()                      ,0.0)),
    M_time_end               (dfile((section + "/time/timeend").data()                      ,1.0)),
    M_time_step              (dfile((section + "/time/timestep").data()                     ,0.1)),
        //! Discretization
    M_mesh_file              (dfile((section + "/discretization/mesh_file").data()          ,"mesh.mesh")),
    M_mesh_dir               (dfile((section + "/discretization/mesh_dir").data()           ,"./")),
    M_x_left                 (dfile((section + "/discretization/x_left").data()             ,0.)),
    M_x_right                (dfile((section + "/discretization/x_right").data()            ,1.)),
    M_nb_elem                (dfile((section + "/discretization/nb_elem").data()            ,10)),
        //! Miscellaneous
    M_post_dir               (dfile((section + "/miscellaneous/post_dir").data()            ,"./")),
    M_post_file              (dfile((section + "/miscellaneous/post_file").data()           ,"sol")),
    M_verbose                (dfile((section + "/miscellaneous/verbose").data()             ,1)),
    M_postProcessTimeStep    (dfile((section + "/miscellaneous/postprocess_timestep").data(),0.01)),
    M_adaptiveMesh           (dfile((section + "/problem/adaptiveMesh").data()              , 0)),
    M_CFL                    (dfile((section + "miscellaneous/showCFL").data(),                        0)),
    M_UW                     (dfile((section + "miscellaneous/alternate_solver").data(),               0)),
    M_inertial_wall          (dfile((section + "miscellaneous/inertial_wall").data(),                  0)),
    M_viscoelastic_wall      (dfile((section + "miscellaneous/viscoelastic_wall").data(),              0)),
    M_linearize_string_model (dfile((section + "miscellaneous/linearize_string_model").data(),         1)),
    M_linearize_equations    (dfile((section + "miscellaneous/linearize_equations").data(),            0)),
    M_longitudinal_wall      (dfile((section + "miscellaneous/longitudinal_wall").data(),              0)),
    M_flux_second_der        (dfile((section + "miscellaneous/compute_flux_second_derivative").data(), 0)),
    M_dP_dt_steps            (dfile((section + "miscellaneous/pressure_derivative_steps").data(),      1)),
    M_firstNode              (dfile((section + "/initialize/firstnode").data(),                        1)),
    M_lastNode               (dfile((section + "/initialize/lastnode").data(),                         2)),
    M_initVar                (dfile((section + "/initialize/var").data(),                            "P")),
    M_initValue              (dfile((section + "/initialize/init_value").data(),                      0.)),
    M_restValue              (dfile((section + "/initialize/rest_value").data(),                      0.)),
    M_multiplier             (dfile((section + "/initialize/multiplier").data(),                      1.)),
    M_width                  (dfile((section + "/initialize/width").data(),                           5.)),
    M_mesh                   ( )
{

    std::cout << "    1d- Mesh setup ... " << std::endl;

    M_mesh.reset(new DataOneDModel::mesh_raw_type);

    M_mesh->setUp( M_x_left, M_x_right, M_nb_elem );
    std::cout << "    ok." << std::endl;
}


Real const& DataOneDModel::timestep() const
{
    return M_time_step;
}
Real const& DataOneDModel::inittime() const
{
    return M_time_beg;
}
Real const& DataOneDModel::endtime() const
{
    return M_time_end;
}

void DataOneDModel::settimestep( const Real& dt )
{
    M_time_step = dt;
}
void DataOneDModel::setinittime( const Real& t0 )
{
    M_time_beg = t0;
}
void DataOneDModel::setendtime( const Real& T )
{
    M_time_end = T;
}

Real DataOneDModel::xLeft() const
{
    return M_x_left;
}
Real DataOneDModel::xRight() const
{
    return M_x_right;
}
Real DataOneDModel::nbElem() const
{
    return M_nb_elem;
}


void DataOneDModel::showMe(std::ostream& c) const
{
    //! Time
    c << "\n*** Values for data [time]\n";
    c << "time_beg = " << M_time_beg << "\n";
    c << "time_end = " << M_time_end << "\n";
    c << "time_step = " << M_time_step << "\n";

    //! Discretization
    c << "\n*** Values for data [discretization]\n";
    c << "mesh_file = " << M_mesh_file << "\n";
    c << "mesh_dir  = " << M_mesh_dir << "\n";
    c << "x_left = " << M_x_left << "\n";
    c << "x_right = " << M_x_right << "\n";
    c << "nb_elem = " << M_nb_elem << "\n";

    //! Miscellaneous
    c << "\n*** Values for data [miscellaneous]\n\n";
    c << "post_dir         = " << M_post_dir << "\n";
    c << "post_file        = " << M_post_file << "\n";
    c << "verbose          = " << M_verbose << std::endl;
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
