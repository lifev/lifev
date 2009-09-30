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
#include <life/lifecore/life.hpp>
#include <life/lifecore/GetPot.hpp>
//#include <life/lifemesh/basicOneDMesh.hpp>
//#include <life/lifemesh/regionMesh2D.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
//! boost shared pointers
#include <boost/shared_ptr.hpp>

#include <life/lifealg/SolverAmesos.hpp>

#include <lifemc/lifemesh/regionMesh1D.hpp>
//#include "life/lifemesh/regionMesh2D.hpp"
//#include "life/lifemesh/regionMesh3D.hpp"

namespace ublas = boost::numeric::ublas;


namespace LifeV
{
/*!
  \class DataOneDModel

  Base class which holds usual data for the Blood flow One Dimensional Model solvers

*/

class DataOneDModel
{
public:

    typedef RegionMesh1D     <LinearLine>        mesh_raw_type;
    typedef boost::shared_ptr<mesh_raw_type>     mesh_type;


    //! 2D vector
    //  typedef std::pair< Real, Real > Vec2D;

    typedef ublas::bounded_array<Real, 2>        Vec2D;

    typedef ublas::vector<double>                ScalVec;
    typedef ublas::vector<ScalVec>               ScalVec_vector;
    typedef ublas::vector_range<ScalVec_vector>  ScalVec_vector_range;

    typedef SolverAmesos                         solver_type;
    typedef solver_type::vector_type             vector_type;

    typedef solver_type::matrix_type             matrix_type;
    typedef boost::shared_ptr<matrix_type>       matrix_ptrtype;

//     //! Constructors
    DataOneDModel(const GetPot& dfile);
//     //DataOneDModel(const DataOneDModel& dataOneDModel);

    //    DataOneDModel(const DataOneDModel& dataOneDModel);


    // return the different time data
    Real const& timestep() const;
    Real const& inittime() const;
    Real const& endtime()  const;

    // set the different time data
    void settimestep( const Real& dt );
    void setinittime( const Real& t0 );
    void setendtime( const Real& T );

    // return the two extremal points
    Real xLeft() const;
    Real xRight() const;
    Real nbElem() const;

    const std::string PostFile     () const {return _M_post_file;}
    const std::string PostDirectory() const {return _M_post_dir;}

    const int verbose()               const {return _M_verbose;}

    mesh_type mesh() {return _M_mesh;}
    //BasicOneDMesh& mesh() {return _M_mesh;}
    //! Output
    void showMe(std::ostream& c=std::cout) const;

protected:

    //
    //! Time
    //
    Real            _M_time_beg;  //! starting time
    Real            _M_time_end; //! finishing time
    Real            _M_time_step; //! time step
    //
    //! Discretization
    //
    std::string     _M_mesh_file;
    std::string     _M_mesh_dir;
    Real            _M_x_left;  //! left coordinate
    Real            _M_x_right; //! right coordinate
    UInt            _M_nb_elem;    //! number of elements
    //
    //! Miscellaneous
    //
    std::string     _M_post_dir; //! full directory name (including path)
    std::string     _M_post_file; //! output file name
    int             _M_verbose;
    /*!
      Write down results to file each _M_postProcessTimeStep seconds
    */
    Real            _M_postProcessTimeStep;
    /*!
      Mesh
     */
    const UInt      _M_adaptiveMesh;
    //BasicOneDMesh   _M_mesh;
    mesh_type       _M_mesh;


};



Real dot(const DataOneDModel::Vec2D& vec1,
         const DataOneDModel::Vec2D& vec2);
// {
//     ASSERT_PRE( vec1.size() == 2 && vec2.size() == 2,
//                 "dot : dot works only for 2D vectors." )
//         return vec1[0]*vec2[0] + vec1[1]*vec2[1];
// }

}
#endif
