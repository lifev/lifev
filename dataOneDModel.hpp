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
    DataOneDModel(const GetPot& dfile, const std::string& section = "");
//     //DataOneDModel(const DataOneDModel& dataOneDModel);

    //    DataOneDModel(const DataOneDModel& dataOneDModel);

    void timeAdvance(){M_time += M_time_step;}

    // return the different time data
    Real const& timestep() const;
    Real const& inittime() const;
    Real const& endtime()  const;
    Real const& time()     const;
    // set the different time data
    void settimestep( const Real& dt );
    void setinittime( const Real& t0 );
    void setendtime ( const Real& T );
    void settime    ( const Real& time);
    // return the two extremal points
    Real xLeft() const;
    Real xRight() const;
    Real nbElem() const;

    const std::string PostFile     ()         const {return M_post_file;}
    const std::string PostDirectory()         const {return M_post_dir;}

    const int         verbose()               const {return M_verbose;}

    const Real        CFL()                   const {return M_CFL;}
    const bool        UW()                    const {return M_UW;}
    const bool        inertialWall()          const {return M_inertial_wall;}
    const bool        viscoelasticWall()      const {return M_viscoelastic_wall;}
    const bool        linearizeStringModel()  const {return M_linearize_string_model;}
    const bool        linearizeEquations()    const {return M_linearize_equations;}
    const bool        longitudinalWall()      const {return M_longitudinal_wall;}

    const bool        fluxSecondDer()         const {return M_flux_second_der;}

    const int         DPdtSteps()             const {return M_dP_dt_steps;}

    const std::string initVar()               const {return M_initVar;}
    const int         firstNode()             const {return M_firstNode;}
    const int         lastNode()              const {return M_lastNode;}
    const Real        restValue()             const {return M_restValue;}
    const Real        initValue()             const {return M_initValue;}
    const Real        multiplier()            const {return M_multiplier;}
    const Real        width()                 const {return M_width;}


    mesh_type mesh() {return M_mesh;}
    //BasicOneDMesh& mesh() {return M_mesh;}
    //! Output
    void showMe(std::ostream& c=std::cout) const;

protected:

    //
    //! Time
    //
    Real            M_time_beg;  //! starting time
    Real            M_time_end; //! finishing time
    Real            M_time_step; //! time step
    Real            M_time;
    //
    //! Discretization
    //
    std::string     M_mesh_file;
    std::string     M_mesh_dir;
    Real            M_x_left;  //! left coordinate
    Real            M_x_right; //! right coordinate

    UInt            M_nb_elem;    //! number of elements
    //
    //! Miscellaneous
    //
    std::string     M_post_dir; //! full directory name (including path)
    std::string     M_post_file; //! output file name
    int             M_verbose;
    //
    Real            M_CFL;
    bool            M_UW;
    //! boolean: activate inertial/ viscoelastic/ longitudinal term in pressure-area relationship?
    bool            M_inertial_wall;
    bool            M_viscoelastic_wall;
    bool            M_linearize_string_model;
    bool            M_linearize_equations;
    bool            M_longitudinal_wall;

    //! boolean: compute second spatial derivative of flux?
    bool            M_flux_second_der;
    //! approximation of pressure temporal derivative: how many time steps?
    int             M_dP_dt_steps;
    //! initialize
    std::string     M_initVar;
    int             M_firstNode;
    int             M_lastNode;
    Real            M_restValue;
    Real            M_initValue;
    Real            M_multiplier;
    Real            M_width;
    /*!
      Write down results to file each M_postProcessTimeStep seconds
    */
    Real            M_postProcessTimeStep;
    /*!
      Mesh
     */
    const UInt      M_adaptiveMesh;
    //BasicOneDMesh   M_mesh;
    mesh_type       M_mesh;

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
