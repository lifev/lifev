/* -*- mode: c++ -*-

This file is part of the LifeV library

Author(s): Daniele Antonio Di Pietro <dipietro@unibg.it>
Date: 26-1-2005

Copyright (C) 2005 Università degli Studi di Bergamo

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
/**
   \file LevelSetSolver.hpp
   \author Daniele Antonio Di Pietro <dipietro@unibg.it>
   \date 1-28-2005
*/

#ifndef _LEVELSETSOLVER_H_
#define _LEVELSETSOLVER_H_

#include <tab.hpp>

#include <HyperbolicSolverIP.hpp>
#include <LevelSetSolverUtils.hpp>

namespace LifeV {
    /*!
      \class Level set solver
      \brief Level set solver class

      \c A level set solver

      @author Daniele Antonio Di Pietro <dipietro@unibg.it>
      @see
    */
    template<typename MeshType>
    class LevelSetSolver 
        : 
        public HyperbolicSolverIP<MeshType> 
    {
    public:
        /**
           @name Subclasses
        */
        //@{
        template<int NUMNODES = 3, typename PointType = geo_point_type>
        class Face{
        public:
            /** @name Typedefs
             */
            //@{
            typedef geo_point_type point_type;
            typedef geo_point_type vector_type;

            typedef boost::numeric::ublas::bounded_vector<Real, NUMNODES> face_container_type;

            typedef std::vector<point_type> point_container_type;
            typedef UInt id_type;
            //@}

            /** @name Constructors, destructor
             */
            //@{
            Face() {
                _M_points.resize(3);
                _M_has_normal = false;
            }

            Face(ID id) : _M_id(id) {
                _M_points.resize(3);
            }
            //@}

            /** @name Accessors
             */
            //@{

            
            /**
               \Return a const reference to the i-th point
            */
            const point_type& point(const id_type i) const {
                return _M_points[i - 1];
            }
            //@}

            /** @name Mutators
             */
            //@{
            /**
               \Return the i-th point
            */
            point_type& point(const id_type i) {
                return _M_points[i - 1];
            }

            /**
               \Return the id
            */
            id_type& id() {
                return _M_id;
            }

            /**
               \Set i-th point
            */
            void setPoint(const id_type i, const point_type& P) {
                _M_points[i - 1] = P;
            }

            //@}

            /** @name Methods
             */
            //@{

            /**
               \Return the normal versor
            */

            inline vector_type normal() {
                vector_type n;

                if(_M_has_normal)
                    n = _M_normal;
                else {
                    vector_type v1 = _M_points[1] - _M_points[0];
                    vector_type v2 = _M_points[2] - _M_points[0];
                    vector_type n = cross_prod(v1, v2);

                    Real d = boost::numeric::ublas::norm_2(n);            
                    n /= d;
                    _M_normal = n;
                    _M_has_normal = true;
                }
                return n;
            }

            /**
               \Compute the distance between point P and the face_id-th face
            */

            inline Real point_to_face_distance(point_type& P) {
                point_type P0 = _M_points[0];
                vector_type n = normal();
                point_type Q;

                bool projection_is_on_face;

                Real d;

                point_projection_on_plane(P, P0, n, Q);
                projection_is_on_face = point_is_on_face(Q, *this);

                if (projection_is_on_face) {
                    d = point_to_point_distance(P, Q);
                } else {
                    d = point_to_point_distance(P, P0);
                    for(int iP = 2; iP <= 3; iP++)
                        d = std::min( d, point_to_point_distance(P, point(iP)) );
                }

                return d;
            }

            /**
               \Display the face
            */

            void showMe() {
                std::cout << "Face " << _M_id << ":" <<std::endl;
                std::cout << " points" << std::endl
                          << " (1) " << _M_points[0] << std::endl
                          << " (2) " << _M_points[1] << std::endl
                          << " (3) " << _M_points[2] << std::endl;
            }
            //@}
        private:
            //! Point container
            point_container_type _M_points;

            //! Face id
            id_type _M_id;

            //! The normal was already computed?
            bool _M_has_normal;

            //! The normal to the face
            vector_type _M_normal;
        };
        //@}

        /**
           @name Friend classes
        */
        //@{
        //friend class LevelSetSolver<MeshType>;
        //@}

        /** @name Typedefs
         */
        //@{
        typedef MeshType mesh_type;
        typedef typename HyperbolicSolverIP<mesh_type>::velocity_type velocity_type;

        typedef geo_point_type point_type;
        typedef geo_point_type vector_type;
        typedef Face<3> face_type;

        typedef std::vector<point_type> point_list_type;
        typedef std::list<face_type> face_list_type;
        typedef typename face_list_type::iterator face_list_iterator;

        typedef std::vector<UInt> orientation_type;

        typedef typename LevelSetSolver<MeshType>::u_type lsfunction_type;

        typedef enum {fluid1 = 1, fluid2 = 2} fluid_type;

        //@}

        /** @name Constructors and destructors
         */
        //@{
        LevelSetSolver(const GetPot& data_file,
                       const std::string& data_section,
                       const RefFE& reffe, 
                       const QuadRule& qr, 
                       const QuadRule& qr_bd, 
                       const BCHandler& bc_h,
                       CurrentFE& fe_velocity,
                       const Dof& dof_velocity,
                       velocity_type& velocity0)
            :
            HyperbolicSolverIP<mesh_type>(data_file,
                                          data_section,
                                          reffe, 
                                          qr, 
                                          qr_bd, 
                                          bc_h,
                                          fe_velocity,
                                          dof_velocity,
                                          velocity0) 
        {
            _M_tol = 1e-6;
        }
        ~LevelSetSolver() {}
        //@}

        /** @name Accessors
         */
        //@{

        /**
           \Return the level set function
        */

        lsfunction_type const & lsfunction() const {
            return u();
        }
        //@}

        /** @name Mutators
         */
        //@{
        void setTol(Real tol) {
            _M_tol = tol;
        }
        //@}

        /** @name Methods
         */
        //@{

        /**
           \Reinitialize the interface using direct method
        */

        void directReinitialization() {
            _M_u = _M_bdf.extrap();

            build_interface();

            for(UInt iP = 1; iP < _mesh.numPoints(); iP++) {

                point_type P;
                convert_point_type(P, _mesh.pointList(iP));

                Real d = _M_face_list.begin()->point_to_face_distance(P);
            
                for(face_list_iterator faces_it = _M_face_list.begin(); faces_it != _M_face_list.end(); faces_it++) {
                    Real d1 = faces_it->point_to_face_distance(P);
                    d = std::min(d, d1);
                }

                _M_u[iP] = d;
            }
        }
  
        /**
           \Compute the mass of i-th fluid
        */

        Real computeMassOfFluid(fluid_type fluid_id) {
            Real mass = 0.;
            Real ls_fun;
            int jg;

            for(UInt i = 1; i <= _mesh.numVolumes(); i++) {
                _M_fe.updateJac( _mesh.volumeList(i) );

                for(int iq = 0; iq < _M_fe.nbQuadPt; iq++) {
                    ls_fun = 0.;
                    for(int j = 0; j < _M_fe.nbNode; j++) {
                        jglo = _M_dof.localToGlobal(i, j + 1) - 1;
                        ls_fun += _M_fe.phi(j, iq) * U(jglo);
                    }
                    if(fluid_id == fluid1)
                        mass += (ls_fun < 0 ? 1. : 0.) * _M_fe.weightDet(iq);
                    else
                        mass += (ls_fun >= 0 ? 1. : 0.) * _M_fe.weightDet(iq);
                }
            }
            return mass;
        }
        //@}

    private:
        //! The face list
        face_list_type _M_face_list;

        //! Tolerance on interface position
        Real _M_tol;

        /**
           \Build the interface given the level set function
        */

        inline void build_interface() {
            for(UInt iV = 1; iV <= _mesh.numVolumes(); iV++) {

                point_list_type locPointList;

                for(UInt ie = 1; ie <= _mesh.numLocalEdges(); ie++) {
                    
                    UInt j1 = mesh_type::ElementShape::eToP(ie, 1);
                    UInt j2 = mesh_type::ElementShape::eToP(ie, 2);

                    UInt jg1 = _M_dof.localToGlobal(iV, j1) - 1;
                    UInt jg2 = _M_dof.localToGlobal(iV, j2) - 1;

                    Real u1 = _M_u[jg1];
                    Real u2 = _M_u[jg2];
                  
                    if(u1 * u2 < 0) {
                        point_type P;
                        point_type P1;
                        point_type P2;
                        
                        P1[0] = _mesh.pointList(jg1 + 1).x();
                        P1[1] = _mesh.pointList(jg1 + 1).y(); 
                        P1[2] = _mesh.pointList(jg1 + 1).z();

                        P2[0] = _mesh.pointList(jg2 + 1).x();
                        P2[1] = _mesh.pointList(jg2 + 1).y(); 
                        P2[2] = _mesh.pointList(jg2 + 1).z();

                        find_zero_on_edge(P1, u1, P2, u2, P);

                        points_coincide<point_type> compare(_M_tol);
                        bool isAlreadyThere = false;

                        for(point_list_type::iterator point_it = locPointList.begin(); point_it != locPointList.end(); point_it++)
                            if(compare(P, *point_it)) {
                                isAlreadyThere = true;
                                break;
                            }

                        if(!isAlreadyThere)
                            locPointList.push_back(P);
                    }
                }

                // Handle several cases according to the number of local 
                // interface points found in the previous step

                switch( locPointList.size() ) {
                case 3: {
                    face_type face;

                    int i = 1;
                    for(point_list_type::iterator point_it = locPointList.begin(); point_it != locPointList.end(); point_it++) {
                        point_type P(*point_it);
                        face.setPoint(i, P);
                        i++;
                    }

                    _M_face_list.push_back(face);

                    break;
                }
                case 4: {
                    face_type face1;
                    face_type face2;

                    std::vector<UInt> orientation(4);
                    estabilish_orientation_4pt(locPointList, orientation);

                    face1.setPoint(1, locPointList[orientation[0]]);
                    face1.setPoint(2, locPointList[orientation[1]]);
                    face1.setPoint(3, locPointList[orientation[2]]);

                    face2.setPoint(1, locPointList[orientation[2]]);
                    face2.setPoint(2, locPointList[orientation[3]]);
                    face2.setPoint(3, locPointList[orientation[0]]);

                    _M_face_list.push_back(face1);
                    _M_face_list.push_back(face2);

                    break;
                }
                }
            }
        }

        /**
           \Given a vector of 4 nodes, it returns a point ordering suitable
           \to obtain a proper quadrilateral. Orientation list is given below:
           \1 : [P0, P1, P2, P3]
           \2 : [P0, P1, P3, P2]
           \3 : [P0, P2, P1, P3]
         */
        inline void estabilish_orientation_4pt(point_list_type& v, orientation_type& orientation) {
            point_type P0 = v[0];
            point_type P1 = v[1];
            point_type P2 = v[2];
            point_type P3 = v[3];

            vector_type v1 = P1 - P0;
            vector_type v2 = P2 - P0;
            vector_type v3 = P3 - P0;
            vector_type n = cross_prod(v1, v2);

            Real a12 = inner_prod(cross_prod(v1, v2), n);
            Real a13 = inner_prod(cross_prod(v1, v3), n);
            Real a23 = inner_prod(cross_prod(v2, v3), n);

            if (a12 * a13 > 0)
                if (a23 > 0) {
                    orientation[0] = 0;
                    orientation[1] = 1;
                    orientation[2] = 2;
                    orientation[3] = 3;
                }
                else {
                    orientation[0] = 0;
                    orientation[1] = 1;
                    orientation[2] = 3;
                    orientation[3] = 2;
                }
            else {
                orientation[0] = 0;
                orientation[1] = 2;
                orientation[2] = 1;
                orientation[3] = 3;
         
            }
        }
    };
}
#endif
