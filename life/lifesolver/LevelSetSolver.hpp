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
/*!
  \file LevelSetSolver.hpp
  \author Daniele Antonio Di Pietro <dipietro@unibg.it>
  \date 1-28-2005
*/

#ifndef _LEVELSETSOLVER_HPP_
#define _LEVELSETSOLVER_HPP_

#define LSS_DEBUG_REINI 1

#include <algorithm>

#include <life/lifearray/tab.hpp>

#include <life/lifesolver/HyperbolicSolverIP.hpp>
#include <life/lifesolver/LevelSetSolverUtils.hpp>
#include <life/lifefem/subelements.hpp>

namespace LifeV {
    /*
      A type for points and vectors. Points are defined by their position
      vectors.
    */
    typedef boost::numeric::ublas::bounded_vector<Real, 3> geo_point_type;

    /*!
      \class LevelSetSolver
      \brief Level set solver class

      \author Daniele Antonio Di Pietro <dipietro@unibg.it>
    */
    template<typename MeshType>
    class LevelSetSolver
        :
        public HyperbolicSolverIP<MeshType>
    {
    public:
        /*! \name Subclasses */
        //@{
        template<int NUMNODES = 3, typename PointType = geo_point_type>
        class Face{
        public:
            /*! \name Typedefs */
            //@{
            typedef geo_point_type point_type;
            typedef geo_point_type vector_type;

            typedef boost::numeric::ublas::bounded_vector<Real, NUMNODES> face_container_type;

            typedef std::vector<point_type> point_container_type;
            typedef UInt id_type;
            //@}

            /*! \name Constructors, destructor */
            //@{
            Face() {
                _M_points.resize(NUMNODES);
                _M_has_normal = false;
            }

            Face(ID id) : _M_id(id) {
                _M_points.resize(NUMNODES);
            }
            //@}

            /*! \name Accessors */
            //@{

            /*! \return the number of local nodes */
            const int numLocalNodes() const {
                return NUMNODES;
            }

            /*! \return a const reference to the i-th point */
            const point_type& point(const id_type i) const {
                return _M_points[i - 1];
            }
            //@}

            /*! \name Mutators */
            //@{

            /*! \return the i-th point */
            point_type& point(const id_type i) {
                return _M_points[i - 1];
            }

            /*! \return the id */
            id_type& id() {
                return _M_id;
            }

            //! Set i-th point
            void setPoint(const id_type i, const point_type& P) {
                _M_points[i - 1] = P;
            }

            //@}

            /*! \name Methods */
            //@{

            /*! \return the normal versor */
            inline vector_type normal() {
                vector_type n;

                if(_M_has_normal)
                    n = _M_normal;
                else {
                    vector_type v1 = _M_points[1] - _M_points[0];
                    vector_type v2 = _M_points[2] - _M_points[0];
                    vector_type n = crossProd(v1, v2);

                    Real d = boost::numeric::ublas::norm_2(n);
                    n /= d;
                    _M_normal = n;
                    _M_has_normal = true;
                }
                return n;
            }

            //! Compute the distance between point P and the face
            inline Real pointToFaceDistance(point_type& P) {
                point_type P0 = _M_points[0];
                vector_type n = normal();
                point_type Q;

                bool projection_is_on_face;

                Real d;

                pointProjectionOnPlane(P, P0, n, Q, d);
                //projection_is_on_face = pointIsOnFace(Q);
                projection_is_on_face = false;

                if (! projection_is_on_face) {
                    d = pointToPointDistance(P, point(1));
                    for(int iP = 2; iP <= 3; iP++)
                        d = std::min<Real>( d, pointToPointDistance(P, point(iP)) );
                }

                return d;
            }

            /*!
              Determine whether the point P belongs to the face.
              It is assumed that point P belongs to the plane the face is
              a subset of.
            */
            bool pointIsOnFace(const point_type& P) {
                Real D =
                    _M_points[0][0] * (_M_points[2][1] - _M_points[1][1]) +
                    _M_points[1][0] * (_M_points[0][1] - _M_points[2][1]) +
                    _M_points[2][0] * (_M_points[1][1] - _M_points[0][1]);

                Real xi =
                    (P[0]-_M_points[0][0])*(_M_points[0][1]-_M_points[2][1]) +
                    (P[1]-_M_points[0][1])*(_M_points[2][0]-_M_points[2][0]);
                xi /= D;

                Real eta =
                    (P[0]-_M_points[0][0])*(_M_points[1][1]-_M_points[0][1]) +
                    (P[1]-_M_points[0][1])*(_M_points[0][0]-_M_points[1][0]);
                eta /= D;

                return (xi >= 0 && xi <= 1 && eta >= 0 && eta <= 1 - xi);
            }

            //! Display the face
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

        /*! \name Typedefs */
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

        /*! \name Constructors and destructors */
        //@{
        LevelSetSolver(mesh_type& mesh,
                       const GetPot& data_file,
                       const std::string& data_section,
                       const RefFE& reffe,
                       const QuadRule& qr,
                       const QuadRule& qr_bd,
                       const BCHandler& bc_h,
                       CurrentFE& fe_velocity,
                       const Dof& dof_velocity,
                       velocity_type& velocity0)
            :
            HyperbolicSolverIP<mesh_type>(mesh,
                                          data_file,
                                          data_section,
                                          reffe,
                                          qr,
                                          qr_bd,
                                          bc_h,
                                          fe_velocity,
                                          dof_velocity,
                                          velocity0),
            _M_tol(1e-10) {}
        //@}

        /*! \name Accessors */
        //@{

        /*! \return the level set function */
        lsfunction_type const & lsfunction() const {
            return this->u();
        }

        //@}

        /*! \name Mutators */
        //@{

        //! set the tolerance on interface position
        void setTol(Real tol) {
            _M_tol = tol;
        }
        //@}

        /*! \name Methods */
        //@{

        //! reinitialize the interface using direct method
        void directReinitialization() {
            build_interface();
            if ( this->verboseMode() ) {
                std::cout << "** LSS ** " << _M_face_list.size() << " faces"
                          << std::endl;
                std::cout << "** LSS ** Reinitializing the interface"
                          << std::endl;
            }
#if LSS_DEBUG_REINI
            exportToMatlab("./results/before.m");
#endif
            for(UInt iP = 1; iP <= this->mesh().numPoints(); iP++) {
                point_type P;
                convertPointType(P, this->mesh().pointList(iP));

                Real d = _M_face_list.begin()->pointToFaceDistance(P);

                for(face_list_iterator faces_it = _M_face_list.begin();
                    faces_it != _M_face_list.end(); faces_it++)
                    d = std::min<Real>(d, faces_it->pointToFaceDistance(P));


                this->u()[iP - 1] = signum(this->u()[iP - 1]) * d;
            }
#if LSS_DEBUG_REINI
            build_interface();
            exportToMatlab("./results/after.m");
#endif
        }

        //! Compute the mass of i-th fluid
        Real computeMass(fluid_type __fluid_id) {
            Real __mass = 0.;
            Real __ls_fun;
            UInt __jglo;

            for(UInt i = 1; i <= this->mesh().numVolumes(); i++) {
                this->fe().updateJac( this->mesh().volumeList(i) );

                for(int iq = 0; iq < this->fe().nbQuadPt; iq++) {
                    __ls_fun = 0.;
                    for(int j = 0; j < this->fe().nbNode; j++) {
                        __jglo = this->dof().localToGlobal(i, j + 1) - 1;
                        __ls_fun += this->fe().phi(j, iq) * this->u()(__jglo);
                    }
                    __mass += (__fluid_id == fluid1) ?
                        (__ls_fun <  0 ? 1. : 0.) * this->fe().weightDet(iq) :
                        (__ls_fun >= 0 ? 1. : 0.) * this->fe().weightDet(iq);
                }
            }
            return __mass;
        }

        /*!
          Export the interface to a Matlab script for debugging/visualization
          purposes
        */
        void exportToMatlab(std::string __file_name) {
            std::ofstream ofile(__file_name.c_str());
            ASSERT( ofile, "Error exporting interface to Matlab format: Output file cannot be open" );

            // POINTS

            ofile << "% ========================================" << std::endl;
            ofile << "% POINTS" << std::endl;
            ofile << "% ----------------------------------------" << std::endl;
            ofile << "% P(i, j) = j-th coordinate of point i" << std::endl;
            ofile << "% ========================================" << std::endl;
            ofile << "P = [..." << std::endl;
            for(face_list_iterator faces_it = _M_face_list.begin();
                faces_it != _M_face_list.end(); faces_it++) {
                for(int node_id = 1;
                    node_id <= faces_it->numLocalNodes(); node_id++)
                    ofile << faces_it->point(node_id)[0] << ", "
                          << faces_it->point(node_id)[1] << ", "
                          << faces_it->point(node_id)[2] << ";... "
                          << std::endl;
            }
            ofile << "];" << std::endl;

            // CONNECTIVITIES

            ofile << "% ========================================" << std::endl;
            ofile << "% CONNECTIVITIES" << std::endl;
            ofile << "% ----------------------------------------" << std::endl;
            ofile << "% C(i, j) = global id of j-th local node  " << std::endl;
            ofile << "% of face i" << std::endl;
            ofile << "% ========================================" << std::endl;
            ofile << " C = [..." << std::endl;
            int nP = 1;
            for(face_list_iterator faces_it = _M_face_list.begin();
                faces_it != _M_face_list.end(); faces_it++) {
                for(int node_id = 1;
                    node_id <= faces_it->numLocalNodes(); node_id++, nP++)
                    ofile << nP << " ";
                ofile << ";..." << std::endl;
            }
            ofile << "];" << std::endl;

            // SOLUTION

            ofile << "% ========================================" << std::endl;
            ofile << "% SOLUTION" << std::endl;
            ofile << "% ========================================" << std::endl;
            ofile << "U = [..." << std::endl;
            for(UInt i = 0; i < this->u().size(); i++)
                ofile << this->u()[i] << ";..." << std::endl;
            ofile << "];" << std::endl;

            ofile.close();
        }
        //@}

    private:
        //! The face list
        face_list_type _M_face_list;

        //! Tolerance on interface position
        Real _M_tol;

        //! Build the interface given the level set function
        inline void build_interface() {
            // Clear the face list

            _M_face_list.clear();

            Subelements se( this->fe() );

            // A loop on the volumes to find the local intersection with the
            // interface
            for(UInt iV = 1; iV <= this->mesh().numVolumes(); iV++) {

                this->fe().updateJac( this->mesh().volumeList( iV ) );

                for(UInt iSe = 0; iSe < se.numSubelements(); ++iSe ) {

                    point_list_type locPointList;

                    for(UInt ie = 1; ie <= se.numLocalEdges(); ie++) {

                        UInt j1 = se.edge2Point( iSe, ie, 1 );
                        UInt j2 = se.edge2Point( iSe, ie, 2 );

                        UInt jg1 = this->dof().localToGlobal(iV, j1) - 1;
                        UInt jg2 = this->dof().localToGlobal(iV, j2) - 1;

                        Real u1 = this->u()[jg1];
                        Real u2 = this->u()[jg2];

                        // If the solution is zero somewhere on the segment
                        if(u1 * u2 < 0) {
                            point_type P;
                            point_type P1;
                            point_type P2;

                            se.coord( P1[0], P1[1], P1[2], j1 );
                            se.coord( P2[0], P2[1], P2[2], j2 );

                            //this->fe().coorMap( P1[0], P1[1], P1[2], se.xi(j1),
                            //               se.eta(j1), se.zeta(j1) );
                            //this->fe().coorMap( P2[0], P2[1], P2[2], se.xi(j2),
                            //               se.eta(j2), se.zeta(j2) );

                            //P1[0] = this->mesh().pointList(jg1 + 1).x();
                            //P1[1] = this->mesh().pointList(jg1 + 1).y();
                            //P1[2] = this->mesh().pointList(jg1 + 1).z();

                            //P2[0] = this->mesh().pointList(jg2 + 1).x();
                            //P2[1] = this->mesh().pointList(jg2 + 1).y();
                            //P2[2] = this->mesh().pointList(jg2 + 1).z();

                            findZeroOnEdge(P1, u1, P2, u2, P);

                            pointsCoincide<point_type> compare(_M_tol);
                            bool isAlreadyThere = false;

                            for(point_list_type::iterator point_it =
                                    locPointList.begin();
                                point_it != locPointList.end(); point_it++)
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
                            for(point_list_type::iterator point_it =
                                    locPointList.begin();
                                point_it != locPointList.end(); point_it++) {
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
                            establish_orientation_4pt(locPointList,
                                                      orientation);

                            face1.setPoint(1, locPointList[orientation[0]]);
                            face1.setPoint(2, locPointList[orientation[1]]);
                            face1.setPoint(3, locPointList[orientation[2]]);

                            face2.setPoint(1, locPointList[orientation[2]]);
                            face2.setPoint(2, locPointList[orientation[3]]);
                            face2.setPoint(3, locPointList[orientation[0]]);

                            _M_face_list.push_back(face1);
                            _M_face_list.push_back(face2);

                            break;
                        } // case 4:
                    } // switch
                } // loop over local edges
            } // loop over subelements
        } // loop over mesh volumes

        /*!
           Given a vector of 4 nodes, it returns a point ordering suitable
           to obtain a proper quadrilateral. Orientation list is given below:<br>
           1 : [P0, P1, P2, P3]<br>
           2 : [P0, P1, P3, P2]<br>
           3 : [P0, P2, P1, P3]<br>
         */
        inline void establish_orientation_4pt(point_list_type& v,
                                              orientation_type& orientation) {
            point_type P0 = v[0];
            point_type P1 = v[1];
            point_type P2 = v[2];
            point_type P3 = v[3];

            vector_type v1 = P1 - P0;
            vector_type v2 = P2 - P0;
            vector_type v3 = P3 - P0;
            vector_type n = crossProd(v1, v2);

            Real a12 = inner_prod(crossProd(v1, v2), n);
            Real a13 = inner_prod(crossProd(v1, v3), n);
            Real a23 = inner_prod(crossProd(v2, v3), n);

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

    }; // class LevelSetSolver

} // namespace LifeV

#endif /* _LEVELSETSOLVER_HPP_ */
