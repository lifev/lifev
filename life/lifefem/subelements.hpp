/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
      Date: 2005-02-14

 Copyright (C) 2004 EPFL

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
   \file subelement.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2005-02-14
*/
#ifndef SUBELEMENT_HPP
#define SUBELEMENT_HPP 1

#include <vector>
#include <basisElSh.hpp> // for eToP
#include <currentFE.hpp>

namespace LifeV
{

/*!
  \class Subelements
  \brief P1 subdivision of tetrahedra for level set interface reconstruction.
  
  P1 is subdivided into 1 P1 tetrahedron, P1-bubble is subdivided into 4
  tetrahedra (cross-grid), P2 is subdivided into 8 tetrahedra.
  \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
*/
class Subelements
{
public:
    /*!
      Constructor.
      \param fe the macro element to subdivide
    */
    Subelements( const CurrentFE& fe )
        : M_fe( fe )
        {
            M_xi.resize(fe.nbNode);
            M_eta.resize(fe.nbNode);
            M_zeta.resize(fe.nbNode);
            M_xi[0] = 0.; M_eta[0] = 0.; M_zeta[0] = 0.;
            M_xi[1] = 1.; M_eta[1] = 0.; M_zeta[1] = 0.;
            M_xi[2] = 0.; M_eta[2] = 1.; M_zeta[2] = 0.;
            M_xi[3] = 0.; M_eta[3] = 0.; M_zeta[3] = 1.;
            switch ( fe.nbNode )
            {
                case 4:
                    {
                        M_numSubelements = 1;
                        M_numLocalEdges = 6;
                        M_nodes.resize( 4 * M_numSubelements );
                        M_nodes[0] = 1; M_nodes[1] = 2;
                        M_nodes[2] = 3; M_nodes[3] = 4;
                        break;
                    }
                case 5:
                    {
                        M_numSubelements = 4;
                        M_numLocalEdges = 6;
                        M_xi[4] = 1./3.; M_eta[4] = 1./3.; M_zeta[4] = 1./3.;
                        M_nodes.resize( 4 * M_numSubelements );
                        M_nodes[ 0] = 1; M_nodes[ 1] = 2;
                        M_nodes[ 2] = 3; M_nodes[ 3] = 5;
                        M_nodes[ 4] = 1; M_nodes[ 5] = 2;
                        M_nodes[ 6] = 4; M_nodes[ 7] = 5;
                        M_nodes[ 8] = 2; M_nodes[ 9] = 3;
                        M_nodes[10] = 4; M_nodes[11] = 5;
                        M_nodes[12] = 3; M_nodes[13] = 1;
                        M_nodes[14] = 4; M_nodes[15] = 5;
                        break;
                    }
                case 10:
                    {
                        M_numSubelements = 8;
                        M_numLocalEdges = 6;
                        M_xi[4] = 0.5; M_eta[4] = 0. ; M_zeta[4] = 0.;
                        M_xi[5] = 0.5; M_eta[5] = 0.5; M_zeta[5] = 0.;
                        M_xi[6] = 0. ; M_eta[6] = 0.5; M_zeta[6] = 0.;
                        M_xi[7] = 0. ; M_eta[7] = 0. ; M_zeta[7] = 0.5;
                        M_xi[8] = 0.5; M_eta[8] = 0. ; M_zeta[8] = 0.5;
                        M_xi[9] = 0. ; M_eta[9] = 0.5; M_zeta[9] = 0.5;
                        M_nodes.resize( 4 * M_numSubelements );
                        M_nodes[ 0] =  1; M_nodes[ 1] =  5;
                        M_nodes[ 2] =  7; M_nodes[ 3] =  8;
                        M_nodes[ 4] =  5; M_nodes[ 5] =  2;
                        M_nodes[ 6] =  6; M_nodes[ 7] =  9;
                        M_nodes[ 8] =  7; M_nodes[ 9] =  6;
                        M_nodes[10] =  3; M_nodes[11] = 10;
                        M_nodes[12] =  5; M_nodes[13] =  6;
                        M_nodes[14] =  7; M_nodes[15] =  9;
                        M_nodes[16] =  5; M_nodes[17] =  9;
                        M_nodes[18] =  7; M_nodes[19] =  8;
                        M_nodes[20] =  6; M_nodes[21] =  7;
                        M_nodes[22] =  9; M_nodes[23] = 10;
                        M_nodes[24] =  7; M_nodes[25] =  8;
                        M_nodes[26] =  9; M_nodes[27] = 10;
                        M_nodes[28] =  8; M_nodes[29] =  9;
                        M_nodes[30] = 10; M_nodes[31] =  4;
                        break;
                    }
                default:
                    {
                        ERROR_MSG( "Subelements implemented only for tetrahedra" );
                        break;
                    }
            }
        }
    
    /*! \return number of subelements */
    UInt numSubelements() const { return M_numSubelements; }
    
    /*! \return number of local edges in each subelement */
    UInt numLocalEdges() const { return M_numLocalEdges; }
    
    /*!
      \param iSubelement number of subelement, 0 ... numSubelements()-1
      \param iEdge number of edge in subelement, 1 or 2
      \param iPoint number of point in edge, 1 ... numLocalEdges()
      \return point index in macro element, starting from 1
    */
    UInt edge2Point( UInt iSubelement, UInt iEdge, UInt iPoint )
        {
            UInt subPoint = LinearTetra::eToP( iEdge, iPoint ) - 1;
            return M_nodes[ 4 * iSubelement + subPoint ];
        }

    /*! real world coordinates of macro-element node
      \param x x coordinate (output)
      \param y y coordinate (output)
      \param z z coordinate (output)
      \param iLocPoint node index in macro element, starting from 1
    */
    void coord( Real& x, Real& y, Real& z, UInt iLocPoint )
        {
            M_fe.coorMap( x, y, z, M_xi[iLocPoint],
                          M_eta[iLocPoint], M_zeta[iLocPoint] );
        }
private:
    //! current finite element
    const CurrentFE& M_fe;
    
    //! number of subelements
    UInt M_numSubelements;
    
    //! number of edges per subelement
    UInt M_numLocalEdges;
    
    //! reference coordinates xi of nodes in macroelement
    std::vector<Real> M_xi;
    
    //! reference coordinates eta of nodes in macroelement
    std::vector<Real> M_eta;
    
    //! reference coordinates zeta of nodes in macroelement
    std::vector<Real> M_zeta;
    
    //! 4 consecutive node numbers in this vector define one subelement
    std::vector<UInt> M_nodes;
    
}; // class Subelements

} // namespace LifeV

#endif /* SUBELEMENT_HPP */
