/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
  \file bdf.h
  \author A. Veneziani
  \date 04/2003
  \version 1.0

  \brief File containing a class for an easy handling of different order time
  discretizations/extrapolations BDF based

*/
#ifndef _BDF_H
#define _BDF_H
#include <string>
#include <iostream>
#include "GetPot.hpp"
#include "lifeV.hpp"
#include "vecUnknown.hpp"

namespace LifeV
{
const uint BDF_MAX_ORDER = 3;

typedef Real ( *Funct ) ( const Real&, const Real&, const Real&, const Real&,
                          const ID& );

/*!
  \class Bdf
  \brief Backward differencing formula time discretization

  A differential equation of the form

  \f$ M u' = A u + f \f$

  is discretized in time as

  \f$ M p'(t_{k+1}) = A u_{k+1} + f_{k+1} \f$

  where p denotes the polynomial of order n in t that interpolates
  (t_i,u_i) for i = k-n+1,...,k+1.

  The approximative time derivative \f$ p'(t_{k+1}) \f$ is a linear
  combination of state vectors u_i:

  \f$ p'(t_{k+1}) = \frac{1}{\Delta t} (\alpha_0 u_{k+1} - \sum_{i=0}^n \alpha_i u_{k+1-i} )\f$

  Thus we have

  \f$ \frac{\alpha_0}{\Delta t} M u_{k+1} = A u_{k+1} + f + M \bar{p} \f$

  with

  \f$ \bar{p} = \frac{1}{\Delta t} \sum_{i=1}^n \alpha_i u_{k+1-i} \f$

  This class stores the n last state vectors in order to be able to
  calculate \f$ \bar{p} \f$. It also provides alpha_i
  and can extrapolate the the new state from the n last states with a
  polynomial of order n-1:

  \f$ u_{k+1} \approx \sum_{i=0}^{n-1} \beta_i u_{k-i} \f$
*/
class Bdf
{
public:
    /*! Constructor
     *  @param n order of the BDF
     */
    Bdf( const UInt n );

    //! Initialize all the entries of the unknown vector to be derived with the
    //! vector u0 (duplicated)
    void initialize_unk( Vector u0 );

    //! Initialize all the entries of the unknown vector to be derived with a
    //! set of vectors uv0
    void initialize_unk( std::vector<Vector> uv0 );

    /*! Initialize all the entries of the unknonwn vectors with a given function
        The array of initial conditions needed by the selected BDF is
        initialized as follows: _unk=[ u0(t0), u0(t0-dt), u0(t0-2*dt), ...]
        For the space dependence of the initial conditions we need informations
        on:
        -# the mesh (coordinates of points)
        -# the reference and the current FE (basis functions)
        -# is it a vector or a scalar problem ? bdf doesn't know it
        -# which is the initial time (t0) and the time step (for solutions
           before the initial instant)
        Based on the NavierStokesHandler::initialize by M. Fernandez
    */
    template <typename Mesh, typename RefFE, typename CurrFE, typename Dof>
    void initialize_unk( const Funct& u0, Mesh& mesh, RefFE& refFE, CurrFE& currFE,
                         Dof& dof, Real t0, Real dt, UInt nbComp );

    /*! Update the vectors of the previous time steps by shifting on the right
     *  the old values.
     *  @param u_curr current (new) value of the state vector
     */
    void shift_right( Vector u_curr );

    //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative
    //! formula
    Vector time_der( Real dt ) const;

    //! Returns the right hand side \f$ \bar{p} \Delta t \f$ of the time
    //! derivative formula. The timestep is taken into account elsewhere,
    //! e. g. in the mass matrix.
    Vector time_der() const;

    //! Compute the polynomial extrapolation approximation of order n-1 of
    //! u^{n+1} defined by the n stored state vectors
    Vector extrap() const;

    //! Return the i-th coefficient of the time derivative alpha_i
    double coeff_der( UInt i ) const;

    //! Return the i-th coefficient of the time extrapolation beta_i
    double coeff_ext( UInt i ) const;

    //! Return a vector with the last n state vectors
    const std::vector<Vector>& unk() const;

    void showMe() const;

    ~Bdf();

private:
    //! Order of the BDF derivative/extrapolation: the time-derivative
    //! coefficients vector has size n+1, the extrapolation vector has size n
    UInt _n;

    //! Size of the unknown vector
    UInt _s;

    //! Coefficients \f$ \alpha_i \f$ of the time bdf discretization
    Vector _alpha;

    //! Coefficients \f$ \beta_i \f$ of the extrapolation
    Vector _beta;

    //! Last n state vectors
    std::vector<Vector> _unk;
};


///
// template implementations
//

template <typename Mesh, typename RefFE, typename CurrFE, typename Dof>
void Bdf::initialize_unk( const Funct& u0, Mesh& mesh, RefFE& refFE,
                          CurrFE& currFE, Dof& dof, Real t0, Real dt,
                          UInt nbComp = 1 )
{
    typedef typename Mesh::VolumeShape GeoShape; // Element shape

    UInt nDofpV = refFE.nbDofPerVertex; // number of Dof per vertex
    UInt nDofpE = refFE.nbDofPerEdge;   // number of Dof per edge
    UInt nDofpF = refFE.nbDofPerFace;   // number of Dof per face
    UInt nDofpEl = refFE.nbDofPerVolume; // number of Dof per Volume

    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF = GeoShape::numFaces;    // Number of element's faces

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element
    UInt nDofElemF = nElemF * nDofpF; // number of face's Dof on a Element

    std::vector< Vector >::iterator iter = _unk.begin();
    std::vector< Vector >::iterator iter_end = _unk.end();


    UInt size_comp = dof.numTotalDof();
    _s = size_comp * nbComp; // Inizialization of the dimension of the vector

    Vector aux( _s );
    aux = 0.0;


    for ( iter = _unk.begin() ; iter != iter_end; iter++ )
    {
        *iter = aux;
    }

    Real x, y, z;
    UInt backtime;
    ID lDof;

    // Loop on elements of the mesh
    for ( ID iElem = 1; iElem <= mesh.numVolumes(); ++iElem )
    {

        currFE.updateJac( mesh.volume( iElem ) );

        // Vertex based Dof
        if ( nDofpV )
        {

            // loop on element vertices
            for ( ID iVe = 1; iVe <= nElemV; ++iVe )
            {

                // Loop number of Dof per vertex
                for ( ID l = 1; l <= nDofpV; ++l )
                {
                    lDof = ( iVe - 1 ) * nDofpV + l; // Local dof in this element

                    // Nodal coordinates
                    currFE.coorMap( x, y, z, currFE.refFE.xi( lDof - 1 ),
                                    currFE.refFE.eta( lDof - 1 ),
                                    currFE.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components ****************
                    backtime = 0;
                    for ( std::vector<Vector>::iterator it = _unk.begin(); it<_unk.end();
                            it++ )
                    {
                        for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                        {
                            // std::cout << "comp " << icmp*size_comp + dof.localToGlobal(iElem,lDof) - 1 << std::endl;
                            ( *it )
                            ( icmp * size_comp + dof.localToGlobal( iElem, lDof ) - 1 ) =
                                u0( t0 - backtime * dt, x, y, z, icmp + 1 );
                            backtime++;
                        }
                    }
                }
            }
        }

        // Edge based Dof
        if ( nDofpE )
        {

            // loop on element edges
            for ( ID iEd = 1; iEd <= nElemE; ++iEd )
            {

                // Loop number of Dof per edge
                for ( ID l = 1; l <= nDofpE; ++l )
                {
                    // Local dof in the adjacent Element
                    lDof = nDofElemV + ( iEd - 1 ) * nDofpE + l;

                    // Nodal coordinates
                    currFE.coorMap( x, y, z, currFE.refFE.xi( lDof - 1 ),
                                    currFE.refFE.eta( lDof - 1 ),
                                    currFE.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    backtime = 0;
                    for ( std::vector<Vector>::iterator it = _unk.begin(); it<_unk.end();
                            it++ )
                    {
                        for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                        {
                            ( *it )
                            ( icmp * size_comp + dof.localToGlobal( iElem, lDof ) - 1 ) =
                                u0( t0 - backtime * dt, x, y, z, icmp + 1 );
                            backtime++;
                        }
                    }
                }
            }
        }

        // Face based Dof
        if ( nDofpF )
        {

            // loop on element faces
            for ( ID iFa = 1; iFa <= nElemF; ++iFa )
            {

                // Loop on number of Dof per face
                for ( ID l = 1; l <= nDofpF; ++l )
                {

                    lDof = nDofElemE + nDofElemV + ( iFa - 1 ) * nDofpF + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    currFE.coorMap( x, y, z, currFE.refFE.xi( lDof - 1 ), currFE.refFE.eta( lDof - 1 ), currFE.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    backtime = 0;
                    for ( std::vector<Vector>::iterator it = _unk.begin();it<_unk.end();it++ )
                    {
                        for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                        {
                            ( *it )
                            ( icmp * size_comp + dof.localToGlobal( iElem, lDof ) - 1 ) = u0( t0 - backtime * dt, x, y, z, icmp + 1 );
                            backtime++;
                        }
                    }
                }
            }
        }

        // Element based Dof
        // Loop on number of Dof per Element
        for ( ID l = 1; l <= nDofpEl; ++l )
        {
            lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

            // Nodal coordinates
            currFE.coorMap( x, y, z, currFE.refFE.xi( lDof - 1 ), currFE.refFE.eta( lDof - 1 ), currFE.refFE.zeta( lDof - 1 ) );

            // Loop on data vector components
            backtime = 0;
            for ( std::vector<Vector>::iterator it = _unk.begin();it<_unk.end();it++ )
            {
                for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                {
                    ( *it )
                    ( icmp * size_comp + dof.localToGlobal( iElem, lDof ) - 1 ) = u0( t0 - backtime * dt, x, y, z, icmp + 1 );
                    backtime++;
                }
            }
        }
    }
    return ;
}
}
#endif
