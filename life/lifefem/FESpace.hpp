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
  \file FESpace.hpp
  \author Gilles Fourestey
  \date 06/2007
  \version 1.0

  \brief This file contains an abstract class for NavierStokes solvers.

*/

#ifndef _FESPACE_H_
#define _FESPACE_H_

#include <sstream>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>


#include <life/lifecore/life.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifefem/geoMap.hpp>
//#include <life/lifefem/bcHandler.hpp>
#include <life/lifearray/pattern.hpp>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <ext/slist>
#include <life/lifearray/SimpleVect.hpp>
#include <utility>
#include <life/lifemesh/partitionMesh.hpp>


using std::pair;

namespace LifeV
{

/*!
  \class FESpace

  Abstract class which defines the general structure of a NavierStokes solver.
  For each new NavierStokes solver  we have to implement the corresponding
  timeAdvance and an iterate methods

*/


inline Real
elemL22( boost::function<double( double, double, double, double, UInt )> fct,
         const CurrentFE& fe, const Real t, const UInt nbcomp )
{
    int ig;
    UInt ic;
    Real s = 0., f, x, y, z;
    for ( ig = 0;ig < fe.nbQuadPt;ig++ )
    {
        fe.coorQuadPt( x, y, z, ig );
        for ( ic = 0; ic < nbcomp; ic++ )
        {
            f = fct( t, x, y, z, ic + 1 );
            s += f * f * fe.weightDet( ig );
        }
    }
    return s;
}



template <typename Mesh, typename Map>
class FESpace
{

public:

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&,
                                 const Real&, const ID& );
    typedef Mesh mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef Map  map_type;

//    typedef typename NSStabilization NSStabilization;
    //! Constructors
    /*!
      \param data_file GetPot data file
      \param refFE_u reference FE for the velocity
      \param refFE_p reference FE for the pressure
      \param Qr_u volumic quadrature rule for the velocity
      \param bdQr_u surface quadrature rule for the velocity
      \param Qr_p volumic quadrature rule for the pressure
      \param bdQr_p surface quadrature rule for the pressure
      \param BCh_fluid boundary conditions for the fluid
      \param ord_bdf order of the bdf time advancing scheme and incremental pressure approach (default: Backward Euler)
    */

    FESpace(  partitionMesh<Mesh>&  mesh,
              const RefFE&     refFE,
              const QuadRule&  Qr,
              const QuadRule&  bdQr,
              const int        fDim,
              Epetra_Comm&     comm
              );

    FESpace(  mesh_ptrtype            mesh,
              const RefFE&     refFE,
              const QuadRule&  Qr,
              const QuadRule&  bdQr,
              const int        fDim,
              Epetra_Comm&     comm
              );

    void initMap( const int    fDim,
                  Epetra_Comm& comm );

//    const mesh_type& mesh()       {return M_mesh;}
//    const mesh_type& mesh() const {return *M_mesh;}
//    mesh_type& mesh() {return *M_mesh;}

    const mesh_ptrtype mesh() const {return M_mesh;}
    mesh_ptrtype mesh() {return M_mesh;}

//    NSStabilization stabilization(){return M_dataType.stabilization();}



    //! Returns the velocity dof
    const Dof& dof() const;
    Dof& dof();

    //! returns the res FE

//    RefFE& refFE()             {return M_refFE;}
    const RefFE& refFE() const {return M_refFE;}

    //! returns the current FE
    CurrentFE&   fe()   {return M_fe;}
    const CurrentFE&   fe() const  {return M_fe;}

    //! returns the current boundary FE
    CurrentBdFE& feBd() {return M_feBd;}

    //! returns the volumic quadratic rule

    const QuadRule& qr() {return M_Qr;}

    //! returns the surfasic quadratic rule

    const QuadRule& bdQr() {return M_bdQr;}

    //! FE space dimension

    const UInt dim()      const {return M_dim;}
    const UInt fieldDim() const {return M_fieldDim;}
    //! map getter

    const map_type& map() const {return M_map;}
    map_type& map() {return M_map;}

    //! Interpolate a given velocity function nodally onto a velocity vector
    template<typename vector_type>
    void interpolate( const Function& fct, vector_type& vect, const Real time = 0. );

    //! calculate L2 pressure error for given exact pressure function
    //! takes into account a possible offset by a constant
    //! \param pexact the exact pressure as a function
    //! \param time the time
    //! \param relError Real* to store the relative error in
//     Real pErrorL2( const Function&            pexact,
//                    const ScalUnknown<Vector>& p,
//                    Real                       time,
//                    Real*                      relError = 0 );

    //! calculate L2 velocity error for given exact velocity function
    //! \param pexact the exact velocity as a function
    //! \param time the time
    //! \param relError Real* to store the relative error in

    template<typename vector_type>
    void L2ScalarProduct( const Function& fct,
                       vector_type& vec,
                       const double t
                       );

        template<typename vector_type>
    double L20Error( const Function& fexact,
                     const vector_type& vec,
                     const Real time,
                     Real* relError = 0 );

    template<typename vector_type>
    double L2Error( const Function&                fexact,
                    const vector_type&             vec,
                    const double                   time,
                    double*                        relError=0 );


    template<typename vector_type>
    double L2Norm( vector_type& vec );


    BasePattern::PatternType patternType();


    //! Do nothing destructor
    virtual ~FESpace()
    {}

private:

    //! copy constructor

    FESpace( const FESpace& fespace );


    //! reference to the mesh
    mesh_ptrtype           M_mesh;

    //! Reference FE for the velocity
    const RefFE&    M_refFE;

    //! The Dof object
    Dof             M_dof;

    //! The number of total dofs
    UInt            M_dim;

    //! dimention of the field variable ( scalar/vector field)
    UInt            M_fieldDim;

    //! Quadrature rule for volumic elementary computations
    const QuadRule& M_Qr;

    //! Quadrature rule for surface elementary computations
    const QuadRule& M_bdQr;


    //! Current FE
    CurrentFE       M_fe;
    CurrentBdFE     M_feBd;

    //! Map

    map_type        M_map;
//     //! MPI communicator
//     Epetra_Comm*    M_comm;




//     BasePattern::PatternType M_patternType;


};



//
// IMPLEMENTATION
//


// Constructors
template <typename Mesh, typename Map>
FESpace<Mesh, Map>::
FESpace( partitionMesh<Mesh>& mesh,
         const RefFE&         refFE,
         const QuadRule&      Qr,
         const QuadRule&      bdQr,
         const int            fDim,
         Epetra_Comm&         comm
         ):
    M_mesh                             ( mesh.mesh() ),
    M_refFE                            ( refFE ),
    M_dof                              ( *M_mesh, refFE ),
    M_dim                              ( M_dof.numTotalDof() ),
    M_fieldDim                         ( fDim ),
    M_Qr                               ( Qr ),
    M_bdQr                             ( bdQr ),
    M_fe                               ( M_refFE,
                                         getGeoMap( *M_mesh ),
                                         M_Qr ),
    M_feBd                             ( refFE.boundaryFE(),
                                         getGeoMap( *M_mesh ).boundaryMap(),
                                         M_bdQr ),
    M_map                              ()
{
    Map map( M_refFE, mesh, comm );
    for (int ii = 0; ii < fDim; ++ii)
        M_map += map;
}


template <typename Mesh, typename Map>
FESpace<Mesh, Map>::
FESpace( mesh_ptrtype         mesh,
         const RefFE&         refFE,
         const QuadRule&      Qr,
         const QuadRule&      bdQr,
         const int            fDim,
         Epetra_Comm&         comm
         ):
    M_mesh                             ( mesh ),
    M_refFE                            ( refFE ),
    M_dof                              ( *M_mesh, refFE ),
    M_dim                              ( M_dof.numTotalDof() ),
    M_fieldDim                         ( fDim ),
    M_Qr                               ( Qr ),
    M_bdQr                             ( bdQr ),
    M_fe                               ( M_refFE,
                                         getGeoMap( *M_mesh ),
                                         M_Qr ),
    M_feBd                             ( refFE.boundaryFE(),
                                         getGeoMap( *M_mesh ).boundaryMap(),
                                         M_bdQr ),
    M_map                              ()

{
    Map map( M_refFE, *M_mesh, comm );
    for (int ii = 0; ii < fDim; ++ii)
        M_map += map;
}


// Returns the velocity Dof
template <typename Mesh, typename Map>
const Dof&
FESpace<Mesh, Map>::dof() const
{
    return M_dof;
}

template <typename Mesh, typename Map>
Dof&
FESpace<Mesh, Map>::dof()
{
    return M_dof;
}






template <typename Mesh, typename Map>
template<typename vector_type>
void
FESpace<Mesh, Map>::interpolate( const Function& fct,
                                 vector_type&    vect,
                                 const Real      time)
{
    typedef typename mesh_type::VolumeShape GeoShape ; // Element shape

    UInt nDofpV    = refFE().nbDofPerVertex; // number of Dof per vertex
    UInt nDofpE    = refFE().nbDofPerEdge;   // number of Dof per edge
    UInt nDofpF    = refFE().nbDofPerFace;   // number of Dof per face
    UInt nDofpEl   = refFE().nbDofPerVolume; // number of Dof per Volume

    UInt nElemV    = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE    = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF    = GeoShape::numFaces;    // Number of element's faces

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element
    UInt nDofElemF = nElemF * nDofpF; // number of face's Dof on a Element

    ID nbComp = M_fieldDim; // Number of components of the mesh velocity

    Real x, y, z;

    ID lDof;

    // Loop on elements of the mesh
    for ( ID iElem = 1; iElem <= mesh()->numVolumes(); ++iElem )
    {

        fe().updateJac( mesh()->volume( iElem ) );
        int elemId = mesh()->volume( iElem ).localId();

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
                    fe().coorMap( x, y, z, refFE().xi( lDof - 1 ), refFE().eta( lDof - 1 ), refFE().zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        int iDof = icmp * dim() + dof().localToGlobal( elemId, lDof  );
                        if (vect.Map().LID( iDof ) >= 0)
                            vect( iDof ) = fct( time, x, y, z, icmp + 1 );
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
                    lDof = nDofElemV + ( iEd - 1 ) * nDofpE + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    fe().coorMap( x, y, z, refFE().xi( lDof - 1 ), refFE().eta( lDof - 1 ), refFE().zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        int iDof = icmp * dim() + dof().localToGlobal( elemId, lDof  );
                        if (vect.Map().LID( iDof ) >= 0)
                            vect( iDof ) = fct( time, x, y, z, icmp + 1 );
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
                    fe().coorMap( x, y, z, refFE().xi( lDof - 1 ), refFE().eta( lDof - 1 ), refFE().zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        vect( icmp * dim() + dof().localToGlobal( elemId, lDof ) - 0 ) = fct( time, x, y, z, icmp + 1 );
                    }
                }
            }
        }
//         // Element based Dof
//         // Loop on number of Dof per Element
//         for ( ID l = 1; l <= nDofpEl; ++l )
//         {
//             lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

//             // Nodal coordinates
//             fe().coorMap( x, y, z, refFE().xi( lDof - 1 ), refFE().eta( lDof - 1 ), refFE().zeta( lDof - 1 ) );

//             // Loop on data vector components
//             for ( UInt icmp = 0; icmp < nbComp; ++icmp )
//             {
// //                vect( icmp * dim() + dof().localToGlobal( elemId, lDof ) - 0 ) = fct( time, x, y, z, icmp + 1 );
//             }
//         }
    }
}






template <typename Mesh, typename Map>
template<typename vector_type>
void
FESpace<Mesh, Map>::L2ScalarProduct( const Function& fct, vector_type& vec, const double t)
{



    for ( UInt iVol = 1; iVol <= this->mesh()->numVolumes(); iVol++ )
    {
        this->fe().updateFirstDeriv( this->mesh()->volumeList( iVol ) );

        Real s = 0., f, x, y, z;

        int i, inod, ig;
        UInt eleID = this->fe().currentLocalId();
        UInt ic;
        Real u_ig;

        for ( ic = 0; ic < M_fieldDim; ic++ )
        {
            for ( ig = 0; ig < this->fe().nbQuadPt; ig++ )
            {
                this->fe().coorQuadPt( x, y, z, ig );
                f = fct( t, x, y, z, ic + 1 );
                u_ig = 0.;
                for ( i = 0;i < this->fe().nbNode; i++ )
                {
                    inod = this->dof().localToGlobal( eleID, i + 1 ) + ic * dim();
                    u_ig = f*this->fe().phi( i, ig );
                    vec.sumIntoGlobalValues(inod,u_ig*this->fe().weightDet( ig ));
                }

            }
        }
    }

    vec.GlobalAssemble();
}



template <typename Mesh, typename Map>
template<typename vector_type>
double
FESpace<Mesh, Map>::L20Error( const Function& fexact,
                              const vector_type& vec,
                              const Real time,
                              Real* relError )
{
    Real sum2      = 0.;
    Real sum1      = 0.;
    Real sum0      = 0.;
    Real sumExact2 = 0.;
    Real sumExact1 = 0.;

    for ( UInt iVol = 1; iVol <= this->mesh()->numVolumes(); iVol++ )
    {
        this->fe().updateFirstDeriv( this->mesh()->volumeList( iVol ) );
        sum2 += elem_L2_diff_2( vec, fexact, this->fe(), this->dof(), time, M_fieldDim );
        sum1 += elem_integral_diff( vec, fexact, this->fe(), this->dof(), time, M_fieldDim );
        sum0 += this->fe().measure();
        if (relError)
        {
            sumExact2 += elem_f_L2_2( fexact, this->fe(), time, 1 );
            sumExact1 += elem_integral( fexact, this->fe(), time, 1 );
        }
    }


    double sendbuff[5] = {sum0, sum1, sum2, sumExact1, sumExact2};
    double recvbuff[5];

    this->map().Comm().SumAll(&sendbuff[0], &recvbuff[0], 5);


    int me = this->map().Comm().MyPID();

//     if (me == 0)
//     {
    sum0      = recvbuff[0];
    sum1      = recvbuff[1];
    sum2      = recvbuff[2];
    sumExact1 = recvbuff[3];
    sumExact2 = recvbuff[4];

    Real absError = sqrt( sum2 - sum1*sum1/sum0 );

    if (relError)
    {
        Real normExact = sqrt( sumExact2 - sumExact1*sumExact1/sum0 );
        *relError = absError / normExact;
    }
    return absError;
}







template <typename Mesh, typename Map>
template<typename vector_type>
double
FESpace<Mesh, Map>::L2Error( const Function&    fexact,
                             const vector_type& vec,
                             const double       time,
                             double*            relError )
{
    Real normU       = 0.;
    Real sumExact    = 0.;

    CurrentFE p1(feTetraP1, getGeoMap( *M_mesh ), M_Qr );

    for ( UInt iVol  = 1; iVol <= this->mesh()->numVolumes(); iVol++ )
    {
        this->fe().updateFirstDeriv( this->mesh()->volumeList( iVol ) );
        //p1.updateFirstDeriv( this->mesh()->volumeList( iVol ) );

        normU += elem_L2_diff_2( vec, fexact,
                                 this->fe(),
                                 //p1,
                                 this->dof(),
                                 time,
                                 M_fieldDim );
        if (relError)
        {
            sumExact += elemL22( fexact,
                                 this->fe(),
                                 //p1,
                                 time,
                                 M_fieldDim );
        }
    }


    double sendbuff[2] = {normU, sumExact};
    double recvbuff[2];

    this->map().Comm().SumAll(&sendbuff[0], &recvbuff[0], 2);


    int me = this->map().Comm().MyPID();

//     if (me == 0)
//     {
    normU    = recvbuff[0];
    sumExact = recvbuff[1];

    if (relError)
    {
        *relError = sqrt( normU / sumExact );
    }
//     }

    return sqrt( normU );
}




// template <typename Mesh, typename Map>
// template<typename vector_type>
// double
// FESpace<Mesh, Map>::L2Norm(vector_type vec)
// {

//     double norm = 0.;
//     double tmp  = 0.;

//    for ( UInt ielem = 1; ielem <= this->mesh()->numVolumes(); ielem++ )
//     {
//         //UInt elem = M_FESpace.mesh()->volumeList( ielem ).id();
//         this->fe().updateFirstDerivQuadPt( M_FESpace.mesh()->volumeList( ielem ) );

//         //int    marker = M_FESpace.mesh()->volumeList( ielem ).marker();
//         double s      = 0;
//         double volume = M_FESpace.fe().detJac(0);

//         for ( int ig = 0; ig < M_FESpace.fe().nbQuadPt; ++ig )
//         {
//             for ( int k = 0; k < M_FESpace.fe().nbNode; ++k )
//             {
//                 for ( int j = 0; j < M_fDim; ++j)
//                 {
//                     int i    = M_FESpace.fe().patternFirst(k);
//                     int idof = M_FESpace.dof().localToGlobal(M_FESpace.fe().currentLocalId(), i + 1);

//                     tmp = this->fe().phi( k, 0 , ig )*
//                         vec[idof + j*this->dim()];

//                     norm += this->fe().weigh(ig)*tmp*tmp;
//                 }
//             }
//         }


//         return norm;

// }





}
#endif
