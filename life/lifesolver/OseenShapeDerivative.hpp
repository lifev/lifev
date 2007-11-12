/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Miguel A. Fernandez <miguel.fernandez@inria.fr>
            Christoph Winkelmann <christoph.winkelmann@epfl.ch>
      Date: 2003-06-09

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
  \file Oseen.hpp
  \author G. Fourestey
  \date 2/2007

  \brief This file contains a Oseen equation solver class with shape derivative
  for FSI
*/



#ifndef _OSEENSD_H_
#define _OSEENSD_H_

#include <life/lifesolver/Oseen.hpp>


namespace LifeV
{
template< typename Mesh,
          typename SolverType = LifeV::Epetra::SolverTrilinos >
class OseenShapeDerivative:
        public Oseen<Mesh, SolverType>
{

public:

    typedef Oseen<Mesh, SolverType>             super;
    typedef typename super::vector_type         vector_type;
    typedef typename super::matrix_type         matrix_type;
    typedef typename super::matrix_ptrtype      matrix_ptrtype;
    typedef typename super::data_type           data_type;

    typedef typename super::bchandler_raw_type bchandler_raw_type;


    OseenShapeDerivative( const data_type&          dataType,
                          FESpace<Mesh, EpetraMap>& uFESpace,
                          FESpace<Mesh, EpetraMap>& pFESpace,
                          BCHandler&                bcHandler,
                          Epetra_Comm&              comm );


    OseenShapeDerivative( const data_type&          dataType,
                          FESpace<Mesh, EpetraMap>& uFESpace,
                          FESpace<Mesh, EpetraMap>& pFESpace,
                          Epetra_Comm&              comm );


    void iterateLin( bchandler_raw_type& bch );
    void updateLinearSystem( double       alpha,
                             vector_type& betaVec,
                             vector_type& w,
                             vector_type& dw,
                             vector_type& disp,
                             vector_type& sourceVec );

private:


    vector_type               M_rhsLinNoBC;
    vector_type               M_rhsLinFull;


    ElemVec                   M_elvec_du; // Elementary right hand side for the linearized velocity
    ElemVec                   M_elvec_dp; // Elementary right hand side for the linearized pressure
    ElemVec                   M_w_loc;    // Elementary mesh velocity
    ElemVec                   M_uk_loc;   // Elementary velocity
    ElemVec                   M_pk_loc;   // Elementary pressure
    ElemVec                   M_elvec;    // Elementary convection velocity
    ElemVec                   M_d_loc;    // Elementary displacement for right hand side
    ElemVec                   M_dw_loc;   // Elementary mesh velocity for right hand side



};


template<typename Mesh, typename SolverType>
OseenShapeDerivative<Mesh, SolverType>::
OseenShapeDerivative( const data_type&          dataType,
                      FESpace<Mesh, EpetraMap>& uFESpace,
                      FESpace<Mesh, EpetraMap>& pFESpace,
                      BCHandler&                BCh_u,
                      Epetra_Comm&              comm ):
    super            (dataType,
                      uFESpace,
                      pFESpace,
                      BCh_u,
                      comm),
    M_rhsLinNoBC     ( this->getMap()),
    M_rhsLinFull     ( this->getMap()),
    M_elvec_du       ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_elvec_dp       ( this->M_pFESpace.fe().nbNode, nDimensions ),
    M_w_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_uk_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_pk_loc         ( this->M_pFESpace.fe().nbNode, 1 ),
    M_elvec          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_d_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_dw_loc         ( this->M_uFESpace.fe().nbNode, nDimensions )
{

}

template<typename Mesh, typename SolverType>
OseenShapeDerivative<Mesh, SolverType>::
OseenShapeDerivative( const data_type&          dataType,
                      FESpace<Mesh, EpetraMap>& uFESpace,
                      FESpace<Mesh, EpetraMap>& pFESpace,
                      Epetra_Comm&              comm ):
    super            (dataType,
                      uFESpace,
                      pFESpace,
                      comm),
    M_rhsLinNoBC     ( this->getMap()),
    M_rhsLinFull     ( this->getMap()),
    M_elvec_du       ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_elvec_dp       ( this->M_pFESpace.fe().nbNode, nDimensions ),
    M_w_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_uk_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_pk_loc         ( this->M_pFESpace.fe().nbNode, 1 ),
    M_elvec          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_d_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_dw_loc         ( this->M_uFESpace.fe().nbNode, nDimensions )
{

}



template<typename Mesh, typename SolverType>
void OseenShapeDerivative<Mesh, SolverType>::iterateLin( bchandler_raw_type& bch )
{
    Chrono chrono;

    // matrix and vector assembling communication

    if (this->M_verbose)
        {
            std::cout << "  f-  Finalizing the matrix and vectors ...    ";
        }

    chrono.start();


    this->M_matrNoBC->GlobalAssemble();
    M_rhsLinNoBC.GlobalAssemble();

    matrix_ptrtype matrFull( new matrix_type(*this->M_matrNoBC) );
    vector_type    rhsFull = M_rhsLinNoBC;

    chrono.stop();

    if (this->M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    // boundary conditions update
    this->M_comm->Barrier();

    if (this->M_verbose) std::cout << "  f-  Applying boundary conditions...          "
              << std::flush;

    chrono.start();
    applyBoundaryConditions( *matrFull, rhsFull, bch );

    chrono.stop();

    this->M_comm->Barrier();

    if (this->M_verbose) std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;

    // solving the system

    solveSystem( matrFull, rhsFull );

    this->M_residual  = M_rhsLinNoBC;
    this->M_residual -= *this->M_matrNoBC*this->M_sol;
}



template<typename Mesh, typename SolverType>
void OseenShapeDerivative<Mesh, SolverType>::updateLinearSystem( double       alpha,
                                                                 vector_type& betaVec,
                                                                 vector_type& w,
                                                                 vector_type& dw,
                                                                 vector_type& disp,
                                                                 vector_type& sourceVec )
{
    std::cout << "  f-  LINEARIZED FLUID SYSTEM\n";

    Chrono chrono;

    if (this->M_verbose)
            std::cout << "  f-  Updating mass term on right hand side... "
                  << std::flush;

    chrono.start();

    UInt velTotalDof   = this->M_uFESpace.dof().numTotalDof();
//    UInt pressTotalDof = M_pFESpace.dof().numTotalDof();

    // Right hand side for the velocity at time

    M_rhsLinNoBC = sourceVec;

    chrono.stop();

    if (this->M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"  << std::flush;


    this->M_updated = false;

//

    if (this->M_recomputeMatrix)
        this->buildSystem();

    if (this->M_verbose)
          std::cout << "  f-  Copying the matrices ...                 "
                    << std::flush;

    chrono.start();

    this->M_matrFull.reset( new matrix_type(this->M_localMap, this->M_matrStokes->getMeanNumEntries() ));

    *(this->M_matrFull) += *(this->M_matrStokes);

    if (alpha != 0. )
    {
        *(this->M_matrFull) += *(this->M_matrMass)*alpha;
    }


    chrono.stop();
    if (this->M_verbose) std::cout << "done in " << chrono.diff() << " s.\n"
                             << std::flush;


    //! managing the convective term

    double normInf;
    betaVec.NormInf(&normInf);





    if (alpha != 0.)
     {
         this->M_matrFull->GlobalAssemble();
     }


    int nbCompU = nDimensions;

    std::cout << "    F-  Updating right hand side... ";

    //
    // RIGHT HAND SIDE FOR THE LINEARIZED ALE SYSTEM
    //
    chrono.start();

    //initialize right hand side
 //    _f_duWithOutBC = ZeroVector( _f_duWithOutBC.size() );
//     _f_p = ZeroVector( _f_p.size() );

    // Loop on elements

    vector_type betaVecRep( betaVec, *this->M_localMap.getRepeatedEpetra_Map() );
    vector_type wRep      ( w,       *this->M_localMap.getRepeatedEpetra_Map() );
    vector_type dwRep     ( dw,      *this->M_localMap.getRepeatedEpetra_Map() );
    vector_type dispRep   ( disp,    *this->M_localMap.getRepeatedEpetra_Map() );


    for ( UInt i = 1; i <= this->mesh().numVolumes(); i++ )
    {

        this->M_pFESpace.update( this->mesh().volumeList( i ) );
        this->M_uFESpace.updateFirstDerivQuadPt( this->mesh().volumeList( i ) );

        // initialization of elementary vectors
        M_elvec_du.zero();
        M_elvec_dp.zero();

        for ( UInt k = 0 ; k < ( UInt ) this->M_uFESpace.fe().nbNode ; k++ )
        {
            UInt iloc = this->M_uFESpace.patternFirst( k );

            for ( UInt ic = 0; ic < nbCompU; ++ic )
            {
                UInt ig = this->M_uFESpace.dof().localToGlobal( i, iloc + 1 ) + ic * this->_dim_u;

                M_elvec( )      [ iloc + ic * this->M_uFESpace.fe().nbNode ] = betaVecRep( ig );             // u^n - w^k local
                M_w_loc.vec( )  [ iloc + ic * this->M_uFESpace.fe().nbNode ] = wRep( ig );  // w^k local
                M_uk_loc.vec( ) [ iloc + ic * this->M_uFESpace.fe().nbNode ] = betaVecRep( ig );      // u^k local
                M_d_loc.vec( )  [ iloc + ic * this->M_uFESpace.fe().nbNode ] = dispRep( ig );  // d local
                M_dw_loc.vec( ) [ iloc + ic * this->M_uFESpace.fe().nbNode ] = dwRep( ig ); // dw local
            }
        }

        for ( UInt k = 0 ; k < ( UInt ) this->M_pFESpace.nbNode ; k++ )
        {
            UInt iloc = this->M_pFESpace.fe().patternFirst( k );
            UInt ig   = this->M_pFESpace.dof().localToGlobal( i, iloc + 1 );

            M_pk_loc[ iloc ] = this->p( ig );  // p^k local
        }

        //
        // Elementary vectors
        //

        //  - \rho ( \grad( u^n-w^k ):[I\div d - (\grad d)^T] u^k + ( u^n-w^k )^T[I\div d - (\grad d)^T] (\grad u^k)^T , v  )
        source_mass1( - this->density(), M_uk_loc, M_elvec, M_d_loc, M_elvec_du, this->M_uFESpace );

        //  + \rho * ( \grad u^k dw, v  )
        source_mass2( this->density(), M_uk_loc, M_dw_loc, M_elvec_du, this->M_uFESpace );

        //  - ( [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] , \grad v  )
        source_stress( -1.0, this->viscosity(), M_uk_loc, M_pk_loc, M_d_loc, M_elvec_du, this->M_uFESpace, this->M_pFESpace );

        // + \mu ( \grad u^k \grad d + [\grad d]^T[\grad u^k]^T : \grad v )
        source_stress2( this->viscosity(), M_uk_loc, M_d_loc, M_elvec_du, this->M_uFESpace );

        //  + ( (\grad u^k):[I\div d - (\grad d)^T] , q  )
        source_press( 1.0, M_uk_loc, M_d_loc, M_elvec_dp, this->M_uFESpace, this->M_pFESpace );

        //
        // Assembling
        //

        // assembling presssure right hand side
        assemb_vec( M_rhsLinNoBC, M_elvec_dp, this->M_pFESpace, this->_dof_p, nbCompU );

        // loop on velocity components
        for ( UInt ic = 0; ic < nbCompU; ic++ )
        {
            // assembling velocity right hand side
            assemb_vec( M_rhsLinNoBC, M_elvec_du, this->M_uFESpace, this->_dof_u, ic );
        }
    }


    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

}

}


#endif
