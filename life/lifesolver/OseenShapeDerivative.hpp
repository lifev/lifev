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
    typedef typename super::prec_type           prec_type;
    typedef typename super::prec_raw_type       prec_raw_type;

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

    OseenShapeDerivative( const data_type&          dataType,
                          FESpace<Mesh, EpetraMap>& uFESpace,
                          FESpace<Mesh, EpetraMap>& pFESpace,
                          Epetra_Comm&              comm,
                          const EpetraMap           bigMap,
                          const UInt                offset=0);

    ~OseenShapeDerivative();

    void setUp( const GetPot& dataFile );

    void iterateLin( bchandler_raw_type& bch );

    void updateLinearSystem( const matrix_type& matrNoBC,
                             double       alpha,
                             const vector_type& un,
                             const vector_type& uk,
                             const vector_type& disp,
                             const vector_type& w,
                             const vector_type& dw,
                             const vector_type& sourceVec);

private:


    vector_type               M_rhsLinNoBC;
    vector_type               M_rhsLinFull;

    vector_type               M_linSol;

    SolverType                M_linearLinSolver;
    prec_type                 M_linPrec;


//    ElemVec                   M_elvec_du; // Elementary right hand side for the linearized velocity
    ElemVec                   M_elvec_du; // Elementary right hand side for the linearized pressure
    ElemVec                   M_elvec_dp; // Elementary right hand side for the linearized pressure
    ElemVec                   M_w_loc;    // Elementary mesh velocity
    ElemVec                   M_uk_loc;   // Elementary velocity
    ElemVec                   M_pk_loc;   // Elementary pressure
    ElemVec                   M_elvec;    // Elementary convection velocity
    ElemVec                   M_d_loc;    // Elementary displacement for right hand side
    ElemVec                   M_dw_loc;   // Elementary mesh velocity for right hand side
    ElemVec                   M_u_loc;


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
    M_linSol         ( this->getMap()),
    M_linearLinSolver( ),
    M_linPrec        ( new prec_raw_type() ),
//    M_elvec_du       ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_elvec_du       ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_elvec_dp       ( this->M_pFESpace.fe().nbNode, 1 ),
    M_w_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_uk_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_pk_loc         ( this->M_pFESpace.fe().nbNode, 1 ),
    M_elvec          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_d_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_dw_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_u_loc          ( this->M_uFESpace.fe().nbNode, nDimensions )
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
    M_linSol         ( this->getMap()),
    M_linearLinSolver( ),
    M_linPrec        ( new prec_raw_type() ),
    M_elvec_du       ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_elvec_dp       ( this->M_pFESpace.fe().nbNode, 1 ),
//    M_elvec_dp       ( this->M_pFESpace.fe().nbNode, nDimensions ),
    M_w_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_uk_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_pk_loc         ( this->M_pFESpace.fe().nbNode, 1 ),
    M_elvec          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_d_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_dw_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_u_loc          ( this->M_uFESpace.fe().nbNode, nDimensions )
{

}

template<typename Mesh, typename SolverType>
OseenShapeDerivative<Mesh, SolverType>::
OseenShapeDerivative( const data_type&          dataType,
                      FESpace<Mesh, EpetraMap>& uFESpace,
                      FESpace<Mesh, EpetraMap>& pFESpace,
                      Epetra_Comm&              comm,
                      const EpetraMap           bigMap,
                      const UInt                offset):
    super            (dataType,
                      uFESpace,
                      pFESpace,
                      comm,
                      bigMap,
                      offset),
    M_rhsLinNoBC     ( this->getMap()),
    M_rhsLinFull     ( this->getMap()),
    M_linSol         ( this->getMap()),
    M_linearLinSolver( ),
    M_linPrec        ( new prec_raw_type() ),
    M_elvec_du       ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_elvec_dp       ( this->M_pFESpace.fe().nbNode, 1 ),
//    M_elvec_dp       ( this->M_pFESpace.fe().nbNode, nDimensions ),
    M_w_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_uk_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_pk_loc         ( this->M_pFESpace.fe().nbNode, 1 ),
    M_elvec          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_d_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_dw_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_u_loc          ( this->M_uFESpace.fe().nbNode, nDimensions )
{

}


template<typename Mesh, typename SolverType>
OseenShapeDerivative<Mesh, SolverType>::
~OseenShapeDerivative()
{

}

template<typename Mesh, typename SolverType>
void OseenShapeDerivative<Mesh, SolverType>::setUp( const GetPot& dataFile )
{
    super::setUp( dataFile );

    M_linearLinSolver.setDataFromGetPot( dataFile, "lin_fluid/solver" );
    M_linPrec->setDataFromGetPot( dataFile, "lin_fluid/prec" );

}



template<typename Mesh, typename SolverType>
void OseenShapeDerivative<Mesh, SolverType>::iterateLin( bchandler_raw_type& bch )
{
    Chrono chrono;

    // matrix and vector assembling communication

    this->leaderPrint("  f-  Finalizing the matrix and vectors ...    ");

    chrono.start();


    this->M_matrNoBC->GlobalAssemble();
    if (this->M_stab)
        this->M_matrStab->GlobalAssemble();

    M_rhsLinNoBC.GlobalAssemble();

    matrix_ptrtype matrFull( new matrix_type( this->M_localMap, this->M_matrNoBC->getMeanNumEntries()));
    *matrFull += *this->M_matrNoBC;
    if (this->M_stab)
        *matrFull += *this->M_matrStab;

    vector_type    rhsFull(M_rhsLinNoBC);

    chrono.stop();

    this->leaderPrintMax("done in " , chrono.diff());

    // boundary conditions update

    this->leaderPrint("  f-  Applying boundary conditions...          ");

//    std::cout << "    norm_inf( rhsFullNoBC )     = " << rhsFull.NormInf() << std::endl;

    chrono.start();

    bch.bdUpdate( *this->M_uFESpace.mesh(), this->M_uFESpace.feBd(), this->M_uFESpace.dof() );
    applyBoundaryConditions( *matrFull, rhsFull, bch );

    chrono.stop();

//    this->M_comm->Barrier();

    this->leaderPrintMax("done in ", chrono.diff() );

    this->leaderPrint("    norm_inf( rhsFull )     = " , rhsFull.NormInf());
    // solving the system

    // using the same preconditioner as for the non linear problem (the matrix changes only in the
    // boundary terms.
    solveSystem( matrFull, rhsFull, M_linSol, M_linearLinSolver, this->M_prec);

    this->M_residual  = M_rhsLinNoBC;
    this->M_residual -= *this->M_matrNoBC*this->M_linSol;

    this->leaderPrintMax( "NormInf Residual Lin = " , this->M_residual.NormInf());
    this->leaderPrintMax( "NormInf Solution Lin = " , this->M_linSol.NormInf());
} // iterateLin



template<typename Mesh, typename SolverType>
void
OseenShapeDerivative<Mesh, SolverType>::updateLinearSystem( const matrix_type& matrNoBC,
                                                            double       /*alpha*/,
                                                            const vector_type& un,
                                                            const vector_type& uk,
                                                            const vector_type& disp,
                                                            const vector_type& w,
                                                            const vector_type& dw,
                                                            const vector_type& sourceVec)
{
    this->leaderPrint("  f-  LINEARIZED FLUID SYSTEM\n");

    Chrono chrono;

     this->M_matrNoBC.reset( new matrix_type(matrNoBC));

    int nbCompU = nDimensions;

    M_rhsLinNoBC = sourceVec;//which is usually zero

    if(this->M_data.useShapeDerivatives())
        {
            this->leaderPrint("  f-  Updating right hand side... ");

            //
            // RIGHT HAND SIDE FOR THE LINEARIZED ALE SYSTEM
            //
            chrono.start();

            // Loop on elements

            vector_type vel(this->M_uFESpace.map());
            vector_type press(this->M_pFESpace.map());

            vel.subset(uk);
            press.subset(uk, this->M_uFESpace.dim()*this->M_uFESpace.fieldDim());

            vector_type unRep  ( un  , Repeated );
            vector_type ukRep  ( uk  , Repeated );
            vector_type dispRep( disp, Repeated );
            vector_type wRep   ( w   , Repeated );
            vector_type dwRep  ( dw  , Repeated );

//     std::cout << wRep.NormInf() << std::endl;
//     std::cout << dwRep.NormInf() << std::endl;
//     std::cout << dispRep.NormInf() << std::endl;

            vector_type rhsLinNoBC( M_rhsLinNoBC, Repeated, Zero); // ignoring non-local entries, Otherwise they are summed up lately



            for ( UInt i = 1; i <= this->M_uFESpace.mesh()->numVolumes(); i++ )
                {

                    this->M_pFESpace.fe().update( this->M_pFESpace.mesh()->volumeList( i ) );
                    this->M_uFESpace.fe().updateFirstDerivQuadPt( this->M_uFESpace.mesh()->volumeList( i ) );

                    // initialization of elementary vectors
                    M_elvec_du.zero();
                    M_elvec_dp.zero();

                    for ( UInt k = 0 ; k < ( UInt ) this->M_uFESpace.fe().nbNode ; k++ )
                        {
                            UInt iloc = this->M_uFESpace.fe().patternFirst( k ); // iloc = k

                            for ( int ic = 0; ic < nbCompU; ++ic )
                                {
                                    UInt ig    = this->M_uFESpace.dof().localToGlobal( i, iloc + 1 ) + ic * this->dim_u();

                                    M_elvec.vec( )  [ iloc + ic*this->M_uFESpace.fe().nbNode ] = unRep(ig) - wRep( ig ); // u^n - w^k local
                                    M_w_loc.vec( )  [ iloc + ic*this->M_uFESpace.fe().nbNode ] = wRep( ig );             // w^k local
                                    M_uk_loc.vec( ) [ iloc + ic*this->M_uFESpace.fe().nbNode ] = ukRep( ig );            // u^k local
                                    M_d_loc.vec( )  [ iloc + ic*this->M_uFESpace.fe().nbNode ] = dispRep( ig );          // d local
                                    M_dw_loc.vec( ) [ iloc + ic*this->M_uFESpace.fe().nbNode ] = dwRep( ig );            // dw local
                                    M_u_loc.vec()   [ iloc + ic*this->M_uFESpace.fe().nbNode ] = unRep( ig );            // un local
                                }
                        }
                    /*                    std::cout << M_elvec.vec() << std::endl;
                    std::cout << M_w_loc.vec() << std::endl;
                    std::cout << M_uk_loc.vec() << std::endl;
                    std::cout << M_d_loc.vec() << std::endl;
                    std::cout << M_dw_loc.vec() << std::endl;
                    std::cout << M_u_loc.vec() << std::endl;*/

                    for ( UInt k = 0 ; k < ( UInt ) this->M_pFESpace.fe().nbNode ; k++ )
                        {
                            UInt iloc = this->M_pFESpace.fe().patternFirst( k ); // iloc = k
                            UInt ig   = this->M_pFESpace.dof().localToGlobal( i, iloc + 1 ) + nbCompU*this->dim_u();
                            M_pk_loc[ iloc ] = ukRep[ ig ];  // p^k local

                        }

        //
        // Elementary vectors
        //
                    //commented the code to print out the elementary data. Useful for debugging.

                    //  - \rho ( \grad( u^n-w^k ):[I\div d - (\grad d)^T] u^k + ( u^n-w^k )^T[I\div d - (\grad d)^T] (\grad u^k)^T , v  )
                    source_mass1( - this->M_data.density(), M_uk_loc, M_w_loc, M_elvec, M_d_loc, M_elvec_du, this->M_uFESpace.fe() );
                    //                    std::cout << "source_mass1 -> norm_inf(M_elvec_du)" << std::endl;
                    //                    M_elvec_du.showMe(std::cout);

                    //  + \rho * ( \grad u^k dw, v  )
                    source_mass2( this->M_data.density(), M_uk_loc, M_dw_loc, M_elvec_du, this->M_uFESpace.fe() );
                    //std::cout << "source_mass2 -> norm_inf(M_elvec_du)" << std::endl;
                    //M_elvec_du.showMe(std::cout);

                    //  - \rho/2 ( \nabla u^n:[2 * I\div d - (\grad d)^T]  u^k , v  )
                    source_mass3( - 0.5*this->M_data.density(), M_u_loc, M_uk_loc, M_d_loc, M_elvec_du, this->M_uFESpace.fe() );
                    //std::cout << "source_mass3 -> norm_inf(M_elvec_du)" << std::endl;
                    //M_elvec_du.showMe(std::cout);

                    //  - ( [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] , \grad v  )
                    source_stress( - 1.0, this->M_data.viscosity(), M_uk_loc, M_pk_loc, M_d_loc, M_elvec_du, this->M_uFESpace.fe(), this->M_pFESpace.fe() );
                    //std::cout << "source_stress -> norm_inf(M_elvec_du)" << std::endl;
                    //M_elvec_du.showMe(std::cout);

                    // + \mu ( \grad u^k \grad d + [\grad d]^T[\grad u^k]^T : \grad v )
                    source_stress2( this->M_data.viscosity(), M_uk_loc, M_d_loc, M_elvec_du, this->M_uFESpace.fe() );
                    //std::cout << "source_stress2 -> norm_inf(M_elvec_du)" << std::endl;
                    //M_elvec_du.showMe(std::cout);

                    //  + ( (\grad u^k):[I\div d - (\grad d)^T] , q  )
                    source_press( 1.0, M_uk_loc, M_d_loc, M_elvec_dp, this->M_uFESpace.fe(), this->M_pFESpace.fe() );
                    //std::cout << "source_press -> norm_inf(M_elvec_du)"  << std::endl;
                    //M_elvec_dp.showMe(std::cout);

                    //
                    // Assembling
                    //

//         std::cout << "debut ====================" << std::endl;
//         M_elvec_dp.showMe(std::cout);
//         M_elvec_du.showMe(std::cout);
//         std::cout << "fin   ====================" << std::endl;

                    // assembling presssure right hand side
                    assembleVector( rhsLinNoBC, M_elvec_dp, this->M_pFESpace.fe(), this->M_pFESpace.dof(), 0, nbCompU*this->dim_u() );

           // loop on velocity components
                    for ( int ic = 0; ic < nbCompU; ic++ )
                        {
                            // assembling velocity right hand side
                            assembleVector( rhsLinNoBC, M_elvec_du, this->M_uFESpace.fe(), this->M_uFESpace.dof(), ic, ic*this->dim_u() );
                        }
                }


            M_rhsLinNoBC = rhsLinNoBC;

            this->leaderPrint( "norm( M_rhsLinNoBC)  = " , M_rhsLinNoBC.NormInf() );

        }

    chrono.stop();
    this->leaderPrintMax("done in ", chrono.diff() );

}

}


#endif
