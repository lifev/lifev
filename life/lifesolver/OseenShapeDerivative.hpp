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
          typename SolverType = LifeV::SolverTrilinos >
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
                          boost::shared_ptr<Epetra_Comm>& comm,
                          const int                 lagrangeMultiplier = 0);

    OseenShapeDerivative( const data_type&          dataType,
                          FESpace<Mesh, EpetraMap>& uFESpace,
                          FESpace<Mesh, EpetraMap>& pFESpace,
                          boost::shared_ptr<Epetra_Comm>& comm,
                          const EpetraMap           bigMap,
                          const UInt                offset=0);

    OseenShapeDerivative( const data_type&          dataType,
                          FESpace<Mesh, EpetraMap>& uFESpace,
                          FESpace<Mesh, EpetraMap>& pFESpace,
                          FESpace<Mesh, EpetraMap>& mmFESpace,
                          boost::shared_ptr<Epetra_Comm>& comm,
                          const EpetraMap           bigMap,
                          const UInt                offset=0);

    ~OseenShapeDerivative();

    void setUp( const GetPot& dataFile );

    void iterateLin( bchandler_raw_type& bch );

    void updateLinearSystem( const matrix_type& matrNoBC,
                             double&       alpha,
                             const vector_type& un,
                             const vector_type& uk,
                             const vector_type& disp,
                             const vector_type& w,
                             const vector_type& dw,
                             const vector_type& sourceVec);
    //getter
    vector_type& rhsLinNoBC() {return M_rhsLinNoBC;}

    Real GetLinearFlux    ( const EntityFlag& flag );
    Real GetLinearPressure( const EntityFlag& flag );

    //! Get the lagrange multiplier related to a flux imposed on a given part of the boundary.
    /*!
     * @param Flag flag of the boundary face associated with the flux and the Lagrange multiplier we want.
     * @param BC BChandler containing the boundary conditions of the problem.
     * @return Lagrange multiplier
     */
    Real LinearLagrangeMultiplier( const EntityFlag& Flag, bchandler_raw_type& BC );

    //! Get the solution of the Shape Derivative problem.
    /*!
     * @return vector containing the solution of the Shape Derivative problem.
     */
    const vector_type& LinearSolution() const { return M_linSol; }

    void updateShapeDerivatives(
                                matrix_type& matrNoBC,
                                double&       alpha,
                                const vector_type& un,
                                const vector_type& uk,
                                //const vector_type& disp,
                                const vector_type& w,
                                UInt offset,
                                FESpace<Mesh, EpetraMap>& dFESpace,
                                bool wImplicit=true,
                                bool convectiveTermDerivative=false);

    void updateRhsLinNoBC( const vector_type& rhs){M_rhsLinNoBC=rhs;}
    bool stab()         { return this->M_stab; }
private:


    vector_type               M_rhsLinNoBC;
    vector_type               M_rhsLinFull;

    vector_type               M_linSol;

    SolverType                M_linearLinSolver;
    prec_type                 M_linPrec;


//    ElemVec                   M_elvec_du; // Elementary right hand side for the linearized velocity
    ElemVec                   M_elvec_du; // Elementary right hand side for the linearized pressure
    //    boost::shared_ptr<ElemMat>                   M_elmat_du;
    //    boost::shared_ptr<ElemMat>                   M_elmat_du_convective;
    ElemVec                   M_elvec_dp; // Elementary right hand side for the linearized pressure
    //    boost::shared_ptr<ElemMat>                   M_elmat_dp;    // Elementary displacement for right hand side
    ElemVec                   M_w_loc;    // Elementary mesh velocity
    ElemVec                   M_uk_loc;   // Elementary velocity
    ElemVec                   M_pk_loc;   // Elementary pressure
    ElemVec                   M_elvec;    // Elementary convection velocity
    ElemVec                   M_d_loc;    // Elementary displacement for right hand side
    ElemVec                   M_dw_loc;   // Elementary mesh velocity for right hand side
    ElemVec                   M_u_loc;
    bool                      M_reusePrecLin;
    FESpace<Mesh, EpetraMap>* M_mmFESpace;
};


template<typename Mesh, typename SolverType>
OseenShapeDerivative<Mesh, SolverType>::
OseenShapeDerivative( const data_type&          dataType,
                      FESpace<Mesh, EpetraMap>& uFESpace,
                      FESpace<Mesh, EpetraMap>& pFESpace,
                      boost::shared_ptr<Epetra_Comm>& comm,
                      const int                 lagrangeMultiplier):
    super            (dataType,
                      uFESpace,
                      pFESpace,
                      comm,
                      lagrangeMultiplier),
    M_rhsLinNoBC     ( this->getMap()),
    M_rhsLinFull     ( this->getMap()),
    M_linSol         ( this->getMap()),
    M_linearLinSolver( comm ),
    M_linPrec        ( ),
    M_elvec_du       ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_elvec_dp       ( this->M_pFESpace.fe().nbNode, 1 ),
//    M_elvec_dp       ( this->M_pFESpace.fe().nbNode, nDimensions ),
    M_w_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_uk_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_pk_loc         ( this->M_pFESpace.fe().nbNode, 1 ),
    M_elvec          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_d_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_dw_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_u_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_reusePrecLin   (true)
{

}

template<typename Mesh, typename SolverType>
OseenShapeDerivative<Mesh, SolverType>::
OseenShapeDerivative( const data_type&          dataType,
                      FESpace<Mesh, EpetraMap>& uFESpace,
                      FESpace<Mesh, EpetraMap>& pFESpace,
                      boost::shared_ptr<Epetra_Comm>& comm,
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
    M_linearLinSolver( comm ),
    M_linPrec        ( ),
    M_elvec_du       ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_elvec_dp       ( this->M_pFESpace.fe().nbNode, 1 ),
//    M_elvec_dp       ( this->M_pFESpace.fe().nbNode, nDimensions ),
    M_w_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_uk_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_pk_loc         ( this->M_pFESpace.fe().nbNode, 1 ),
    M_elvec          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_d_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_dw_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_u_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_reusePrecLin   (true)
{

}

template<typename Mesh, typename SolverType>
OseenShapeDerivative<Mesh, SolverType>::
OseenShapeDerivative( const data_type&          dataType,
                      FESpace<Mesh, EpetraMap>& uFESpace,
                      FESpace<Mesh, EpetraMap>& pFESpace,
                      FESpace<Mesh, EpetraMap>& mmFESpace,
                      boost::shared_ptr<Epetra_Comm>& comm,
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
    M_linPrec        ( ),
    M_elvec_du       ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_elvec_dp       ( this->M_pFESpace.fe().nbNode, 1 ),
//    M_elvec_dp       ( this->M_pFESpace.fe().nbNode, nDimensions ),
    M_w_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_uk_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_pk_loc         ( this->M_pFESpace.fe().nbNode, 1 ),
    M_elvec          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_d_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_dw_loc         ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_u_loc          ( this->M_uFESpace.fe().nbNode, nDimensions ),
    M_reusePrecLin   (true),
    M_mmFESpace  (&mmFESpace)
{

}


template<typename Mesh, typename SolverType>
OseenShapeDerivative<Mesh, SolverType>::
~OseenShapeDerivative()
{

}

template<typename Mesh, typename SolverType>
Real
OseenShapeDerivative<Mesh, SolverType>::GetLinearFlux( const EntityFlag& flag )
{
    return flux( flag, M_linSol );
}

template<typename Mesh, typename SolverType>
Real
OseenShapeDerivative<Mesh, SolverType>::GetLinearPressure( const EntityFlag& flag )
{
    return pressure( flag, M_linSol );
}

template<typename Mesh, typename SolverType>
Real
OseenShapeDerivative<Mesh, SolverType>::LinearLagrangeMultiplier( const EntityFlag& Flag, bchandler_raw_type& BC )
{
  return LagrangeMultiplier( Flag, BC, M_linSol );
}

template<typename Mesh, typename SolverType>
void OseenShapeDerivative<Mesh, SolverType>::setUp( const GetPot& dataFile )
{
    //M_linearLinSolver.setDataFromGetPot( dataFile, "lin_fluid/solver" );
    //    M_linearLinSolver.setAztecooPreconditioner( dataFile, "lin_fluid/solver" );

    super::setUp( dataFile );

    M_reusePrecLin = dataFile( "lin_fluid/prec/reuse", true);

    //std::string precType = dataFile( "lin_fluid/prec/prectype", "Ifpack");

    //    M_linPrec            = prec_type( PRECFactory::instance().createObject( precType ) ); //linPrec is not used
    //    M_linPrec->setDataFromGetPot( dataFile, "lin_fluid/prec" );

}



template<typename Mesh, typename SolverType>
void OseenShapeDerivative<Mesh, SolverType>::iterateLin( bchandler_raw_type& bch )
{
    this->M_Displayer.leaderPrint(" LF-  Finalizing the matrix and vectors ...    ");

    Chrono chrono;
    chrono.start();

    // matrix and vector assembling communication
    this->M_matrNoBC->GlobalAssemble();

    M_rhsLinNoBC.GlobalAssemble();

    matrix_ptrtype matrFull( new matrix_type( this->M_localMap, this->M_matrNoBC->getMeanNumEntries()));

    updateStab(*matrFull);
    getFluidMatrix(*matrFull);

    vector_type    rhsFull(M_rhsLinNoBC);

    chrono.stop();
    this->M_Displayer.leaderPrintMax("done in " , chrono.diff());

    // boundary conditions update
    this->M_Displayer.leaderPrint(" LF-  Applying boundary conditions ...         ");
    chrono.start();

    applyBoundaryConditions( *matrFull, rhsFull, bch );

    chrono.stop();
    this->M_Displayer.leaderPrintMax("done in ", chrono.diff() );

    // solving the system

    // using the same preconditioner as for the non linear problem (the matrix changes only in the
    // boundary terms).
    matrFull->GlobalAssemble();
    this->M_linearSolver.setMatrix(*matrFull);
    this->M_linearSolver.setReusePreconditioner( M_reusePrecLin );
    this->M_linearSolver.solveSystem( rhsFull, M_linSol, matrFull );

    this->M_residual  = M_rhsLinNoBC;
    this->M_residual -= *this->M_matrNoBC*this->M_linSol;

//     if(S_verbose)
//         {
//             this->M_Displayer.leaderPrintMax( "NormInf Residual Lin = " , this->M_residual.NormInf());
//             this->M_Displayer.leaderPrintMax( "NormInf Solution Lin = " , this->M_linSol.NormInf());
//         }
} // iterateLin

template<typename Mesh, typename SolverType>
void
OseenShapeDerivative<Mesh, SolverType>::updateLinearSystem( const matrix_type& /*matrNoBC*/, //Fluid Matrix withoud BC
                                                                  double&      /*alpha*/,    //alpha
                                                            const vector_type& un,           //Beta
                                                            const vector_type& uk,           //Fluid Solution
                                                            const vector_type& disp,         //Mesh deltaX
                                                            const vector_type& w,            //Mesh Velocity
                                                            const vector_type& dw,           //Mesh deltaVelocity
                                                            const vector_type& sourceVec)    //RHS (usually 0 )
{
    this->M_Displayer.leaderPrint(" LF-  Updating the right hand side ...         ");
    Chrono chrono;
    chrono.start();

    int nbCompU = nDimensions;
    M_rhsLinNoBC = sourceVec;//which is usually zero

    if(this->M_data.useShapeDerivatives())
        {
            //
            // RIGHT HAND SIDE FOR THE LINEARIZED ALE SYSTEM
            //

            // Loop on elements
            vector_type unRep  ( un  , Repeated );
            vector_type ukRep  ( uk  , Repeated );
            vector_type dispRep( disp, Repeated );
            vector_type wRep   ( w   , Repeated );
            vector_type dwRep  ( dw  , Repeated );

            //     std::cout << wRep.NormInf() << std::endl;
            //     std::cout << dwRep.NormInf() << std::endl;
            //     std::cout << dispRep.NormInf() << std::endl;

            vector_type rhsLinNoBC( M_rhsLinNoBC.getMap(), Repeated);

            for ( UInt i = 1; i <= this->M_uFESpace.mesh()->numVolumes(); i++ )
                {

                    this->M_pFESpace.fe().update( this->M_pFESpace.mesh()->volumeList( i ) );
                    this->M_uFESpace.fe().updateFirstDerivQuadPt( this->M_uFESpace.mesh()->volumeList( i ) );

                    // initialization of elementary vectors
                    M_elvec_du.zero();
                    M_elvec_dp.zero();

                    for ( UInt k = 0 ; k < this->M_uFESpace.fe().nbNode ; k++ )
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
                    /*
                      std::cout << M_elvec.vec() << std::endl;
                      std::cout << M_w_loc.vec() << std::endl;
                      std::cout << M_uk_loc.vec() << std::endl;
                      std::cout << M_d_loc.vec() << std::endl;
                      std::cout << M_dw_loc.vec() << std::endl;
                      std::cout << M_u_loc.vec() << std::endl;
                    */
                    for ( UInt k = 0 ; k < this->M_pFESpace.fe().nbNode ; k++ )
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
                    /*
                    std::cout << "source_mass1 -> norm_inf(M_elvec_du)" << std::endl;
                    M_elvec_du.showMe(std::cout);
                    */
                    //  + \rho * ( \grad u^k dw, v  )
                    source_mass2( this->M_data.density(), M_uk_loc, M_dw_loc, M_elvec_du, this->M_uFESpace.fe() );
                    /*
                    std::cout << "source_mass2 -> norm_inf(M_elvec_du)" << std::endl;
                    M_elvec_du.showMe(std::cout);
                    */
                    //  - \rho/2 ( \nabla u^n:[2 * I\div d - (\grad d)^T]  u^k , v  )
                    source_mass3( - 0.5*this->M_data.density(), M_u_loc, M_uk_loc, M_d_loc, M_elvec_du, this->M_uFESpace.fe() );
                    /*
                    std::cout << "source_mass3 -> norm_inf(M_elvec_du)" << std::endl;
                    M_elvec_du.showMe(std::cout);
                    */
                    //  - ( [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] , \grad v  )
                    source_stress( - 1.0, this->M_data.viscosity(), M_uk_loc, M_pk_loc, M_d_loc, M_elvec_du, this->M_uFESpace.fe(), this->M_pFESpace.fe() );
                    /*
                    std::cout << "source_stress -> norm_inf(M_elvec_du)" << std::endl;
                    M_elvec_du.showMe(std::cout);
                    */
                    // + \mu ( \grad u^k \grad d + [\grad d]^T[\grad u^k]^T : \grad v )
                    source_stress2( this->M_data.viscosity(), M_uk_loc, M_d_loc, M_elvec_du, this->M_uFESpace.fe() );
                    /*
                    std::cout << "source_stress2 -> norm_inf(M_elvec_du)" << std::endl;
                    M_elvec_du.showMe(std::cout);
                    */
                    //  + ( (\grad u^k):[I\div d - (\grad d)^T] , q  )
                    source_press( -1.0, M_uk_loc, M_d_loc, M_elvec_dp, this->M_uFESpace.fe(), this->M_pFESpace.fe() );
                    /*
                    std::cout << "source_press -> norm_inf(M_elvec_du)"  << std::endl;
                    M_elvec_dp.showMe(std::cout);
                    */
                    //
                    // Assembling
                    //
                    /*
                      std::cout << "debut ====================" << std::endl;
                      M_elvec_dp.showMe(std::cout);
                      M_elvec_du.showMe(std::cout);
                      std::cout << "fin   ====================" << std::endl;
                    */
                    // assembling pressure right hand side
                    assembleVector( rhsLinNoBC, M_elvec_dp, this->M_pFESpace.fe(), this->M_pFESpace.dof(), 0, nbCompU*this->dim_u() );

                    // loop on velocity components
                    for ( int ic = 0; ic < nbCompU; ic++ )
                        {
                            // assembling velocity right hand side
                            assembleVector( rhsLinNoBC, M_elvec_du, this->M_uFESpace.fe(), this->M_uFESpace.dof(), ic, ic*this->dim_u() );
                        }
                }

            rhsLinNoBC.GlobalAssemble();
            M_rhsLinNoBC += rhsLinNoBC;
            //            M_rhsLinNoBC *= -1.;
//             if(S_verbose)
//                 this->M_Displayer.leaderPrint( "norm( M_rhsLinNoBC)  = " , M_rhsLinNoBC.NormInf() );

        }

    chrono.stop();
    this->M_Displayer.leaderPrintMax("done in ", chrono.diff() );
}


//#if UNDEF
template<typename Mesh, typename SolverType>
void
OseenShapeDerivative<Mesh, SolverType>::updateShapeDerivatives( matrix_type& M_matr,
                                                                double&       alpha,
                                                                const vector_type& un,
                                                                const vector_type& uk,
                                                                //const vector_type& disp,
                                                                const vector_type& w,
                                                                UInt offset,
                                                                FESpace<Mesh, EpetraMap>& mmFESpace,
                                                                bool wImplicit,
                                                                bool convectiveTermDerivative
                                                                )
{
    Chrono chrono;

    UInt nbCompU = nDimensions;

    //    M_rhsLinNoBC = sourceVec;//which is usually zero

    if(this->M_data.useShapeDerivatives())
        {
            this->M_Displayer.leaderPrint(" LF-  Updating shape derivative blocks ...     ");

            //
            // RIGHT HAND SIDE FOR THE LINEARIZED ALE SYSTEM
            //
            chrono.start();

            // Loop on elements

            vector_type unRep  ( un  , Repeated );
            vector_type ukRep  ( uk  , Repeated );
            //vector_type dispRep( disp, Repeated );
            vector_type wRep   ( w   , Repeated );
            //            vector_type dwRep  ( dw  , Repeated );

//     std::cout << wRep.NormInf() << std::endl;
//     std::cout << dwRep.NormInf() << std::endl;
//     std::cout << dispRep.NormInf() << std::endl;

//            vector_type rhsLinNoBC( M_rhsLinNoBC.getMap(), Repeated);

            for ( UInt i = 1; i <= this->M_uFESpace.mesh()->numVolumes(); i++ )
                {



                    this->M_pFESpace.fe().update( this->M_pFESpace.mesh()->volumeList( i ) );
                    this->M_uFESpace.fe().updateFirstDerivQuadPt( this->M_uFESpace.mesh()->volumeList( i ) );
                    this->M_pFESpace.fe().updateFirstDerivQuadPt( this->M_uFESpace.mesh()->volumeList( i ) );

                    //this->M_pFESpace.fe().updateFirstDeriv( this->M_uFESpace.mesh()->volumeList( i ) ); // just to provide the id number in the assem_mat_mixed
                    //this->M_uFESpace.fe().updateFirstDeriv( this->M_uFESpace.mesh()->volumeList( i ) ); //as updateFirstDer
                    mmFESpace.fe().updateFirstDerivQuadPt( mmFESpace.mesh()->volumeList( i ) );

                    // initialization of elementary vectors
                    boost::shared_ptr<ElemMat> elmat_dp       ( new ElemMat(this->M_pFESpace.fe().nbNode, 1, 0, mmFESpace.fe().nbNode, 0, nDimensions ));
                    boost::shared_ptr<ElemMat> elmat_du       ( new ElemMat(this->M_uFESpace.fe().nbNode, nDimensions, 0, this->M_uFESpace.fe().nbNode, 0, nDimensions ));
                    boost::shared_ptr<ElemMat> elmat_du_convective;

                    if(convectiveTermDerivative)
                    {
                        elmat_du_convective.reset( new ElemMat(this->M_uFESpace.fe().nbNode, nDimensions, 0, mmFESpace.fe().nbNode, 0, nDimensions ));
                        elmat_du_convective->zero();
                    }



                    elmat_dp->zero();
                    elmat_du->zero();


                    for ( UInt k = 0 ; k < this->M_uFESpace.fe().nbNode ; k++ )
                        {
                            UInt iloc = this->M_uFESpace.fe().patternFirst( k ); // iloc = k

                            for ( UInt ic = 0; ic < nbCompU; ++ic )
                                {
                                    UInt ig    = this->M_uFESpace.dof().localToGlobal( i, iloc + 1 ) + ic * this->dim_u();

                                    //                                    if(!wImplicit)
                                    M_elvec.vec( )  [ iloc + ic*this->M_uFESpace.fe().nbNode ] = unRep(ig) - wRep( ig ); // u^n - w^k local
//                                     else
//                                         M_elvec.vec( )  [ iloc + ic*this->M_uFESpace.fe().nbNode ] = ukRep(ig) - wRep( ig ); // u^n - w^k local
                                    M_w_loc.vec( )  [ iloc + ic*this->M_uFESpace.fe().nbNode ] = wRep( ig );             // w^k local
                                    M_uk_loc.vec( ) [ iloc + ic*this->M_uFESpace.fe().nbNode ] = ukRep( ig );            // u^k local
                                    //M_d_loc.vec( ) [ iloc + ic*this->M_uFESpace.fe().nbNode ] = dispRep( ig );            // dw local
                                    //M_dw_loc.vec( ) [ iloc + ic*this->M_uFESpace.fe().nbNode ] = dwRep( ig );            // dw local
                                    M_u_loc.vec()   [ iloc + ic*this->M_uFESpace.fe().nbNode ] = unRep( ig );            // un local
                                }
                       }
                    /*
                    std::cout << M_elvec.vec() << std::endl;
                    std::cout << M_w_loc.vec() << std::endl;
                    std::cout << M_uk_loc.vec() << std::endl;
                    std::cout << M_d_loc.vec() << std::endl;
                    std::cout << M_dw_loc.vec() << std::endl;
                    std::cout << M_u_loc.vec() << std::endl;
                    */
                    for ( UInt k = 0 ; k < this->M_pFESpace.fe().nbNode ; k++ )
                        {
                            UInt iloc = this->M_pFESpace.fe().patternFirst( k ); // iloc = k
                            UInt ig   = this->M_pFESpace.dof().localToGlobal( i, iloc + 1 ) + nbCompU*this->dim_u();
                            M_pk_loc[ iloc ] = ukRep[ ig ];  // p^k local
                        }


                    shape_terms(
                                //M_d_loc,
                                this->M_data.density(),
                                this->M_data.viscosity(),
                                M_u_loc,
                                M_uk_loc,
                                M_w_loc,
                                M_elvec,
                                M_pk_loc,
                                *elmat_du,
                                this->M_uFESpace.fe(),
                                this->M_pFESpace.fe(),
                                (ID) mmFESpace.fe().nbNode,
                                *elmat_dp,
                                0,
                                wImplicit,
                                alpha//,
                                //elmat_du_convective
                                );

                    //elmat_du->showMe(std::cout);

                    //source_mass2( this->M_data.density(), M_uk_loc, *M_elmat_du_convective, this->M_uFESpace.fe(), alpha );

                    source_press( 1.0, M_uk_loc,*elmat_dp, this->M_uFESpace.fe(), this->M_pFESpace.fe(), (ID) mmFESpace.fe().nbNode);


                    if(convectiveTermDerivative)//derivative of the convective term
                        mass_gradu(this->M_data.density(), M_uk_loc, *elmat_du_convective, this->M_uFESpace.fe());
                      /*
                        std::cout << "source_press -> norm_inf(M_elvec_du)"  << std::endl;
                    M_elvec_dp.showMe(std::cout);
                    */
                    //
                    // Assembling
                    //
                    /*
         std::cout << "debut ====================" << std::endl;
         M_elvec_dp.showMe(std::cout);
         M_elvec_du.showMe(std::cout);
         std::cout << "fin   ====================" << std::endl;
                    */
                    UInt const velTotalDof (this->M_uFESpace.dof().numTotalDof());
                    for(UInt icomp=0; icomp<nbCompU; ++icomp)
                        {
                            for(UInt jcomp=0; jcomp<nbCompU; ++jcomp)
                            {
                                assembleMatrix( M_matr,
                                                *elmat_du,
                                                this->M_uFESpace.fe(),
                                                mmFESpace.fe(),
                                                this->M_uFESpace.dof(),
                                                mmFESpace.dof(),
                                                icomp, jcomp,
                                                icomp*velTotalDof, offset+jcomp*velTotalDof
                                                );
                                if(convectiveTermDerivative)//assembling the derivative of the convective term
                                    assembleMatrix( M_matr,
                                                    *elmat_du_convective,
                                                    this->M_uFESpace.fe(),
                                                    this->M_uFESpace.fe(),
                                                    this->M_uFESpace.dof(),
                                                    this->M_uFESpace.dof(),
                                                    icomp, jcomp,
                                                    icomp*velTotalDof, jcomp*velTotalDof
                                                    );
                            }
                            assembleMatrix( M_matr,
                                            *elmat_dp,
                                            this->M_pFESpace.fe(),
                                            mmFESpace.fe(),
                                            this->M_pFESpace.dof(),
                                            mmFESpace.dof(),
                                            (UInt)0, icomp,
                                            (UInt)nbCompU*velTotalDof, offset+icomp*velTotalDof
                                            );
                        }
                }

        }
    chrono.stop();
    this->M_Displayer.leaderPrintMax("done in ", chrono.diff() );
}

}

#endif
