/* -*- mode: c++ -*-

 This file is part of the LifeV library

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

#ifndef _NS_ASSEMBLER_H_
#define _NS_ASSEMBLER_H_

#include <life/lifealg/EpetraMap.hpp>

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>
//
#include <life/lifecore/chrono.hpp>

#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/postProc.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifesolver/nsipterms.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>

#include <life/lifearray/MatrixContainer.hpp>

//
#include <boost/shared_ptr.hpp>


#include <list>

namespace LifeV
{

/*!
  This file contains an Oseen equation solver class.
  The resulting linear systems are solved by GMRES on the full
  matrix ( u and p coupled ).
*/

template< typename Mesh >
class NS_Assembler
{

public:

    typedef Mesh                                  mesh_type;
    typedef DataNavierStokes                      data_type;
    typedef boost::shared_ptr<data_type>          data_ptr;

    typedef BCHandler                             bchandler_raw_type;
    typedef boost::shared_ptr<bchandler_raw_type> bchandler_type;

    typedef EpetraMatrix<double>      matrix_type;
    typedef boost::shared_ptr<matrix_type>        matrix_ptrtype;
    typedef EpetraVector      vector_type;
    typedef boost::shared_ptr<vector_type>        vector_ptrtype;

    //! Constructor
    /*!
      \param dataType
      \param velocity FE space
      \param pressure FE space
      \param communicator
    */
    NS_Assembler( const data_ptr&          dataType,
                  FESpace<Mesh, EpetraMap>& uFESpace,
                  FESpace<Mesh, EpetraMap>& pFESpace,
                  boost::shared_ptr<Epetra_Comm>& comm
                );

    //! virtual destructor
    virtual ~NS_Assembler();

    virtual void setUp        ( const GetPot& dataFile );

    virtual void updateMatrix(MatrixContainer<std::string> & matrixContainer,
                              const double       alpha,
                              const vector_type& betaVec,
                              bchandler_raw_type& BCh,
                              Real time = 0);

    virtual void updateRHS(vector_type & rhs,
                           bchandler_raw_type& BCh,
                           Real time = 0);


    //! returns the FeSpaces
    FESpace<Mesh, EpetraMap>& velFESpace()   {return M_uFESpace;}
    FESpace<Mesh, EpetraMap>& pressFESpace() {return M_pFESpace;}

    // return the density and the viscosity of the fluid
    Real density()   const { return M_data->density(); }
    Real viscosity() const { return M_data->viscosity(); }

    boost::shared_ptr<Epetra_Comm> const& comm()      const {return M_Displayer.comm();}

    Displayer   const& getDisplayer() const { return M_Displayer; }

protected:

    UInt dim_u() const           { return M_uFESpace.dim(); }
    UInt dim_p() const           { return M_pFESpace.dim(); }

    void buildSystem(MatrixContainer<std::string> & matrixContainer);
    void linearConvectiveTerm(const vector_type& betaVec, matrix_ptrtype A);
    void applyBoundaryConditions(  matrix_type&        C,
                                   matrix_type&        Bt,
                                   bchandler_raw_type& BCh,
                                   Real time);

    void applyBoundaryConditions(  vector_type&        rhs,
                                   bchandler_raw_type& BCh,
                                   Real time);


    //private members

    //! data for NS solvers
    data_ptr                       M_data;

    // FE spaces
    FESpace<Mesh, EpetraMap>&      M_uFESpace;
    FESpace<Mesh, EpetraMap>&      M_pFESpace;

    //! MPI communicator
    Displayer                      M_Displayer;

    bool                           M_stiffStrain;
    bool 						   M_recomputeConstantMatrices;

    UInt                           M_count;

    //private:

    //! Elementary matrices and vectors
    ElemMat                        M_elmatStiff;      // velocity Stokes
    ElemMat                        M_elmatMass;       // velocity mass
    ElemMat                        M_elmatP;          // (p,q) bloc for preconditioners
    ElemMat                        M_elmatDiv;
    ElemMat                        M_elmatGrad;
    ElemVec                        M_elvec;           // Elementary right hand side
    ElemVec                        M_uLoc;

    //! big value
    Real M_tgv;
}; // class Oseen



//
// IMPLEMENTATION
//
template<typename Mesh>
NS_Assembler<Mesh>::
NS_Assembler( const data_ptr&          dataType,
              FESpace<Mesh, EpetraMap>& uFESpace,
              FESpace<Mesh, EpetraMap>& pFESpace,
              boost::shared_ptr<Epetra_Comm>&              comm):
        M_data                   ( dataType ),
        M_uFESpace               ( uFESpace ),
        M_pFESpace               ( pFESpace ),
        M_Displayer              ( comm ),
        M_stiffStrain            ( false),
        M_recomputeConstantMatrices(true),
        M_elmatStiff             ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
        M_elmatMass              ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
        M_elmatP                 ( M_pFESpace.fe().nbNode, 1, 1 ),
        M_elmatDiv               ( M_pFESpace.fe().nbNode, 1, 0, M_uFESpace.fe().nbNode, 0, nDimensions ),
        M_elmatGrad              ( M_uFESpace.fe().nbNode, nDimensions, 0, M_pFESpace.fe().nbNode, 0, 1 ),
        M_elvec                  ( M_uFESpace.fe().nbNode, nDimensions ),
        M_uLoc                   ( M_uFESpace.fe().nbNode, nDimensions ),
        M_tgv                    (1e0)
{}


template<typename Mesh>
NS_Assembler<Mesh>::
~NS_Assembler()
{}

template<typename Mesh>
void NS_Assembler<Mesh>::setUp( const GetPot& dataFile )
{
    //FIXME: This variable should be in DataNavierStokes!!!
    M_stiffStrain = dataFile( "fluid/space_discretization/stiff_strain",false); // Enable grad( u )^T in stress tensor
}

template<typename Mesh>
void NS_Assembler<Mesh>::
updateMatrix(MatrixContainer<std::string> & matrixContainer,
             const double       alpha,
             const vector_type& betaVec,
             bchandler_raw_type& BCh,
             Real time
            )
{

    matrixContainer.setParameter("sigma", alpha);
    matrixContainer.setParameter("recomputedConstantMatrix", M_recomputeConstantMatrices);

    if (M_recomputeConstantMatrices)
        buildSystem(matrixContainer);

    M_recomputeConstantMatrices = false;

    matrix_ptrtype A, C;

    A.reset(new matrix_type(*matrixContainer.getMatrix("K")));
    M_Displayer.leaderPrint("A has been constructed\n");
    linearConvectiveTerm( betaVec, A);

    C.reset(new matrix_type(*A));
    C->add(alpha, *matrixContainer.getMatrix("M"));

    applyBoundaryConditions(*C, *matrixContainer.getMatrix("Bt"), BCh, time);

    matrixContainer.set("C", C);
    matrixContainer.set("A", A);

}

template<typename Mesh>
void NS_Assembler<Mesh>::updateRHS(vector_type & rhs,
                                   bchandler_raw_type& BCh,
                                   Real time)
{
    applyBoundaryConditions(rhs, BCh, time);
}


//===================================================================================//
// Private Functions                                                                 //
//===================================================================================//

template<typename Mesh>
void NS_Assembler<Mesh>::buildSystem(MatrixContainer<std::string> & matrixContainer)
{

    matrix_ptrtype M,K,Bt,B;

    M.reset(new matrix_type(M_uFESpace.map()));
    K.reset(new matrix_type(M_uFESpace.map()));
    Bt.reset(new matrix_type(M_uFESpace.map()));
    B.reset(new matrix_type(M_pFESpace.map()) );
    M_Displayer.leaderPrint("  F-  Computing constant matrices ...          ");

    Chrono chrono;

    Chrono chronoDer;
    Chrono chronoStiff;
    Chrono chronoMass;
    Chrono chronoGrad;

    Chrono chronoStiffAssemble;
    Chrono chronoMassAssemble;
    Chrono chronoGradAssemble;
    Chrono chronoDivAssemble;
    Chrono chronoStab;
    Chrono chronoZero;

    // Number of velocity components
    UInt nbCompU = nDimensions;

    // Elementary computation and matrix assembling
    // Loop on elements

    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();
//    UInt pressTotalDof = M_pFESpace.dof().numTotalDof();

    chrono.start();

    for ( UInt iVol = 1; iVol <= M_uFESpace.mesh()->numVolumes(); iVol++ )
    {
        chronoDer.start();
        M_pFESpace.fe().update( M_uFESpace.mesh()->volumeList( iVol ) ); // just to provide the id number in the assem_mat_mixed
        M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) );

        chronoDer.stop();

        chronoZero.start();
        M_elmatStiff.zero();
        M_elmatMass.zero();
        M_elmatP.zero();
        M_elmatDiv.zero();
        M_elmatGrad.zero();
        chronoZero.stop();

        // stiffness
        chronoStiff.start();
        if ( M_stiffStrain )
            stiff_strain( 2.0*M_data->viscosity(), M_elmatStiff, M_uFESpace.fe() );
        else
            stiff( M_data->viscosity(), M_elmatStiff,  M_uFESpace.fe(), 0, 0, nDimensions );
        chronoStiff.stop();

        // mass
        chronoMass.start();
        mass( M_data->density(), M_elmatMass, M_uFESpace.fe(), 0, 0, nDimensions );
        chronoMass.stop();

        for ( UInt iComp = 0; iComp < nbCompU; iComp++ )
        {
            // stiffness
            chronoStiffAssemble.start();
            if ( M_stiffStrain ) // sigma = 0.5 * mu (grad( u ) + grad ( u )^T)
            {
                for ( UInt jComp = 0; jComp < nbCompU; jComp++ )
                    assembleMatrix( *K,
                                    M_elmatStiff,
                                    M_uFESpace.fe(),
                                    M_uFESpace.fe(),
                                    M_uFESpace.dof(),
                                    M_uFESpace.dof(),
                                    iComp, jComp,
                                    iComp*velTotalDof, jComp*velTotalDof);
            }
            else // sigma = mu grad( u )
            {
                assembleMatrix( *K,
                                M_elmatStiff,
                                M_uFESpace.fe(),
                                M_uFESpace.fe(),
                                M_uFESpace.dof(),
                                M_uFESpace.dof(),
                                iComp, iComp,
                                iComp*velTotalDof, iComp*velTotalDof);
            }
            chronoStiffAssemble.stop();

            // mass
            chronoMassAssemble.start();
            assembleMatrix( *M,
                            M_elmatMass,
                            M_uFESpace.fe(),
                            M_uFESpace.fe(),
                            M_uFESpace.dof(),
                            M_uFESpace.dof(),
                            iComp, iComp,
                            iComp*velTotalDof, iComp*velTotalDof);
            chronoMassAssemble.stop();

            // div
            chronoGrad.start();
            grad( iComp,  1.0, M_elmatGrad, M_uFESpace.fe(), M_pFESpace.fe(), iComp,     0 );
            chronoGrad.stop();

            chronoGradAssemble.start();
            assembleMatrix( *Bt,
                            M_elmatGrad,
                            M_uFESpace.fe(),
                            M_pFESpace.fe(),
                            M_uFESpace.dof(),
                            M_pFESpace.dof(),
                            iComp, 0,
                            iComp*velTotalDof, 0
                          );
            chronoGradAssemble.stop();

            chronoDivAssemble.start();
            assembleTransposeMatrix( *B,
                                     1.,
                                     M_elmatGrad,
                                     M_pFESpace.fe(),
                                     M_uFESpace.fe(),
                                     M_pFESpace.dof(),
                                     M_uFESpace.dof(),
                                     0 , iComp,
                                     0, iComp*velTotalDof
                                   );
            chronoDivAssemble.stop();
        }
    }

    comm()->Barrier();

    chrono.stop();
    M_Displayer.leaderPrintMax( "done in " , chrono.diff());


    M_Displayer.leaderPrint( "  F-  Finalizing the matrices ...              ");

    chrono.start();

    K->GlobalAssemble();
    M->GlobalAssemble();
    B->GlobalAssemble(M_uFESpace.mapPtr(), M_pFESpace.mapPtr());
    Bt->GlobalAssemble(M_pFESpace.mapPtr(), M_uFESpace.mapPtr());

    matrixContainer.set("M", M);
    matrixContainer.set("B", B);
    matrixContainer.set("Bt", Bt);
    matrixContainer.set("K", K);



    chrono.stop();
    M_Displayer.leaderPrintMax("done in " , chrono.diff() );

    if (false)
        std::cout << "partial times:  \n"
                  << " Der            " << chronoDer.diffCumul() << " s.\n"
                  << " Zero           " << chronoZero.diffCumul() << " s.\n"
                  << " Stiff          " << chronoStiff.diffCumul() << " s.\n"
                  << " Stiff Assemble " << chronoStiffAssemble.diffCumul() << " s.\n"
                  << " Mass           " << chronoMass.diffCumul() << " s.\n"
                  << " Mass Assemble  " << chronoMassAssemble.diffCumul() << " s.\n"
                  << " Grad           " << chronoGrad.diffCumul() << " s.\n"
                  << " Grad Assemble  " << chronoGradAssemble.diffCumul() << " s.\n"
                  << " Div Assemble   " << chronoDivAssemble.diffCumul() << " s.\n"
                  << std::endl;

}



template<typename Mesh>
void NS_Assembler<Mesh>::
linearConvectiveTerm(const vector_type& betaVec, matrix_ptrtype A)
{
    Chrono chrono;

    chrono.start();

    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();
//    UInt pressTotalDof = M_pFESpace.dof().numTotalDof();

    // Right hand side for the velocity at time

    chrono.stop();

    M_Displayer.leaderPrintMax("done in ", chrono.diff());

    UInt nbCompU       = nDimensions;

    //! managing the convective term

    double normInf;
    betaVec.NormInf(&normInf);

    if (normInf != 0.)
    {
        EpetraMatrix<double> N(M_uFESpace.map());
        M_Displayer.leaderPrint("  F-  Sharing convective term ...              ");
        chrono.start();

        // vector with repeated nodes over the processors

        vector_type betaVecRep(betaVec, Repeated );

        chrono.stop();

        M_Displayer.leaderPrintMax( "done in " , chrono.diff() );
        M_Displayer.leaderPrint("  F-  Updating the convective terms ...        ");
        chrono.start();

        for ( UInt iVol = 1; iVol<= M_uFESpace.mesh()->numVolumes(); ++iVol )
        {

            M_pFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) ); // just to provide the id number in the assem_mat_mixed
            M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) ); //as updateFirstDer

            M_elmatStiff.zero();

            UInt eleID = M_uFESpace.fe().currentLocalId();
            // Non linear term, Semi-implicit approach
            // M_elvec contains the velocity values in the nodes
            for ( UInt iNode = 0 ; iNode < M_uFESpace.fe().nbNode ; iNode++ )
            {
                UInt  iloc = M_uFESpace.fe().patternFirst( iNode );
                for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
                {
                    UInt ig = M_uFESpace.dof().localToGlobal( eleID, iloc + 1 ) + iComp*dim_u();
                    M_elvec.vec()[ iloc + iComp*M_uFESpace.fe().nbNode ] = betaVecRep[ig]; // BASEINDEX + 1
                }
            }

            // compute local convective terms
            advection( M_data->density(), M_elvec, M_elmatStiff, M_uFESpace.fe(), 0, 0, nbCompU );

            //stiff_sd( 1.0, M_elvec, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.fe(), 0, 0, nbCompU );


            // loop on components
            for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
            {
                // compute local convective term and assembling

                assembleMatrix( N,
                                M_elmatStiff,
                                M_uFESpace.fe(),
                                M_uFESpace.fe(),
                                M_uFESpace.dof(),
                                M_uFESpace.dof(),
                                iComp, iComp,
                                iComp*velTotalDof, iComp*velTotalDof
                              );

            }
        }
        N.GlobalAssemble();
        *A += N;
    }
    chrono.stop();
    M_Displayer.leaderPrintMax( "done in " , chrono.diff() );

}

template<typename Mesh>
void NS_Assembler<Mesh>::applyBoundaryConditions( matrix_type& C,
                                                  matrix_type& Bt,
                                                  bchandler_raw_type& BCh,
                                                  Real time)
{
    // BC manage for the velocity
    if ( !BCh.bdUpdateDone() )
        BCh.bdUpdate( *M_uFESpace.mesh(), M_uFESpace.feBd(), M_uFESpace.dof() );

    bcManageMatrix( C, *M_uFESpace.mesh(), M_uFESpace.dof(), BCh, M_uFESpace.feBd(),
                    M_tgv, time );


    for (bchandler_raw_type::BCBase_Iterator it = BCh.begin(); it != BCh.end(); ++it)
    {
        if (it->type() == Essential)
            bcEssentialManageMatrix( Bt, M_uFESpace.dof(), *it, 0, 0 );
    }
} // applyBoundaryCondition

template<typename Mesh>
void NS_Assembler<Mesh>::applyBoundaryConditions( vector_type & rhs,
                                                  bchandler_raw_type& BCh,
                                                  Real time)
{
    // BC manage for the velocity
    if ( !BCh.bdUpdateDone() )
        BCh.bdUpdate( *M_uFESpace.mesh(), M_uFESpace.feBd(), M_uFESpace.dof() );


    bcManageVector( rhs, *M_uFESpace.mesh(), M_uFESpace.dof(), BCh, M_uFESpace.feBd(), time, M_tgv);
} // applyBoundaryCondition



} // namespace LifeV


#endif //_NS_ASSEMBLER_H_
