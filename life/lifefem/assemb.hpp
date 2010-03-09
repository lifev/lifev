/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
/* --------------------------------------------------------------------------*
/                                                                            /
/      ...                                                                   /
/                                                                            /
/                                                                            /
/ ASSEMBLY CLASS                                                             /
/                                                                            /
/ #Version 0.1 Experimental:   31/01/2001 Alessandro Veneziani               /
/ #Version 0.2 Added features for vector problems: VBR or global CSR formats /
/                              02/01/2002 Alessandro Veneziani       /=> la specializzazione per il caso simmetrico e' ancora da fare
/ #Version 0.3 6/10/2004 Added features for discontinuous finite elements
/ #Purpose: Container for the assembly routines                              /
/                                                                            /
/                                                                            /
/---------------------------------------------------------------------------*/
#ifndef _ASSEMBLY
#define _ASSEMBLY

#include <life/lifefem/currentFE.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/localDofPattern.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>
#include <vector>

namespace LifeV
{

// For the moment, the selection of the components of the differential operators to be computed is carried out
// inside the differential operator. This is simpler, but not the best: indeed an "if" is processed
// for every integral computation. It would be better to set up a table of simple operators
// and then to plug the correct diff oper component in the by-components loop
// in "assemble". It will be done (I have to understand "how") **** AV January 2002

template <typename Oper, typename DOF, typename RegionMesh,
typename UsrSourceFct>
void
assemble( Oper oper, const RegionMesh& mesh, CurrentFE& fe,
          const DOF& dof, const UsrSourceFct& source_fct, EpetraMatrix<double>& A, EpetraVector& b )
{
    UInt i, ic, jc;
    UInt nc = b.size() / dof.numTotalDof();
    ElemMat elmat( fe.nbNode, nc, nc );
    ElemVec elvec( fe.nbNode, nc );
    //
    for ( i = 1; i <= mesh.numElements(); ++i )
    {
        fe.updateFirstDeriv( mesh.element( i ) ); // as updateFirstDer
        //
        elmat.zero();
        //
        elvec.zero();
        for ( ic = 0; ic < nc; ++ic )
        {
            for ( jc = 0; jc < nc; ++jc )
            {
                compute_mat( elmat, oper, fe, ic, jc ); // compute local matrix
                // the previous line would become:  compute_mat(elmat,oper.comp(ic,jc),fe); ****
                assemb_mat( A, elmat, fe, dof, ic, jc ); // assemble local matrix into global one
            }
#ifdef STABILIZED
            compute_vec_stab( source_fct, elvec, fe, ( int ) ic ); // compute local vector
#else

            compute_vec( source_fct, elvec, fe, ( int ) ic ); // compute local vector
#endif
            // the previous line would become:  compute_vec(source_fct.comp(ic),elvec,fe); ****
            assemb_vec( b, elvec, fe, dof, ic ); // assebmle local vector into global one
        }
    }
}
// version with source term depending on time t
template <typename Oper, typename DOF, typename RegionMesh,
typename UsrSourceFct>
void
assemble( Oper oper, const RegionMesh& mesh, CurrentFE& fe, const DOF& dof,
          const UsrSourceFct& source_fct, EpetraMatrix<double>& A, EpetraVector& b, Real const t )
{
    UInt i, ic, jc;
    UInt nc = b.size() / dof.numTotalDof(); // it should be F.Size() if F is for instance a Physical Vector (like the unknown....)
    ElemMat elmat( fe.nbNode, nc, nc );
    ElemVec elvec( fe.nbNode, nc );
    //
    for ( i = 1; i <= mesh.numElements(); ++i )
    {
        fe.updateFirstDeriv( mesh.element( i ) ); // as updateFirstDer
        //
        elmat.zero();
        //
        elvec.zero();
        for ( ic = 0; ic < nc; ++ic )
        {
            for ( jc = 0; jc < nc; ++jc )
            {
                compute_mat( elmat, oper, fe, ic, jc ); // compute local matrix
                // the previous line would become:  compute_mat(elmat,oper.comp(ic,jc),fe); ****
                assemb_mat( A, elmat, fe, dof, ic, jc ); // assemble local matrix into global one
            }
#ifdef STABILIZED
            compute_vec_stab( source_fct, elvec, fe, t, ic ); // compute local vector
#else

            compute_vec( source_fct, elvec, fe, t, ic ); // compute local vector
#endif
            // the previous line would become:  compute_vec(source_fct.comp(ic),elvec,fe); ****
            assemb_vec( b, elvec, fe, dof, ic ); // assemble local vector into global one
        }
    }
}
//
//////////////////////////////////
//
//! assembling of operator matrix for mixed FE.
//! Alain, 31/01/02.
//! I have added the number of component of the unknown corresponding
//! to FE1 and FE2: respectively nc1 and nc2. Is it possible to know it
//! through the other parameters ?
template <typename Oper, typename DOF1, typename DOF2, typename RegionMesh,
typename UsrSourceFct>
void
assemble_mixed( Oper oper, const RegionMesh& mesh, CurrentFE& fe1,
                CurrentFE& fe2, const DOF1& dof1, const DOF2& dof2,
                UInt const nc1, UInt const nc2, const UsrSourceFct& source_fct,
                EpetraMatrix<double>& A, EpetraVector& b )
{
    UInt i, ic, jc;
    ElemMat elmat( fe1.nbNode, nc1, 0, fe2.nbNode, 0, nc2 );
    ElemVec elvec( fe1.nbNode, nc1 );
    for ( i = 1; i <= mesh.numElements(); ++i )
    {
        fe1.updateFirstDeriv( mesh.element( i ) ); // as updateFirstDer
        fe2.updateFirstDeriv( mesh.element( i ) ); // as updateFirstDer
        //
        elmat.zero();
        //
        elvec.zero();

        for ( ic = 0; ic < nc1; ++ic )
        {
            for ( jc = 0; jc < nc2; ++jc )
            {

                compute_mat_mixed( elmat, oper, fe1, fe2, ic, jc ); // compute local matrix

                assemb_mat_mixed( A, elmat, fe1, fe2, dof1, dof2, ic, jc ); // assemble local matrix into global one
            }
            compute_vec( source_fct, elvec, fe1, ic ); // compute local vector

            assemb_vec( b, elvec, fe1, dof1, ic ); // assemble local vector into global one
        }
    }
}
// version with source term depending on time
template <typename Oper, typename DOF1, typename DOF2, typename RegionMesh,
typename UsrSourceFct>
void
assemble_mixed( Oper oper, const RegionMesh& mesh, CurrentFE& fe1,
                CurrentFE& fe2, const DOF1& dof1, const DOF2& dof2,
                UInt const nc1, UInt const nc2, const UsrSourceFct& source_fct,
                EpetraMatrix<double>& A, EpetraVector& b, Real const t )
{
    UInt i, ic, jc;
    ElemMat elmat( fe1.nbNode, nc1, 0, fe2.nbNode, 0, nc2 );
    ElemVec elvec( fe1.nbNode, nc1 );
    for ( i = 1; i <= mesh.numElements(); ++i )
    {
        fe1.updateFirstDeriv( mesh.element( i ) ); // as updateFirstDer
        fe2.updateFirstDeriv( mesh.element( i ) ); // as updateFirstDer
        //
        elmat.zero();
        //
        elvec.zero();

        for ( ic = 0; ic < nc1; ++ic )
        {
            for ( jc = 0; jc < nc2; ++jc )
            {

                compute_mat_mixed( elmat, oper, fe1, fe2, ic, jc ); // compute local matrix

                assemb_mat_mixed( A, elmat, fe1, fe2, dof1, dof2, ic, jc ); // assemble local matrix into global one
            }
            compute_vec( source_fct, elvec, fe1, t, ic ); // compute local vector

            assemb_vec( b, elvec, fe1, dof1, ic ); // assemble local vector into global one
        }
    }
}
///////////////////
// version with source term depending on time t
template <typename Oper, typename DOF, typename RegionMesh,
typename UsrSourceFct>
void
assemble_symm( Oper oper, const RegionMesh& mesh, CurrentFE& fe, const DOF& dof,
               const UsrSourceFct& source_fct, EpetraMatrix<double>& A, EpetraVector& b, Real const t )
{
    UInt i, ic, jc;
    UInt nc = b.size() / dof.numTotalDof(); // it should be F.Size() if F is for instance a Physical Vector (like the unknown....)
    ElemMat elmat( fe.nbNode, nc, nc );
    ElemVec elvec( fe.nbNode, nc );
    //
    for ( i = 1; i <= mesh.numElements(); ++i )
    {
        fe.updateFirstDeriv( mesh.element( i ) ); // as updateFirstDer
        //
        elmat.zero();
        //
        elvec.zero();
        for ( ic = 0; ic < nc; ++ic )
        {
            for ( jc = 0; jc < nc; ++jc )
            {
                compute_mat_symm( elmat, oper, fe, ic, jc ); // compute local matrix expoiting symmetry of the operator
                // the previous line would become:  compute_mat(elmat,oper.comp(ic,jc),fe); ****
                assemb_mat( A, elmat, fe, dof, ic, jc ); // assemble local matrix into global one
            }
            compute_vec( source_fct, elvec, fe, t, ic ); // compute local vector
            // the previous line would become:  compute_vec(source_fct.comp(ic),elvec,fe); ****
            assemb_vec( b, elvec, fe, dof, ic ); // assebmle local vector into global one
        }
    }
}
//
//////////////////////////////////
///////////////////
// BLOCK DIAGONAL OPERATOR ASSEMBLING
// It assembles the same scalar operator over the (3) diagonal blocks
template <typename Oper, typename DOF, typename RegionMesh,
typename UsrSourceFct>
void
assemble_symm_block_diagonal( Oper oper, const RegionMesh& mesh, CurrentFE& fe, const DOF& dof,
                              const UsrSourceFct& source_fct, EpetraMatrix<double>& A, EpetraVector& b, Real const t )
{
    UInt i, ic = 0;
    UInt nc = b.size() / dof.numTotalDof(); // it should be F.Size() if F is for instance a Physical Vector (like the unknown....)
    ElemMat elmat( fe.nbNode, nc, nc );
    ElemVec elvec( fe.nbNode, nc );
    //
    for ( i = 1; i <= mesh.numElements(); ++i )
    {
        fe.updateFirstDeriv( mesh.element( i ) ); // as updateFirstDer
        //
        elmat.zero();
        //
        elvec.zero();
        compute_mat_symm( elmat, oper, fe ); // compute local matrix expoiting symmetry of the operator

        ic = 0;
        assemb_mat( A, elmat, fe, dof, ic, ic ); // assemble local matrix into global one
        compute_vec( source_fct, elvec, fe, t, ic ); // compute local vector
        assemb_vec( b, elvec, fe, dof, ic ); // assebmle local vector into global one

        for ( ic = 1; ic < nc; ++ic )
        {
            elmat.block( ic, ic ) = elmat.block( 0, 0 );
            assemb_mat( A, elmat, fe, dof, ic, ic ); // assemble local matrix into global one

            compute_vec( source_fct, elvec, fe, t, ic ); // compute local vector
            // the previous line would become:  compute_vec(source_fct.comp(ic),elvec,fe); ****
            assemb_vec( b, elvec, fe, dof, ic ); // assebmle local vector into global one
        }
    }
}


////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////
//
template <typename Oper>
void compute_mat( ElemMat& elmat, Oper& oper,
                  const CurrentFE& fe, int iblock = 0, int jblock = 0 )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt ig;
    unsigned int i, j;
    double s;
    Real x = 0, y = 0, z = 0;
    for ( i = 0;i < ( UInt ) fe.nbNode;i++ )
    {
        for ( j = 0;j < ( UInt ) fe.nbNode;j++ )
        {
            s = 0;
            for ( ig = 0;ig < ( UInt ) fe.nbQuadPt;ig++ )
            {
#ifdef SPACE_DEP_OPERATOR
                fe.coorQuadPt( x, y, z, ig );
#endif

                s += oper( i, j, ig, x, y, z, iblock, jblock ) * fe.weightDet( ig );
                // the previous line would become:  oper(i,j,ig,x,y,z)*fe.weightDet(ig); ****
            }
            mat( ( int ) i, ( int ) j ) += s;
        }
    }
}



////////////////////////////////////////////////////////////////////////////////
//
//! compute local matrix element for mixed FE
//Alain 1/02/02.
template <typename Oper>
void compute_mat_mixed( ElemMat& elmat, Oper& oper,
                        const CurrentFE& fe1, const CurrentFE& fe2, int iblock = 0, int jblock = 0 )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int ig;
    int i, j;
    double s;
    Real x, y, z;

    for ( i = 0;i < fe1.nbNode;i++ )
    {
        for ( j = 0;j < fe2.nbNode;j++ )
        {
            s = 0;
            // Alain, 08/02/02:
            // REMARK : the choice of quadrature points depends on the definition
            //          of the type FE1 and FE2 in mainAZ.h.
            //          Therefore the same quadrature rule must be chosen
            //          for both FE1 and FE2.
            for ( ig = 0;ig < fe1.nbQuadPt;ig++ )
            {
                fe1.coorQuadPt( x, y, z, ig );
                s += oper( i, j, ig, x, y, z, iblock, jblock ) * fe1.weightDet( ig );
            }
            mat( ( int ) i, ( int ) j ) += s;
        }
    }
}
//
template <typename Oper>
void compute_mat_symm( ElemMat& elmat, Oper& oper,
                       const CurrentFE& fe, int iblock = 0, int jblock = 0 )
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt ig;
    unsigned int i, iloc, jloc;
    double s;
    Real x = 0, y = 0, z = 0;
    //
    // diagonal
    //
    for ( i = 0;i < ( UInt ) fe.nbDiag;i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < ( UInt ) fe.nbQuadPt;ig++ )
        {
#ifdef SPACE_DEP_OPERATOR
            fe.coorQuadPt( x, y, z, ig );
#endif

            s += oper( iloc, iloc, ig, x, y, z, iblock, jblock ) * fe.weightDet( ig );
        }
        mat( ( int ) iloc, ( int ) iloc ) += s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag;i < ( UInt ) fe.nbDiag + fe.nbUpper;i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < ( UInt ) fe.nbQuadPt;ig++ )
        {
#ifdef SPACE_DEP_OPERATOR
            fe.coorQuadPt( x, y, z, ig );
#endif

            s += oper( iloc, jloc, ig, x, y, z, iblock, jblock ) * fe.weightDet( ig );
        }

        mat( ( int ) iloc, ( int ) jloc ) += s;
        mat( ( int ) jloc, ( int ) iloc ) += s;
    }
}

//
/////////////////////////////////
//




template <typename DOF>
void
assembleVector( EpetraVector&    vec,
                ElemVec&         elvec,
                const CurrentFE& fe,
                const DOF&       dof,
                int              iblock,
                int              ipos = 0 )

{
//    elmat.showMe();
    ElemVec::vector_view vecView = elvec.block( iblock );

//    UInt totdof = dof.numTotalDof();

    int i;
    UInt ig;

    UInt eleID = fe.currentLocalId();

    bool verbose = false;

    for ( i = 0 ; i < fe.nbNode ; i++ )
    {
        ig = dof.localToGlobal( eleID, i + 1 ) + ipos;/*+ iblock*totdof*/  // damned 1-base vs 0-base !
            vec.sumIntoGlobalValues( ig, vecView( i ) );
    }

    if (verbose)
        std::cout << "ok." << std::endl;
}




template <typename DOF>
void
assembleMatrix( EpetraMatrix<double>& M,
                ElemMat&              elmat,
                const CurrentFE&      fe,
                const DOF&            dof,
                int                   iblock,
                int                   jblock,
                int                   ipos,
                int                   jpos)

{
    return assembleMatrix( M, elmat, fe, fe, dof, dof,
                           iblock, jblock, ipos, jpos);

    /*
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );

    int i, j, k;
    UInt ig, jg;

    UInt eleID = fe.currentLocalId();

//    std::cout << "fe.nbPattern = " << fe.nbPattern << std::endl;

    for ( k = 0 ; k < fe.nbPattern ; k++ )
    {
        i = fe.patternFirst ( k );
        j = fe.patternSecond( k );

        if (mat(i,j) != 0.)
            {
                ig = dof.localToGlobal( eleID, i + 1 ) - 1 + ipos ; //iblock*totdof1;  // damned 1-base vs 0-base !
                jg = dof.localToGlobal( eleID, j + 1 ) - 1 + jpos ; //jblock*totdof2;  // damned 1-base vs 0-base !
                M.set_mat_inc( ig, jg, mat( i, j ) );
            }
    }
    */

}

template <typename DOF1, typename DOF2, typename SUBMAT>
void
assembleMatrix( EpetraMatrix<double>& M,
                UInt const&           eleID1,
                UInt const&           eleID2,
                SUBMAT&              mat,
                const CurrentFE& fe1,
                const CurrentFE& fe2,
                const DOF1&      dof1,
                const DOF2&      dof2,
                int              ipos,
                int              jpos,
                bool             /*verbose = false*/)

{
    int i, j, k1, k2;

    std::vector<int> ilist(fe1.nbNode);
    std::vector<int> jlist(fe2.nbNode);

    double* matPtr[fe2.nbNode];


    for ( k1 = 0 ; k1 < fe1.nbNode ; k1++ )
    {
        i = k1;//fe1.patternFirst( k1 );
        ilist[k1] = dof1.localToGlobal( eleID1, i + 1 ) - 1 + ipos ; //iblock*totdof1;  // damned 1-base vs 0-base !
    }

    for ( k2 = 0 ; k2 < fe2.nbNode ; k2++ )
    {
        j = k2;//fe2.patternFirst( k2 );
        jlist[k2]  = dof2.localToGlobal( eleID2, j + 1 ) - 1 + jpos ; //iblock*totdof1;  // damned 1-base vs 0-base !
        matPtr[k2] = &(mat(0,k2));
    }

    // coded a version to insert the little matrix directly.
    // This needs that mat has the shape checked by the following line:
    assert(mat.indexij( int (1), int(0) ) == 1);

    M.set_mat_inc( fe1.nbNode, fe2.nbNode, ilist, jlist, matPtr, Epetra_FECrsMatrix::COLUMN_MAJOR );

#ifdef ONLY_FOR_DEBUGGING

    for ( k1 = 0 ; k1 < fe1.nbNode ; k1++ )
    {
        i = fe1.patternFirst( k1 );
        for ( k2 = 0 ; k2 < fe2.nbNode ; k2++ )
        {
            j  = fe2.patternSecond( k2 );
            //            if (mat(i,j) != 0.)
                {
//                     ig = dof1.localToGlobal( eleID1, i + 1 ) - 1 + ipos ; //iblock*totdof1;  // damned 1-base vs 0-base !
//                     jg = dof2.localToGlobal( eleID2, j + 1 ) - 1 + jpos ; //jblock*totdof2;  // damned 1-base vs 0-base !
                    matPtr[k2] = &(mat(0,j));

                    assert(matPtr[k2][k1] ==  mat( i, j ));

//                     std::cout << "ig, jg, mat( i, j ) = "
//                               << ig << " " <<  jg << " "<<  mat( i, j ) << " "
//                               << matPtr[k1] - mat( i, j )
//                               << std::endl
//                               << "                      "
//                               << ilist[k1] << " " <<  jg << " "<<  matPtr[k1]
//                               << std::endl;

                    //M.set_mat_inc( ig, jg, mat( i, j ) );
                }
        }
    }
#endif
}

template <typename DOF1, typename DOF2>
void
assembleMatrix( EpetraMatrix<double>&          M,
                ElemMat&         elmat,
                const CurrentFE& fe1,
                const CurrentFE& fe2,
                const DOF1&      dof1,
                const DOF2&      dof2,
                int              iblock,
                int              jblock,
                int              ipos,
                int              jpos,
                bool             verbose = false)

{
//    elmat.showMe();
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );

    UInt eleID1 = fe1.currentLocalId();
    UInt eleID2 = fe2.currentLocalId();

    assembleMatrix( M, eleID1, eleID2,
                    mat, fe1, fe2,
                    dof1,  dof2, ipos, jpos, verbose);

    return;

}

template <typename DOF1, typename DOF2>
void
assembleTransposeMatrix( EpetraMatrix<double>& M,
                         Real                  val,
                         ElemMat&              elmat,
                         const CurrentFE&      fe1,
                         const CurrentFE&      fe2,
                         const DOF1&           dof1,
                         const DOF2&           dof2,
                         int                   iblock,
                         int                   jblock,
                         int                   ipos ,
                         int                   jpos )

{
    ElemMat::matrix_type mat(elmat.block( jblock, iblock ));
    mat *= val;

    int i, j, k1, k2;

    UInt eleID1 = fe1.currentLocalId();
    UInt eleID2 = fe2.currentLocalId();

    std::vector<int> ilist(fe1.nbNode);
    std::vector<int> jlist(fe2.nbNode);

    double* matPtr[fe1.nbNode];


    for ( k1 = 0 ; k1 < fe1.nbNode ; k1++ )
    {
        i =  k1;
        ilist[k1] = dof1.localToGlobal( eleID1, i + 1 ) - 1 + ipos ; //iblock*totdof1;  // damned 1-base vs 0-base !
        matPtr[k1] = &(mat(0,i));
    }

    for ( k2 = 0 ; k2 < fe2.nbNode ; k2++ )
    {
        j = k2;
        jlist[k2]  = dof2.localToGlobal( eleID2, j + 1 ) - 1 + jpos ; //iblock*totdof1;  // damned 1-base vs 0-base !
    }

    // coded a version to insert the little matrix directly.
    // This needs that mat has the shape checked by the following line:
    assert(mat.indexij( int (1), int(0) ) == 1);

    M.set_mat_inc( fe1.nbNode, fe2.nbNode, ilist, jlist, matPtr, Epetra_FECrsMatrix::ROW_MAJOR );

}



//
/////////////////////////////////
//
template < typename DOF>
void
assemb_mat( EpetraMatrix<double>& M, ElemMat& elmat, const CurrentFE& fe, const DOF& dof, int iblock = 0, int jblock = 0 , int offset =0)
{

    UInt totdof = dof.numTotalDof();
    return assembleMatrix(  M, elmat, fe, dof, iblock, jblock, iblock * totdof, jblock * totdof);

    // Warning: in assemble Matrix
    // UInt eleId = fe.currentId();
    //is replaced by
    // UInt eleID = fe.currentLocalId();


    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int i, j, k;
    UInt ig, jg;
    UInt eleId = fe.currentLocalId();
    for ( k = 0 ; k < fe.nbPattern ; k++ )
    {
        i = fe.patternFirst( k );
        j = fe.patternSecond( k );
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof + offset;  // damned 1-base vs 0-base !
        jg = dof.localToGlobal( eleId, j + 1 ) - 1 + jblock * totdof + offset;  // damned 1-base vs 0-base !
        M.set_mat_inc( ig, jg, mat( i, j ) );
//        std::cout << ig << " " << jg << " " << mat(i, j) << std::endl;

    }
}


//
/////////////////////////////////
//
// Miguel 01/2004: matrix assembling of interior penalty terms
//
template < typename DOF>
void
assemb_mat( EpetraMatrix<double>& M, ElemMat& elmat, const CurrentFE& fe1, const CurrentFE& fe2, const DOF& dof, int iblock = 0, int jblock = 0 )
{
    UInt totdof = dof.numTotalDof();
    return assembleMatrix(  M, elmat, fe1, fe2, dof, dof,
                            iblock, jblock, iblock * totdof, jblock * totdof);

    // Warning: in assemble Matrix
    // UInt eleId = fe.currentId();
    //is replaced by
    // UInt eleID = fe.currentLocalId();

    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int i, j, k;
    UInt ig, jg;
    UInt eleId1 = fe1.currentId();
    UInt eleId2 = fe2.currentId();
    for ( k = 0 ; k < fe1.nbPattern ; k++ )
    {
        i = fe1.patternFirst( k );
        j = fe2.patternSecond( k );
        ig = dof.localToGlobal( eleId1, i + 1 ) - 1 + iblock * totdof;
        jg = dof.localToGlobal( eleId2, j + 1 ) - 1 + jblock * totdof;
        M.set_mat_inc( ig, jg, mat( i, j ) );
    }
}

//
/////////////////////////////////
//
//! Added by V. Martin 09/2002 (slightly different from the previous function...)
//!    Works with a reference hybrid element, and a given number of the current geo element.
template <typename DOF>
void assemb_mat( EpetraMatrix<double>& M, ElemMat& elmat, const LocalDofPattern& fe, const DOF& dof,
                 const UInt feId, int iblock = 0, int jblock = 0 )
{


    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt totdof = dof.numTotalDof();
    return assembleMatrix( M, feId, feId, mat, fe, fe, dof, dof, iblock * totdof, jblock * totdof);

    //  if(elmat.nBlockRow()!=1 || elmat.nBlockCol() != 1){
    //    std::cout << "assemble for vector elem mat not yet implemented\n";
    //    exit(1);
    //  }
    int i, j, k;
    UInt ig, jg;
    UInt eleId = feId; //! direct use of the fe identity number. (different from the other assemb_mat)

    for ( k = 0 ; k < fe.nbPattern() ; k++ )
    { //! instead of currentFE::nbPattern  ...
        i  = fe.patternFirst( k );
        j  = fe.patternSecond( k );
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof;  // damned 1-base vs 0-base !
        jg = dof.localToGlobal( eleId, j + 1 ) - 1 + jblock * totdof;  // damned 1-base vs 0-base !
        M.set_mat_inc( ig, jg, mat( i, j ) );
    }
}


////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////
//
//! assembling of the elementary matrix with a symmetric form: use of the LOWER part
//!
//! How to make it more general? For the upper?? (without too much efficiency loss)
//! CAUTION: TESTED ONLY WITH THE STANDARD PATTERN...
//! V. Martin.
template <typename DOF>
void assemb_mat_symm_lower( EpetraMatrix<double>& M, ElemMat& elmat, const LocalDofPattern& fe,
                            const DOF& dof, const UInt feId, int iblock = 0, int jblock = 0 )
{

    assert(false); // Not tested, not coded

    //  if(elmat.nBlockRow()!=1 || elmat.nBlockCol() != 1){
    //    std::cout << "assemble for vector elem mat not yet implemented\n";
    //    exit(1);
    //  }
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt totdof = dof.numTotalDof();
    int i, j, k;
    UInt ig, jg;
    UInt eleId = feId; //! direct use of the fe identity number.

    //  std::cout << " Lower case " << std::endl;
    //  std::cout << mat << std::endl;

    //
    // diagonal
    //
    for ( k = 0 ; k < fe.nbDiag() ; k ++ )
    {  //! instead of currentFE::nbDiag
        i = fe.patternFirst( k );
        j = fe.patternSecond( k );
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof;  // damned 1-base vs 0-base !
        jg = dof.localToGlobal( eleId, j + 1 ) - 1 + jblock * totdof;  // damned 1-base vs 0-base !
        M.set_mat_inc( ig, jg, mat( i, j ) );
        //    std::cout << mat(i,j)  << "\t";
    }
    //
    // extra diagonal : Lower part.
    //
    //  std::cout << "\nextra diag \nnbupper " <<  fe.nbUpper()  << std::endl;
    for ( k = fe.nbDiag() ; k < fe.nbDiag() + fe.nbUpper() ; k ++ )
    {   //! instead of currentFE::nbUpper
        j = fe.patternFirst( k );   //! transpose of upper -> lower CAUTION: TESTED ONLY WITH THE STANDARD PATTERN...
        i = fe.patternSecond( k );  //! transpose of upper -> lower
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof;  // damned 1-base vs 0-base !
        jg = dof.localToGlobal( eleId, j + 1 ) - 1 + jblock * totdof;  // damned 1-base vs 0-base !
        M.set_mat_inc( ig, jg, mat( i, j ) );
        M.set_mat_inc( jg, ig, mat( i, j ) );
        //    std::cout << mat(i,j)  << "\t";
    }
}

//
/////////////////////////////////
//
//! assembling of the elementary matrix with a symmetric form: use of the UPPER part
//!
//! How to make it more general? For the upper?? (without too much efficiency loss)
//! CAUTION: TESTED ONLY WITH THE STANDARD PATTERN...
//! V. Martin.
template <typename DOF>
void assemb_mat_symm_upper( EpetraMatrix<double>& M, ElemMat& elmat, const LocalDofPattern& fe,
                            const DOF& dof, const UInt feId, int iblock = 0, int jblock = 0 )
{

    assert(false); // Not tested, not coded

    //  if(elmat.nBlockRow()!=1 || elmat.nBlockCol() != 1){
    //    std::cout << "assemble for vector elem mat not yet implemented\n";
    //    exit(1);
    //  }
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    UInt totdof = dof.numTotalDof();
    int i, j, k;
    UInt ig, jg;
    UInt eleId = feId; //! direct use of the fe identity number.

    //  std::cout << " Upper case " << std::endl;

    //
    // diagonal
    //
    for ( k = 0 ; k < fe.nbDiag() ; k ++ )
    {  //! instead of currentFE::nbDiag
        i = fe.patternFirst( k );
        j = fe.patternSecond( k );
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof;  // damned 1-base vs 0-base !
        jg = dof.localToGlobal( eleId, j + 1 ) - 1 + jblock * totdof;  // damned 1-base vs 0-base !
        M.set_mat_inc( ig, jg, mat( i, j ) );
    }
    //
    // extra diagonal : Upper part.
    //
    for ( k = fe.nbDiag() ; k < fe.nbDiag() + fe.nbUpper() ; k ++ )
    {   //! instead of currentFE::nbUpper
        i = fe.patternFirst( k );   //! upper. CAUTION: TESTED ONLY WITH THE STANDARD PATTERN...
        j = fe.patternSecond( k );  //! upper
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof;  // damned 1-base vs 0-base !
        jg = dof.localToGlobal( eleId, j + 1 ) - 1 + jblock * totdof;  // damned 1-base vs 0-base !
        M.set_mat_inc( ig, jg, mat( i, j ) );
        M.set_mat_inc( jg, ig, mat( i, j ) );
    }
}

//
//////////////////
//
//! assembling of the elementary matrix into the global one for
//! mixed FE
//Alain 1/02/02.
template <typename DOF1, typename DOF2>
void
assemb_mat_mixed( EpetraMatrix<double>& M, ElemMat& elmat, const CurrentFE& fe1, const CurrentFE& fe2,
                  const DOF1& dof1, const DOF2& dof2, int iblock = 0, int jblock = 0 )
{

    UInt totdof1 = dof1.numTotalDof();
    UInt totdof2 = dof2.numTotalDof();
    return assembleMatrix( M, elmat, fe1, fe2,
                           dof1, dof2, iblock, jblock,
                           iblock * totdof1, jblock * totdof2);


    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int i, j, k1, k2;
    UInt ig, jg;
    UInt eleID1 = fe1.currentId();
    UInt eleID2 = fe2.currentId();

    for ( k1 = 0 ; k1 < fe1.nbNode ; k1++ )
    {
        i = fe1.patternFirst( k1 );
        for ( k2 = 0 ; k2 < fe2.nbNode ; k2++ )
        {
            j = fe2.patternSecond( k2 );
            ig = dof1.localToGlobal( eleID1, i + 1 ) - 1 + iblock * totdof1;  // damned 1-base vs 0-base !
            jg = dof2.localToGlobal( eleID2, j + 1 ) - 1 + jblock * totdof2;  // damned 1-base vs 0-base !
            M.set_mat_inc( ig, jg, mat( i, j ) );
        }
    }
}

/*!
  Assemble the transposed of an element matrix (multiplied by a
  constant "mulfac") JFG 03/09/02.
*/
template <typename DOF1, typename DOF2>
void
assemb_tr_mat_mixed( Real mulfac, EpetraMatrix<double>& M, ElemMat& elmat,
                     const CurrentFE& fe1, const CurrentFE& fe2,
                     const DOF1& dof1, const DOF2& dof2,
                     int iblock = 0, int jblock = 0 )
{
    UInt totdof1 = dof1.numTotalDof();
    UInt totdof2 = dof2.numTotalDof();
    return assembleTransposeMatrix( M, mulfac, elmat, fe1, fe2,
                                    dof1, dof2, iblock, jblock,
                                    iblock * totdof1, jblock * totdof2);


    ElemMat::matrix_view mat = elmat.block( jblock, iblock );
    int i, j, k1, k2;
    UInt ig, jg;
    UInt eleID1 = fe1.currentId();
    UInt eleID2 = fe2.currentId();

    for ( k1 = 0 ; k1 < fe1.nbNode ; k1++ )
    {
        i = fe1.patternFirst( k1 );
        for ( k2 = 0 ; k2 < fe2.nbNode ; k2++ )
        {
            j = fe2.patternSecond( k2 );
            ig = dof1.localToGlobal( eleID1, i + 1 ) - 1 + iblock * totdof1;  // damned 1-base vs 0-base !
            jg = dof2.localToGlobal( eleID2, j + 1 ) - 1 + jblock * totdof2;  // damned 1-base vs 0-base !
            M.set_mat_inc( ig, jg, mulfac * mat( j, i ) ); // note the mat(j,i) which is the
            // only difference with assemb_mat_mixed
        }
    }
}
//
//////////////////
//
void compute_vec( Real constant, ElemVec& elvec, const CurrentFE& fe, int iblock = 0 );
//////////////////
//
// Suggestion by Luca F.
/*
  Why dont we let CurrentFE to store the time as well?
  we may choose then a common layout for the user functions that always has
  the time as entry point.
 */
//
//
//////////////////
//
template <typename UsrFct>
void compute_vec( const UsrFct& fct, ElemVec& elvec, const CurrentFE& fe, int iblock = 0 )
{
    int i, ig;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s, x, y, z;
    for ( i = 0;i < fe.nbNode;i++ )
    {
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            fe.coorQuadPt( x, y, z, ig );
            s += fe.phi( i, ig ) * fct( x, y, z, iblock + 1 ) * fe.weightDet( ig );
            //the previous line would become: s += fe.phi(i,ig) * fct(x,y,z) * fe.weightDet(ig); ****
        }
        vec( i ) += s;
    }
}
//
template <typename UsrFct>
void compute_vec( const UsrFct& fct, ElemVec& elvec, const CurrentFE& fe )
{
    compute_vec( fct, elvec, fe, 0 );
}

// ! versione per il caso stabilizzato (o comunque son termine forzante wrappato)
template <typename OperFct>
void compute_vec_stab( OperFct& fct, ElemVec& elvec, const CurrentFE& fe,
                       int iblock )
{
    int i, ig;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s, x, y, z;
    for ( i = 0;i < fe.nbNode;i++ )
    {
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            fe.coorQuadPt( x, y, z, ig );
            s += fct( i, ig, x, y, z, iblock + 1 ) * fe.weightDet( ig );
        }
        vec( i ) += s;
    }
}

// version with source term depending on time
template <typename UsrFct>
void compute_vec( const UsrFct& fct, ElemVec& elvec, const CurrentFE& fe, const Real& t, int iblock = 0 )
{
    int i, ig;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s, x, y, z;
    for ( i = 0;i < fe.nbNode;i++ )
    {
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            x = fe.quadPt( ig, 0 );
            y = fe.quadPt( ig, 1 );
            z = fe.quadPt( ig, 2 );
            s += fe.phi( i, ig ) * fct( t, x, y, z, iblock + 1 ) * fe.weightDet( ig );
            //the previous line would become: s += fe.phi(i,ig) * fct(x,y,z,t) * fe.weightDet(ig); ****
        }
        vec( i ) += s;
    }
}
// ! Stabilized case with time dependence
template <typename OperFct>
void compute_vec_stab( OperFct& fct, ElemVec& elvec, const CurrentFE& fe, Real t,
                       int iblock )
{
    int i, ig;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s, x, y, z;
    for ( i = 0;i < fe.nbNode;i++ )
    {
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
        {
            fe.coorQuadPt( x, y, z, ig );
            s += fct( i, ig, x, y, z, t, iblock + 1 ) * fe.weightDet( ig );
            //the previous line would become: s += fe.phi(i,ig) * fct(x,y,z) * fe.weightDet(ig); ****
        }
        vec( i ) += s;
    }
}
//


template <typename DOF, typename ElemVec>
void
assemb_vec( EpetraVector& V, ElemVec& elvec, const CurrentFE& fe, const DOF& dof, int iblock=0 )
{
    UInt totdof = dof.numTotalDof();
    typename ElemVec::vector_view vec = elvec.block( iblock );
    int i;
    UInt ig;
    UInt eleId = fe.currentId();
    for ( i = 0 ; i < fe.nbNode ; i++ )
    {
        ig = dof.localToGlobal( eleId, i + 1 ) - 1 + iblock * totdof;
        V[ ig ] += vec( i );
    }
}
///

////////////////////////////////////////////////////////////////////////////////

///
//////////////////
/// V. Martin  09/2002
//! version of assemb_vec that works with a LocalDofPattern (and also a RefHybridFE)...
template <typename DOF, typename ElemVec>
void
assemb_vec( EpetraVector& V, ElemVec& elvec, const LocalDofPattern& fe, const DOF& dof,
            const UInt feId, int iblock )
{
    //  if(elvec.nBlockRow()!=1){
    //    std::cout << "assemble for vector elem vec not yet implemented\n";
    //    exit(1);
    //  }
    UInt totdof = dof.numTotalDof();
    typename ElemVec::vector_view vec = elvec.block( iblock );
    int i;
    //  std::cout << "in assemb_vec" << std::endl;
    UInt ig;
    //  UInt eleId = feId;  //simplify.
    for ( i = 0 ; i < fe.nbLocalDof() ; i++ )
    {    //! instead of CurrentFE::nbNode
        ig = dof.localToGlobal( feId, i + 1 ) - 1 + iblock * totdof;
        //    std::cout << "i= " << i << std::endl;
        //    std::cout << "ig= " << ig << std::endl;
        V[ ig ] += vec( i );
    }
}
///

///
//////////////////
/// V. Martin  09/2002
//! \function extract_vec
//! \brief from a global vector, extract a vector corresponding to the element feId.
//! works with a LocalDofPattern (and also a RefHybridFE)...
template <typename DOF, typename ElemVec>
void
extract_vec( EpetraVector& V, ElemVec& elvec, const LocalDofPattern& fe, const DOF& dof,
             const UInt feId, int iblock )
{
    //  if(elvec.nBlockRow()!=1){
    //    std::cout << "assemble for vector elem vec not yet implemented\n";
    //    exit(1);
    //  }

    UInt totdof = dof.numTotalDof();
    typename ElemVec::vector_view vec = elvec.block( iblock );
    int i;
    //  std::cout << "in assemb_vec" << std::endl;
    UInt ig;
    for ( i = 0 ; i < fe.nbLocalDof() ; i++ )
    {
        ig = dof.localToGlobal( feId, i + 1 ) - 1 + iblock * totdof;
        //    std::cout << "i= " << i << std::endl;
        //    std::cout << "ig= " << ig << std::endl;
        vec( i ) = V[ ig ];
    }
}
}

#endif
