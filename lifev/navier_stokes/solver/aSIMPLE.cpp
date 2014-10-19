/*
 * HighOrderYosida.cpp
 *
 *  Created on: Aug 26, 2010
 *      Author: uvilla
 */

#include <lifev/navier_stokes/solver/aSIMPLE.hpp>

namespace LifeV
{
namespace Operators
{
//===========================================================================//
// Constructors
//===========================================================================//

aSIMPLE::aSIMPLE():
M_label("aSIMPLE"),
M_sigma(0),
M_useTranspose(false)
{
    
}

aSIMPLE::~aSIMPLE()
{
    
}

void aSIMPLE::SetUp(const matrixPtr_Type & Bt,
                    const lumpedMatrixPtr_Type & ml,
                    const matrixPtr_Type & A,
                    const operatorPtr_Type & DL,
                    int p)
{
    // Set the communicator
    //	M_comm = M_DL->getComm_ptr();
    
    // Set the matrices
    M_Bt = Bt;
    M_ml = ml;
    M_A  = A;
    M_DL = DL;
    
    M_p = p;
}
//Set Sigma and matrix A. This function will become handy when solving time adaptive NS
void aSIMPLE::setSigma(Real sigma)
{
    ASSERT_PRE(sigma > 0, "Sigma must be strictly positive");
    M_sigma=sigma;
    // The inverse lumped mass matrix
    M_h.reset(new lumpedMatrix_Type(*M_ml));
    M_h->Scale(M_sigma);
    M_h->Reciprocal(*M_h);
    
}

void aSIMPLE::setA(const matrixPtr_Type & A)
{
    M_A = A;
}

//Change the order of the splitting
void aSIMPLE::setSplittingOrder(int p)
{
    M_p = p;
}

//Solve the Approximate Shur-Complement S = 1/sigma*BHB'
//Contains a trick to deal with rank deficient problems:
//if the problem is only Dirichlet, it returns a pressure with null mean.
void aSIMPLE::solveS(const vector_Type & X, vector_Type & Y) const
{
    
    ASSERT_PRE0(M_sigma!=0, "_sigma has not been set yet\n");
    M_DL->ApplyInverse(X,Y);
    Y.Scale(-M_sigma);
    
}

//Solve Yosida Pressure Correction with order p.
//If the problem is only Dirichlet, it returns a pressure with null mean.
int aSIMPLE::ApplyInverse(const vector_Type &X, vector_Type &Y) const
{
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(),
               "Multivector X and Y don't have the same number of vectors");
    ASSERT_PRE(X.Map().PointSameAs(this->OperatorDomainMap()),
               "Multivector X has a map not compatible with this HighOrderYosida operator");
    ASSERT_PRE(Y.Map().PointSameAs(this->OperatorRangeMap()),
               "MultiVector Y has a map not compatible with this HighOrderYosida operator");
    
    // vector of pressure corrections z[i] = O(dt^i)
    Zdata z(M_p+1);
    // Smart data structure for implementing Yosida corrections
    ZZdata zz(M_p, M_p), bzz(M_p, M_p);
    
    for(Zdata::iterator it=z.begin(); it!= z.end(); ++it)
        (*it).reset(new vector_Type(Y.Map(), Y.NumVectors()) );
    
    for(ZZdata::iterator1 it1=zz.begin1(); it1!=zz.end1(); ++it1)
        for(ZZdata::iterator2 it2=it1.begin(); it2!=it1.end(); ++it2)
            (*it2).reset(new vector_Type(M_Bt->OperatorRangeMap(), Y.NumVectors()));
    
    for(ZZdata::iterator1 it1=bzz.begin1(); it1!=bzz.end1(); ++it1)
        for(ZZdata::iterator2 it2=it1.begin(); it2!=it1.end(); ++it2)
            (*it2).reset(new vector_Type(Y.Map(), Y.NumVectors()));
    
    
    // Temp vectors
    vector_Type tmp1(M_Bt->OperatorRangeMap(), Y.NumVectors());
    vector_Type tmp2(M_Bt->OperatorRangeMap(), Y.NumVectors());
    vector_Type cc(Y.Map(),  Y.NumVectors());
    
    Y.PutScalar(0.0);
    solveS(X, *z[0] );
    Y.Update(1.0, *z[0], 1.0);
    
    for(UInt i=0; i<M_p; ++i)
    {
        //ZZ(i,0) = - h.*(A*h.*B'*z(i))
        EPETRA_CHK_ERR( M_Bt->Multiply(false, *z[i], tmp1) );
        tmp1.Scale(-1.0);
        EPETRA_CHK_ERR( tmp2.Multiply(1.0, tmp1, *M_h, 0.0) );
        EPETRA_CHK_ERR( M_A->Multiply(M_useTranspose, tmp2, tmp1) );
        zz(i,0)->Multiply(1.0, tmp1, *M_h, 0.0);
        //BZZ(i,0) = a.B*ZZ(i,0);
        EPETRA_CHK_ERR( M_Bt->Multiply(true, *zz(i,0), *bzz(i, 0) ) );
        cc = *bzz(i, 0);
        for(UInt j=1; j<i+1;++j)
        {
            //ZZ(i-j,j)  = -h.*(A*ZZ(i-j,j-1));
            EPETRA_CHK_ERR( M_A->Multiply(M_useTranspose, *zz(i-j,j-1), tmp1) );
            tmp1.Scale(-1.0);
            EPETRA_CHK_ERR( zz(i-j,j)->Multiply(1.0, tmp1, *M_h, 0.0) );
            //BZZ(i-j,j) = a.B*ZZ(i-j,j);
            EPETRA_CHK_ERR(  M_Bt->Multiply(true, *zz(i-j,j), *bzz(i-j,j) ) );
            cc.Update(1.0, *bzz(i-j,j), 1.0);
        }
        solveS(cc, *z[i+1]);
        EPETRA_CHK_ERR( Y.Update(1.0, *z[i+1], 1.0) );
    }
    
    
    
    return 0;
}

// show information about the class
void aSIMPLE::showMe(){
    std::cout<<"Dimension u: "<< M_Bt->NumGlobalRows()<<
    ", Dimension p: "<<M_Bt->NumGlobalCols()<<std::endl;
    std::cout<<"Pressure correction order "<<M_p<<std::endl;
}
    
} /* end namespace Operators */
} /*end namespace */
