/*
 * VerifySolutions.cpp
 *
 *  Created on: May 10, 2014
 *      Author: Simone Deparis
 */

#include "VerifySolutions.hpp"
#include <vector>
#include <iostream>

namespace LifeV {

VerifySolutions::VerifySolutions() :
                M_NormMean(0),
                M_CorrelationMatrix(),
                M_VectorList(0)
{
}

VerifySolutions::~VerifySolutions()
{
}

bool VerifySolutions::Check( Epetra_SerialDenseMatrix const& refM, Real tol) const
{
    if ( refM.M() != M_CorrelationMatrix.M() || refM.N() != M_CorrelationMatrix.N() )
    {
        return false;
    }

    for (int i(0); i < refM.M(); ++i)
    {
        for (int j(0); j < refM.M(); ++j)
        {
            if ( abs(refM(i,j) - M_CorrelationMatrix(i,j) ) > tol )
            {
                std::cout << "Correlation matrix changed" << std::endl;
                return false;
            }
        }
    }

    std::cout << "Correlation matrix unchanged" << std::endl;

    return true;
}

bool VerifySolutions::Check( Real referenceMean, Real tol) const
{
    if ( abs(referenceMean - M_NormMean) > tol )
    {
        std::cout << "Norm of the mean vector changed" << std::endl;
        return false;
    }

    std::cout << "Norm of the mean vector unchanged" << std::endl;

    return true;
}

void VerifySolutions::ComputeCorrelation()
{
    if (M_VectorList.size() == 0)
    {
        M_NormMean = 0;
        M_CorrelationMatrix.Shape(0,0);
        return;
    }

    // Compute Mean Vector
    std::list<VectorEpetra>::iterator it=M_VectorList.begin();

    VectorEpetra MeanVector(*it);
    ++it;
    for ( ; it != M_VectorList.end(); ++it)
    {
        MeanVector += *it;
    }

    MeanVector /= M_VectorList.size();
    M_NormMean = MeanVector.norm2();

    // Compute u_i-u, inside a vector
    std::vector<VectorEpetra> corrVectorsList(M_VectorList.begin(),M_VectorList.end());
    std::vector<VectorEpetra>::iterator cit;
    for (cit = corrVectorsList.begin(); cit != corrVectorsList.end(); ++cit)
    {
        *cit -= MeanVector;
    }


    // Compute Correlation matrix
    M_CorrelationMatrix.Shape(corrVectorsList.size(),corrVectorsList.size());

    for (int i(0); i < M_VectorList.size(); ++i)
    {
        for (int j(i); j < M_VectorList.size(); ++j)
        {
            M_CorrelationMatrix(i,j) = corrVectorsList[i].dot(corrVectorsList[j]);
            M_CorrelationMatrix(j,i) = M_CorrelationMatrix(i,j);
        }
    }
}

void VerifySolutions::PushBack (VectorEpetra const& newVector)
{

    M_VectorList.push_back(newVector);
}


void VerifySolutions::Print () const
{

    std::cout.precision(15);

    std::cout << "Real referenceMeanNorm = " << M_NormMean << ";" << std::endl;

    std::cout << "Epetra_SerialDenseMatrix refM("
              << M_CorrelationMatrix.M() << "," << M_CorrelationMatrix.N() << ");"
              << std::endl;

    for (int i(0); i < M_CorrelationMatrix.M(); ++i)
    {
        for (int j(0); j < M_CorrelationMatrix.M(); ++j)
        {
            std::cout << "refM(" << i << "," << j << ") = " << M_CorrelationMatrix(i,j) << "; ";
        }
        std::cout << std::endl;
    }


}

} /* namespace LifeV */
