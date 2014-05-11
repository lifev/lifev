/*
 * VerifySolutions.hpp
 *
 *  Created on: May 10, 2014
 *      Author: Simone Deparis
 */

#ifndef VERIFYSOLUTIONS_HPP_
#define VERIFYSOLUTIONS_HPP_

#include <Epetra_SerialDenseMatrix.h>
#include <lifev/core/LifeV.hpp>
#include <list>

#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV {

/** This class helps the testsuites in verifying that the solutions
 * Do not change from one execution to the other.
 * We want to check if the solutions are correct between versions.
 * we want to produce a small number of scalars and check they are always the same.
 *
 * Given: u_0, u_1, ...
 * mean: u = mean u_j
 * scalars: (u_i-u , u_j-u) i,j=0,1,... (correlation matrix)
 *
 * Store the correlation matrix and (u,u)_X
 * each entry of the matrix shall be equial to the stored one up to a tolerance.
 *
 */
class VerifySolutions
{
public:

    //! Constructor and destructor

    VerifySolutions();
    virtual ~VerifySolutions();

    //! Methods
    bool Check( Epetra_SerialDenseMatrix const& refM, Real tol) const;

    bool Check( Real referenceMean, Real tol) const;

    void ComputeCorrelation();

    void PushBack (VectorEpetra const& newVector);

private:
    Real                     M_NormMean;
    Epetra_SerialDenseMatrix M_CorrelationMatrix;

    std::list<VectorEpetra> M_VectorList;

};

} /* namespace LifeV */
#endif /* VERIFYSOLUTIONS_HPP_ */
