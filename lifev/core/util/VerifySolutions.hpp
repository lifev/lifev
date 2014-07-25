//@HEADER
/*
*******************************************************************************

    Copyright (C) 2014 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file VerifySolutions.hpp
    @brief This files contains the description of the VerifySolutions class, which is a helper class for checking results by
           Computing the mean norm and correlation matrix of vectors

    @author Simone Deparis <simone.deparis@epfl.ch>
    @date 12-05-2014
    @mantainer Simone Deparis <simone.deparis@epfl.ch>

 */

#ifndef VERIFYSOLUTIONS_HPP_
#define VERIFYSOLUTIONS_HPP_

#include <Epetra_SerialDenseMatrix.h>
#include <lifev/core/LifeV.hpp>
#include <list>

#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{

/** This class helps the testsuites in verifying that the solutions
 * do not change from one execution to the other.
 * We want to check if the solutions are correct between versions.
 * we want to produce a small number of scalars and check they are always the same.
 *
 * Given: u_0, u_1, ...
 *
 * mean: u = mean u_j
 *
 * correlation matrix: (u_i-u , u_j-u) i,j=0,1,...
 *
 * Store the correlation matrix and (u,u)
 * each entry of the matrix shall be equal to the stored one up to a tolerance.
 *
 * Example of use taken from \ref testsuite/verify_solution/main.cpp
 *
 * \snippet testsuite/verify_solution/main.cpp Example of use of VerifySolutions
 */
class VerifySolutions
{
public:

    /*! @name Constructors and destructor
      */
    //@{

    VerifySolutions();
    virtual ~VerifySolutions();

    //@}

    /*! @name  Methods
     */
    //@{

    //! Add a vector to the list of stored vectors.
    /*!
      *  The vector is copied in a list container.
      * \param newVector vector to be copied and stored
      */
    void PushBack (VectorEpetra const& newVector);

    //! Compute the correlation matrix and norm of the mean.
    /*!
      *  This has to be called after all vectors have been PushBack and befor the Checks
      */
    void ComputeCorrelation();

    //!  Checks that the current Correlation matrix is the same as the one provided by the user.
    /*!
      *  Usually the reference matrix is computed with a solution set which is known to be correct.
      * \param refM precomputed correlation matrix
      * \param tol Tolerance for checking same results
      */
    bool Check ( Epetra_SerialDenseMatrix const& refM, Real tol) const;

    //!  Checks that the current norm of the mean of the vectors is the same as the one provided by the user.
    /*!
      *  Usually the norm of the mean is computed with a solution set which is known to be correct.
      * \param referenceMean precomputed norm of the mean of the vectors
      * \param tol Tolerance for checking same results
      */
    bool Check ( Real referenceMean, Real tol) const;

    //! Print Mean and Correlation matrix as c++ lines ready to be inserted into the test.
    void Print () const;

private:
    Real                     M_NormMean;
    Epetra_SerialDenseMatrix M_CorrelationMatrix;

    std::list<VectorEpetra> M_VectorList;

};

} /* namespace LifeV */
#endif /* VERIFYSOLUTIONS_HPP_ */
