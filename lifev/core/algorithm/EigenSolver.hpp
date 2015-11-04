/* -*- mode: c++ -*-*/
//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

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

/**
   \file
   \author Paolo Crosetto <paolo.crosetto@epfl.ch>
   \date 2009-03-25
   \brief Class handling an eigensolver
 */


// #ifndef HAVE_TRILINOS_ANASAZI
// #warning: you should use ANASAZI
// #else

// #endif

#ifdef  HAVE_TRILINOS_ANASAZI

#ifndef EIGENSOLVER_HPP
#define EIGENSOLVER_HPP


#include <Epetra_MultiVector.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>

#include <Teuchos_RefCountPtrDecl.hpp>

#include <AnasaziBlockDavidson.hpp>
#include <AnasaziLOBPCG.hpp>
#include <AnasaziBasicOutputManager.hpp>
#include <AnasaziBasicSort.hpp>
#include <AnasaziConfigDefs.hpp>
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <AnasaziLOBPCGSolMgr.hpp>
#include <AnasaziEpetraAdapter.hpp>


#include <lifev/core/LifeV.hpp>

#include <lifev/core/filter/GetPot.hpp>

namespace LifeV
{

class UNDEF_EIGENSOLVER_EXCEPTION;

// namespace Epetra
// {

//  template <typename DataType, typename solver_Type, typename vector_Type>
class EigenSolver
{
    /**
       @class
     * class handling the Anasazi BlockKrylovSchur eigensolver. It works for any Epetra_Operator:
     * in particular in Lifev it shuld be used with the ComposedOperator, that inherits direcly
     * from the Epetra_Operator. Since the class ComposedOperator is templated it can contain
     *  both matrices or preconditioners. See Monolithic.cpp for an example using the ComposedOperator.
     */

public:

    //!@name Public Types
    //@{
    typedef double data_Type;
    typedef Epetra_Operator solver_Type;
    typedef Epetra_MultiVector vector_Type;

    typedef Anasazi::BasicEigenproblem<data_Type, vector_Type, solver_Type>                     eigenpb_Type;
    typedef Teuchos::RCP<Anasazi::BasicEigenproblem<data_Type, vector_Type, solver_Type> >      eigenpbPtr_Type;
    typedef Anasazi::BlockKrylovSchurSolMgr <data_Type, vector_Type, solver_Type>                 eigensolver_Type;
    typedef std::shared_ptr<eigensolver_Type>                                  eigensolverPtr_Type;

    //@}
    //!@name Constructor and Destructor
    //@{

    EigenSolver (std::shared_ptr<solver_Type> const matrix, Epetra_BlockMap const& block_map, long unsigned int numvec);

    virtual ~EigenSolver()
    {}

    //@}
    //!@name Public Methods
    //@{

    /**sets the parameters for the eigenproblem:
     */
    void setDataFromGetPot ( GetPot const& dataFile, std::string const& section/*="eigensolver"*/ );

    /** set to true if the operator is symmetric*/
    void setHermitian (bool flag)
    {
        MyProblem->setHermitian (flag);
    }

    /** fills the input vectors with the real and imaginary part of the eigenvalues*/
    void eigenvalues (std::vector< data_Type>& realPart, std::vector< data_Type>& imgPart);

    /** solves the eigenproblem*/
    int solve();
    //@}

private :

    //!@name Private Members
    //@{
    Teuchos::RCP<vector_Type> M_eigenVectors;
    eigenpbPtr_Type MyProblem;
    Teuchos::ParameterList MyPL;
    eigensolverPtr_Type MySolver;
    //@}
};

class UNDEF_EIGENSOLVER_EXCEPTION
{
public:
    UNDEF_EIGENSOLVER_EXCEPTION() {}
    virtual ~UNDEF_EIGENSOLVER_EXCEPTION() {}
};

}
#endif
#endif
