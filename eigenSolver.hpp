/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Paolo Crosetto <paolo.crosetto@epfl.ch>
       Date: 2009-03-25

  Copyright (C) 2009 EPFL

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
   \file eigenSolver.hpp
   \author Paolo Crosetto <paolo.crosetto@epfl.ch>
   \date 2009-03-25
 */


// #ifndef HAVE_TRILINOS_ANASAZI
// #warning: you should use ANASAZI
// #else

// #endif

#include <lifeconfig.h>

#ifdef  HAVE_TRILINOS_ANASAZI


#ifndef EIGENSOLVER_HPP
#define EIGENSOLVER_HPP

#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_RefCountPtrDecl.hpp"

#include "AnasaziBlockDavidson.hpp"
#include "AnasaziLOBPCG.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include <cstdlib>
#include <boost/shared_ptr.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/GetPot.hpp>

namespace LifeV
{
/**
 * class handling the Anasazi BlockKrylovSchur eigensolver. It works for any Epetra_Operator:
 * in particular in Lifev it shuld be used with the ComposedOperator, that inherits direcly
 * from the Epetra_Operator. Since the class ComposedOperator is templated it can contain
 *  both matrices or preconditioners. See Monolithic.cpp for an example using the ComposedOperator.
 */

class UNDEF_EIGENSOLVER_EXCEPTION;

// namespace Epetra
// {

//  template <typename DataType, typename Solver, typename Vector>
  class EigenSolver
  {
  public:
      typedef double DataType;
      typedef Epetra_Operator Solver;
      typedef Epetra_MultiVector Vector;

      typedef Anasazi::BasicEigenproblem<DataType, Vector, Solver>                     eigenpb_raw_type;
      typedef Teuchos::RCP<Anasazi::BasicEigenproblem<DataType, Vector, Solver> >      eigenpb_type;
      typedef Anasazi::BlockKrylovSchurSolMgr <DataType,Vector,Solver>                 eigensolver_raw_type;
      typedef boost::shared_ptr<eigensolver_raw_type>                                  eigensolver_type;

      EigenSolver(boost::shared_ptr<Solver> const matrix, Epetra_BlockMap const& block_map, long unsigned int numvec);

      virtual ~EigenSolver()
      {}

      /**sets the parameters for the eigenproblem:
       */
      void setDataFromGetPot( GetPot const& dataFile, std::string const& section/*="eigensolver"*/ );

      /** set to true if the operator is symmetric*/
      void setHermitian(bool flag){ MyProblem->setHermitian(flag);}

      /** fills the input vectors with the real and imaginary part of the eigenvalues*/
      void eigenvalues(std::vector< DataType>& realPart, std::vector< DataType>& imgPart);

      /** solves the eigenproblem*/
      int solve();


  private :

      Teuchos::RCP<Vector> M_eigenVectors;
      eigenpb_type MyProblem;
      Teuchos::ParameterList MyPL;
      eigensolver_type MySolver;
  };

class UNDEF_EIGENSOLVER_EXCEPTION{
public:
    UNDEF_EIGENSOLVER_EXCEPTION(){}
    virtual ~UNDEF_EIGENSOLVER_EXCEPTION(){}
};

}
#endif
#endif
