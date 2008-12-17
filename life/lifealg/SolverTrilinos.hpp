/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Simone Deparis   <simone.deparis@epfl.ch>
            Gilles Fourestey <gilles.fourestey@epfl.ch>
      Date: 2006-11-08

 Copyright (C) 2006 EPFL

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
   \file SolverTrilinos.hpp
   \author Simone Deparis   <simone.deparis@epfl.ch>
   \author Gilles Fourestey <gilles.fourestey@epfl.ch>
   \date 2006-11-08
*/

#ifndef __SolverTrilinos_H
#define __SolverTrilinos_H 1

#include <boost/shared_ptr.hpp>

#include "AztecOO_config.h"
#include "AztecOO.h"

#include "Teuchos_ParameterList.hpp"
#include "life/lifearray/EpetraVector.hpp"
#include "life/lifearray/EpetraMatrix.hpp"

#include "life/lifealg/EpetraPreconditioner.hpp"
#include "life/lifealg/IfpackPreconditioner.hpp"

class GetPot;

namespace LifeV
{
// namespace Epetra
// {
/*!
  \class SolverTrilinos
  \brief wrap trilinos linear solvers

  By default the solver is gmres and the preconditioner is ilu.

  \author Simone Deparis   <simone.deparis@epfl.ch>
  \author Gilles Fourestey <gilles.fourestey@epfl.ch>
  @see
*/

class SolverTrilinos
{
public:


    /** @name Typedefs
     */
    //@{

    typedef double                           value_type;

    typedef SolverTrilinos                   solver_type;

    typedef EpetraMatrix<double>             matrix_type;
    typedef EpetraVector                     vector_type;

    typedef EpetraPreconditioner             prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type> prec_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    //! @param filename GetPot data file containing options for solver in
    //!        section "aztec"
    SolverTrilinos();

    //@}

    /** @name Accessors
     */
    //@{

    //! Total number of iterations
    int    NumIters();

    //! True Residual
    double TrueResidual();

    //@}

    /** @name  Mutators
     */
    //@{

    //! set matrix from EpetraMatrix
    void setMatrix(matrix_type& m);

    void setOperator(Epetra_Operator& op);

    //! set Epetra_Operator preconditioner
    void setPreconditioner( prec_type _prec );

    void setDataFromGetPot( const GetPot& dfile, const std::string& section );
    void SetParameters(bool cerr_warning_if_unused = false);

    void useGMRES(const int restart=300);
    void useCG();
    void useBICGSTAB();
    void useJacobyPrec();
    void useNoPrec();
    void useDDPrec();
    void useNeumannPrec();
    void setReuse();
    void setTolMaxiter(const double tol, const int maxiter=-1);


    //! describes the verobisty level
    enum VerboseLevel { NONE, SUMMARY, LAST };
    //! sets verbosity level
    void SetVerbose(const VerboseLevel verb=NONE);


    //@}


    /** @name  Methods
     */
    //@{

    /*
      solve the problem \f$ A x = b \f$

      \c A has been entered via \c setMatrix .

      return the number of iterations, M_maxIter+1 if solve failed

    */
    int solve( vector_type& x, vector_type& b );

    double computeResidual( vector_type& __X, vector_type& __B );

    bool precSet() const {return (M_prec.get() !=0 && M_prec->getPrec() != 0);}
    void precReset() { M_prec->precReset(); }
//@}


private:

    prec_ptr               M_prec;

    AztecOO                M_solver;

    Teuchos::ParameterList M_TrilinosParameterList;

    int                    M_maxIter;
    double                 M_tol;

};

// } // namespace Epetra
} // namespace LifeV

#endif /* __SolverTrilinos_H */

    /*
    //! set the level of recursion in the solution of a linear system
    void setRecursionLevel( int newLevel )
    {
        M_options[ AZ_recursion_level ] = newLevel;
    }
    */

