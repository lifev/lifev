/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
      Date: 2004-09-22

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
   \file SolverAztec.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2004-09-22
*/
#ifndef SolverAztec_H
#define SolverAztec_H 1

#include <dataAztec.hpp>
#include <vecUnknown.hpp>

class GetPot;

namespace LifeV
{

/*!
  \class SolverAztec
  \brief wrap aztec linear solvers

  @author Christoph Winkelmann
*/
class SolverAztec
{
public:

    /** @name Typedefs
     */
    //@{

    typedef double value_type;
    typedef SolverAztec solver_type;
    //typedef St::Array::SArray<double,1>::type array_type;
    typedef Vector array_type;

    //@}

    //! describes if/how the preconditioner shall be recalculated
    enum PreCalc { SAME_PRECONDITIONER,
                   SAME_NONZERO_PATTERN,
                   DIFFERENT_NONZERO_PATTERN };
    
    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    //! @param filename GetPot data file containing options for solver in
    //!        section "aztec"
    SolverAztec( std::string filename="" );

    //! create a new instance
    static SolverAztec* New();

    //! destructor
    ~SolverAztec();

    //@}

    /** @name Accessors
     */
    //@{

    double residualNorm() const;

    int iterations() const;

    bool converged() const;

    //@}

    /** @name  Mutators
     */
    //@{

    //! set matrix from MSRMatr
    void setMatrix( MSRMatr<value_type> const& newMatrix );

    /** set matrix from CSRMatr
     *
     *  Warning: The matrix is converted to MSR. This method provides ease of
     *  use, possibly for the sake of efficiency.
     */
    void setMatrix( CSRMatr<CSRPatt, value_type> const& newMatrix );

    /** set matrix free data
     *  @param nEq number of equations
     *  @data data for matvec
     *  @matvec user provided function for matrix vector product
     *  @see Aztec documentation, http://www.cs.sandia.gov/CRF/aztec1.html
     */
    void setMatrixFree( int nEq, void* data,
                        void ( *matvec ) ( double*, double*, AZ_MATRIX_STRUCT*,
                                           int* ) );

    void setTolerance( double newTolerance )
    {
        M_params[ AZ_tol ] = newTolerance;
    }

    //! set the level of recursion in the solution of a linear system
    void setRecursionLevel( int newLevel )
    {
        M_options[ AZ_recursion_level ] = newLevel;
    }

    //@}

    /** @name  Methods
     */
    //@{

    /*
      solve the problem \f$ A x = b \f$

      \c A has been entered via \c setMatrix .

      \c preCalc can have the following values :

      -# \c SAME_PRECONDITIONER -
      Pmat is identical during successive linear solves.
      This option is intended for folks who are using
      different Amat and Pmat matrices and want to reuse the
      same preconditioner matrix.  For example, this option
      saves work by not recomputing incomplete factorization
      for ILU/ICC preconditioners.

      -# \c SAME_NONZERO_PATTERN :
      Pmat has the same nonzero structure during
      successive linear solves.

      -# \c DIFFERENT_NONZERO_PATTERN -
      Pmat does not have the same nonzero structure.
    */
    void solve( array_type& x,
                array_type const& b,
                PreCalc preCalc = DIFFERENT_NONZERO_PATTERN );

    //@}

    /*! Sets options from data file for this solver.
     *  @param dataFile GetPot object containing the options from the data file
     *  @param section section in the GetPot object containing the Aztec stuff
     *
     *  @author Christoph Winkelmann
     */
    void setOptionsFromGetPot( GetPot const& dataFile,
                               std::string section = "aztec" );

private:
    //! private method containing code shared by public setMatrix methods
    void F_setMatrix( MSRMatr<value_type> const& newMatrix );

    //! data organisation for C
    int M_dataOrg[ AZ_COMM_SIZE ];

    //! processor information
    int M_procConfig[ AZ_PROC_SIZE ];

    //! array used to select solver options
    int M_options[ AZ_OPTIONS_SIZE ];

    //! user selected solver parameters
    double M_params[ AZ_PARAMS_SIZE ];

    //! information returned from AZ_solve() indicating success or failure
    double M_status[ AZ_STATUS_SIZE ];

    //! aztec matrix
    AZ_MATRIX* M_matrix;

    //! aztec preconditioner
    AZ_PRECOND* M_precond;

    //! number of the next solver to be created
    static UInt S_solverNumber;

    //! MSRPatt converted from CSRPatt if matrix given as CSRMatr
    std::auto_ptr<MSRPatt> M_tempPattern;

    //! MSRMatr converted from CSRMatr if given as such
    std::auto_ptr<MSRMatr<value_type> > M_tempMatrix;
    
};

}
#endif /* SolverAztec_H */
