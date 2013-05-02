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

/*!
    @file
    @brief Preconditioner Solver Belos

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 17-09-2011
 */

#ifndef PRECONDITIONERSOLVERBELOS_HPP
#define PRECONDITIONERSOLVERBELOS_HPP 1

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"
#include <Epetra_Operator.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/filter/GetPot.hpp>

class GetPot;

namespace LifeV
{

//! PreconditionerLinearSolver - Class to wrap linear solver
/*!
  @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
*/
class PreconditionerLinearSolver
    : public Preconditioner
{
public:

    //! @name Public Types
    //@{

    typedef Epetra_Operator                        prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>       prec_type;

    typedef LinearSolver                           solver_Type;
    typedef boost::shared_ptr<solver_Type>         solverPtr_Type;

    typedef MatrixEpetra<Real>                     operator_raw_type;
    typedef boost::shared_ptr<operator_raw_type>   operator_type;

    typedef Preconditioner                         preconditioner_Type;
    typedef boost::shared_ptr<preconditioner_Type> preconditionerPtr_Type;

    typedef Displayer::comm_Type                   comm_Type;
    typedef Displayer::commPtr_Type                commPtr_Type;

    typedef Teuchos::ParameterList                 list_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Default constructor
    /*!
     * @param comm The communicator.
     */
#ifdef HAVE_MPI
    PreconditionerLinearSolver ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) ) );
#else
    PreconditionerLinearSolver ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_SerialComm ) );
#endif

    //! Destructor
    virtual ~PreconditionerLinearSolver();

    //@}

    //! @name Methods
    //@{

    //! Create the list of parameters of the preconditioner
    /*!
      @param list A Parameter list to be filled
      @param dataFile A GetPot object containing the data about the preconditioner
      @param section The section in "dataFile" where to find data about the preconditioner
      @param subSection The subsection in "dataFile" where to find data about the preconditioner
     */
    virtual void createParametersList ( list_Type& list,
                                        const GetPot& dataFile,
                                        const std::string& section,
                                        const std::string& subSection );

    //! Create the list of parameters of the preconditioner
    /*!
      @param list A Parameter list to be filled
      @param dataFile A GetPot object containing the data about the preconditioner
      @param section The section in "dataFile" where to find data about the preconditioner
      @param subSection The subsection in "dataFile" where to find data about the preconditioner
     */
    static void createLinearSolverList ( list_Type&         list,
                                         const GetPot&      dataFile,
                                         const std::string& section,
                                         const std::string& subSection = "SolverLinear",
                                         const bool&        verbose = true );

    //! Build a preconditioner based on the given matrix
    /*!
      @param matrix Matrix upon which construct the preconditioner
     */
    virtual Int buildPreconditioner ( operator_type& matrix );

    //! Reset the preconditioner
    virtual void resetPreconditioner();

    //! Return An estimation of the condition number of the preconditioner
    virtual Real condest();

    //! Show informations about the preconditioner
    virtual void showMe ( std::ostream& output = std::cout ) const;

    //@}


    //! @name Epetra Operator Interface Methods
    //@{

    //! Set the matrix to be used transposed (or not)
    /*!
      @param useTranspose If true the preconditioner is transposed
     */
    Int SetUseTranspose ( const bool useTranspose = false );

    //! Return true if the preconditioner is transposed
    bool UseTranspose();

    //! Apply the inverse of the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    Int Apply ( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const;

    //! Apply the inverse of the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    Int ApplyInverse ( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const;


    //! Return the Range map of the operator
    const Epetra_Map& OperatorRangeMap() const;

    //! Return the Domain map of the operator
    const Epetra_Map& OperatorDomainMap() const;
    //@}


    //! @name Set Methods
    //@{

    //! Set the data of the preconditioner using a GetPot object
    /*!
      @param dataFile A GetPot object containing the data about the preconditioner
      @param section The section in "dataFile" where to find data about the preconditioner
     */
    void setDataFromGetPot ( const GetPot& dataFile, const std::string& section );

    //! Set the internal solver
    /*!
      Note: the argument is unused
      @param solver SolverAztecOO
     */
    void setSolver ( SolverAztecOO& /*solver*/ );

    //@}


    //! @name Get Methods
    //@{

    //! Return true if the preconditioner is set
    bool isPreconditionerSet() const;

    //! Return a raw pointer on the preconditioner
    prec_raw_type* preconditioner();

    //! Return a shared pointer on the preconditioner
    prec_type preconditionerPtr();

    //! Return a shared pointer on the solver
    solverPtr_Type solverPtr();

    //! Return the type of preconditioner
    std::string preconditionerType();

    //@}

private:

    solverPtr_Type                 M_solver;
    preconditionerPtr_Type         M_preconditioner;
    GetPot                         M_dataFile;

    bool                           M_printSubiterationCount;
    std::string                    M_precName;
    std::string                    M_precDataSection;

};

inline Preconditioner* createLinearSolverPreconditioner()
{
    return new PreconditionerLinearSolver();
}
namespace
{
static bool registerSB = PRECFactory::instance().registerProduct ( "LinearSolver", &createLinearSolverPreconditioner );
}

} // namespace LifeV

#endif /* PRECONDITIONERAMESOS_HPP */
