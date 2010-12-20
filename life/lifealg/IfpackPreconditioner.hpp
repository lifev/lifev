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
    @brief Ifpack preconditioner

    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 09-11-2006
 */

#ifndef _IFPACKPRECONDITIONER_HPP_
#define _IFPACKPRECONDITIONER_HPP_

#include <boost/shared_ptr.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Ifpack_ConfigDefs.h>
#include <Ifpack.h>
#include <Ifpack_Preconditioner.h>
#include <Ifpack_AdditiveSchwarz.h>
#include <Ifpack_Amesos.h>
#include <Ifpack_ILU.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifefilters/GetPot.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraPreconditioner.hpp>

namespace LifeV
{

//! IfpackPreconditioner - Class implementing overlapping Schwarz preconditioner
/*!
  @author Simone Deparis   <simone.deparis@epfl.ch>
*/
class IfpackPreconditioner:
        public EpetraPreconditioner
{
public:

    //! @name Public Types
    //@{

    typedef EpetraPreconditioner             super;

    typedef Ifpack_Preconditioner            prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type> prec_type;

    typedef super::operator_raw_type         operator_raw_type;
    typedef super::operator_type             operator_type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor.
    IfpackPreconditioner();

    //! Destructor
    ~IfpackPreconditioner();

    //@}


    //! @name Methods
    //@{

    //! Build a preconditioner based on the given matrix
    /*!
      @param matrix Matrix upon which construct the preconditioner
     */
    Int buildPreconditioner( operator_type& matrix );

    //! Reset the preconditioner
    void precReset();

    //! Create the list of parameters of the preconditioner
    /*!
      @param list A Parameter list to be filled
      @param dataFile A GetPot object containing the data about the preconditioner
      @param section The section in "dataFile" where to find data about the preconditioner
      @param subSection The subsection in "dataFile" where to find data about the preconditioner
     */
    virtual void createList( list_Type&         list,
                             const GetPot&      dataFile,
                             const std::string& section,
                             const std::string& subSection );

    //! Create the list of parameters of the preconditioner
    /*!
      @param list A Parameter list to be filled
      @param dataFile A GetPot object containing the data about the preconditioner
      @param section The section in "dataFile" where to find data about the preconditioner
      @param subSection The subsection in "dataFile" where to find data about the preconditioner
     */
    static void createIfpackList( list_Type&         list,
                                  const GetPot&      dataFile,
                                  const std::string& section,
                                  const std::string& subSection = "ifpack" );

    //! Apply the inverse of the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    virtual Int ApplyInverse( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const;

    //! Apply the inverse of the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    virtual Int Apply( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const;

    //! Show informations about the preconditioner
    virtual void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the data of the preconditioner using a GetPot object
    /*!
      @param dataFile A GetPot object containing the data about the preconditioner
      @param section The section in "dataFile" where to find data about the preconditioner
     */
    void setDataFromGetPot ( const GetPot&      dataFile,
                             const std::string& section );

    //! Set the matrix to be used transposed (or not)
    /*!
      @param useTranspose If true the preconditioner is transposed
     */
    Int SetUseTranspose( bool useTranspose = false );

    //@}


    //! @name Get Methods
    //@{

    //! Return true if the preconditioner is set
    bool set() const;

    //! Return An estimation of the condition number of the preconditioner
    Real Condest ();

    //! Return a raw pointer on the preconditioner
    super::prec_raw_type* getPrec();

    //! Return a shared pointer on the preconditioner
    super::prec_type getPrecPtr();

    //! Return the type of preconditioner
    std::string precType();

    //! Return the overlap level
    const Int& getOverlapLevel() const;

    //! Return true if the preconditioner is transposed
    bool UseTranspose();

    //! Return the Range map of the operator
    const Epetra_Map & OperatorRangeMap() const;

    //! Return the Domain map of the operator
    const Epetra_Map & OperatorDomainMap() const;

    //@}

protected:

    prec_type M_preconditioner;

private:

    Int M_overlapLevel;
    operator_raw_type::matrix_ptrtype M_operator;

};


inline EpetraPreconditioner* createIfpack() { return new IfpackPreconditioner(); }
namespace
{
static bool registerIF = PRECFactory::instance().registerProduct( "Ifpack", &createIfpack );
}

} // namespace LifeV

#endif
