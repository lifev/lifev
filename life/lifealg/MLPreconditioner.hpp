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
    @brief ML preconditioner

    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 09-11-0006
 */

#ifndef _MLPRECONDITIONER_HPP_
#define _MLPRECONDITIONER_HPP_

#include <boost/shared_ptr.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <ml_MultiLevelPreconditioner.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/GetPot.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraPreconditioner.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>

namespace LifeV
{

//! MLPreconditioner - Class of multilevels preconditioner
/*!
  @author Simone Deparis   <simone.deparis@epfl.ch>
*/
class MLPreconditioner:
        public EpetraPreconditioner
{
public:

    //! @name Public Types
    //@{

    typedef EpetraPreconditioner                 super;

    typedef ML_Epetra::MultiLevelPreconditioner  prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>     prec_type;

    typedef super::operator_raw_type             operator_raw_type;
    typedef super::operator_type                 operator_type;

    //@}


    //! @name Constructors & Destructor
    //@{
    //! Empty constructor.
    MLPreconditioner();

    //! destructor.
    ~MLPreconditioner();

    //! Constructor from a matrix
    /*!
      @param matrix Matrix upon which construct the preconditioner
    */
    MLPreconditioner( operator_type& matrix );

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

    //! ?
    void testSmoothers( operator_type& matrix );

    //! Returns true if the preconditioner is set
    bool  set() const {return M_preconditioner;}

    //! Create the list of parameters of the preconditioner
    /*!
      @param list A Parameter list to be filled
      @param dataFile A GetPot object containing the data about the preconditioner
      @param section The section in "dataFile" where to find data about the preconditioner
      @param subSection The subsection in "dataFile" where to find data about the preconditioner
     */
    virtual void createList( list_Type& list,
                             const GetPot& dataFile,
                             const std::string& section,
                             const std::string& subSection ) { createMLList( list, dataFile, section, subSection ); }

    //! Create the list of parameters of the preconditioner
    /*!
      @param list A Parameter list to be filled
      @param dataFile A GetPot object containing the data about the preconditioner
      @param section The section in "dataFile" where to find data about the preconditioner
      @param subSection The subsection in "dataFile" where to find data about the preconditioner
     */
    static void createMLList( list_Type& list,
                              const GetPot& dataFile,
                              const std::string& section,
                              const std::string& subSection = "ML" );

    //! Apply the inverse of the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    virtual Int ApplyInverse( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const
    {
        return M_preconditioner->ApplyInverse( vector1, vector2 );
    }

    //! Apply the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    virtual Int Apply( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const
    {
        return M_preconditioner->Apply( vector1, vector2 );
    }

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
    Int SetUseTranspose( bool useTranspose=false ) { return M_preconditioner->SetUseTranspose(useTranspose); }

    //@}


    //! @name Get Methods
    //@{

    //! Return An estimation of the condition number of the preconditioner
    Real Condest ();

    //! Return a raw pointer on the preconditioner
    super::prec_raw_type* getPrec();

    //! Return a shared pointer on the preconditioner
    super::prec_type getPrecPtr() { return M_preconditioner; }

    //! Return the type of preconditioner
    std::string precType() { return M_precType; }

    //! Return true if the preconditioner is transposed
    bool UseTranspose() { return M_preconditioner->UseTranspose(); }

    //! Return the Range map of the operator
    const Epetra_Map & OperatorRangeMap() const
    { return M_preconditioner->OperatorRangeMap(); }

    //! Return the Domain map of the operator
    const Epetra_Map & OperatorDomainMap() const
    { return M_preconditioner->OperatorDomainMap(); }

    //@}

protected:

    list_Type  M_IfpackSubList;

private:

    operator_raw_type::matrix_ptrtype M_operator;

    prec_type               M_preconditioner;

    bool                    M_analyze;

};


inline EpetraPreconditioner* createML() {return new MLPreconditioner(); }
namespace
{
static bool registerML = PRECFactory::instance().registerProduct( "ML", &createML );
}


} // namespace LifeV

#endif
