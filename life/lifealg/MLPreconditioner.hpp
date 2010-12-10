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

// #include <Ifpack_config.h>
// #include <Ifpack.h>
// #include <Ifpack_Preconditioner.h>
// #include <Ifpack_AdditiveSchwarz.h>
// #include <Ifpack_Amesos.h>
// #include <Ifpack_ILU.h>
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
    //! default constructor.
    MLPreconditioner();

    //! destructor.
    ~MLPreconditioner();

    //! constructor from matrix A.
    //! @param A EpetraMatrix<Real> matrix upon which construct the preconditioner
    MLPreconditioner( operator_type& matrix );

    //@}

    //! @name Methods
    //@{

    Int buildPreconditioner( operator_type& matrix );

    void precReset();

    void testSmoothers( operator_type& matrix );

    //! returns true if prec exists
    /*const*/
    bool  set() const {return M_preconditioner;}

    virtual void createList( list_Type& list,
                             const GetPot&      dataFile,
                             const std::string& section,
                             const std::string& subSection ) { createMLList( list, dataFile, section, subSection ); }

    static void createMLList( list_Type&   list,
                              const GetPot&      dataFile,
                              const std::string& section,
                              const std::string& subSection = "ML" );

    virtual Int ApplyInverse( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const
    {
        return M_preconditioner->ApplyInverse( vector1, vector2 );
    }

    virtual Int Apply( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const
    {
        return M_preconditioner->Apply( vector1, vector2 );
    }

    virtual void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    void setDataFromGetPot ( const GetPot&      dataFile,
                             const std::string& section );

    Int SetUseTranspose( bool useTranspose=false ) { return M_preconditioner->SetUseTranspose(useTranspose); }

    //@}


    //! @name Get Methods
    //@{

    Real Condest ();

    super::prec_raw_type* getPrec();

    super::prec_type getPrecPtr() { return M_preconditioner; }

    std::string precType() { return M_precType; }

    bool UseTranspose(  ) { return M_preconditioner->UseTranspose(); }

    const Epetra_Map & OperatorRangeMap() const
    { return M_preconditioner->OperatorRangeMap(); }

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
