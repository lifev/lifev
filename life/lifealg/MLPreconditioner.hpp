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

// #include <Ifpack_config.h>
// #include <Ifpack.h>
// #include <Ifpack_Preconditioner.h>
// #include <Ifpack_AdditiveSchwarz.h>
// #include <Ifpack_Amesos.h>
// #include <Ifpack_ILU.h>

#include <ml_MultiLevelPreconditioner.h>

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
    //! @param A EpetraMatrix<double> matrix upon which construct the preconditioner
    MLPreconditioner(operator_type& A);

    //@}

    //! @name Methods
    //@{

    int                     buildPreconditioner(operator_type& A);

    void                    precReset();

    void                    testSmoothers(operator_type& A);

    //! returns true if prec exists
    /*const*/
    bool  set() const {return M_Prec;}

    virtual void createList(       list_Type&   list,
                                   const GetPot&      dataFile,
                                   const std::string& section,
                                   const std::string& subSection ) {createMLList( list, dataFile, section, subSection);}

    static void createMLList(       list_Type&   list,
                                    const GetPot&      dataFile,
                                    const std::string& section,
                                    const std::string& subSection = "ML" );

    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
        return M_Prec->ApplyInverse(X, Y);
    }

    virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
        return M_Prec->Apply(X, Y);
    }

    //@}


    //! @name Set Methods
    //@{

    void                    setDataFromGetPot ( const GetPot&      dataFile,
                                                const std::string& section );

    int            SetUseTranspose( bool useTranspose=false ) {return M_Prec->SetUseTranspose(useTranspose);}

    //@}


    //! @name Get Methods
    //@{

    double                  Condest ();

    super::prec_raw_type*   getPrec();

    super::prec_type   getPrecPtr() {return M_Prec;}

    std::string             precType() {return M_precType;}

    bool            UseTranspose(  ) {return M_Prec->UseTranspose();}

    const Epetra_Map & OperatorRangeMap() const
    {return M_Prec->OperatorRangeMap();}

    const Epetra_Map & OperatorDomainMap() const
    {return M_Prec->OperatorDomainMap();}

    //@}

protected:

    list_Type  M_IFPACKSubList;

private:

    operator_raw_type::matrix_ptrtype   M_Oper;

    prec_type               M_Prec;

    bool                    M_analyze;

};


inline EpetraPreconditioner* createML() {return new MLPreconditioner(); }
namespace
{
static bool registerML = PRECFactory::instance().registerProduct( "ML", &createML );
}


} // namespace LifeV

#endif
