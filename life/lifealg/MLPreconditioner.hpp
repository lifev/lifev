/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2006-11-09

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
   \file EpetraPreconditioner.hpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2006-11-09
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

namespace LifeV{

class MLPreconditioner:
        public EpetraPreconditioner
{
public:

    /** @name Typedefs
     */
    //@{
    //$$    typedef Ifpack_Preconditioner prec_raw_type;
    typedef EpetraPreconditioner                 super;

    typedef ML_Epetra::MultiLevelPreconditioner  prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>     prec_type;

    typedef super::operator_raw_type             operator_raw_type;
    typedef super::operator_type                 operator_type;
    //@}


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    MLPreconditioner();

    //@{
    //! destructor.

    ~MLPreconditioner();

    //! constructor from matrix A.
    //! @param A EpetraMatrix<double> matrix upon which construct the preconditioner
    MLPreconditioner(operator_type& A);

    //@}


    /** @name  Methods
     */

    void                    setDataFromGetPot ( const GetPot&      dataFile,
                                                const std::string& section );

    double                  Condest ();

    super::prec_raw_type*   getPrec();

    super::prec_type   getPrecPtr(){return M_Prec;}

    std::string             precType(){return M_precType;}

    int                     buildPreconditioner(operator_type& A);

    void                    precReset();

    void                    testSmoothers(operator_type& A);

    //! returns true if prec exists
    /*const*/ bool  set() const {return M_Prec;}

    virtual void createList(       list_Type&   list,
                             const GetPot&      dataFile,
                             const std::string& section,
                             const std::string& subSection ) {createMLList( list, dataFile, section, subSection);}

    static void createMLList(       list_Type&   list,
                              const GetPot&      dataFile,
                              const std::string& section,
                              const std::string& subSection = "ML" );

    int            SetUseTranspose( bool useTranspose=false ) {return M_Prec->SetUseTranspose(useTranspose);}
    bool            UseTranspose(  ) {return M_Prec->UseTranspose();}

    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
        return M_Prec->ApplyInverse(X, Y);
    }

    virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
        return M_Prec->Apply(X, Y);
    }

    const Epetra_Map & OperatorRangeMap() const
    {return M_Prec->OperatorRangeMap();}

    const Epetra_Map & OperatorDomainMap() const
    {return M_Prec->OperatorDomainMap();}

protected:

    list_Type  M_IFPACKSubList;

private:

    operator_raw_type::matrix_ptrtype   M_Oper;

    prec_type               M_Prec;

    bool                    M_analyze;

};


inline EpetraPreconditioner* createML(){return new MLPreconditioner(); }
namespace
{
	static bool registerML = PRECFactory::instance().registerProduct( "ML", &createML );
}


} // namespace LifeV

#endif
