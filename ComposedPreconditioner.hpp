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
   \file ComposedPreconditioner.hpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2009-05-07
 */


#ifndef _IFPACKCOMPOSEDPREC_HPP_
#define _IFPACKCOMPOSEDPREC_HPP_

#include <boost/shared_ptr.hpp>

#include <Ifpack_config.h>
#include <Ifpack.h>
#include <Ifpack_Preconditioner.h>
#include <ml_MultiLevelPreconditioner.h>
#include <Ifpack_AdditiveSchwarz.h>
#include <Ifpack_Amesos.h>
#include <Ifpack_ILU.h>

#include <life/lifecore/GetPot.hpp>
#include <life/lifearray/EpetraMatrix.hpp>

#include <life/lifealg/EpetraPreconditioner.hpp>
#include "ComposedOperator.hpp"

namespace LifeV
{

class ComposedPreconditioner:
        public EpetraPreconditioner
{
public:

    /** @name Typedefs
     */
    //@{
    typedef EpetraPreconditioner                 super;

//     typedef super::prec_raw_type                 prec_raw_type;
//     typedef super::prec_type                     prec_type;

    typedef ComposedOperator<EpetraPreconditioner>    prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>          prec_type;

    typedef boost::shared_ptr<EpetraPreconditioner>                epetra_prec_type;
    typedef boost::shared_ptr<ML_Epetra::MultiLevelPreconditioner>  ml_prec_type;

    typedef super::operator_raw_type              operator_raw_type;
    typedef boost::shared_ptr<operator_raw_type>  operator_type;
    //@}


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    ComposedPreconditioner( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm>() );

    ComposedPreconditioner( ComposedPreconditioner& P );

    //! constructor from matrix A.
    //! @param A EpetraMatrix<double> matrix upon which construct the preconditioner
    //    ComposedPreconditioner(operator_type& A);

    //! default destructor

    ~ComposedPreconditioner();

    //@}


    /** @name  Methods
     */

    void                   setDataFromGetPot ( const GetPot&      dataFile,
                                               const std::string& section );

    void                   createList( list_Type& /*list*/, const GetPot& dataFile, const std::string& section, const std::string& subSection );
    double                 Condest ();

    super::prec_raw_type*  getPrec ();

    const UInt getNumber() const {return M_Prec->getNumber();}

    super::prec_type              getPrecPtr() {return M_Prec;}

    std::string            precType() {return "composedPreconditioner";}

    int                    buildPreconditioner(operator_type& A);
    int                    buildPreconditioner(operator_type& A,
                                               const bool useInverse,
                                               const bool useTranspose=false);
    //! Build a preconditioner based on A and push it back in the composedPreconditioner.
    int                    push_back          (operator_type& A,
                                               const bool useInverse=false,
                                               const bool useTranspose=false
                                              );

    //! Build a preconditioner based on A and replace it in the composedPreconditioner.
    int                    replace            (operator_type& A,
                                               const UInt index,
                                               const bool useInverse=false,
                                               const bool useTranspose=false);

    void                   precReset();

    //! returns true if prec exists
    /*const*/
    bool                   set() const {return M_Prec;}



    const Epetra_Comm& Comm() {return getPrec()->Comm(); }

    const Epetra_Map& OperatorDomainMap() { return  M_Prec->OperatorDomainMap(); }
    const Epetra_Map& OperatorRangeMap() { return  M_Prec->OperatorRangeMap(); }
    std::vector<operator_type>& getOperVector() {return M_OperVector;}

    int            SetUseTranspose( const bool useTranspose=false )
    {
        return M_Prec->SetUseTranspose(useTranspose);
    }

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

    static EpetraPreconditioner* createComposedPreconditioner()
    {
        return new ComposedPreconditioner();
    }

private:

    int createPrec (operator_type& oper,
                    boost::shared_ptr<EpetraPreconditioner>& prec);

//     int createPrec (operator_type& oper,
//                            boost::shared_ptr<MultiLevelPreconditioner>& prec);

    prec_type                  M_Prec;
    std::vector<operator_type> M_OperVector; // we need to keep track of all the operators.

    //std::string            M_precType;
    static bool registerComposed;
};

} // namespace LifeV

#endif
