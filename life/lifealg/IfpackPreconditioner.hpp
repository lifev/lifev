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

#include <Ifpack_config.h>
#include <Ifpack.h>
#include <Ifpack_Preconditioner.h>
#include <Ifpack_AdditiveSchwarz.h>
#include <Ifpack_Amesos.h>
#include <Ifpack_ILU.h>

#include <life/lifecore/GetPot.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraPreconditioner.hpp>

namespace LifeV
{

class IfpackPreconditioner:
        public EpetraPreconditioner
{
public:

    //! @name Public Types
    //@{

    typedef EpetraPreconditioner                 super;

    typedef Ifpack_Preconditioner                prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>     prec_type;

    typedef super::operator_raw_type             operator_raw_type;
    typedef super::operator_type                 operator_type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! default constructor.
    IfpackPreconditioner();

    //! constructor from matrix A.
    //! @param A EpetraMatrix<double> matrix upon which construct the preconditioner
    //    IfpackPreconditioner(operator_type& A);

    //! default destructor
    ~IfpackPreconditioner();

    //@}


    //! @name Methods
    //@{

    int                    buildPreconditioner(operator_type& A);

    void                   precReset();



    virtual void createList( list_Type&         list,
                             const GetPot&      dataFile,
                             const std::string& section,
                             const std::string& subSection );

    static void createIfpackList(       list_Type&   list,
                                        const GetPot&      dataFile,
                                        const std::string& section,
                                        const std::string& subSection = "ifpack" );

    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //@}


    //! @name Set Methods
    //@{

    void                   setDataFromGetPot ( const GetPot&      dataFile,
                                               const std::string& section );

    int            SetUseTranspose( bool useTranspose=false );

    //@}


    //! @name Get Methods
    //@{

    bool                   set() const;

    double                 Condest ();

    super::prec_raw_type*  getPrec();

    super::prec_type  getPrecPtr();

    std::string            precType();

    const int& getOverlapLevel() const;

    bool            UseTranspose(  );

    const Epetra_Map & OperatorRangeMap() const;

    const Epetra_Map & OperatorDomainMap() const;

    //@}

protected:

    prec_type              M_Prec;

private:

    int                                 M_overlapLevel;
    operator_raw_type::matrix_ptrtype   M_Oper;

};


inline EpetraPreconditioner* createIfpack() { return new IfpackPreconditioner(); }
namespace
{
static bool registerIF = PRECFactory::instance().registerProduct( "Ifpack", &createIfpack );
}

} // namespace LifeV

#endif
