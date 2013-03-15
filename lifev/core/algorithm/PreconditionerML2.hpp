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
    @brief ML2 preconditioner

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 21-10-2011
 */

#ifndef _PRECONDITIONERML2_HPP_
#define _PRECONDITIONERML2_HPP_

#include <boost/shared_ptr.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"

#include <ml_MultiLevelPreconditioner.h>
#include "ml_LevelWrap.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wextra"

#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

namespace LifeV
{

//! PreconditionerML2 - Class of multilevels preconditioner
/*!
  @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
  @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
*/
class PreconditionerML2:
    public Preconditioner
{
public:

    //! @name Public Types
    //@{

    typedef Preconditioner                       super;

    typedef Epetra_Operator                      prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>     prec_type;

    typedef super::operator_raw_type             operator_raw_type;
    typedef super::operator_type                 operator_type;

    typedef RegionMesh<LinearTetra>              mesh_Type;
    typedef MapEpetra                            map_Type;
    typedef boost::shared_ptr<FESpace<mesh_Type, map_Type> >  FESpacePtr_Type;
    typedef MatrixEpetra<Real>                   matrix_Type;
    typedef boost::shared_ptr<matrix_Type>       matrixPtr_Type;
    //@}


    //! @name Constructors & Destructor
    //@{
    //! Empty constructor
#ifdef HAVE_MPI
    PreconditionerML2 ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) ) );
#else
    PreconditionerML2 ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_SerialComm ) );
#endif

    //! destructor.
    virtual ~PreconditionerML2();

    //@}

    //! @name Methods
    //@{
    void createParametersList ( list_Type&         list,
                                const GetPot&      dataFile,
                                const std::string& section,
                                const std::string& subSection );

    //! Build a preconditioner based on the given matrix
    /*!
      @param matrix Matrix upon which construct the preconditioner
     */
    Int buildPreconditioner ( operator_type& matrix );

    //! Reset the preconditioner
    void resetPreconditioner();

    //! Returns true if the preconditioner is set
    bool  isPreconditionerSet() const
    {
        return M_preconditioner;
    }

    //! Apply the inverse of the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    virtual Int ApplyInverse ( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const
    {
        return M_preconditioner->ApplyInverse ( vector1, vector2 );
    }

    //! Apply the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    virtual Int Apply ( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const
    {
        return M_preconditioner->Apply ( vector1, vector2 );
    }

    //! Show informations about the preconditioner
    virtual void showMe ( std::ostream& output = std::cout ) const;

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

    //! Set parameter from a Teuchos parameters list
    /*!
      @param list ParameterList from Teuchos
     */
    void setParameters ( Teuchos::ParameterList& list );

    //! Set the matrix to be used transposed (or not)
    /*!
      @param useTranspose If true the preconditioner is transposed
     */
    Int SetUseTranspose ( bool useTranspose = false )
    {
        return M_preconditioner->SetUseTranspose (useTranspose);
    }

    //! Setter for the FESpace
    /*!
        This method set the pointer for the FESpaces needed
        for the construction of the operators Ap, Fp and Mp.
        @param uFESpace Boost::shared_ptr on the FESpace for the velocity
        @param pFESpace Boost::shared_ptr on the FESpace for the pressure
     */
    void setFESpace ( FESpacePtr_Type uFESpace, FESpacePtr_Type pFESpace );

    //@}


    //! @name Get Methods
    //@{

    //! Return An estimation of the condition number of the preconditioner
    Real condest ();

    //! Return a raw pointer on the preconditioner
    super::prec_raw_type* preconditioner();

    //! Return a shared pointer on the preconditioner
    super::prec_type preconditionerPtr()
    {
        return M_preconditioner;
    }

    //! Return the type of preconditioner
    std::string preconditionerType()
    {
        return M_precType;
    }

    //! Return true if the preconditioner is transposed
    bool UseTranspose()
    {
        return M_preconditioner->UseTranspose();
    }

    //! Return the Range map of the operator
    const Epetra_Map& OperatorRangeMap() const
    {
        return M_preconditioner->OperatorRangeMap();
    }

    //! Return the Domain map of the operator
    const Epetra_Map& OperatorDomainMap() const
    {
        return M_preconditioner->OperatorDomainMap();
    }

    //@}

protected:

    list_Type                         M_MLList;
    operator_raw_type::matrix_ptrtype M_operator;
    prec_type                         M_preconditioner;
    matrixPtr_Type                    M_restriction;
    matrixPtr_Type                    M_prolongation;
    bool                              M_useP2ToP1Wrapper;
    FESpacePtr_Type                   M_uFESpace;
    FESpacePtr_Type                   M_pFESpace;
    boost::shared_ptr<Epetra_Comm>    M_comm;

};


inline Preconditioner* createML2()
{
    return new PreconditionerML2();
}
namespace
{
static bool registerML2 = PRECFactory::instance().registerProduct ( "ML2", &createML2 );
}


} // namespace LifeV

#endif // _PRECONDITIONERML2_
