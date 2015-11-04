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
    @brief Epetra preconditioner

    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 09-11-2006
 */

#ifndef _EPETRAPRECONDITIONER_HPP_
#define _EPETRAPRECONDITIONER_HPP_



#include <Teuchos_ParameterList.hpp>


#include <lifev/core/LifeV.hpp>

#include <lifev/core/util/Factory.hpp>
#include <lifev/core/util/FactorySingleton.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>

namespace LifeV
{

// Forward declaration
class SolverAztecOO;

//! Preconditioner - Abstract preconditioner class
/*!
  @author Simone Deparis   <simone.deparis@epfl.ch>
*/
class Preconditioner
{
public:

    //! @name Public Types
    //@{

    typedef Epetra_Operator                      prec_raw_type;
    typedef std::shared_ptr<prec_raw_type>     prec_type;

    typedef MatrixEpetra<Real>                   operator_raw_type;
    typedef std::shared_ptr<operator_raw_type> operator_type;

    typedef Displayer::comm_Type                 comm_Type;
    typedef Displayer::commPtr_Type              commPtr_Type;

    typedef Teuchos::ParameterList               list_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
      @param comm Comminicator
     */
    Preconditioner ( const commPtr_Type& comm = commPtr_Type() );

    //! Copy constructor
    /*!
      @param preconditioner Preconditioner
      @param comm Comminicator
     */
    Preconditioner ( const Preconditioner& preconditioner, const commPtr_Type& comm = commPtr_Type() );

    //! Destructor
    virtual ~Preconditioner();

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
                                        const std::string& subSection ) = 0;

    //! Build a preconditioner based on the given matrix
    /*!
      @param matrix Matrix upon which construct the preconditioner
     */
    virtual Int buildPreconditioner ( operator_type& matrix ) = 0;

    //! Reset the preconditioner
    virtual void resetPreconditioner() = 0;

    //! Return An estimation of the condition number of the preconditioner
    virtual Real condest() = 0;

    //! Show informations about the preconditioner
    virtual void showMe ( std::ostream& output = std::cout ) const;

    //@}


    //! @name Epetra Operator Interface Methods
    //@{

    //! Set the matrix to be used transposed (or not)
    /*!
      @param useTranspose If true the preconditioner is transposed
     */
    virtual Int SetUseTranspose ( const bool useTranspose = false );

    //! Return true if the preconditioner is transposed
    virtual bool UseTranspose();

    //! Apply the inverse of the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    virtual Int Apply ( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const;

    //! Apply the inverse of the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    virtual Int ApplyInverse ( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const;


    //! Return the Range map of the operator
    virtual const Epetra_Map& OperatorRangeMap() const;

    //! Return the Domain map of the operator
    virtual const Epetra_Map& OperatorDomainMap() const;
    //@}


    //! @name Set Methods
    //@{

    //! The the internal list
    /*!
      @param list List to be set into the preconditioner
     */
    void setParametersList ( const list_Type& list );

    //! Set the data of the preconditioner using a GetPot object
    /*!
      @param dataFile A GetPot object containing the data about the preconditioner
      @param section The section in "dataFile" where to find data about the preconditioner
     */
    virtual void setDataFromGetPot ( const GetPot& dataFile, const std::string& section ) = 0;

    //! Set the internal solver
    /*!
      Note: the argument is unused
      @param solver SolverAztecOO
     */
    virtual void setSolver ( SolverAztecOO& /*solver*/ );

    //@}


    //! @name Get Methods
    //@{

    //! Return true if the preconditioner has been created
    const bool& preconditionerCreated();

    //! Return a raw pointer on the preconditioner
    virtual prec_raw_type* preconditioner() = 0;

    //! Return a shared pointer on the preconditioner
    virtual prec_type preconditionerPtr() = 0;

    //! Return the type of preconditioner
    virtual std::string preconditionerType() = 0;

    //! Return the parameters list
    const list_Type& parametersList() const;

    //! Return the parameters list
    list_Type& parametersList();

    //@}

protected:

    std::string M_precType;
    Displayer   M_displayer;
    list_Type   M_list;
    bool        M_preconditionerCreated;

};

typedef FactorySingleton<Factory<Preconditioner, std::string> > PRECFactory;

} // namespace LifeV

#endif
