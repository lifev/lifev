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
    @brief NavierStokesPreconditionerOperator - Abstract interface for preconditioners for Navier-Stokes

    @author Davide Forti <davide.forti@epfl.ch>
    @date 08-12-2014

    @maintainer Davide Forti <davide.Forti@epfl.ch>
 */

#ifndef NSPRECONDITIONEROPERATOR_HPP
#define NSPRECONDITIONEROPERATOR_HPP 1

#include <lifev/core/linear_algebra/LinearOperator.hpp>
#include <lifev/core/util/Factory.hpp>
#include <lifev/core/util/FactorySingleton.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{
namespace Operators
{

class NavierStokesPreconditionerOperator : public LinearOperator
{
public:

    typedef MatrixEpetra<Real> matrixEpetra_Type;

    typedef boost::shared_ptr<matrixEpetra_Type> matrixEpetraPtr_Type;

    typedef Epetra_Comm comm_Type;

    typedef boost::shared_ptr<comm_Type> commPtr_Type;

    typedef Epetra_Map map_Type;

    typedef boost::shared_ptr<map_Type> mapPtr_Type;

    typedef boost::shared_ptr<const map_Type> constMapPtr_Type;

    typedef Epetra_Operator operator_Type;

    typedef boost::shared_ptr<operator_Type> operatorPtr_Type;

    typedef  boost::shared_ptr<Teuchos::ParameterList> parameterListPtr_Type;

    NavierStokesPreconditionerOperator();

    ~NavierStokesPreconditionerOperator();

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface
        does not support transpose use, this method should return a value of -1.

        \param In
           UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

        \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
     */
    virtual int SetUseTranspose(bool UseTranspose){};

    //! @name Mathematical functions
    //@{

    //! Returns the result of a raw_operator applied to a raw_vector X in Y.
    /*!
        \param In
           X - A raw_vector of dimension NumVectors to multiply with matrix.
        \param Out
           Y -A raw_vector of dimension NumVectors containing result.

        \return Integer error code, set to 0 if successful.
     */
    virtual int Apply(const vector_Type& X, vector_Type& Y) const {};

    //! Returns the result of a raw_operator inverse applied to an raw_vector X in Y.
    /*!
        \param In
           X - A raw_vector of dimension NumVectors to solve for.
        \param Out
           Y -A raw_vector of dimension NumVectors containing result.

        \return Integer error code, set to 0 if successful.

        \warning In order to work with AztecOO, any implementation of this method must
                  support the case where X and Y are the same object.
     */
    virtual int ApplyInverse(const vector_Type& X, vector_Type& Y) const {};

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
           \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

           \warning This method must not be called unless HasNormInf() returns true.
     */
    virtual double NormInf() const {};

    //! Returns a character string describing the operator
    virtual const char * Label() const {};

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const {};

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const {};

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const comm_Type & Comm() const {};

    //! Returns the raw_map object associated with the domain of this operator.
    virtual const map_Type & OperatorDomainMap() const {};

    //! Returns the raw_map object associated with the range of this operator.
    virtual const map_Type & OperatorRangeMap() const {};
    //@}

    //! @name Methods used by the preconditioners that derive from here
    //@{

    virtual void setUp ( const matrixEpetraPtr_Type & F,
    		   	   	     const matrixEpetraPtr_Type & B,
    		   	   	     const matrixEpetraPtr_Type & Btranspose ){};

    virtual void setUp ( const matrixEpetraPtr_Type & F,
    					 const matrixEpetraPtr_Type & B,
    					 const matrixEpetraPtr_Type & Btranspose,
    					 const matrixEpetraPtr_Type & D ){};

    virtual void setOptions ( const Teuchos::ParameterList& solversOptions){};

    virtual void setDomainMap ( const boost::shared_ptr<BlockEpetra_Map> & domainMap){};

    virtual void setRangeMap ( const boost::shared_ptr<BlockEpetra_Map> & rangeMap){};

    virtual void updateApproximatedMomentumOperator ( ){};

    virtual void updateApproximatedSchurComplementOperator ( ){};

    virtual void updateApproximatedPressureMassOperator ( ){};

    virtual void setUp ( const matrixEpetraPtr_Type & F, const matrixEpetraPtr_Type & B, const matrixEpetraPtr_Type & Btranspose,
               	   	   	 const matrixEpetraPtr_Type & Fp, const matrixEpetraPtr_Type & Mp, const matrixEpetraPtr_Type & Mu){};

    virtual void setMomentumOptions(const parameterListPtr_Type & _oList) {};

    virtual void setSchurOptions(const parameterListPtr_Type & _oList) {};

    virtual void setPressureMassOptions(const parameterListPtr_Type & _oList) {};

    virtual matrixEpetraPtr_Type const& F() const {};

    virtual matrixEpetraPtr_Type const& B() const {};

    virtual matrixEpetraPtr_Type const& Btranspose() const {};

    //! Returns the High Order Yosida approximation of the inverse pressure Schur Complement applied to \c (Xu, Xp).
    virtual int ApplyInverse( VectorEpetra const& X_velocity,
    						  VectorEpetra const& X_pressure,
    						  VectorEpetra & Y_velocity,
    						  VectorEpetra & Y_pressure) const {};


private:

};

inline NavierStokesPreconditionerOperator::NavierStokesPreconditionerOperator()
{
}

inline NavierStokesPreconditionerOperator::~NavierStokesPreconditionerOperator()
{
}

typedef FactorySingleton<Factory<NavierStokesPreconditionerOperator, std::string> >  NSPreconditionerFactory;

} // namespace Operators

} // namespace LifeV


#endif // NSPRECONDITIONEROPERATOR_HPP
