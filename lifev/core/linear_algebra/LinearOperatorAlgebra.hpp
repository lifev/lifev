/*
 * LinearOperator.hpp
 *
 *  Created on: Sep 3, 2010
 *      Author: uvilla
 */

#ifndef LINEAR_OPERATOR_ALGEBRA_HPP_
#define LINEAR_OPERATOR_ALGEBRA_HPP_

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_MultiVector.h>

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

//Forward declaration (save a lot of re-compiling time when VectorEpetra or MatrixEpetra are modified)
class VectorEpetra;

namespace Operators
{
//! @class LinearOperator
/*! @brief Abstract class which defines the interface of a Linear Operator.
 *
 * LinearOperator is an abstract which inherits from Epetra_Operator.
 * LinearOperator should be the base class for all the LifeV class which implements a linear operator.
 *
 * LinearOperator ensures perfect compatibility with all the Trilinos solvers,
 * plus it supports directly the LifeV::VectorEpetra data.
 *
 * Two concrete methods are implemented in LinearOperator
 * int apply(const VectorEpetra &X, VectorEpetra & Y) const ;
 * int applyInverse(const VectorEpetra &X, VectorEpetra &Y) const.
 *
 *
 * Such methods extract a raw Epetra_MultiVector from the VectorEpetra and then call the virtual methods
 * Int Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
 * or
 * Int ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
 * respectively.
 *
 *
 */
class LinearOperatorAlgebra : public Epetra_Operator
{
public:

    //! @name Public Types
    //@{
    typedef Epetra_Comm comm_Type;
    typedef std::shared_ptr<comm_Type> commPtr_Type;
    typedef Epetra_Map map_Type;
    typedef std::shared_ptr<map_Type> mapPtr_Type;
    typedef std::shared_ptr<const map_Type> constMapPtr_Type;
    typedef Epetra_Operator operator_Type;
    typedef std::shared_ptr<operator_Type> operatorPtr_Type;
    typedef Epetra_MultiVector vector_Type;
    typedef std::shared_ptr<vector_Type> vectorPtr_Type;
    //@}

    //! @name Destructor
    //@{
    //! Destructor
    virtual ~LinearOperatorAlgebra() {};
    //@}

    //! @name Attribute set methods
    //@{

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface
    does not support transpose use, this method should return a value of -1.

    \param In
       UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
    */
    virtual int SetUseTranspose(bool UseTranspose) = 0;
    //@}

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
    virtual int Apply(const vector_Type& X, vector_Type& Y) const = 0;

    //! Returns the result of a LinearOperator applied to a VectorEpetra X in Y.
    /*!
    \param In
       X - A VectorEpetra to multiply with matrix.
    \param Out
       Y -A VectorEpetra containing result.

    \return Integer error code, set to 0 if successful.
    */
    int apply(const VectorEpetra & X, VectorEpetra & Y) const;

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
    virtual int ApplyInverse(const vector_Type& X, vector_Type& Y) const = 0;

    //! Returns the result of a LinearOperator inverse applied to an VectorEpetra X in Y.
    /*!
    \param In
       X - A VectorEpetra to solve for.
    \param Out
       Y -A VectorEpetra containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method must
              support the case where X and Y are the same object.
    */

    int applyInverse(const VectorEpetra & X, VectorEpetra & Y);

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

       \warning This method must not be called unless HasNormInf() returns true.
    */
    virtual double NormInf() const = 0;
    //@}

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char * Label() const = 0;

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const = 0;

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const = 0;

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const comm_Type & Comm() const = 0;

    //! Returns the raw_map object associated with the domain of this operator.
    virtual const map_Type & OperatorDomainMap() const = 0;

    //! Returns the raw_map object associated with the range of this operator.
    virtual const map_Type & OperatorRangeMap() const = 0;
    //@}

};

//! @class IdentityOperator
/*! @brief Identity operator x = I*x. */
class IdentityOperator : public LinearOperatorAlgebra
{
public:

    typedef LinearOperatorAlgebra super;
    typedef super::map_Type map_Type;
    typedef super::mapPtr_Type mapPtr_Type;
    typedef super::vector_Type vector_Type;

    IdentityOperator():M_name("identity"), M_useTranspose(false) {};
    void setUp(const mapPtr_Type & map)
    {
        M_map = map;
    }
    int SetUseTranspose(bool useTranspose)
    {
        M_useTranspose = useTranspose;
        return 0;
    }
    int Apply(const vector_Type & X, vector_Type & Y) const
    {
        Y = X;
        return 0;
    }
    int ApplyInverse(const vector_Type & X, vector_Type & Y) const
    {
        Y = X;
        return 0;
    }
    double NormInf() const
    {
        return 1.0;
    }
    const char * Label() const
    {
        return M_name.c_str();
    }
    bool UseTranspose() const
    {
        return M_useTranspose;
    }
    bool HasNormInf() const
    {
        return true;
    }
    const comm_Type & Comm() const
    {
        return M_map->Comm();
    }
    const map_Type & OperatorDomainMap() const
    {
        return *M_map;
    }
    const map_Type & OperatorRangeMap()  const
    {
        return *M_map;
    }
private:
    std::string M_name;
    mapPtr_Type M_map;
    bool M_useTranspose;
};

//! @class NullOperator
/*! @brief Null operator 0 = Z*x. */
class NullOperator : public LinearOperatorAlgebra
{
public:

    typedef LinearOperatorAlgebra super;
    typedef super::map_Type map_Type;
    typedef super::mapPtr_Type mapPtr_Type;
    typedef super::vector_Type vector_Type;

    NullOperator():M_name("Null Operator"), M_useTranspose(false) {};
    void setUp(const mapPtr_Type & domainMap,
               const mapPtr_Type &  rangeMap)
    {
        M_domainMap = domainMap;
        M_rangeMap = rangeMap;
    }
    int SetUseTranspose(bool useTranspose)
    {
        M_useTranspose = useTranspose;
        return 0;
    }
    int Apply(const vector_Type & /*X*/, vector_Type & Y) const
    {
        Y.PutScalar(0.0);
        return 0;
    }
    int ApplyInverse(const vector_Type & /*X*/, vector_Type & Y) const
    {
        Y.PutScalar(1.0/0.0);
        return -1;
    }
    double NormInf() const
    {
        return 0.0;
    }
    const char * Label() const
    {
        return M_name.c_str();
    }
    bool UseTranspose() const
    {
        return M_useTranspose;
    }
    bool HasNormInf() const
    {
        return true;
    }
    const comm_Type & Comm() const
    {
        return M_rangeMap->Comm();
    }
    const map_Type & OperatorDomainMap() const
    {
        return *M_domainMap;
    }
    const map_Type & OperatorRangeMap()  const
    {
        return *M_rangeMap;
    }
private:
    std::string M_name;
    mapPtr_Type M_domainMap;
    mapPtr_Type M_rangeMap;
    bool M_useTranspose;
};

} /*end namespace Operators*/
} /*end namespace */
#endif /* LINEAR_OPERATOR_ALGEBRA_HPP_ */
