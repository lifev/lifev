/*
 * OP_LinearOperator.hpp
 *
 *  Created on: Sep 3, 2010
 *      Author: uvilla
 */

#ifndef LINEAROPERATOR_HPP_
#define LINEAROPERATOR_HPP_


#include <life/lifecore/life.hpp>

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_MultiVector.h>

#include <life/lifearray/EpetraVector.hpp>

namespace LifeV
{
namespace Operators
{
//! @class LinearOperator
/*! @brief Abstract class which defines the interface of a Linear Operator.
 *
 * LinearOperator is an abstract which inherits from Epetra_Operator.
 * LinearOperator should be the base class for all the LifeV class which implements a linear operator.
 *
 * LinearOperator ensures perfect compatibility with all the Trilinos solvers,
 * plus it supports directly the LifeV::EpetraVector data.
 *
 * Two concrete methods are implemented in LinearOperator
 * int apply(const EpetraVector &X, EpetraVector & Y) const ;
 * int applyInverse(const EpetraVector &X, EpetraVector &Y) const.
 *
 *
 * Such methods extract a raw Epetra_MultiVector from the EpetraVector and then call the virtual methods
 * Int Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
 * or
 * Int ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
 * respectively.
 *
 *
 */
class LinearOperator : public Epetra_Operator
{
public:

	   //! @name Destructor
	  //@{
	    //! Destructor
	    virtual ~LinearOperator() {};
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

	    //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
	    /*!
	    \param In
		   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
	    \param Out
		   Y -A Epetra_MultiVector of dimension NumVectors containing result.

	    \return Integer error code, set to 0 if successful.
	  */
	    virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

	    //! Returns the result of a LinearOperator applied to a EpetraVector X in Y.
	    /*!
	    \param In
		   X - A EpetraVector to multiply with matrix.
	    \param Out
		   Y -A EpetraVector containing result.

	    \return Integer error code, set to 0 if successful.
	  */
	    inline int apply(const EpetraVector & X, EpetraVector & Y) const
	    		{
	    	return Apply(X.getEpetraVector(), Y.getEpetraVector());
	    		}

	    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
	    /*!
	    \param In
		   X - A Epetra_MultiVector of dimension NumVectors to solve for.
	    \param Out
		   Y -A Epetra_MultiVector of dimension NumVectors containing result.

	    \return Integer error code, set to 0 if successful.

	    \warning In order to work with AztecOO, any implementation of this method must
	              support the case where X and Y are the same object.
	  */
	    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

	    //! Returns the result of a LinearOperator inverse applied to an EpetraVector X in Y.
	    /*!
	    \param In
		   X - A EpetraVector to solve for.
	    \param Out
		   Y -A EpetraVector containing result.

	    \return Integer error code, set to 0 if successful.

	    \warning In order to work with AztecOO, any implementation of this method must
	              support the case where X and Y are the same object.
	  */

	    inline int applyInverse(const EpetraVector & X, EpetraVector Y)
	    {
	    	return ApplyInverse(X.getEpetraVector(), Y.getEpetraVector());
	    }

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
	    virtual const Epetra_Comm & Comm() const = 0;

	    //! Returns the Epetra_Map object associated with the domain of this operator.
	    virtual const Epetra_Map & OperatorDomainMap() const = 0;

	    //! Returns the Epetra_Map object associated with the range of this operator.
	    virtual const Epetra_Map & OperatorRangeMap() const = 0;
	  //@}

	};

//! @class IdentityOperator
/*! @brief Identity operator x = I*x. */
class IdentityOperator : public LinearOperator
{
public:
	IdentityOperator():M_name("identity"), M_useTranspose(false){};
	void setUp(const boost::shared_ptr<Epetra_Map> & map){M_map = map;}
	int SetUseTranspose(bool useTranspose){M_useTranspose = useTranspose; return 0; }
	int Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const {Y = X; return 0;}
	int ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const {Y=X; return 0;}
	double NormInf() const {return 1.0;}
	const char * Label() const {return M_name.c_str();}
	bool UseTranspose() const {return M_useTranspose;}
	bool HasNormInf() const {return true;}
	const Epetra_Comm & Comm() const {return M_map->Comm();}
	const Epetra_Map & OperatorDomainMap() const {return *M_map;}
	const Epetra_Map & OperatorRangeMap()  const {return *M_map;}
private:
	std::string M_name;
	boost::shared_ptr<Epetra_Map> M_map;
	bool M_useTranspose;
};

//! @class NullOperator
/*! @brief Null operator 0 = Z*x. */
class NullOperator : public LinearOperator
{
public:
	NullOperator():M_name("Null Operator"), M_useTranspose(false){};
	void setUp(const boost::shared_ptr<Epetra_Map> & domainMap,
			   const boost::shared_ptr<Epetra_Map> &  rangeMap)
	{
				M_domainMap = domainMap;
				M_rangeMap = rangeMap;
	}
	int SetUseTranspose(bool useTranspose){M_useTranspose = useTranspose; return 0; }
	int Apply(const Epetra_MultiVector & /*X*/, Epetra_MultiVector & Y) const {Y.PutScalar(0.0); return 0;}
	int ApplyInverse(const Epetra_MultiVector & /*X*/, Epetra_MultiVector & Y) const {Y.PutScalar(1.0/0.0); return -1;}
	double NormInf() const {return 0.0;}
	const char * Label() const {return M_name.c_str();}
	bool UseTranspose() const {return M_useTranspose;}
	bool HasNormInf() const {return true;}
	const Epetra_Comm & Comm() const {return M_rangeMap->Comm();}
	const Epetra_Map & OperatorDomainMap() const {return *M_domainMap;}
	const Epetra_Map & OperatorRangeMap()  const {return *M_rangeMap;}
private:
	std::string M_name;
	boost::shared_ptr<Epetra_Map> M_domainMap;
	boost::shared_ptr<Epetra_Map> M_rangeMap;
	bool M_useTranspose;
};

} /*end namespace Operators*/
} /*end namespace */
#endif /* LINEAROPERATOR_HPP_ */
