/*
 * MatrixContainer.hpp
 *
 *  Created on: Sep 30, 2010
 *      Author: uvilla
 */

#ifndef MATRIXCONTAINER_HPP_
#define MATRIXCONTAINER_HPP_

#include <boost/shared_ptr.hpp>
#include <map>

#include <life/lifearray/EpetraMatrix.hpp>
#include <Teuchos_ParameterList.hpp>

namespace LifeV
{

//! @class MatrixContainer
/*!
 * This class is a container for matrices and other parameters (stored in a Teuchos::ParameterList).
 */

template<typename KEYTYPE = std::string>
class MatrixContainer
{
public:
	//! @name Public Typedef
	//@{

	typedef KEYTYPE KeyType;

	typedef boost::shared_ptr<EpetraMatrix<double> > MatrixType_ptr;
	typedef boost::shared_ptr<Epetra_CrsMatrix> MatrixRawType_ptr;
	typedef typename std::map<KeyType, MatrixType_ptr > Container;
	typedef typename Container::iterator Iterator;
	typedef typename Container::const_iterator CIterator;
	typedef typename Container::value_type ValuePair;
	//@}

	//! Empty constructor
	MatrixContainer():M_pList(new Teuchos::ParameterList){};

	//! @name setters
	//@{
	void  set(const KeyType & _name, const MatrixType_ptr & _matrix);
	template<typename T>
	void setParameter(const std::string & _name, const T & par);
	//@}
	//! @name getters
	//@{
	MatrixType_ptr    getMatrix(const KeyType & _name) const;
	MatrixRawType_ptr get(const KeyType & _name) const;
	template<typename T>
	T getParameter(const std::string & _name, const T & dpar) const;
	//@}

private:
	//! map containing the matrices
	Container M_container;
	//! list for additional parameters
	boost::shared_ptr<Teuchos::ParameterList> M_pList;
};

template<typename KEYTYPE>
void MatrixContainer<KEYTYPE>::set(const KeyType & _name, const MatrixType_ptr & _matrix)
{
	ASSERT_PRE(_matrix.get() != 0, "Setting a null pointer to matrix");
	ValuePair valuePair(_name, _matrix);
	M_container.insert(valuePair);
}

template<typename KEYTYPE>
boost::shared_ptr<EpetraMatrix<double> > MatrixContainer<KEYTYPE>::getMatrix(const KeyType & _name) const
{
	CIterator it(M_container.find(_name));
	ASSERT_POS(it != M_container.end(), "Matrix not found");
	return it->second;
}

template<typename KEYTYPE>
boost::shared_ptr<Epetra_CrsMatrix> MatrixContainer<KEYTYPE>::get(const KeyType & _name) const
{
	CIterator it(M_container.find(_name));
	ASSERT_POS(it != M_container.end(), "Matrix not found");
	return boost::shared_dynamic_cast<Epetra_CrsMatrix>(it->second->getMatrixPtr());
}


template<typename KEYTYPE>
template<typename T>
void MatrixContainer<KEYTYPE>::setParameter(const std::string & _name, const T & par)
{
	M_pList->set(_name, par);
}

template<typename KEYTYPE>
template<typename T>
T MatrixContainer<KEYTYPE>::getParameter(const std::string & _name, const T & dpar) const
{
	// get should be a const method but it is not.
	//This is the reason why I'm using a shared_ptr and not an object.
	return M_pList->get(_name, dpar);
}


} /*end namespace */

#endif /* MATRIXCONTAINER_HPP_ */
