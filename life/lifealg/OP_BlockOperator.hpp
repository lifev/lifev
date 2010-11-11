/*
 * OP_BlockOperator.hpp
 *
 *  Created on: Sep 20, 2010
 *      Author: uvilla
 */

#ifndef BLOCKOPERATOR_HPP_
#define BLOCKOPERATOR_HPP_

#include <life/lifealg/OP_LinearOperator.hpp>
#include <Epetra_Import.h>

namespace LifeV
{
namespace Operators
{
//! @class BlockOperator
/*! @brief A abstract class for handling n-by-m block operators
 * This class inherits from LifeV::LinearOperator.
 *
 * The Transpose is not supported yet.
 */

class BlockOperator: public LinearOperator
{
public:

	//! @name Public Typedefs
	//@{
	typedef LinearOperator super;
	typedef Epetra_Operator raw_operator;
	typedef boost::shared_ptr<raw_operator> operator_ptr;
	typedef boost::numeric::ublas::matrix<operator_ptr> BlockOper;

	typedef Epetra_MultiVector raw_vector;
	typedef boost::shared_ptr<raw_vector> vector_ptr;
	typedef std::vector<vector_ptr> vector_container;

	enum Structure{Diagonal, LowerTriangular, UpperTriangular, NoStructure, Rectangular};
	//@}

	//! Empty Constructor
	BlockOperator();

	//! @name Set Methods
	//@{
	//! SetUp for a "square operator"
	/*!
	 * @param nBlocks:   number of blocks in a square n-by-n matrix
	 * @param domainMap: the map of a vector in the domain of \e this
	 *                   is obtained by concatenating the block maps in
	 *                   domainMap.
	 *                   rangeMap is assumed to be the same of domainMap.
	 * @param comm:      the communicator.
	 */
	void setUp(UInt nBlocks,
			   const std::vector<boost::shared_ptr<Epetra_Map> > domainMap,
			   const boost::shared_ptr<Epetra_Comm> & comm);

	//! SetUp for a "rectangular operator"
	/*!
	 * @param nRowBlocks, nColBlocks:  \e this is a nRowBlocks-by-nColBlocks
	 *                                  block operator
	 * @param domainMap: the map of a vector in the domain of \e this
	 *                   is obtained by concatenating the block maps in
	 *                   domainMap.
	 * @param rangeMap:  the map of a vector in the range of \e this
	 *                   is obtained by concatenating the block maps in
	 *                   rangeMap.
	 * @param comm:      the communicator.
	 */
	void setUp(UInt nRowBlocks, UInt nColBlocks,
			const std::vector<boost::shared_ptr<Epetra_Map> > & domainMap,
			const std::vector<boost::shared_ptr<Epetra_Map> > & rangeMap,
			const boost::shared_ptr<Epetra_Comm> & comm
			   );

	//! SetUp when the operator is given like a boost::matrix
	/*!
	 * @param blockOper: a dense matrix to describe the block operator
	 * @param comm:      the communicator
	 */
	void setUp(const BlockOper & blockOper, const boost::shared_ptr<Epetra_Comm> & comm);

	//! set a component of the block operator
	/*!
	 * @param iblock, jblock: The position of the block is (iblock, jblock).
	 * @param operBlock     : an operator_ptr representing the block
	 */
	void setBlock(UInt iblock, UInt jblock, const operator_ptr & operBlock);

	//! Complete the block matrix with null operators
	void fillComplete();

	//! If true the transpose of the operator will be computed.
	/*
	 * Not Supported yet
	 */
	int SetUseTranspose(bool useTranspose){M_useTranspose = useTranspose; return -1;}

	//! Compute Y = Op*X;
	virtual int Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const;
	//! Compute Y = Op\X;
	/*!
	 * ApplyInverse is implemented for the matrices with block diagonal, lowerTriangular, upperTriangular form
	 */
	virtual int ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const;

	//! Compute the Inf norm of the operator
	/*!
	 * Not implemented yet.
	 */
    double NormInf() const {return -1;}

    //! Returns a character string describing the operator
    virtual const char * Label() const {return M_name.c_str();}

    //! Returns the current UseTranspose setting.
    bool UseTranspose() const {return M_useTranspose;}

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    bool HasNormInf() const {return false;}

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    const Epetra_Comm & Comm() const {return *M_comm;}

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map & OperatorDomainMap() const {return *M_domainMap;}
    //! Returns the Epetra_Map object associated with the domain of this operator as a pointer
    const boost::shared_ptr<Epetra_Map> & OperatorDomainMap_ptr() const {return M_domainMap;}

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map & OperatorRangeMap() const {return *M_rangeMap;}
    //! Returns the Epetra_Map object associated with the range of this operator as a pointer
    const boost::shared_ptr<Epetra_Map> & OperatorRangeMap_ptr() const {return M_rangeMap;}

    //! Merge two vectors using the domain map
    int merge(const Epetra_MultiVector & vBlock, Epetra_MultiVector & vMono, UInt jblock) const;
    //! Extract vectors using the range map
    int extract(Epetra_MultiVector & vBlock, const Epetra_MultiVector & vMono, UInt jblock) const;

	int split(const Epetra_MultiVector & up,
					vector_container & vi) const;

	int merge(      Epetra_MultiVector & up,
			   const vector_container & vi) const;


protected:
    //! Y = diag(block(i,i)^-1)*X
    int blockJacobi(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const;
    //! Y = backwardsubstitution(X)
    int blockUpperTriangularSolve(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const;
    //! Y = forwardsubstitution(X)
    int blockLowerTriangularSolve(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const;
    //! Change the name of the operator, (avaible for derivate classes).
    void setName(const std::string & name){M_name =name;}
private:
    //! Construct the maps and the importers
    void buildImporter(UInt nblocks,
					   std::vector<boost::shared_ptr<Epetra_Map> > & blockMaps,
					   boost::shared_ptr<Epetra_Map> & fullMap,
					   std::vector< boost::shared_ptr<Epetra_Import> > & block2mono,
					   std::vector< boost::shared_ptr<Epetra_Import> > & mono2block );

    //! Number of blocks in each row
    UInt M_nBlockRows;
    //! Number of blocks in each column
    UInt M_nBlockCols;

    //! Name of the object
    std::string M_name;
    //! Communicator
    boost::shared_ptr<Epetra_Comm> M_comm;

    //! @name Maps
    //@{
    //! Domain Map
    boost::shared_ptr<Epetra_Map> M_domainMap;
    //! Range Map
    boost::shared_ptr<Epetra_Map> M_rangeMap;

    //! Domain Map of each block
	std::vector< boost::shared_ptr<Epetra_Map> > M_domainBlockMaps;
	//! Range Map of each block
	std::vector< boost::shared_ptr<Epetra_Map> > M_rangeBlockMaps;

	//! Shifted domain Map of each block
	std::vector< boost::shared_ptr<Epetra_Map> > M_domainBlockMapsShift;
	//! Shifted domain Map of each block
	std::vector< boost::shared_ptr<Epetra_Map> > M_rangeBlockMapsShift;
	//@}

	//! block operator represented like a dense matrix of pointers to Operators
	BlockOper M_oper;

	//! @name Importers
	//@{
	//! merge block domain vectors in the monolithic domain vector
	std::vector< boost::shared_ptr<Epetra_Import> > M_block2monoDomain;
	//! split the monolithic domain vector in the block domain vectors
	std::vector< boost::shared_ptr<Epetra_Import> > M_mono2blockDomain;
	//! merge block range vectors in the monolithic range vector
	std::vector< boost::shared_ptr<Epetra_Import> > M_block2monoRange;
	//! split the monolithic domain vector in the block domain vectors
	std::vector< boost::shared_ptr<Epetra_Import> > M_mono2blockRange;
	//@}

	//! whenever transpose should be used
	bool M_useTranspose;

	//! structure of the block operator (Diagonal, LowerDiagonal, UpperDiagonal, NoStructure)
	Structure M_structure;
};

} /*end namespace Operators*/
} /*end namespace */
#endif /* BLOCKOPERATOR_HPP_ */
