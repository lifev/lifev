/*
 * BlockEpetra_MultiVector.hpp
 *
 *  Created on: Aug 14, 2011
 *      Author: uvilla
 */
//@ HEADER

/*!
 * \file BlockEpetra_MultiVector.hpp
 * \author Umberto Villa
 * \date 2011-08-13
 * A derived class from \c Epetra_MultiVector specialized to handle parallel block structured Vectors.
 */

#ifndef BLOCKEPETRA_MULTIVECTOR_HPP_
#define BLOCKEPETRA_MULTIVECTOR_HPP_

#include <lifev/core/LifeV.hpp>

// Tell the compiler to ignore specific kind of warnings
LIFEV_SUPPRESS_WARNINGS

#include <Epetra_MultiVector.h>

// Tell the compiler to restore the warnings
LIFEV_RESTORE_WARNINGS

#include <lifev/operator/linear_algebra/BlockEpetra_Map.hpp>

namespace LifeV
{
//! @class
/*!
 * @brief A derived class from \c Epetra_MultiVector specialized to handle parallel block structured Vectors.
 *
 * \c BlockEpetra_MultiVector provides monolithic or block access to an Epetra_MultiVector.
 * A \c BlockEpetra_MultiVector object can be constructed or by striding \c Epetra_MultiVector or from a
 * \c BlockEpetra_Map.
 *
 * This class can be also used to overly a block structure to a given Epetra_MultiVector. In this case, to
 * maximize performances and reduce memory footprint the original object and it's block view will share the same
 * data in memory. If a \c BlockEpetra_MultiVector object is a view of another Epetra_MultiVector each modification on the
 * data value of one vector will also affect the other vector and viceversa.
 */
class BlockEpetra_MultiVector : public Epetra_MultiVector
{

public:

    //@name Typedef
    //@{
    typedef Epetra_MultiVector vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;
    typedef std::vector<vectorPtr_Type>    vectorPtrContainer_Type;
    //@}

    //@name Constructors
    //@{
    //! Generate a BlockEpetra_MultiVector from a BlockEpetra_Map
    BlockEpetra_MultiVector(const BlockEpetra_Map& map, int numVectors, bool zeroOut = true);
    //! Copy Constructor
    BlockEpetra_MultiVector(const BlockEpetra_MultiVector& Source);
    //! Overlay a block structure to a source.
    /*!
     * If the DataAccess type is Copy, then this object will make a deep copy of all the data structure of source
     * If the DataAccess type is View, then this object will only copy the pointers to the datas in source. (shallow copy).
     * \warning View access is faster but the user must ensure that \c source will not be deallocated before \c this
     */
    BlockEpetra_MultiVector(Epetra_DataAccess CV, const vector_Type & source, const BlockEpetra_Map& map);
    //@}

    //@name Getters
    //@{
    //! Return a Epetra_MultiVector (view) object for block iblock
    /*!
     * \warning The returned object doesn't own the internal data structure. User must ensure correct scoping.
     */
    Epetra_MultiVector & block(UInt iblock);
    //! const version
    const Epetra_MultiVector & block(UInt iblock) const;

    //! retrieve the \c BlockEpetra_Map
    const BlockEpetra_Map & blockEpetraMap() const;
    //@}
private:

    void createBlockViews();

    UInt M_nBlocks;
    BlockEpetra_Map M_blockMap;
    std::vector<UInt> M_myLocalOffsets;
    vectorPtrContainer_Type M_blocks;
};

//! Generate a BlockEpetra_MultiVector from a list of Epetra_MultiVector.
BlockEpetra_MultiVector * stride(std::vector<const BlockEpetra_MultiVector::vector_Type *> vector);
//! Generate a BlockEpetra_MultiVector from two Epetra_MultiVectors
BlockEpetra_MultiVector * stride(const BlockEpetra_MultiVector::vector_Type & v1,
                                 const BlockEpetra_MultiVector::vector_Type & v2);
//! Generate a BlockEpetra_MultiVector from three Epetra_MultiVectors
BlockEpetra_MultiVector * stride(const BlockEpetra_MultiVector::vector_Type & v1,
                                 const BlockEpetra_MultiVector::vector_Type & v2,
                                 const BlockEpetra_MultiVector::vector_Type & v3);
//! Generate a BlockEpetra_MultiVector from three Epetra_MultiVectors
BlockEpetra_MultiVector * stride(const BlockEpetra_MultiVector::vector_Type & v1,
                                 const BlockEpetra_MultiVector::vector_Type & v2,
                                 const BlockEpetra_MultiVector::vector_Type & v3,
                                 const BlockEpetra_MultiVector::vector_Type & v4);
/*!
 *  Generate a BlockEpetra_MultiVector from a Epetra_MultiVector and a BlockEpetra_Map
 *  Note the \c BlockEpetra_MultiVector object and the \c source vector will share the same internal data structure.
 *  Each modification in the entry of any of the two vector will affect also the other vector.
 *  \warning The user has the responsibility to keep the scope of the \c BlockEpetra_MultiVector inside the scope of the vector \c source.
 *  If the vector source is distructed, then the internal data structure of the returned object will contain
 *   dangling pointers.
 */
BlockEpetra_MultiVector * createBlockView(const BlockEpetra_MultiVector::vector_Type & source, const BlockEpetra_Map& map);

} /* end namespace */


#endif /* BLOCKEPETRA_MULTIVECTOR_HPP_ */
