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
 *  @file
 *  @brief The file contains two classes implementing a wrap up
           of Standard Library vector class to allow indexing from one.
 *
 *  @date 30-08-1999
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *
 *  @contributor Laura Cattaneo
 *  @mantainer Laura Cattaneo
 */

#ifndef _SIMPLEVECT_HH_
#define _SIMPLEVECT_HH_

#include <cstdlib>
#include <vector>
#include <life/lifecore/LifeV.hpp>
#include <algorithm>
#include <iterator>
#include <boost/function.hpp>

namespace LifeV
{
//! A simple predicate to test the boolean flag on a mesh entity
/*!
   @prerequisite MeshEntity must have a method flag_Type flag();

    The ComparisonPolicy passed as (possible) second template parameter must be a functor
    capable of being constructed from a bool (*)(flag_Type const&, flag_Type const &)
    and so that bool operator()(flag_Type const&, flag_Type const &) is defined.
    By default it is a boost function and the class uses the testOneSet policy by default
    Usage: if you want a predicate that tests if a boolean flag is equal to a given flag MYFLAG
    you create a object of type

    EntityFlagInterrogator(MYFLAG,EntityFlagInterrogator::comparisonPolicy_Type(testAllSet))

    which can now be used on all std algorithms operating on containers of mesh entities

 @author Luca Formaggia
 */
template <typename MeshEntity,
typename ComparisonPolicy=boost::function2<bool,flag_Type,flag_Type> >
class EntityFlagInterrogator{
public:
    typedef ComparisonPolicy comparisonPolicy_Type;
    EntityFlagInterrogator(flag_Type flag,
                                    ComparisonPolicy p=ComparisonPolicy(Flag::testOneSet) ):
                                    M_flag(flag),M_policy(p){}

    bool operator()(const MeshEntity & entity){
        return M_policy(entity.flag(),M_flag);
    }
private:
    flag_Type M_flag;
    ComparisonPolicy M_policy;
};

// VectorSimple
/*!
    @author Luca Formaggia

    The class is a wrap up of Standard Library vector class.
    Its role is to held meshEntities.
    Its old name VectorSimple has been changed to MeshEntityContainer
    Any other use of class is deprecated and it should be replaced by std::vector<T>

 */

template <typename DataType, class Allocator = std::allocator<DataType> >
class VectorSimple : public std::vector<DataType,Allocator>
{
public:

    //! @name Public Types
    //@{

    typedef DataType                                data_Type;
    typedef std::vector<DataType,Allocator>         vector_Type;
    typedef typename vector_Type::size_type         size_type;
    typedef typename vector_Type::reference         reference;
    typedef typename vector_Type::const_reference   const_reference;
    typedef typename vector_Type::iterator          iterator;
    typedef typename vector_Type::const_iterator    const_iterator;
    typedef typename vector_Type::const_reverse_iterator const_reverse_iterator;
    typedef typename vector_Type::pointer           pointer;
    typedef typename vector_Type::const_pointer     const_pointer;
    typedef typename vector_Type::reverse_iterator  reverse_iterator;
    typedef typename vector_Type::allocator_type    allocator_type;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    VectorSimple() : vector_Type() {}

    //! Constructor
    /*!
        @param vectorSize size of the vector
     */
    explicit VectorSimple( size_type vectorSize ) : vector_Type( vectorSize ){};

    //! Copy constructor
    /*!
        @param vector VectorSimple vector to copy
     */
    VectorSimple( const VectorSimple<DataType,Allocator> & vector );

    //! Constructor
    /*!
        Construct by copying a Standard Library vector
        @param vector Standard Library vector
     */
    explicit VectorSimple( const vector_Type & vector );

    //! Destructor
    ~VectorSimple() {}

    //@}


    //! @name Operators
    //@{

    //! Assignement operator
    /*!
        Copies source VectorSimple vector into "this"
        @param vector VectorSimple vector
        @return Reference to a new VectorSimple vector with the same
                content of VectorSimple vector
     */
    VectorSimple<DataType, Allocator> & operator=( const VectorSimple<DataType,Allocator> & vector );

    //! Access operator
    /*!
        Example: a(10)=90; // a[10] will contain 90.0
        @param i index of the element of the VectorSimple vector
        @return a vector reference
     */
    reference operator() ( size_type const i )
    {
        return ( this->operator[] ( i ) );
    }

    //! Const access operator
    /*!
        Example: a(10)=90; // a[10] will contain 90.0
        @param i index of the element of the VectorSimple vector
        @return a vector const reference
     */
    const_reference operator() ( size_type const i ) const
    {
        return ( this->operator[] ( i ) );
    }

    //@}


    //! @name Methods
    //@{

    //! Completely clear out the container, returning memory to the system
    inline void clearVector();

    //! Check if the VectorSimple contains an element with index i
    /*!
        @param i index
        @return boolean
     */
    bool checkIndex( size_type const i ) const
    {
        return i >= 0 && i < this->size() ;
    }

    //! Returns the number of elements in the container
    UInt numItems()const
    {
        return this->size();
    }

    //! Returns the capacity of the container
    UInt maxNumItems()const
    {
        return this->capacity();
    }

    //! It sets the capacity of the container
    //! It returns a bool to allow to test whether the container data pool has changed
    //! since in this case pointers to the container elements are invalidated
    //! It does not change the value of the elements currently stored, nor the container size
    //! @param size the new capacity of the container
    //! @return a bool which is true if the data pool has changed
    bool setMaxNumItems(UInt size);
    /*! @brief It extracts all elements that satisfy a condition on the flag
     *
     *  It operates only container elements  where the method
     *  flag_Type flag()  is defined.
     *  Examples:
     *
     *  #include "MeshEntity.hpp"
     *  v=sv.extractElementsWithFlag(PHYSICAL_BOUNDARY,testOneSet)
     *  will extracts all elements in the boundary
     *  v=sv.extractElementsWithFlag(PHYSICAL_BOUNDARY|INTERNAL_INTERFACE, testOneSet)
     *  will extracts all elements on the boundary or on an internal interface
     *
     *  @param refFlag the flag against which the test is made
     *  @param policy. A functor/function pointer which implements
     *  bool policy(const flag_Type & inputFlag, const flag_Type & refFlag)
     *  Available policies testOneSet and testAllSet (defined in Lifev.hpp)
     *  @return a VectorSimple<DataType> containing the elements whose at least one
     *  flag is set as that of refFlag according to the policy
     */
    template<typename Policy>
    VectorSimple<DataType,Allocator> extractElementsWithFlag(const flag_Type & refFlag, Policy const & policy=&Flag::testOneSet);

    /*! @brief It counts all elements that satisfy a condition on the flag
     *  It only operates on container elements where the method
     *  flag_Type flag() is defined.
     *  Examples:
     *
     *  #include "MeshEntity.hpp"
     *  Uint i=sv.countElementsWithFlag(PHYSICAL_BOUNDARY,testOneSet)
     *  will count all elements on the boundary
     *
     *  @param refFlag the flag against which the test is made
     *  @param policy. A functor/function pointer which implements
     *  bool policy(const flag_Type & inputFlag, const flag_Type & refFlag)
     *  Available policies testOneSet and testAllSet (defined in Lifev.hpp)
     *  @return an unsigned integer
     */
    template<typename Policy>
    UInt countElementsWithFlag(const flag_Type & refFlag, Policy const & policy);
    /*!
     * @brief Set element with flag set first
     *
     * It reorders the container starting from the given offset up to the end of the container
     * so that elements with the given flag set are first. A policy must be passed.
     * Typically either testOneSet or testAllSet. A policy object must have a method
     * bool operator()(DataType const & i, DataType const & r)
     * that compares the boolean flag i with the reference flag r
     * Using the algoritms of the std it is possible to do more complex manipulations
     */
    template<typename Policy>
    UInt reorderAccordingToFlag(const flag_Type & refFlag, UInt offset=0, Policy const & policy=
                    &Flag::testOneSet);
    //@}
};


//============================================================================
//                      IMPLEMENTATIONS

// Constructors
//============================================================================
template <typename DataType, class Allocator>
VectorSimple<DataType,Allocator>::VectorSimple( const VectorSimple<DataType,Allocator> & vector )
:vector_Type( vector )
{}


//============================================================================
// Operators
//============================================================================
template <typename DataType,class Allocator>
VectorSimple<DataType,Allocator> &
VectorSimple<DataType,Allocator>::operator=( const VectorSimple<DataType,Allocator> & vector )
{
    vector_Type::operator=( vector );
    return *this;
}


//============================================================================
// Methods
//============================================================================
template <typename DataType, class Allocator>
void
VectorSimple<DataType, Allocator>::clearVector()
{
    this->clear();
    this->vector_Type::swap(vector_Type());
}

template<typename DataType, class Allocator>
bool VectorSimple<DataType,Allocator>::setMaxNumItems(UInt size)
{
    bool _check=(this->capacity()) < size;
    this->reserve(size);
    return _check;
}

template<typename DataType, class Allocator>
template<typename Policy>
VectorSimple<DataType,Allocator>
VectorSimple<DataType,Allocator>::extractElementsWithFlag(
                const flag_Type & refFlag,
                const Policy & policy)
                {
    VectorSimple<DataType,Allocator> M_tmp;
    // There is no copy_if std algorithm (sic), see Stroustoup for
    // the explanation... they just forgot about it. So I do
    // the for loop directly. No need of EntityFlagInterrogator

    for(const_iterator i=this->begin();i<this->end();++i)
        if ( policy(i->flag(),refFlag) ) M_tmp.push_back(*i);

    // To save memory we make sure that capacity()==size()
    // by returning a copy
    return VectorSimple<DataType,Allocator>(M_tmp);
                }

template<typename DataType, class Allocator>
template<typename Policy>
UInt VectorSimple<DataType,Allocator>::countElementsWithFlag(
                const flag_Type & refFlag,
                const Policy & policy=&Flag::testOneSet)
                {
    return std::count_if(this->begin(),this->end(),
                         EntityFlagInterrogator<DataType> (refFlag,
                                                          EntityFlagInterrogator<DataType>::comparisonPolicy_Type(
                                                                          refFlag, policy))
                                                                          );
                }

template<typename DataType, class Allocator>
template<typename Policy>
inline UInt VectorSimple<DataType,Allocator>::reorderAccordingToFlag(const flag_Type & refFlag,
                                                                     UInt offset,
                                                                     Policy const & policy)
{
    iterator last= std::stable_partition(this->begin(),this->begin()+offset,
                                         EntityFlagInterrogator<DataType>(
                                                                  refFlag,
                                                                  EntityFlagInterrogator<DataType>::comparisonPolicy_Type(
                                                                                          refFlag,
                                                                                          policy)
                                                                          )
    );
    return std::distance(this->begin(),last);
}
}
// Namespace LifeV

#endif /* _SIMPLEVECT_HH_ */

