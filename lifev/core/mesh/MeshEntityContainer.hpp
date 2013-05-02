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

#ifndef _MESHENTITYCONTAINER_HH_
#define _MESHENTITYCONTAINER_HH_

#include <algorithm>
#include <iterator>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/Marker.hpp>
#include <lifev/core/mesh/MeshEntity.hpp>

namespace LifeV
{
//! @todo Move to aspecific file meshentityutility

namespace Comparers
{

/** @defgroup ComparisonOperators
 * They define a comparison operators for mesh entities according to some of
 * their attributes
 * The comparison operator is of the form
 * @code
 * bool operator(MeshEntity const & a, MeshEntity const & b)
 * @endcode
 * and it should generate a well posed comparison.
 * @{
 */
/** Compare according to local ID.
 *  It compares according to the ID (local ID) of a MeshEntity and it relies on std comparison operators
 *  It defaults to std::less<ID>.
 *  We rely on the fact that std::less<ID>(ID a, ID b) is defined (otherwise the user must supply it).
 */
template < typename MeshEntity,
         typename Policy = boost::function2<bool, ID, ID> >
class CompareAccordingToLocalId
{
public:
    //! Constructor
    /**
     * The constructor can receive anything which is convertible the chosen
     * policy type. This allows great flexibility (maybe too much!)
     */
    CompareAccordingToLocalId ( Policy const& p = std::less<ID>() ) : M_policy ( p ) {}

    bool operator() ( MeshEntity const& a, MeshEntity const& b )
    {
        return M_policy ( a.localId(), b.localId() );
    }

private:
    const Policy M_policy;
};

/** Compare according to Marker ID.
 *  It compares according to the Marker of a MeshEntity and it relies on std comparison operators
 *  It defaults to std::less<markerID_Type>. We rely on the fact that less<ID>(ID a, ID b) (otherwise the user must supply it).
 */
template <typename MeshEntity, typename Policy = boost::function2< bool, markerID_Type, markerID_Type> >
class CompareAccordingToMarker
{
public:
    //! Constructor
    /**
     * The constructor can receive anything which is convertible the chosen
     * policy type. This allows great flexibility (maybe too much!)
     */
    CompareAccordingToMarker ( Policy const& p = std::less<markerID_Type>() ) : M_policy ( p ) {}

    bool operator() ( MeshEntity const& a, MeshEntity const& b )
    {
        return M_policy ( a.markerID(), b.markerID() );
    }

private:
    const Policy M_policy;
};

/** @}*/

}

namespace Predicates
{
/** @defgroup Predicates
 *
 * Predicates are functors that take MeshEntity as template argument and
 * implement
 *
 * @code
 *  bool operator()(const MeshEntity & entity)const
 * @endcode
 * @{
 */
//! A simple predicate to test the boolean flag on a mesh entity
/*!
   @prerequisite MeshEntity must have a method flag_Type flag();

    The ComparisonPolicy passed as (possible) second template parameter must be a functor
    capable of being constructed from a bool (*)(flag_Type const&, flag_Type const &)
    and so that bool operator()(flag_Type const&, flag_Type const &) is defined.
    By default it is a boost function and the class uses the testOneSet policy by default
    Usage: if you want a predicate that tests if a boolean flag is equal to a given flag MYFLAG
    you create a object of type

@code
    EntityFlagInterrogator<faceType> interrogator(MYFLAG,Flags::testAllSet)

    interrogator(myFace); // true if flags of myface all all equal to MYFLAG
@endcode

    which can now be used on all std algorithms operating on containers of mesh entities

 @author Luca Formaggia
 */
template < typename MeshEntity,
         typename ComparisonPolicy = boost::function2<bool, flag_Type, flag_Type> >
class EntityFlagInterrogator
{
public:
    typedef ComparisonPolicy comparisonPolicy_Type;

    EntityFlagInterrogator ( flag_Type flag,
                             ComparisonPolicy const& p = ComparisonPolicy ( &Flag::testOneSet ) ) :
        M_flag ( flag ), M_policy ( p ) {}

    bool operator() ( const MeshEntity& entity ) const
    {
        return M_policy ( entity.flag(), M_flag );
    }

private:
    const flag_Type M_flag;
    const ComparisonPolicy M_policy;
};

//! A simple predicate to test the marker ID flag
/*!
   @prerequisite MeshEntity must have a method markerID_Type markerID();

    The ComparisonPolicy passed as (possible) second template parameter must be a functor
    capable of being constructed from a bool (*)(markerID_Type const&, entityflag_Type const &)
    and so that bool operator()(markerID_Type const&, markerID_Type const &) is defined.
    By default we use the std::equal_to functor

    Usage: if you want a predicate that tests if an markerID is equal to a given flag MYFLAG
    you create a object of type

@code
    EntityMarkerIDInterrogator<face_Type,ComparisonPolicy>(MYFLAG,mycomparisonpolicy())
@endcode

    which can now be used on all std algorithms operating on containers of mesh entities

 @author Luca Formaggia
 */
template < typename MeshEntity,
         typename ComparisonPolicy = boost::function2< bool, markerID_Type, markerID_Type> >
class EntityMarkerIDInterrogator
{
public:
    typedef ComparisonPolicy comparisonPolicy_Type;

    EntityMarkerIDInterrogator ( markerID_Type flag, ComparisonPolicy const& p = std::equal_to<markerID_Type>() ) :
        M_flag ( flag ), M_policy ( p ) {}

    bool operator() ( const MeshEntity& entity ) const
    {
        return M_policy ( entity.markerID(), M_flag );
    }

private:
    markerID_Type const M_flag;
    ComparisonPolicy const M_policy;
};
/**@} */
} // End Predicates

// MeshEntityContainer
/*!
    @author Luca Formaggia

    The class is a wrap up of Standard Library vector class.
    Its role is to held meshEntities.
    Its old name VectorSimple has been changed to MeshEntityContainer
    Any other use of class is deprecated and it should be replaced by std::vector<T>

 */

template <typename DataType, class Allocator = std::allocator<DataType> >
class MeshEntityContainer : public std::vector<DataType, Allocator>
{
public:

    //! @name Public Types
    //@{

    typedef DataType                                data_Type;
    typedef std::vector<DataType, Allocator>        vector_Type;
    typedef typename vector_Type::size_type         size_type;
    typedef typename vector_Type::value_type        value_type;
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
    MeshEntityContainer() : vector_Type() {}

    //! Constructor
    /*!
        @param vectorSize size of the vector
     */
    explicit MeshEntityContainer ( size_type vectorSize ) : vector_Type ( vectorSize ) {}

    //! Copy constructor
    /*!
        @param vector MeshEntityContainer vector to copy
     */
    MeshEntityContainer ( const MeshEntityContainer<DataType, Allocator>& vector );

    //! Constructor
    /*!
        Construct by copying a Standard Library vector
        @param vector Standard Library vector
     */
    explicit MeshEntityContainer ( const vector_Type& vector );

    //! Destructor
    ~MeshEntityContainer() {}

    //@}


    //! @name Operators
    //@{

    //! Assignment operator
    /*!
        Copies source MeshEntityContainer vector into "this"
        @param vector MeshEntityContainer vector
        @return Reference to a new MeshEntityContainer vector with the same
                content of MeshEntityContainer vector
     */
    MeshEntityContainer<DataType, Allocator>& operator= ( const MeshEntityContainer<DataType, Allocator>& vector );

    //! Access operator
    /*!
        Example: a(10)=90; // a[10] will contain 90.0
        @param i index of the element of the MeshEntityContainer vector
        @return a vector reference
        @note deprecated
        @todo Eliminate!
     */
    reference operator() ( size_type const i )
    {
        return ( this->operator[] ( i ) );
    }

    //! Const access operator
    /*!
        Example: a(10)=90; // a[10] will contain 90.0
        @param i index of the element of the MeshEntityContainer vector
        @return a vector const reference
        @note deprecated
        @todo Eliminate!
     */
    const_reference operator() ( size_type const i ) const
    {
        return ( this->operator[] ( i ) );
    }

    //@}


    //! @name Methods
    //@{

    //! Completely clear out the container, returning memory to the system
    inline void destroyData();

    //! trims the container so that capacity almost equals size
    void trim()
    {
        MeshEntityContainer<DataType, Allocator> (*this).swap ( *this );
    }

    //! It sets the capacity of the container
    //! It returns a bool to allow to test whether the container data pool
    //! has changed
    //! since in this case pointers to the container elements are invalidated
    //! It does not change the value of the elements currently stored, nor the container size
    //! @param size the new capacity of the container
    //! @return a bool which is true if the data pool has changed
    bool setMaxNumItems ( UInt size );

    /** @name Extractors
     *  Utilities to extract from the SimpleVector. They rely on std::algorithms
     */
    //@{
    /** General extractor.
     *  It extracts the a vector of pointers to the stored entities for which a
     *  predicate: i.e. a pointer to function with signature
     *  @code
     *  bool predicate(DataType const &)
     *  @endcode
     *  or functor object with
     *  @code
     *  bool operator()(DataType const &)
     *  @endcode
     *  returns true
     *
     *  The template parameter is the Predicate type (automatically deduced)
     *  @code
     *  aPredicate p; // A certain predicate
     *  vector<face_Type *> theList=mesh.faceList.extractIdAccordingToPredicate(p);
     *  @endcode
     *
     *  @param p The predicate
     *  @return A vector of pointers to constant mesh entities of the same type
     *  of those in the container
     */
    template<typename Predicate>
    std::vector<DataType const*> extractAccordingToPredicate ( Predicate const& p ) const;
    /** Entity Counter.
     *  It returns the number of stored entities for which a predicate
     *  returns true
     *  The template parameter is the Predicate type (automatically deduced): a function object or
     *  a pointer to function taking a constant reference to a mesh entity and returning a bool
     *
     *  @param p  The predicate
     *  @return   Number of entities in the container satisfying the predicate
     */

    template<typename Predicate>
    UInt countAccordingToPredicate ( Predicate const& p ) const;
    /** @brief It extracts all elements that satisfy a condition on the flag
     *
     *  It operates only container elements  where the method
     *  flag_Type flag()  is defined.
     *  Examples:
     *
     *  @code
     *  #include "MeshEntity.hpp"
     *  v=sv.extractElementsWithFlag(PHYSICAL_BOUNDARY,Flag::testOneSet)
     *  @endcode
     *  will extracts all elements in the boundary
     *  @code
     *  v=sv.extractElementsWithFlag(PHYSICAL_BOUNDARY|INTERNAL_INTERFACE, testOneSet)
     *  @endcode
     *  will extracts all elements on the boundary or on an internal interface
     *
     *  @param refFlag the flag against which the test is made
     *  @param policy. A functor/function pointer which implements
     *  bool policy(const flag_Type & inputFlag, const flag_Type & refFlag)
     *  Available policies are Flag::testOneSet and Flag::testAllSet (defined in Lifev.hpp)
     *  @return a vector containing the pointers to constant mesh entities in the container
     *  whose flag is set as that of refFlag according to the given policy
     */
    template<typename Policy>
    std::vector<DataType const*>
    extractElementsWithFlag ( const flag_Type& refFlag, Policy const& policy = &Flag::testOneSet ) const;

    /*! @brief It counts all elements that satisfy a condition on the flag
     *
     *  It only operates on container elements where the method
     *  flag_Type flag() is defined.
     *  Examples:
     *
     *  @code
     *  #include "MeshEntity.hpp"
     *  Uint i=sv.countElementsWithFlag(PHYSICAL_BOUNDARY,Flag::testOneSet)
     *  @endcode
     *  will count all elements on the boundary
     *
     *  @param refFlag the flag against which the test is made
     *  @param policy. A functor/function pointer which implements
     *  bool policy(const flag_Type & inputFlag, const flag_Type & refFlag)
     *  Available policies testOneSet and testAllSet (defined in Lifev.hpp)
     *  @return an unsigned integer
     */
    template<typename Policy>
    UInt countElementsWithFlag ( const flag_Type& refFlag,
                                 Policy const& policy = &Flag::testOneSet ) const;
    //@}

    /*! @brief It counts all elements with a markerID
     *
     *  @param markerID the markerID against which the test is made
     *  @param policy. A functor/function pointer which implements
     *  bool policy(const markerID_Type & inputMarkerID, const markerID_Type & refMarkerID)
     *  By default it uses std::equal_to
     *  @return an unsigned integer
     */
    template<typename Policy>
    UInt countElementsWithMarkerID ( const markerID_Type& markerID,
                                     const Policy& policy = std::equal_to<markerID_Type>() ) const;
    //@}

    /** @name Changers
      * Utilities that change the elements according to policies
      */
    //@{
    /** Set element with certain flag set first
     *
     * It reorders the container starting from the given offset up to the end of the container
     * so that elements with the given flag set are first. A policy must be passed.
     * Typically either testOneSet or testAllSet. A policy object must have a method
     * @code
     * bool operator()(DataType const & i, DataType const & r)
     * @endcode
     * which compares the boolean flag i with the reference flag r.
     * It fixes the id of the stored entities and it returns the newToOld array
     * in case we need to fix some related id
     *
     * @param refFlag the reference flag against which comparison is made.
     * @param offset We reorder only elements from offset to the end of the container.
     * @param policy The policy we use for comparing the flag of the mesh entity with the
     * reference flag.
     * @return A pair with the new offset, i.e the starting point of the elements for which predicate
     * was false, and a  NewToOld vector with the new->old position. Entity now in position
     * i was before in  NewToOld[i].
     *
     */
    template<typename Policy>
    std::pair<UInt, std::vector<ID> > reorderAccordingToFlag ( const flag_Type& refFlag,
                                                               Policy const& policy = &Flag::testOneSet,
                                                               UInt offset = 0 );

    //! Changes content according to a given functor
    /*!
     * This method is here for people which do not remember how std::for_each works.
     * It takes as argument a functor that must have the method
     * @code
     * void operator()(DataType & d)
     * @endcode
     * and which changes d according to the user will.
     * @param fun The functor which implements the change
     *
     */
    template<typename Functor>
    void changeAccordingToFunctor ( Functor const& fun )
    {
        std::for_each ( this->begin(), this->end(), fun );
    }

    //! Resets the id of the mesh entities so that it matches the position in the container
    /**
     * After its call the id of each mesh entity matches its position in the container.
     * The mesh entity must have the method ID & id()
     * @return a vector containing the association new->old. The entity that is now in position id
     * was originally in position newtoold[id]
     */
    std::vector<UInt> resetId()
    {
        std::vector<ID> newToOld;
        newToOld.reserve ( this->size() );
        iterator a ( this->begin() ) ;
        for ( UInt i = 0; i < this->size(); ++i )
        {
            ID old = a->localId();
            (a++)->setLocalId ( i );
            newToOld.push_back ( old );
        }
        return newToOld;
    }
    //@}
};

namespace Utilities
{
/** @defgroup Utilities
 * Some general Utilities on MeshEntityContainer
 * @{
 */

/** Reorder according to a permutation vector.
 *  We can reorder a MeshEntityContainer by passing a permutation vector
 *  std::vector<ID> newToOld;
 *
 *  The entity originally in position newtoold[id] will be moved to position id
 *  All ids (local id!) will be renumbered to reflect the new position in the
 *  container.
 *
 *  @pre The permutation vector must be of the right size
 *  @pre The permutation vector must be a good permutation vector
 */
template <typename EntityContainer >
void reorderAccordingToIdPermutation ( EntityContainer& container, std::vector<ID> const& newToOld )
{
    ASSERT_BD ( newToOld.size() >= container.size() );
    typedef typename EntityContainer::iterator it;
    typedef typename EntityContainer::value_type meshEntity_Type;
    std::vector<ID>::const_iterator start = newToOld.begin();
    // Change the id's
    for ( it i = container.begin(); i < container.end(); ++i)
    {
        i->setLocalId ( * (start++) );
    }
    // Fix the ordering
    std::sort ( container.begin(),
                container.end(),
                Comparers::CompareAccordingToLocalId<meshEntity_Type, std::less<ID> >() );
}

/** Fix pointers after permutation
 *  If a mesh entity contains pointers to other mesh entities (typically Points), after the renumbering of the
 *  referenced mesh entity (for instance using reorderAccordingToPermutation on Points) the address stored in
 *  the pointers will be wrong since it refers to the old numbering! This routine fixes it.
 *  It has to be called AFTER the reordering of the references mesh entities (i.e. the points) and NOT before.
 *
 *  Example:
 *  @code
 *  reorderAccordingToIdPermutation(mesh.points(),newToOld);
 *  // Now all entities storing pointers to points are invalid!!
 *  fixAfterPermutation(mesh.faces(),mesh.points(),newToOld);
 *  //FIxed!

 *  @endcode
 */
template <typename EntityContainer, typename RefEntityContainer >
void fixAfterPermutation ( EntityContainer& container,
                           RefEntityContainer const& refcontainer,
                           std::vector<ID> const& newToOld )
{
    // I need to build the inverse permutation
    std::vector<ID> oldToNew ( newToOld.size() );
    for ( UInt i = 0; i < newToOld.size(); ++i )
    {
        oldToNew[ newToOld[ i ] ] = i;
    }
    typedef typename EntityContainer::iterator it;
    typedef typename EntityContainer::value_type meshEntity_Type;
    const UInt numPoints = meshEntity_Type::geoShape_Type::S_numPoints;
    for ( it i = container.begin(); i < container.end(); ++i )
    {
        for ( UInt j = 0; j < numPoints; ++j )
        {
            ID oldaddress = i->point ( j ).localId();
            ID newaddress = oldToNew[ oldaddress ];
            i->setPoint ( j, & ( refcontainer[ newaddress ] ) );
        }
    }
}

///! Fix pointers after shallow copy
/*!
 *  If a mesh entity contains pointers to other mesh entities (typically Points), after a shallow copy, for instance
 *  by the automatic copy constructor, the pointers may still refer to a wrong list of Points!
 *  This utility allows to change the situation and make the pointers to point to the new list of Points.
 *  By this utility the deep copy of a RegionMesh is now possible!
 *
 *  Example:
 *  @code
 *  PointList newPoints(mesh.pointList); // deep copy since list stores Point(s)
 *  FaceList  newFaces(mesh.faceList); // This is a shallow copy since we store Point*
 *  fixAfterShallowCopy(newFaces,newPoints); now the pointers point to newPoint
 *  @endcode
 *
 *  @param container the entity container with the wrong pointers to Point
 *  @param newPointContainer the container with the list new Points
 *  @pre container must contain pointers to valid points
 *  @note For efficiency reason no consistency checks are made
 */
template <typename EntityContainer, typename PointContainer >
void fixAfterShallowCopy ( EntityContainer& container, PointContainer const& newPointContainer )
{
    typedef typename EntityContainer::iterator it;
    typedef typename EntityContainer::value_type meshEntity_Type;
    const UInt numPoints = meshEntity_Type::geoShape_Type::S_numPoints;

    for ( it i = container.begin(); i < container.end(); ++i )
        for ( UInt j = 0; j < numPoints; ++j )
        {
            i->setPoint ( j, & ( newPointContainer[ i->point ( j ).localId() ] ) );
        }
}
/* @}*/
}// End namespace Utilities



//============================================================================
//                      IMPLEMENTATIONS

// Constructors
//============================================================================
template <typename DataType, class Allocator>
MeshEntityContainer<DataType, Allocator>::MeshEntityContainer ( const MeshEntityContainer<DataType, Allocator>& vector )
    : vector_Type ( vector )
{}


//============================================================================
// Operators
//============================================================================
template <typename DataType, class Allocator>
MeshEntityContainer<DataType, Allocator>&
MeshEntityContainer<DataType, Allocator>::operator= ( const MeshEntityContainer<DataType, Allocator>& vector )
{
    vector_Type::operator= ( vector );
    return *this;
}


//============================================================================
// Methods
//============================================================================
template <typename DataType, class Allocator>
void
MeshEntityContainer<DataType, Allocator>::destroyData()
{
    this->clear();
    this->vector_Type::swap ( vector_Type() );
}

template<typename DataType, class Allocator>
bool MeshEntityContainer<DataType, Allocator>::setMaxNumItems ( UInt size )
{
    bool _check = ( this->capacity() ) < size;
    this->reserve ( size );
    return _check;
}

template<typename DataType, class Allocator>
template<typename Predicate>
std::vector<DataType const*>
MeshEntityContainer<DataType, Allocator>::extractAccordingToPredicate ( Predicate const& p ) const
{
    UInt howmany = this->countAccordingToPredicate ( p );
    std::vector<DataType const*> tmp;
    tmp.reserve ( howmany );
    for ( const_iterator i = this->begin(); i < this->end(); ++i )
        if ( p (*i) )
        {
            tmp.push_back ( & (*i) );
        }
    return tmp;
}

template<typename DataType, class Allocator>
template<typename Predicate>
UInt MeshEntityContainer<DataType, Allocator>::countAccordingToPredicate ( Predicate const& p ) const
{
    return std::count_if ( this->begin(), this->end(), p );
}

template<typename DataType, class Allocator>
template<typename Policy>
std::vector<DataType const*>
MeshEntityContainer<DataType, Allocator>::extractElementsWithFlag ( const flag_Type& refFlag,
                                                                    const Policy& policy ) const
{
    Predicates::EntityFlagInterrogator<DataType> interrogator ( refFlag, policy );
    return this->extractAccordingToPredicate ( interrogator );
}

template<typename DataType, class Allocator>
template<typename Policy>
UInt MeshEntityContainer<DataType, Allocator>::countElementsWithFlag ( const flag_Type& refFlag,
                                                                       const Policy& policy ) const
{
    Predicates::EntityFlagInterrogator<DataType> interrogator ( refFlag, policy );
    return this->countAccordingToPredicate ( interrogator );
}

template<typename DataType, class Allocator>
template<typename Policy>
UInt MeshEntityContainer<DataType, Allocator>::countElementsWithMarkerID ( const markerID_Type& markerID,
                                                                           const Policy& policy ) const
{
    Predicates::EntityMarkerIDInterrogator<DataType> interrogator ( markerID, policy );
    return this->countAccordingToPredicate ( interrogator );
}

template<typename DataType, class Allocator>
template<typename Policy>
std::pair<UInt, std::vector<ID> >
MeshEntityContainer<DataType, Allocator>::reorderAccordingToFlag ( const flag_Type& refFlag,
                                                                   Policy const& policy,
                                                                   UInt offset )
{
    iterator last = std::stable_partition ( this->begin() + offset,
                                            this->end(),
                                            Predicates::EntityFlagInterrogator<DataType> ( refFlag, policy ) );
    std::vector<ID> newToOld = this->resetId();
    return std::make_pair ( std::distance ( this->begin(), last), newToOld );
}
}/// end of LifeV namespace

#endif /* _MESHENTITYCONTAINER_HH_ */
