/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.
  
  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/* -------------------------------------------------------------------------*/
/*!
  \file vecUnknown.h

  Vector classes

  #purpose: provides vector classes handler useful in case of solving
   vector problem.
   Provides Vector and VectorBlock class necessary for IML++ library
   Alain Gauthier.
*/
#ifndef _VEC_UNKNOWN_HH
#define _VEC_UNKNOWN_HH

#include <vector>
#include "tab.hpp"

/*! \class PhysVectUnknown
  vector unknown which has the dimension of the physical domain
   The type VectorType could be a Vector, or a VectorBlock or a
   vector<double> depending on the choice of the linear system solver and
   on the case of scalar or vectorial problem 
*/
template<typename VectorType>
class PhysVectUnknown
    :
    public VectorType

{
    UInt _size;
    static const UInt _nbcomp= nDimensions;
public:

    typedef VectorType super;

    //  PhysVectUnknown(){}
    PhysVectUnknown(UInt const Ndof);
    PhysVectUnknown(PhysVectUnknown<VectorType> const &RhPhysVectUnknown);

    PhysVectUnknown& operator=( double const __val )
	{
	    super & __super = (*this) ;
	    __super = __val;
	    return *this;
	}

    PhysVectUnknown& operator=( VectorType const& __v )
	{
	    if ( this == &__v )
		return *this;
	    super & __super = (*this) ;
	    __super = __v;
	    return *this;
	}

    //! gives the front of the vector
    inline Real * giveVec() {return &((*this)[0]);}
    inline UInt size() const {return _size;}
    inline UInt nbcomp() const {return _nbcomp;}
};

//---------------------------------------------------------------//

/*! \class ScalUnknown
   scalar unknown of dimension=1
   The type VectorType could be a Vector, or a VectorBlock or a
   vector<double> depending on the choice of the linear system solver and
   on the case of scalar or vectorial problem 
*/
template<typename VectorType>
class ScalUnknown
    :
    public VectorType
{
    UInt _size;
    static const UInt _nbcomp= 1;
public:

    typedef VectorType super;

    ScalUnknown(){};
    ScalUnknown(UInt const Ndof);
    ScalUnknown(const ScalUnknown<VectorType> &RhScalUnknown);

    ScalUnknown& operator=( double const __val )
	{
	    super & __super = (*this) ;
		__super = __val;
	    return *this;
	}
    ScalUnknown& operator=( VectorType const& __v )
	{
	    if ( this == &__v )
		return *this;
	    super & __super = (*this) ;
	    __super = __v;
	    return *this;
	}

    //! gives the front of the vector
    inline Real * giveVec(){return &((*this)[0]);}

    inline UInt size() const {return _size;}
    inline UInt nbcomp() const {return _nbcomp;}
};

//---------------------------------------------------------------//

/*! \class GenericVecHdl
   vector problem handler
   The type VectorType could be a Vector, or a VectorBlock or a
   vector<double> depending on the choice of the linear system solver and
   on the case of scalar or vectorial problem 
*/
template<typename VectorType>
class GenericVecHdl
    :
    public VectorType
{
    UInt _size;
    UInt _nbcomp;
public:

    typedef VectorType super;

    GenericVecHdl(){}

    GenericVecHdl(UInt ex_size, UInt ex_nbcomp):
	super(ex_size),_size(ex_size),_nbcomp(ex_nbcomp){};

    GenericVecHdl(const GenericVecHdl<VectorType> &RhGenVec);
    //! construction from a scalar unknown:
    GenericVecHdl(const ScalUnknown<VectorType> &RhScalUnknown);
    //! construction from two scalar unknowns:
    GenericVecHdl(const ScalUnknown<VectorType> &RhScU1,
		  const ScalUnknown<VectorType> &RhScU2);
    //! construction from a physical vectorial unknown:
    GenericVecHdl(const PhysVectUnknown<VectorType> &RhPhVU);
    //! construction from a physical vectorial unknown and a scalar unknown:
    GenericVecHdl(const PhysVectUnknown<VectorType> &RhPhVU,
		  const ScalUnknown<VectorType> &RhScU);
    GenericVecHdl(const ScalUnknown<VectorType> &RhScU,
		  const PhysVectUnknown<VectorType> &RhPhVU);
    //! construction from a physical vectorial unknown and a generic unknown:
    GenericVecHdl(const PhysVectUnknown<VectorType> &RhPhVU,
		  const GenericVecHdl<VectorType> &RhGenVec);
    GenericVecHdl(const GenericVecHdl<VectorType> &RhGenVec,
		  const PhysVectUnknown<VectorType> &RhPhVU);

    //! gives the front of the vector
    inline Real * giveVec(){return &((*this)[0]);}
    inline UInt size() const {return _size;}
    //! gives the size of one block
    inline UInt nbcomp() const {return _nbcomp;}

  
};

//---------------------------------------------------------------//

/////////////////////////////////////////////////////////////
/*! \class Vector
    Class Vector necessary for IML++ library
    Modification of RNM vector: It adds the default
    constructor.
    15/11/01
*/ 
///////////////////////////////////////////////////////////

class Vector
{
  Tab1d _v;
 public:
  //ALAIN : problem with assignment and default constructor...
  // assignment of vector of different sizes !
  Vector():_v(0){}   //!< Default constructor we need for IML++
  Vector(Vector const &ex_v):_v(ex_v.v()){}
  Vector(UInt n):_v((int)n){_v= 0.0;}
  inline Tab1dView v() const {return _v;}
  inline UInt size() const {return _v.size();}
  //operators
  inline Vector operator=(Vector const &ex_v)
    {
      if (&ex_v != this)
	{
	  if (_v.size()==0)
	    {
	      this->~Vector();
	      new (this) Vector(ex_v);
	    }
	  else
	    _v=ex_v.v();
	};
      return *this;
    }
  inline Vector operator=(double const val)
    {
      for (int i=0; i<_v.N(); i++)
	_v(i)=val;
      return *this;
    }
  inline Vector operator-=(Vector const &ex_v)
    {
      _v-=ex_v.v();
      return *this;
    }
  inline Vector operator+=(Vector const &ex_v)
    {
      _v+=ex_v.v();
      return *this;
    }
  inline double & operator()(int i) const
    {
      return _v(i);
    }
  inline double & operator[](int i) const
    {
      return _v(i);
    }

  friend Vector operator+(Vector const &ex_v1, Vector const &ex_v2);
  friend Vector operator-(Vector const &ex_v1, Vector const &ex_v2);
  friend Vector operator*(Vector const &ex_v, double const val);
  friend Vector operator*(double const val, Vector const &ex_v);

  //! vector inner product
  friend double dot(Vector const &ex_v1, Vector const &ex_v2);
  //! norm derived from dot:
  friend inline double norm(Vector const &ex_v)
    {return sqrt(dot(ex_v,ex_v));}
};

//---------------------------------------------------------------//

/*!\class VectorBlock
  Version for block matrices handling.
  19/11/01
*/
class VectorBlock
{
    std::vector<Tab1d> _v; //!< container of block vector
 public:
  // Default constructor we need for IML++
  VectorBlock(){}
  VectorBlock(int n, int blsize=1) //!< default size block =1 !
    {
      Tab1d bloc(blsize);
      bloc=0.0;
      _v.resize(n,bloc);
    }
  VectorBlock(VectorBlock const &ex_v){*this = ex_v;}
  inline vector<Tab1d> v() const {return _v;}
  inline UInt size() const {return _v.size();}
  //operators
  inline VectorBlock & operator=(VectorBlock const &ex_v)
    {
      if (&ex_v != this)
	{
	  _v.clear();
	  _v.reserve(ex_v.v().size());
	  for (UInt ir=0; ir < ex_v.v().size(); ir++)
	    _v.push_back(ex_v.v()[ir]);
	}
      return *this;
    }
  inline VectorBlock operator=(double const val)
    {
      for (vector<Tab1d>::iterator ip=_v.begin(); ip < _v.end(); ip++)
	*ip=val;
      return *this;
    }
  inline VectorBlock operator-=(VectorBlock const &ex_v)
    {
      for (UInt ir=0; ir < ex_v.v().size(); ir++)
	_v[ir]-=ex_v.numBlock(ir);
      return *this;
    }
  inline VectorBlock operator+=(VectorBlock const &ex_v)
    {
      for (UInt ir=0; ir < ex_v.v().size(); ir++)
	_v[ir]+=ex_v.numBlock(ir);
      return *this;
    }
  inline double & operator()(int i) const
    {
      int blsize=_v[0].N();
      int bloc_i= i/blsize;
      int loc_i=i%blsize;
      return _v[bloc_i](loc_i);
    }
  inline double & operator[](int i) const
    {
      int blsize=_v[0].N();
      int bloc_i= i/blsize;
      int loc_i=i%blsize;
      return _v[bloc_i](loc_i);
    }
  inline Tab1d & numBlock(int i)
    {return _v[i];}
  inline const Tab1d & numBlock(int i) const
    {return _v[i];}

  friend VectorBlock operator+(VectorBlock const &ex_v1, VectorBlock
			       const &ex_v2);
  friend VectorBlock operator-(VectorBlock const &ex_v1, VectorBlock
			       const &ex_v2);
  friend VectorBlock operator*(VectorBlock const &ex_v, double const val);
  friend VectorBlock operator*(double const val, VectorBlock const &ex_v);

  //! vector inner product
  friend double dot(VectorBlock const &ex_v1, VectorBlock const &ex_v2);
  //! norm derived from dot:
  friend inline double norm(VectorBlock const &ex_v)
    {return sqrt(dot(ex_v,ex_v));}

};

//---------------------------------------------------------------//
// Useful function

//! assign the values of a pointer to an existing vector
//! No check of the dimensions !!!
template<typename VectorType>
void
point2Vector(double const * point, VectorType & v)
{
  VectorType ans(v.size());
  for (UInt i=0; i < ans.size(); i++)
    {
      ans[i]= *point;
      point++;
    }
  v=ans;
}


/*-------------------------------------------------------/
  IMPLEMENTATION
/-------------------------------------------------------*/

//////////////////////////////
// class PhysVectUnknown
/////////////////////////////

template<typename VectorType>
PhysVectUnknown<VectorType>::PhysVectUnknown(UInt const Ndof)
    :
    super(nDimensions*Ndof),
    //@}
    _size(nDimensions*Ndof)
{
}

//! the case of VectorBlock type
PhysVectUnknown<VectorBlock>::PhysVectUnknown(UInt const Ndof);

template<typename VectorType>
PhysVectUnknown<VectorType>::PhysVectUnknown(PhysVectUnknown<VectorType> const &RhPhysVectUnknown)
    :
    super(RhPhysVectUnknown),
    _size(RhPhysVectUnknown.size())
{
}

//////////////////////////////
// class ScalUnknown
/////////////////////////////

template<typename VectorType>
ScalUnknown<VectorType>::ScalUnknown(UInt const Ndof)
    :
    super(Ndof),
    _size(Ndof)
{}

template<typename VectorType>
ScalUnknown<VectorType>::ScalUnknown(const ScalUnknown<VectorType> &RhScalUnknown)
    :
    super(RhScalUnknown), 
    _size(RhScalUnknown.size())
{}

//////////////////////////////
// class GenericVecHdl
/////////////////////////////

template<typename VectorType>
GenericVecHdl<VectorType>::GenericVecHdl(const GenericVecHdl<VectorType> &RhGenVec)
    :
    super(RhGenVec), 
    _size(RhGenVec.size()), 
    _nbcomp(RhGenVec.nbcomp()){}

//! construction from a scalar unknown:
template<typename VectorType>
GenericVecHdl<VectorType>::
GenericVecHdl(const ScalUnknown<VectorType> &RhScalUnknown)
    :
    super(RhScalUnknown),
    _size(RhScalUnknown.size()),
    _nbcomp(RhScalUnknown.nbcomp())
{}

//! construction from two scalar unknowns:
template<typename VectorType>
GenericVecHdl<VectorType>::GenericVecHdl(const ScalUnknown<VectorType> &RhScU1, const ScalUnknown<VectorType> &RhScU2 )
    :
    super(RhScU1.size()+RhScU2.size()),
    _size(RhScU1.size()+RhScU2.size()),
    _nbcomp(RhScU1.nbcomp()+RhScU2.nbcomp())
{
  for (UInt i=0; i<RhScU1.size(); ++i) (*this)[i]=RhScU1[i];
  for (UInt i=0; i<RhScU2.size(); ++i) (*this)[i+RhScU1.size()]=RhScU2[i];
}
//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::GenericVecHdl(const ScalUnknown<VectorBlock> &RhScU1,
					  const ScalUnknown<VectorBlock> &RhScU2);


//! construction from a physical vectorial unknown:
template<typename VectorType>
GenericVecHdl<VectorType>::GenericVecHdl(const PhysVectUnknown<VectorType> &RhPhVU)
    :
    super(RhPhVU),
    _size(RhPhVU.size()),
    _nbcomp(RhPhVU.nbcomp())
{}

//! construction from a physical vectorial unknown and a scalar unknown:
template<typename VectorType>
GenericVecHdl<VectorType>::GenericVecHdl(const PhysVectUnknown<VectorType> &RhPhVU,
					 const ScalUnknown<VectorType> &RhScU)
    :
    super(RhPhVU.size()+RhScU.size()),
    _size(RhPhVU.size()+RhScU.size()),
    _nbcomp(RhPhVU.nbcomp()+RhScU.nbcomp())
{
  for (UInt i=0; i<RhPhVU.size(); ++i) (*this)[i]=RhPhVU[i];
  for (UInt i=0; i<RhScU.size(); ++i) (*this)[RhPhVU.size()+i]=RhScU[i];
}
//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::
GenericVecHdl(const PhysVectUnknown<VectorBlock> &RhPhVU,
	      const ScalUnknown<VectorBlock> &RhScU);


template<typename VectorType>
GenericVecHdl<VectorType>::GenericVecHdl(const ScalUnknown<VectorType> &RhScU,
					 const PhysVectUnknown<VectorType> &RhPhVU)
    :
    super(RhPhVU.size()+RhScU.size()),
    _size(RhPhVU.size()+RhScU.size()),
    _nbcomp(RhPhVU.nbcomp()+RhScU.nbcomp())
{
  for (UInt i=0; i<RhScU.size(); ++i) (*this)[i]=RhScU[i];
  for (UInt i=0; i<RhPhVU.size(); ++i) (*this)[RhScU.size()+i]=RhPhVU[i];
}
//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::
GenericVecHdl(const ScalUnknown<VectorBlock> &RhScU, const PhysVectUnknown<VectorBlock> &RhPhVU);


//! construction from a physical vectorial unknown and a generic unknown:
template<typename VectorType>
GenericVecHdl<VectorType>::GenericVecHdl(const PhysVectUnknown<VectorType> &RhPhVU,
					 const GenericVecHdl<VectorType> &RhGenVec)
    :
    super(RhPhVU.size()+RhGenVec.size()),
    _size(RhPhVU.size()+RhGenVec.size()),
    _nbcomp(RhPhVU.nbcomp()+RhGenVec.nbcomp())
{
  for (UInt i=0;i<RhPhVU.size();++i) (*this)[i]=RhPhVU[i];
  for (UInt i=0;i<RhGenVec.size();++i) (*this)[RhPhVU.size()+i]=RhGenVec[i];
}
//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::GenericVecHdl(const PhysVectUnknown<VectorBlock> &RhPhVU,
					  const GenericVecHdl<VectorBlock> &RhGenVec);


template<typename VectorType>
GenericVecHdl<VectorType>::
GenericVecHdl(const GenericVecHdl<VectorType> &RhGenVec,
	      const PhysVectUnknown<VectorType> &RhPhVU):
  super(RhPhVU.size()+RhGenVec.size()),
  _size(RhPhVU.size()+RhGenVec.size()),
  _nbcomp(RhPhVU.nbcomp()+RhGenVec.nbcomp())
{
  for (UInt i=0;i<RhGenVec.size();++i) (*this)[i]=RhGenVec[i];
  for (UInt i=0;i<RhPhVU.size();++i) (*this)[i+RhGenVec.size()]=RhPhVU[i];
}
//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::
GenericVecHdl(const GenericVecHdl<VectorBlock> &RhGenVec,
	      const PhysVectUnknown<VectorBlock> &RhPhVU);


//////////////////////////////
// class Vector
/////////////////////////////

inline Vector operator+(Vector const &ex_v1, Vector const &ex_v2)
{
  ASSERT(ex_v1._v.N()==ex_v2._v.N(), " wrong dimension of vector");
  Vector ans(ex_v1._v.N());
  ans =  ex_v1;
  ans += ex_v2;
  return ans;
}
inline Vector operator-(Vector const &ex_v1, Vector const &ex_v2)
{
  ASSERT(ex_v1._v.N()==ex_v2._v.N(), " wrong dimension of vector");
  Vector ans(ex_v1._v.N());
  ans =  ex_v1;
  ans -= ex_v2;
  return ans;
}
inline Vector operator*(Vector const &ex_v, double const val)
{
  Vector ans(ex_v._v.N());
  ans._v = ex_v._v;
  ans._v *= val;
  return ans;
}
inline Vector operator*(double const val, Vector const &ex_v)
{
  return ex_v * val;
}
inline double dot(Vector const &ex_v1, Vector const &ex_v2)
{
  return (ex_v1._v,ex_v2._v);
}

//////////////////////////////
// class VectorBlock
/////////////////////////////

inline VectorBlock operator+(VectorBlock const &ex_v1, VectorBlock
			     const &ex_v2)
{
  ASSERT(ex_v1._v.size()==ex_v2._v.size(), " wrong dimension of vector");
  VectorBlock ans(ex_v1._v.size(),ex_v1.numBlock(0).N());
  ans=0.0;
  ans =  ex_v1;
  ans += ex_v2;
  return ans;
}
inline VectorBlock operator-(VectorBlock const &ex_v1, VectorBlock
			     const &ex_v2)
{
  ASSERT(ex_v1._v.size()==ex_v2._v.size(), " wrong dimension of vector");
  VectorBlock ans(ex_v1._v.size(),ex_v1.numBlock(0).N());
  ans=0.0;
  ans =  ex_v1;
  ans -= ex_v2;
  return ans;
}
inline VectorBlock operator*(VectorBlock const &ex_v, double const val)
{
  VectorBlock ans;
  ans = ex_v;
  for (UInt i=0; i < ex_v._v.size(); i++)
    ans.numBlock(i) *= val;
  return ans;
}
inline VectorBlock operator*(double const val, VectorBlock const &ex_v)
{
  return ex_v * val;
}
inline double dot(VectorBlock const &ex_v1, VectorBlock const &ex_v2)
{
  double ans=0;
  for (UInt i=0; i< ex_v1._v.size(); i++)
    ans+=(ex_v1.numBlock(i),ex_v2.numBlock(i));
  return ans;
}

#endif
