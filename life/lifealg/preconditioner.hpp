
/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2005-04-07

  Copyright (C) 2005 EPFL

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
/**
   \file preconditioner.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-04-07
 */
#ifndef __Preconditioner_H
#define __Preconditioner_H 1


/*!
   \class Preconditioner

   The Preconditioner object performs a preconditioning operation
   based on vector x and stores the result in vector y. The
   \c transSolve() method only need be defined when the preconditioner is
   used with an iterative solver that requires it.

   A preconditioner should define some types using typedef
   directive like this:

   \verbatim
   typedef typename ArrayType::T_numtype T_numtype
   \endverbatim

   Some each preconditioner should provide at least these typedefs:

   \arg numerical type: \c T_numtype
   \arg array type: \c T_arraytype
   \arg preconditioner type: \c T_prectype

   @author Christophe Prud'homme <Christophe.Prudhomme@ann.jussieu.fr>
   @see SolverIterative
  */
template<typename VectorX, typename VectorY=VectorX>
class Preconditioner
{
public:

  /** @name Constructors, destructors and methods
   */
  //@{
  /** base virtual destructor */
  virtual ~Preconditioner() {}

  /** \f[$y \leftarrow M^{-1} x\f$ */
  virtual void solve(const VectorX& x, VectorY& y) const = 0;

  /** \\f$y \leftarrow M^{-T} x\f$ */
  virtual void transSolve(const VectorX& x, VectorY& y) const = 0;
  //@}
};

/*!
  \class PreconditionerAdaptor

  Base class for preconditioner adaptor objects. This is the class used to enable the Adaptor
  Pattern.

  @author Christophe Prud'homme <Christophe.Prudhomme@ann.jussieu.fr>

  @see Gamma et~al.(1995),
  Design Patterns,
  Publisher: Addison Wesley Professional Computing Series. Addison Wesley, 1995.
*/
template <class VecX, class VecY>
class PreconditionerAdaptor
   :
   public Preconditioner<VecX, VecY>
{
public:


};

template<typename VectorX, typename VectorY=VectorX>
class PreconditionerIdentity
    :
    Preconditioner<VectorX, VectorY>
{
public:


    /** @name Constructors, destructors and methods
     */
    //@{

    PreconditionerIdentity()
        :
        Preconditioner<VectorX, VectorY>()
        {}

    /** base virtual destructor */
    ~PreconditionerIdentity()
        {}

    /** \f[$y \leftarrow M^{-1} x\f$ */
    virtual void solve(const VectorX& x, VectorY& y) const
        {
            y.assign( x );
        }

    /** \\f$y \leftarrow M^{-T} x\f$ */
    virtual void transSolve(const VectorX& x, VectorY& y) const
        {
            y.assign( x );
        }
    //@}
};


#endif /* __Preconditioner_H */



