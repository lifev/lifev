/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-13

  Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano

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
   \file typeInfo.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-13
 */
#ifndef __TypeInfo_H
#define __TypeInfo_H 1

#include <typeinfo>

namespace LifeV
{
/*!
  \class TypeInfo
  \brief wrapper for std::type_info

  \sa SFactory, type_info

  @author Christophe Prud'homme
  @version $Id: typeInfo.hpp,v 1.1 2004-10-13 10:18:52 prudhomm Exp $
*/
class TypeInfo
{
public:


    /** @name Typedefs
     */
    //@{

    //@}

    /** @name Constructors, destructor
     */
    //@{

    TypeInfo();
    TypeInfo(const std::type_info&); // non-explicit
    TypeInfo( TypeInfo const & );
    ~TypeInfo();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    //! Access for the wrapped \c std::type_info
    const std::type_info& typeInfo() const;

    //! get the \c name()
    const char* name() const;

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //! Compatibility functions
    bool before(const TypeInfo& rhs) const;


    //@}

private:

    const std::type_info* _M_info;
};

inline bool operator==(const TypeInfo& lhs, const TypeInfo& rhs)
{ return lhs.typeInfo() == rhs.typeInfo(); }

inline bool operator<(const TypeInfo& lhs, const TypeInfo& rhs)
{ return lhs.before(rhs); }

inline bool operator!=(const TypeInfo& lhs, const TypeInfo& rhs)
{ return !(lhs == rhs); }

inline bool operator>(const TypeInfo& lhs, const TypeInfo& rhs)
{ return rhs < lhs; }

inline bool operator<=(const TypeInfo& lhs, const TypeInfo& rhs)
{ return !(lhs > rhs); }

inline bool operator>=(const TypeInfo& lhs, const TypeInfo& rhs)
{ return !(lhs < rhs); }

}// end namespace LifeV

#endif /* __TypeInfo_H */
