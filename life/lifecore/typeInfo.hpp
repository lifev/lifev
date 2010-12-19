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
  @file
  @brief Type information class

  @date 13-10-2004
  @author  Christophe Prud'homme <christophe.prudhomme@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef TYPEINFO_H
#define TYPEINFO_H 1

#include <typeinfo>

namespace LifeV
{

class TypeInfo
{
public:j4 check
    //! @name Constructors, destructor
    //@{
    TypeInfo();
    TypeInfo(const std::type_info&); // non-explicit
    TypeInfo( TypeInfo const & );
    virtual ~TypeInfo();
    //@}

    //! @name Methods
    //@{
    //! Compatibility functions
    bool before(const TypeInfo& rhs) const;
    //@}

    //! @name Get methods
    //@{
    //! Access for the wrapped \c std::type_info
    const std::type_info& typeInfo() const;
    //! get the \c name()
    const char* name() const;
    //@}
private:
    const std::type_info* M_info;
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

} // Namespace LifeV

#endif // TYPEINFO_H
