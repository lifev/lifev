//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file
     @brief This file contains mother classes of the QRAdapter classes.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef QR_ADAPTER_BASE_HPP
#define QR_ADAPTER_BASE_HPP


#include <lifev/core/LifeV.hpp>

namespace LifeV
{

/*!
  ETCurrenteFE is a template class. If fieldDim the general
  case is treated as representing a vectorial FE (only the case
  where fieldDim represents a scalar FE, but this is a partial
  specialization of this class).

*/
template < typename ImplementationType>
class QRAdapterBase
{
public:

    //! @name Public Types
    //@{

    typedef ImplementationType implementation_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty and only constructor
    QRAdapterBase() {}

    //! Destructor
    virtual ~QRAdapterBase() {}

    //@}


    //! @name Methods
    //@{

    //! Method for accessing the actual implementation contained.
    const implementation_Type& implementation() const
    {
        return static_cast<const implementation_Type&> (*this);
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No copy (avoid slicing)
    QRAdapterBase ( const QRAdapterBase<ImplementationType>&);

    //! No assignement (avoid slicing)
    void operator= (const QRAdapterBase<ImplementationType>&);

    //@}
};


} // Namespace LifeV


#endif /* QR_ADAPTER_BASE_HPP */
