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
   \file givens.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-04-07
 */
#ifndef __GIVENS_H
#define __GIVENS_H 1

namespace LifeV
{
/**
   \class GivensRotation
   \brief Givens Plane Rotation

   Input \c a and \c b to the constructor to create a givens plane rotation
   object. Then apply the rotation to two vectors. There is a
   specialization of the givens rotation for complex numbers.

   \verbatim
   [  c  s ] [ a ] = [ r ]
   [ -s  c ] [ b ]   [ 0 ]
   \endverbatim

   -#) the addition operator must be defined for \cT
   -#) the multiplication operator must be defined for \cT
   -#) the division operator must be defined for \c T
   -#) the \c std::abs() function must be defined for \c T

   @author Christophe Prud'homme <Christophe.Prudhomme@epfl.ch>
   @see
*/
template <typename T>
class GivensRotation
{
public:


    /** @name Typedefs
     */
    //@{

    typedef T value_type;
    //@}


    /** @name Constructor
     */
    //@{
    //! Default constructor
    GivensRotation() : _M_a(0), _M_b(0), _M_cos(0), _M_sin(0)
        {
        }

    //! Givens Plane Rotation Constructor
    GivensRotation(value_type _M_ain, value_type _M_bin)
        {
            value_type roe;
            if (std::abs(_M_ain) > std::abs(_M_bin))
                roe = _M_ain;
            else
                roe = _M_bin;

            value_type scal = std::abs(_M_ain) + std::abs(_M_bin);
            value_type r, z;
            if (scal != value_type(0))
            {
                value_type _M_ascl = _M_ain / scal;
                value_type _M_bscl = _M_bin / scal;
                r = scal * sqrt(_M_ascl * _M_ascl + _M_bscl * _M_bscl);
                if (roe < value_type(0))
                    r *= -1;
                _M_cos = _M_ain / r;
                _M_sin = _M_bin / r;
                z = 1;
                if (std::abs(_M_ain) > std::abs(_M_bin))
                    z = _M_sin;
                else if (std::abs(_M_bin) >= std::abs(_M_ain) && _M_cos != value_type(0))
                    z = value_type(1) / _M_cos;
            }
            else
            {
                _M_cos = 1;
                _M_sin = 0;
                r = 0;
                z = 0;
            }
            _M_a = r;
            _M_b = z;
        }
    //@}


    void set_cs(value_type cin, value_type sin)
        {
            _M_cos = cin;
            _M_sin = sin;
        }

    //! Apply plane rotation to two vectors.
    template <
        class VecX,
        class VecY
        >
    void apply(VecX x, VecY y)
        {
            LIFEV_ASSERT(x.size() <= y.size()).error("invalid size");

            typename VecX::iterator xi = x.begin();
            typename VecX::iterator xend = x.end();
            typename VecY::iterator yi = y.begin();

            while (xi != xend)
            {
                apply(*xi, *yi);
                ++xi;
                ++yi;
            }
        }
    //! Apply plane rotation to two real scalars.
    void apply(value_type& x, value_type& y)
        {
            value_type tmp = _M_cos * x + _M_sin * y;
            y = _M_cos * y - _M_sin * x;
            x = tmp;
        }

    value_type a() { return _M_a; }
    value_type b() { return _M_b; }
    value_type c() { return _M_cos; }
    value_type s() { return _M_sin; }

protected:
    value_type _M_a;
    value_type _M_b;
    value_type _M_cos;
    value_type _M_sin;
};

}

#endif /* __GIVENS_H */
