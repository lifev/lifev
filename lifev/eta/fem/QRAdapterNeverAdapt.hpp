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
    @brief Fake quadrature rule adapter for a constant quadrature rule.

    This class has a unique functionality: is returns the quadrature rule
    that it owns.

    This class has been designed in order to facilitate the implementation
    of the ETA framework (avoid duplicating the Integrate** classes). Remark
    that we define all the methods "inline" to allow the compiler to better
    optimize out the effect of this class.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 07 Aug 2012

 */

#ifndef QR_ADAPTER_NEVER_ADAPT_H
#define QR_ADAPTER_NEVER_ADAPT_H 1

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/fem/QRAdapterBase.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>


namespace LifeV
{


class QRAdapterNeverAdapt : public QRAdapterBase< QRAdapterNeverAdapt >
{
public:

    typedef QRAdapterBase< QRAdapterNeverAdapt > base_Type;

    //! @name Constructor & Destructor
    //@{

    //! Constructor with a quadrature rule.
    QRAdapterNeverAdapt (const QuadratureRule& qr) : base_Type(), M_qr (qr) {}

    //! Copy constructor
    QRAdapterNeverAdapt (const QRAdapterNeverAdapt& qrAdapter) : base_Type(), M_qr (qrAdapter.M_qr) {}

    //! Simple destructor
    ~QRAdapterNeverAdapt() {}

    //@}


    //! @name Methods
    //@{

    //! Update this structure with the current element
    /*!
      Usually, it does some computations, but here, nothing!
      This method is static, even if in general, this is not
      the case... We just give more information to the compiler
      so that it can better optimize the code.
     */
    static void update (UInt /*elementID*/) {}

    //@}


    //! @name Set Methods
    //@{

    //@}


    //! @name Get Methods
    //@{

    //! Do we need the adapted quadrature rule?
    /*!
      We never need it! Remark that we put here this method
      as static, while this is not what one want in general.
      In this way, we simply give stronger hints to the
      compiler.
     */
    static bool isAdaptedElement()
    {
        return false;
    }

    //! Getter for the non-adapted quadrature
    const QuadratureRule& standardQR() const
    {
        return M_qr;
    }

    //! Getter for the adapted quadrature
    /*!
      In this case, it should never happen!
     */
    const QuadratureRule& adaptedQR() const
    {
        ERROR_MSG ("No adapted quadrature! Internal error...")
        return M_qr;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No default constructor
    QRAdapterNeverAdapt();

    //@}

    // Quadrature rule
    QuadratureRule M_qr;
};

} // Namespace LifeV

#endif /* QR_ADAPTER_NEVER_ADAPT_H */
