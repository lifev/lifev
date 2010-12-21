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
    @brief File containing a class for an easy handling of different order time
            discretizations/extrapolations BDF based specific for the Navier-Stokes problem

    @date 01-04-2003
    @author Alessandro Veneziani <ale@mathcs.emory.edu>

    @contributor Laura Cattaneo
    @mantainer Laura Cattaneo

 */

#ifndef _BDF_NS_TEMPLATE_H
#define _BDF_NS_TEMPLATE_H

#include <string>
#include <iostream>
#include <algorithm>

#include <life/lifefilters/GetPot.hpp>
#include <life/lifefem/bdf_template.hpp>


namespace LifeV
{

//! \class BdfTNS
/*!
    @author Alessandro Veneziani
    @see Van Kan, Prohl, Guermond

    This class implements an easy handling of different order time discretizations/extrapolations
    Bdf based specific for the Navier-Stokes problem.

    The idea is to couple a Bdf of order q with a pressure incremental approach of order q-1.

    If q=1, we still have an incremental pressure treatment.

    At the moment, the couple BDF of order q + extrapolation of order q seems unstable

 */

template<typename VectorType = EpetraVector >
class BdfTNS
{
public:


    //! @name Public Types
    //@{

    typedef VectorType                          vector_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructor
    BdfTNS();

    //! Destructor
    ~BdfTNS() { }

    //@}

    //! @name Set Methods
    //@{

    //! Setup the bdf velocity and the bdf pressure
    /*!
        @param order Order q of the bdf
     */
    void setup( const UInt order );

    //@}

    //! @name Methods
    //@{

    //! Show the contents of the object
    void showMe() const;

    //! The method returns the Bdf velocity
    /*!
        @return Reference to a new Bdf template which holds the velocity field
     */
    BdfT<VectorType>& bdfVelocity() { return M_bdfVelocity; }

    //! The method returns the Bdf pressure
    /*!
        @return Reference to a new Bdf template which holds the pressure
     */
    BdfT<VectorType>& bdfPressure() { return M_bdfPressure; }

    //@}

private:

    //! Bdf velocity
    BdfT<VectorType> M_bdfVelocity;

    //! Bdf pressure
    BdfT<VectorType> M_bdfPressure;
};


// ===================================================
// Constructors
// ===================================================
template<typename VectorType>
BdfTNS<VectorType>::BdfTNS()
        :
        M_bdfVelocity( ),
        M_bdfPressure( )
{}


template<typename VectorType>
void
BdfTNS<VectorType>::setup( const UInt order )
{
    M_bdfVelocity.setup( order, 1 );
    M_bdfPressure.setup( std::max( UInt( 1 ), order - 1 ), 1 );
}

// ===================================================
// Methods
// ===================================================
template<typename VectorType>
void
BdfTNS<VectorType>::showMe() const
{
    std::cout << " *** Bdf velocity: " << std::endl;
    M_bdfVelocity.showMe();
    std::cout << " *** Bdf pressure: " << std::endl;
    M_bdfPressure.showMe();
}


}// Namespace LifeV

#endif /* _BDF_NS_TEMPLATE_H */
