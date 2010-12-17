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
 *  @brief File containing the MultiScale Explicit Algorithm
 *
 *  @date 26-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleAlgorithmExplicit_H
#define MultiscaleAlgorithmExplicit_H 1

#include <lifemc/lifesolver/MultiscaleAlgorithm.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleAlgorithmExplicit - The MultiScale Algorithm implementation of Explicit
/*!
 *  @author Cristiano Malossi
 *
 *  The MultiscaleAlgorithmExplicit is an implementation of multiscaleAlgorithm_Type
 *  which implements the Explicit method.
 */
class MultiscaleAlgorithmExplicit : public virtual multiscaleAlgorithm_Type
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleAlgorithmExplicit();

    //! Destructor
    virtual ~MultiscaleAlgorithmExplicit() {}

    //@}


    //! @name MultiScale Algorithm Virtual Methods
    //@{

    //! Perform sub-iteration on the coupling variables
    void subIterate() { toleranceSatisfied(); }

    //! Update coupling variables for the next time step.
    void updateCouplingVariables() { M_multiscale->initializeCouplingVariables(); }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleAlgorithmExplicit( const MultiscaleAlgorithmExplicit& algorithm );

    MultiscaleAlgorithmExplicit& operator=( const MultiscaleAlgorithmExplicit& algorithm );

    //@}

};

//! Factory create function
inline multiscaleAlgorithm_Type* createMultiscaleAlgorithmExplicit()
{
    return new MultiscaleAlgorithmExplicit();
}

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleAlgorithmExplicit_H */
