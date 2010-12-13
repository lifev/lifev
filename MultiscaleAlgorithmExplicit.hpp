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
namespace multiscale
{

//! MultiscaleAlgorithmExplicit - The MultiScale Algorithm implementation of Explicit
/*!
 *  @author Cristiano Malossi
 *
 *  The MultiscaleAlgorithmExplicit is an implementation of MS_Algorithm_Type
 *  which implements the Explicit method.
 */
class MultiscaleAlgorithmExplicit : public virtual MS_Algorithm_Type
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MultiscaleAlgorithmExplicit();

    //! Destructor
    virtual ~MultiscaleAlgorithmExplicit() {}

    //@}


    //! @name MultiScale Algorithm Virtual Methods
    //@{

    //! Setup the data of the algorithm using a data file
    /*!
     * @param FileName Name of the data file.
     */
    void setupData( const std::string& fileName );

    //! Perform sub-iteration on the coupling variables
    void subIterate();

    //! Update coupling variables for the next time step.
    void updateCouplingVariables();

    //! Display some information about the algorithm
    void showMe();

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleAlgorithmExplicit( const MultiscaleAlgorithmExplicit& algorithm );

    MultiscaleAlgorithmExplicit& operator=( const MultiscaleAlgorithmExplicit& algorithm );

    //@}

};

//! Factory create function
inline MS_Algorithm_Type* createMultiscaleAlgorithmExplicit()
{
    return new MultiscaleAlgorithmExplicit();
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleAlgorithmExplicit_H */
