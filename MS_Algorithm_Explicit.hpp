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

#ifndef MS_Algorithm_Explicit_H
#define MS_Algorithm_Explicit_H 1

#include <lifemc/lifesolver/MS_Algorithm.hpp>

namespace LifeV
{

//! MS_Algorithm_Explicit - The MultiScale Algorithm implementation of Explicit
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_Algorithm_Explicit is an implementation of MS_Algorithm
 *  which implements the Explicit method.
 */
class MS_Algorithm_Explicit : public virtual MS_Algorithm
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Algorithm_Explicit();

    //! Destructor
    ~MS_Algorithm_Explicit() {}

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

};

//! Factory create function
inline MS_Algorithm* MS_createExplicit()
{
    return new MS_Algorithm_Explicit();
}

} // Namespace LifeV

#endif /* MS_Algorithm_Explicit_H */
