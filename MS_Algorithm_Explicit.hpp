//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief MultiScale Explicit Algorithm
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 26-10-2009
 */

#ifndef MS_Algorithm_Explicit_H
#define MS_Algorithm_Explicit_H 1

#include <lifemc/lifesolver/MS_Algorithm.hpp>

namespace LifeV {

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

    typedef MS_Algorithm                          super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Algorithm_Explicit();

    //! Copy constructor
    /*!
     * @param algorithm MS_Algorithm_Explicit
     */
    MS_Algorithm_Explicit( const MS_Algorithm_Explicit& algorithm );

    //! Destructor
    ~MS_Algorithm_Explicit() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param algorithm MS_Algorithm
     * @return reference to a copy of the class
     */
    MS_Algorithm_Explicit& operator=( const MS_Algorithm_Explicit& algorithm );

    //@}


    //! @name MultiScale Algorithm Virtual Methods
    //@{

    //! Setup the data of the algorithm using a data file
    /*!
     * @param FileName Name of the data file.
     */
    void SetupData( const std::string& FileName );

    //! Perform sub-iteration on the coupling variables
    void SubIterate();

    //! Update coupling variables for the next time step.
    void UpdateCouplingVariables();

    //! Display some information about the algorithm
    void ShowMe();

    //@}

};

//! Factory create function
inline MS_Algorithm* MS_createExplicit()
{
    return new MS_Algorithm_Explicit();
}

} // Namespace LifeV

#endif /* MS_Algorithm_Explicit_H */
