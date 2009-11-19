//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

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
 *  @brief MultiScale Newton Algorithm
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 26-10-2009
 */

#ifndef MS_Algorithm_Newton_H
#define MS_Algorithm_Newton_H 1

#include <life/lifealg/SolverTrilinos.hpp>

#include <life/lifealg/EpetraPreconditioner.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>

#include <lifemc/lifesolver/MS_Algorithm.hpp>

namespace LifeV {

//! MS_Algorithm_Newton - The MultiScale Algorithm implementation of Newton
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_Algorithm_Newton is an implementation of MS_Algorithm
 *  which implements the Newton method.
 */
class MS_Algorithm_Newton : public virtual MS_Algorithm
{
public:

    typedef MS_Algorithm                          super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Algorithm_Newton();

    //! Copy constructor
    /*!
     * @param algorithm MS_Algorithm_Newton
     */
    MS_Algorithm_Newton( const MS_Algorithm_Newton& algorithm );

    //! Destructor
    ~MS_Algorithm_Newton() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param algorithm MS_Algorithm
     * @return reference to a copy of the class
     */
    MS_Algorithm_Newton& operator=( const MS_Algorithm_Newton& algorithm );

    //@}


    //! @name MultiScale Algorithm Virtual Methods
    //@{

    //! Setup the data of the algorithm using a data file
    /*!
     * @param DataFile GetPot data file containing the settings
     */
    void SetupData( const GetPot& DataFile );

    //! Perform sub-iteration on the coupling variables
    void SubIterate();

    //! Display some information about the algorithm
    void ShowMe();

    //@}

protected:

    //! @name Protected Methods
    //@{

    //@}

    SolverTrilinos                           M_solver;
    Matrix_ptrType                           M_Jacobian;
};

//! Factory create function
inline MS_Algorithm* createNewton()
{
    return new MS_Algorithm_Newton();
}

} // Namespace LifeV

#endif /* MS_Algorithm_Newton_H */
