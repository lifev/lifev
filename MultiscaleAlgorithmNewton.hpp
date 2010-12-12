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
 *  @brief File containing the MultiScale Newton Algorithm
 *
 *  @date 26-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MS_Algorithm_Newton_H
#define MS_Algorithm_Newton_H 1

#include <life/lifealg/SolverTrilinos.hpp>

#include <life/lifealg/EpetraPreconditioner.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>

#include <lifemc/lifesolver/MultiscaleAlgorithm.hpp>

namespace LifeV
{

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

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Algorithm_Newton();

    //! Destructor
    virtual ~MS_Algorithm_Newton() {}

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

    //! Display some information about the algorithm
    void showMe();

    //@}

protected:

    //! @name Protected Methods
    //@{

    //@}

    SolverTrilinos                           M_solver;
    MS_Matrix_PtrType                        M_jacobian;

private:

    //! @name Unimplemented Methods
    //@{

    MS_Algorithm_Newton( const MS_Algorithm_Newton& algorithm );

    MS_Algorithm_Newton& operator=( const MS_Algorithm_Newton& algorithm );

    //@}
};

//! Factory create function
inline MS_Algorithm* createMultiscaleAlgorithmNewton()
{
    return new MS_Algorithm_Newton();
}

} // Namespace LifeV

#endif /* MS_Algorithm_Newton_H */
