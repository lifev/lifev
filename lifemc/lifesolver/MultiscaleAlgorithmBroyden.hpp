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
 *  @brief File containing the MultiScale Broyden Algorithm
 *
 *  @date 26-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleAlgorithmBroyden_H
#define MultiscaleAlgorithmBroyden_H 1

#include <life/lifealg/SolverTrilinos.hpp>

#include <life/lifealg/EpetraPreconditioner.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>

#include <lifemc/lifesolver/MultiscaleAlgorithm.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleAlgorithmBroyden - The MultiScale Algorithm implementation of Broyden
/*!
 *  @author Cristiano Malossi
 *
 *  The MultiscaleAlgorithmBroyden is an implementation of multiscaleAlgorithm_Type
 *  which implements the Broyden method.
 */
class MultiscaleAlgorithmBroyden : public virtual multiscaleAlgorithm_Type
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleAlgorithmBroyden();

    //! Destructor
    virtual ~MultiscaleAlgorithmBroyden() {}

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

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleAlgorithmBroyden( const MultiscaleAlgorithmBroyden& algorithm );

    MultiscaleAlgorithmBroyden& operator=( const MultiscaleAlgorithmBroyden& algorithm );

    //@}


    //! @name Private Methods
    //@{

    void assembleJacobianMatrix();

    void broydenJacobianUpdate( const multiscaleVector_Type& delta, const multiscaleVector_Type& minusCouplingResidual );

    //@}

    SolverTrilinos                           M_solver;
    multiscaleMatrixPtr_Type                 M_jacobian;
};

//! Factory create function
inline multiscaleAlgorithm_Type* createMultiscaleAlgorithmBroyden()
{
    return new MultiscaleAlgorithmBroyden();
}

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleAlgorithmBroyden_H */
