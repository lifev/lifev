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
 *  @brief File containing the Multiscale Newton Algorithm
 *
 *  @date 26-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleAlgorithmNewton_H
#define MultiscaleAlgorithmNewton_H 1

#include <lifev/core/algorithm/LinearSolver.hpp>

#include <lifev/multiscale/solver/MultiscaleAlgorithm.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleAlgorithmNewton - The Multiscale Algorithm implementation of Newton
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleAlgorithmNewton is an implementation of multiscaleAlgorithm_Type
 *  which implements the Newton method.
 */
class MultiscaleAlgorithmNewton : public virtual multiscaleAlgorithm_Type
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleAlgorithmNewton();

    //! Destructor
    virtual ~MultiscaleAlgorithmNewton() {}

    //@}


    //! @name Multiscale Algorithm Virtual Methods
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

    MultiscaleAlgorithmNewton( const MultiscaleAlgorithmNewton& algorithm );

    MultiscaleAlgorithmNewton& operator=( const MultiscaleAlgorithmNewton& algorithm );

    //@}


    //! @name Private Methods
    //@{

    void assembleJacobianMatrix();

    //@}

    LinearSolver                             M_solver;
    multiscaleMatrixPtr_Type                 M_jacobian;
};

//! Factory create function
inline multiscaleAlgorithm_Type* createMultiscaleAlgorithmNewton()
{
    return new MultiscaleAlgorithmNewton();
}

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleAlgorithmNewton_H */
