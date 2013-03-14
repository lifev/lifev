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
 *  @brief File containing the Multiscale Broyden Algorithm
 *
 *  @date 26-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleAlgorithmBroyden_H
#define MultiscaleAlgorithmBroyden_H 1

// If the following macro is defined, when a simulation restarts the corresponding saved
// Broyden matrix is loaded from data file.
#define BROYDEN_IMPORTJACOBIAN

#include <lifev/core/algorithm/LinearSolver.hpp>

#include <lifev/multiscale/algorithms/MultiscaleAlgorithm.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleAlgorithmBroyden - The Multiscale Algorithm implementation of Broyden
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
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


    //! @name Multiscale Algorithm Virtual Methods
    //@{

    //! Setup the data of the algorithm using a data file
    /*!
     * @param FileName Name of the data file.
     */
    void setupData ( const std::string& fileName );

    //! Setup coupling variables and other quantities of the algorithm
    void setupAlgorithm();

    //! Perform sub-iteration on the coupling variables
    void subIterate();

    //! Display some information about the algorithm
    void showMe();

    //@}


    //! @name Set Methods
    //@{

    //! Set the the main parameters of the algorithm (tolerance, maximum number of subiterations, etc.)
    /*!
     * @param parameterList teuchos list of parameters
     */
    void setAlgorithmParameters ( const multiscaleParameterList_Type& parameterList );

    //@}

private:

    //! @name Private Types
    //@{

    typedef std::list< multiscaleVector_Type >                    container_Type;
    typedef container_Type::const_iterator                        containerIterator_Type;

    //@}


    //! @name Unimplemented Methods
    //@{

    MultiscaleAlgorithmBroyden ( const MultiscaleAlgorithmBroyden& algorithm );

    MultiscaleAlgorithmBroyden& operator= ( const MultiscaleAlgorithmBroyden& algorithm );

    //@}


    //! @name Private Methods
    //@{

    void assembleJacobianMatrix();

    void broydenJacobianUpdate ( const multiscaleVector_Type& delta );

    void orthogonalizationUpdate ( const multiscaleVector_Type& delta );

#ifdef HAVE_HDF5
    //! Export Jacobian matrix to HDF5 file
    void exportJacobianToHDF5();

    //! Import Jacobian matrix from an HDF5 file
    void importJacobianFromHDF5();
#endif

    //@}

    LinearSolver                             M_solver;
    multiscaleMatrixPtr_Type                 M_jacobian;

    bool                                     M_initializeAsIdentityMatrix;
    bool                                     M_iterationsLimitReached;
    UInt                                     M_iterationsLimitForReset;
    bool                                     M_orthogonalization;
    UInt                                     M_orthogonalizationSize;
    container_Type                           M_orthogonalizationContainer;
    bool                                     M_truncate;
};

//! Factory create function
inline multiscaleAlgorithm_Type* createMultiscaleAlgorithmBroyden()
{
    return new MultiscaleAlgorithmBroyden();
}

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleAlgorithmBroyden_H */
