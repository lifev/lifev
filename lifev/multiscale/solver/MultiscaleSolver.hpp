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
 *  @brief File containing the Multiscale Solver
 *
 *  @date 28-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleSolver_H
#define MultiscaleSolver_H 1

#include <lifev/multiscale/solver/MultiscaleDefinitions.hpp>
#include <lifev/multiscale/solver/MultiscaleModelMultiscale.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleSolver - The Multiscale solver
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleSolver class provides a series of methods to create and
 *  solve a general Geometrical Multiscale problem.
 */
class MultiscaleSolver
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleSolver();

    //! Destructor
    virtual ~MultiscaleSolver() {}

    //@}


    //! @name Methods
    //@{

    //! Setup the problem
    /*!
     * @param fileName Name of the data file.
     * @param problemName the name of the problem (used to save data in a specific folder).
     * @param coresPerNode number of cores for each node (this is mandatory when running on clusters for a correct distribution of the models among the nodes).
     */
    void setupProblem ( const std::string& fileName, const std::string& problemName, const UInt& coresPerNode );

    //! Run the time-loop to solve the Multiscale problem
    /*!
     * If the provided reference solution is positive, the solver make also a check on the last computed solution.
     * @param referenceSolution the reference coupling variables norm 2.
     * @param tolerance the tolerance to check the reference solution with respect to the computed one.
     * @return 0: EXIT_SUCCESS, 1: EXIT_FAILURE
     */
    bool solveProblem ( const Real& referenceSolution = -1., const Real& tolerance = 1e-8 );

    //! Display some information about the Multiscale problem (should be called after setupProblem)
    void showMe() const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the epetra communicator for the Multiscale problem
    /*!
     * @param comm Epetra communicator
     */
    void setCommunicator ( const multiscaleCommPtr_Type& comm )
    {
        M_comm = comm;
    }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleSolver ( const MultiscaleSolver& solver );

    MultiscaleSolver& operator= ( const MultiscaleSolver& solver );

    //@}


    //! @name Private methods
    //@{

    //! Save CPU time at each time step
    /*!
     * @param totalCPUTime total CPU time of the iteration
     * @param buildUpdateCPUTime CPU time to build/update the problem
     * @param solveCPUTime CPU time to solve the problem
     * @param updateSolutionCPUTime CPU time to update the solution of the problem
     * @param saveCPUTime CPU time to save the solution
     */
    void saveCPUTime ( const Real& totalCPUTime, const Real& buildUpdateCPUTime, const Real& solveCPUTime,
                       const Real& updateSolutionCPUTime, const Real& saveCPUTime ) const;

    //! Import iteration number from the CPU file
    void importIterationNumber();

    //@}


    // The main model (can be a specific model or a Multiscale model)
    multiscaleModelPtr_Type          M_model;

    // Container of global data
    multiscaleDataPtr_Type           M_globalData;

    // Communicator
    multiscaleCommPtr_Type           M_comm;
};

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleSolver_H */
