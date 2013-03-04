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
 *  @brief File containing the Multiscale Model 0D
 *
 *  @version 1.0
 *  @date 01-10-2011
 *  @author Mahmoud Jafargholi <mahmoud.jafargholi@epfl.ch>
 *
 *  @contributors Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @mantainer    Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleModel0D_H
#define MultiscaleModel0D_H 1

#include <lifev/bc_interface/fem/BCInterface0D.hpp>
#include <lifev/zero_dimensional/solver/ZeroDimensionalData.hpp>
#include <lifev/zero_dimensional/solver/ZeroDimensionalSolver.hpp>

#include <lifev/multiscale/solver/MultiscaleModel.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleModel0D - Multiscale model for 0D simulations
/*!
 *  @author Mahmoud Jafargholi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleModel0D class is an implementation of the multiscaleModel_Type
 *  for 0D problems.
 */

class MultiscaleModel0D: public virtual multiscaleModel_Type
{
public:

    //! @name Type definitions
    //@{

    typedef ZeroDimensionalBCHandler                                        bc_Type;
    typedef boost::shared_ptr< bc_Type >                                    bcPtr_Type;
    typedef BCInterface0D< bc_Type, MultiscaleGlobalData >                  bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >                           bcInterfacePtr_Type;

    typedef ZeroDimensionalData                                             data_Type;
    typedef boost::shared_ptr< data_Type >                                  dataPtr_Type;
    typedef ZeroDimensionalSolver                                           solver_Type;
    typedef boost::shared_ptr< solver_Type >                                solverPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleModel0D();

    //! Destructor
    virtual ~MultiscaleModel0D() {}

    //@}


    //! @name MultiscaleModel Methods
    //@{

    //! Setup the data of the model.
    /*!
     * @param fileName Name of data file.
     */
    void setupData( const std::string& fileName );

    //! Setup the model.
    void setupModel();

    //! Build the initial model.
    void buildModel();

    //! Update the model.
    void updateModel();

    //! Solve the model.
    void solveModel();

    //! Update the solution.
    void updateSolution();

    //! Save the solution
    void saveSolution();

    //! Display some information about the model.
    void showMe();

    //! Return a specific scalar quantity to be used for a comparison with a reference value.
    /*!
     * This method is meant to be used for night checks.
     * @return reference quantity.
     */
    Real checkSolution() const;

    //@}


    //! @name MultiscaleInterface Methods
    //@{

    //TODO These methods should be implemented after deriving from MultiscaleInterface

    //@}


    //! @name Get Methods
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    bcInterface_Type& bcInterface() { return *M_bc; }

    //! Get the data container of the model
    /*!
     * @return data container
     */
    data_Type& data() const { return *M_data; }

    //! Get the solver of the model
    /*!
     * @return solver
     */
    solver_Type& solver() const { return *M_solver; }

    //@}

private:

    //! @name Private Methods
    //@{

    //! Setup the global data of the model.
    /*!
     * In particular, it replaces the default local values with the ones in the global container.
     * If a value is already specified in the data file, do not perform the replacement.
     *
     * @param fileName File name of the specific model.
     */
    void setupGlobalData( const std::string& fileName );

    //@}


    dataPtr_Type                           M_data;
    solverPtr_Type                         M_solver;
    bcInterfacePtr_Type                    M_bc;
};

//! Factory create function
inline multiscaleModel_Type* createMultiscaleModelZeroDimensional()
{
    return new MultiscaleModel0D();
}

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleModel0D_H */
