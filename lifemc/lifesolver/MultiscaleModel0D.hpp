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
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleModel0D_H
#define MultiscaleModel0D_H 1


// Mathcard includes

#include <lifemc/lifesolver/ZeroDimensionalData.hpp>
#include <lifemc/lifesolver/ZeroDimensionalSolver.hpp>
#include <lifemc/lifesolver/MultiscaleModel.hpp>
#include <lifemc/lifesolver/BCInterface0D.hpp>

// LIFEV
namespace LifeV
{
namespace Multiscale
{

//! MultiscaleModel0D - Multiscale model for 0D  simulations
/*!
 *  @author Mahmoud Jafargholi
 *
 */

class MultiscaleModel0D: public virtual multiscaleModel_Type

{
public:

    //! @name Type definitions
    //@{

    typedef ZeroDimensionalBCHandler                                        bc_Type;

    typedef boost::shared_ptr< bc_Type >                                    bcPtr_Type;

    typedef BCInterface0D< bc_Type, MultiscaleData >                        bcInterface_Type;

    typedef boost::shared_ptr< bcInterface_Type >                           bcInterfacePtr_Type;

    typedef ZeroDimensionalData                                             data_Type;

    typedef ZeroDimensionalSolver                                           solver_Type;
    //@}
        //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleModel0D();

    //! Destructor
    virtual ~MultiscaleModel0D();

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
    void updateSolution() {};

    //! Save the solution
    void saveSolution();

    //! Display some information about the model.
    void showMe();

    //! show some more deteil data
    void showMeVariables();

    //! this method is empty
    Real checkSolution() const;

    //! get model data container
    data_Type& data() const { return *M_data; }

    //! get 0D solver
    solver_Type& solver() const {return *M_solver;}

private:

    boost::shared_ptr< data_Type >         M_data;

    boost::shared_ptr< solver_Type >       M_solver;

    Real                                   M_Tn;

    Real                                   M_TnPlus;

    bcInterfacePtr_Type                    M_bc;

};
} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleModel0D_H */
