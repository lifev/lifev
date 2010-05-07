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
 *  @brief MultiScale Model 1D
 *
 *  @version 1.0
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *  @date 02-26-2010
 *
 *  @version 1.1
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 14-04-2010
 */

#ifndef MS_Model_1D_H
#define MS_Model_1D_H 1

// Mathcard includes
#include <lifemc/lifealg/AztecOOPreconditioner.hpp>
#include <lifemc/lifefem/BCInterface.hpp>
#include <lifemc/lifesolver/MS_PhysicalModel.hpp>

// LifeV includes
#include <life/lifefem/FESpace.hpp>
#include <life/lifefilters/exporter.hpp>

#ifdef HAVE_HDF5
    #include <life/lifefilters/hdf5exporter.hpp>
#endif

#include <lifemc/lifefem/OneDimensionalModel_BCHandler.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Physics_Linear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Physics_NonLinear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Flux_Linear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Flux_NonLinear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Source_Linear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Source_NonLinear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Solver.hpp>

namespace LifeV {

//! MS_Model_1D - MultiScale model for 1D Fluid simulations
/*!
 *  @author Gilles Fourestey, Cristiano Malossi
 *
 *  The MS_Model_1D class is an implementation of the MS_PhysicalModel
 *  for 1D Fluid problem.
 */
class MS_Model_1D: public virtual MS_PhysicalModel
{
public:

    typedef MS_PhysicalModel                                       super;

    typedef OneDimensionalModel_Physics                            Physics_Type;
    typedef boost::shared_ptr< Physics_Type >                      Physics_PtrType;

    typedef OneDimensionalModel_Flux                               Flux_Type;
    typedef boost::shared_ptr< Flux_Type >                         Flux_PtrType;

    typedef OneDimensionalModel_Source                             Source_Type;
    typedef boost::shared_ptr< Source_Type >                       Source_PtrType;

    typedef OneDimensionalModel_Solver                             Solver_Type;
    typedef boost::shared_ptr< Solver_Type >                       Solver_PtrType;

    typedef Solver_Type::Data_Type                                 Data_Type;
    typedef Solver_Type::Mesh_Type                                 Mesh_Type;
    typedef Solver_Type::Vector_Type                               Vector_Type;
    typedef Solver_Type::Vector_PtrType                            Vector_PtrType;
    typedef Solver_Type::Solution_Type                             Solution_Type;
    typedef Solver_Type::Solution_PtrType                          Solution_PtrType;
    typedef Solver_Type::LinearSolver_Type                         LinearSolver_Type;

    typedef Solver_Type::FESpace_Type                              FESpace_Type;
    typedef Solver_Type::FESpace_PtrType                           FESpace_PtrType;

    typedef OneDimensionalModel_BCHandler                          BC_Type;

#ifdef HAVE_HDF5
    typedef Hdf5exporter< Mesh_Type >                              IOFile_Type;
#endif

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Model_1D();

    //! Copy constructor
    /*!
     * @param 1D MS_Model_1D
     */
    MS_Model_1D( const MS_Model_1D& OneD );

    //! Destructor
    ~MS_Model_1D() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param OneD MS_Model_1D
     * @return reference to a copy of the class
     */
    MS_Model_1D& operator = ( const MS_Model_1D& OneD );

    //@}


    //! @name MultiScale PhysicalModel Virtual Methods
    //@{

    //! Setup the data of the model.
    /*!
     * In particular it does the following operations:
     * <ol>
     *     <li> read data from files;
     *     <li> set global parameter for the MS simulation (viscosity, time, ...);
     *     <li> perform preliminary operations which don't depend on the couplings.
     * </ol>
     *
     * @param FileName Name of data file.
     */
    void SetupData( const std::string& FileName );

    //! Setup the model.
    void SetupModel();

    //! Build the initial system (matrix and vectors).
    void BuildSystem();

    //! Update the system for (matrix and vectors).
    void UpdateSystem();

    //! Solve the problem.
    void SolveSystem();

    //! Save the solution
    void SaveSolution();

    //! Display some information about the model.
    void ShowMe();

    //@}


    //! @name Methods
    //@{

    //! Setup the data of the linear model
    /*!
     * @param FileName Name of data file.
     */
    void SetupLinearData( const std::string& FileName );

    //! Setup the linear model
    void SetupLinearModel();

    //! Update the linear system matrix and vectors
    void UpdateLinearModel();

    //! Solve the linear problem
    void SolveLinearModel( bool& SolveLinearSystem );

    //@}


    //! @name Get Methods
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    BC_Type&      GetBC() const;

    //! Get the data container of the 1D model.
    /*!
     * @return 1D Model data container.
     */
    Data_Type&    GetData() const;

    //! Get the Physics of the 1D model.
    /*!
     * @return 1D Model physics.
     */
    Physics_PtrType GetPhysics() const;

    //! Get the Flux of the 1D model.
    /*!
     * @return 1D Model Flux.
     */
    Flux_PtrType GetFlux() const;

    //! Get the Source of the 1D model.
    /*!
     * @return 1D Model Source.
     */
    Source_PtrType GetSource() const;

    //! Get the FESpace of the 1D model.
    /*!
     * @return 1D model FESpace
     */
    FESpace_PtrType GetFESpace() const;

    //! Get the Solver of the 1D model.
    /*!
     * @return 1D model solver.
     */
    Solver_PtrType  GetSolver() const;

    //! Get the solution of the 1D model.
    /*!
     * @return 1D model solution.
     */
    const Solution_PtrType  GetSolution() const;

    //@}


    //! @name Set Methods
    //@{

    void SetBC( boost::shared_ptr<BC_Type>& BC );

    //@}

private:

    //! @name Private Methods
    //@{

    //! Setup the FE space for pressure and velocity
    void SetupFESpace();

    //@}

#ifdef HAVE_HDF5
    boost::shared_ptr< IOFile_Type >       M_Exporter;
    boost::shared_ptr< Mesh_Type >         M_ExporterMesh;
#endif

    // 1D problem
    boost::shared_ptr< Data_Type >         M_Data;
    Physics_PtrType                        M_Physics;
    Flux_PtrType                           M_Flux;
    Source_PtrType                         M_Source;
    Solver_PtrType                         M_Solver;
    boost::shared_ptr< BC_Type >           M_BC;
    std::vector < Vector_PtrType >         M_Solution;

    // FE spaces
    boost::shared_ptr< FESpace_Type >      M_FESpace;

    // Linear solver
    boost::shared_ptr< LinearSolver_Type > M_LinearSolver;
};

//! Factory create function
inline MS_PhysicalModel* createFluid1D()
{
    return new MS_Model_1D();
}

} // Namespace LifeV

#endif /* MS_Model_1D_H */
