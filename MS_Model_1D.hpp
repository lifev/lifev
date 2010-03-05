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
 *  @brief MultiScale Model 1D
 *
 *  @gilles fourestey <gilles.fourestey@epfl.ch>
 *  @date 02-26--2010
 */


#ifndef MS_Model_1D_H
#define MS_Model_1D_H 1

// Mathcard includes
#include <lifemc/lifealg/AztecOOPreconditioner.hpp>
#include <lifemc/lifefem/BCInterface.hpp>
#include <lifemc/lifesolver/MS_PhysicalModel.hpp>

// LifeV includes
//#include <life/lifefilters/ensight.hpp>
//#include <life/lifemesh/partitionMesh.hpp>
//#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifefem/FESpace.hpp>
//#include <life/lifefem/bdfNS_template.hpp>
#include <lifemc/lifesolver/oneDModelSolver.hpp>
#include <lifemc/lifesolver/oneDBCHandler.hpp>
//
#include <life/lifefilters/exporter.hpp>

#ifdef HAVE_HDF5
    #include <life/lifefilters/hdf5exporter.hpp>
#else
    #include <life/lifefilters/ensight.hpp>
#endif

//#include <life/lifesolver/OseenShapeDerivative.hpp>

namespace LifeV {

//! MS_Model_1D - MultiScale model for 3D Fluid simulations
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_Model_1D class is an implementation of the MS_PhysicalModel
 *  for 3D Fluid problem (in particular Oseen with Shape Derivatives).
 */

class MS_Model_1D: public virtual MS_PhysicalModel
{
public:

    typedef MS_PhysicalModel                  super;

    //typedef RegionMesh1D< LinearTetra >       MeshType;



    typedef NonLinearFluxFun1D                    Flux;
    typedef NonLinearSourceFun1D                  Source;
    typedef OneDNonLinModelParam                  Params;

    typedef OneDModelSolver<Params, Flux, Source> SolverType;

    typedef SolverType::data_type                 DataType;
    typedef SolverType::Mesh                      MeshType;
    typedef SolverType::vector_type               VectorType;

    typedef  OneDBCHandler<Flux>                  BCType;

#ifdef HAVE_HDF5
    typedef Hdf5exporter< MeshType >          OutputType;
#else
    typedef Ensight< MeshType >               OutputType;
#endif


//     typedef BCInterface< FluidType >          BCType;
//     typedef BdfTNS< FluidVectorType >         FluidBDFType;
//     typedef DataNavierStokes< MeshType >      FluidDataType;
//     typedef partitionMesh< MeshType >         FluidMeshType;

    typedef FESpace< MeshType, EpetraMap >                  FESpaceType;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Model_1D();

    //! Copy constructor
    /*!
     * @param 1D MS_Model_1D
     */
    MS_Model_1D( const MS_Model_1D& oneD );

    //! Destructor
    ~MS_Model_1D(){};

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param fluid3D MS_Model_1D
     * @return reference to a copy of the class
     */
    MS_Model_1D& operator = ( const MS_Model_1D& fluid3D );

    //@}


    //! @name MultiScale PhysicalModel Virtual Methods
    //@{

    //! Setup the data of the model.
    void SetupData();
    void SetupData(const GetPot& dataFile, std::string section);

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
    void SetupLinearData();

    //! Setup the linear model
    void SetupLinearModel();

    //! Update the linear system matrix and vectors
    void UpdateLinearModel();

    //! Solve the linear problem
    void SolveLinearModel( bool& SolveLinearSystem );

    //@}


    //! @name Get functions
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    BCType&      GetBC();

    SolverType&  GetSolver(){return *M_solver;}
    FESpaceType& GetFESpace();


    //@}

    //! @name Set functions
    //@{

    void SetBC(boost::shared_ptr<BCType> BC)
    {
        M_BC = BC;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! Setup the FE space for pressure and velocity
    void setUpFESpace();

    //@}

    boost::shared_ptr< OutputType >       M_output;

    // 1d physical solver
    boost::shared_ptr< SolverType >       M_solver;

    // linear solver
    boost::shared_ptr<solver_type>        M_linearSolver;

    //
    boost::shared_ptr< BCType >           M_BC;

    //
    boost::shared_ptr< DataType >         M_data;

    //
    boost::shared_ptr< Params >           M_params;
    boost::shared_ptr< Flux>              M_flux;
    boost::shared_ptr< Source >           M_source;

    boost::shared_ptr< VectorType >       M_solution;

    // FE spaces
    boost::shared_ptr< FESpaceType >      M_FESpace;
};

// //! Factory create function
// inline MS_PhysicalModel* create1D()
// {
//     return new MS_Model_1D();
// }

} // Namespace LifeV

#endif /* MS_Model_1D_H */
