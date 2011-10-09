/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
             Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2011-03-09

  Copyright (C) 2010 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/*!
    @file ethiersteiman.hpp
    @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 2011-03-08
 */

#ifndef ETHIERSTEINMAN_H
#define ETHIERSTEINMAN_H 1

#include <life/lifesolver/OseenSolver.hpp>
#include <life/lifemesh/RegionMesh3D.hpp>
#include <life/lifemesh/ElementShapes.hpp>

#include "life/lifefunctions/RossEthierSteinmanDec.hpp"
#include "life/lifefunctions/RossEthierSteinmanInc.hpp"

/*! @enum TimeScheme
    Order of the BDF
 */
enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

/*!
    @class Ethiersteinman
    @brief Ethiersteinman Simulation class

    @author Christophe Prud'homme
    @author Gwenol Grandperrin
    @see C.R. Ethier and D.A. Steinman. Exact fully 3D Navier-Stokes solutions for benchmarking. Int. J. Numer. Methods Fluids, 19(5):369-375, 1994

    This class tests many settings, we propose here a list of the available options.
    <ul>
        <li> RossEthierSteinman/test (none/accuracy/space_convergence)
        <li> RossEthierSteinman/initialization (interpolation/projection)
        <li> RossEthierSteinman/export_norms
        <li> RossEthierSteinman/export_exact_solutions
        <li> RossEthierSteinman/mesh_source (regular_mesh/file)
        <li> exporter/type (ensight/hdf5)
        <li> fluid/problem/Re
        <li> fluid/physics/viscosity
        <li> fluid/physics/density
    </ul>

    In the case of a space convergence test the following options are available
    <ul>
        <li> RossEthierSteinman/space_convergence_tolerance
        <li> fluid/space_discretization/mesh_number
        <li> fluid/space_discretization/mesh_size
        <li> fluid/space_discretization/FE_number
        <li> fluid/space_discretization/vel_order
        <li> fluid/space_discretization/vel_conv_order_order
        <li> fluid/space_discretization/press_order
        <li> fluid/space_discretization/press_conv_order
    </ul>

    In the case of an accuracy test the following options are available
    <ul>
        <li> RossEthierSteinman/accuracy_tolerance
    </ul>

    If the mesh source is a file, the file must be specified in
    <ul>
        <li> fluid/space_discretization
    </ul>
 */

class Ethiersteinman
{
public:
    typedef LifeV::RegionMesh<LifeV::LinearTetra>       mesh_Type;
    typedef LifeV::FESpace< mesh_Type, LifeV::MapEpetra > feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>               feSpacePtr_Type;
    typedef LifeV::OseenSolver< mesh_Type >               fluid_Type;
    typedef fluid_Type::vector_Type                       vector_Type;
    typedef boost::shared_ptr<vector_Type>                vectorPtr_Type;
    typedef fluid_Type::matrix_Type                       matrix_Type;
    typedef LifeV::RossEthierSteinmanUnsteadyDec          problem_Type;

    /** @name Constructors, destructor
     */
    //@{

    //! Constructor
    /*!
        @param argc number of parameter passed through the command line
        @param argv char passed through the command line
     */
    Ethiersteinman( int argc,
                    char** argv );

    //! Destructor
    ~Ethiersteinman()
    {}

    //@}

    /** @name  Methods
     */
    //@{

    //! Launches the simulation
    void run();

    //@}


private:
    /*! @enum TestType
        Order of the BDF
     */
    enum TestType{None, Accuracy,SpaceConvergence};

    /*! @enum InitializationType
        Type of initialization. "Interpolation" just interpolates the value of the exact solution to the DoFs.
        "Projection" solves an Oseen problem where alpha=0, the convective term is linearized by using the exact solution for beta,
        and the time derivative is passed to the right hand side and computed from the exact solution.
     */
    enum InitializationType{Projection,Interpolation};

    /*! @enum MeshSourceType
        Type of mesh source. It can be a "File" or a "RegularMesh" which is generated during the simulation
     */
    enum MeshSourceType{File,RegularMesh};

    struct RESULT_CHANGED_EXCEPTION
    {
    public:
        RESULT_CHANGED_EXCEPTION()
        {
            std::cout << "Some modifications led to changes in the l2 norm of the solution" << std::endl;
        }
    };

    //! Computes the L2 error and the relative L2 error with the exact solution
    /*!
        @param velocityAndPressureSolution number of parameter passed through the command line
        @param uL2Error Variable to store the L2 error for the veloctity
        @param uRelError Variable to store the Relative L2 error for the velocity
        @param uFESpace Variable to store the FE space for the velocity
        @param pL2Error Variable to store the L2 error the pressure
        @param pRelError Variable to store the Relative L2 error for the pressure
        @param pFESpace Variable to store the FE space for the pressure
        @param time Actual timestep of the simulation
     */
    void computeErrors(const vector_Type& velocityAndPressureSolution,
                       LifeV::Real& uL2Error, LifeV::Real& uRelError, feSpacePtr_Type& uFESpace,
                       LifeV::Real& pL2Error, LifeV::Real& pRelError, feSpacePtr_Type& pFESpace,
                       LifeV::Real time);

    //! Method to check the convergence rate of the solution
    /*!
        For each type of finite elements the method uses the computed errors obtained with the
        different meshes to check if the order of convergence follows the theory predictions
        @param uFELabels Vector containing the FE names for the velocity (e.g. P2, P1Bubble, P1)
        @param uL2Error Vector containing the computed errors for the velocity
        @param uConvergenceOrder Vector containing the convergence order corresponding to uFELabel
        @param pFELabels Vector containing the FE names for the pressure (e.g. P2, P1Bubble, P1)
        @param pL2Error Vector containing the computed errors for the pressure
        @param pConvergenceOrder Vector containing the convergence order corresponding to pFELabel
        @param meshDiscretization Vector containing the subdivisions values used to generate the meshes
        @param convTolerance Tolerance for the test. The test is passed if (observed convergence)>convTolerance*(theory error prediction)
     */
    bool checkConvergenceRate(const std::vector<std::string>& uFELabels,
                              const std::vector<std::vector<LifeV::Real> >& uL2Error,
                              const std::vector<LifeV::UInt>& uConvergenceOrder,
                              const std::vector<std::string>& pFELabels,
                              const std::vector<std::vector<LifeV::Real> > pL2Error,
                              const std::vector<LifeV::UInt>& pConvergenceOrder,
                              const std::vector<LifeV::UInt>& meshDiscretizations,
                              LifeV::Real convTolerance);


    struct Private;
    boost::shared_ptr<Private> M_data;

    std::vector<LifeV::UInt>   M_meshDiscretization;
    std::vector<std::string>   M_uFELabels;
    std::vector<std::string>   M_pFELabels;
    std::vector<LifeV::UInt>   M_uConvergenceOrder;
    std::vector<LifeV::UInt>   M_pConvergenceOrder;

    // Test to be performed (accuracy or convergence in space)
    TestType                   M_test;
    LifeV::Real                M_convTol; // Tolerance of the test (should be <1)
                                          // Actually for convTol=1, the test failed
                                          // if the improvement of accuracy is less
                                          // than predicted by the theory.
                                          // convTol lower down the theoretical bounds
    LifeV::Real                M_accuracyTol;

    // Data related to norm export
    bool                       M_exportNorms;
    std::ofstream              M_outNorm;

    // Data related to solution export
    bool                       M_exportExactSolutions;

    // Initialization method
    InitializationType         M_initMethod;

    // Mesh source
    MeshSourceType             M_meshSource;
};

#endif /* ETHIERSTEINMAN_H */
