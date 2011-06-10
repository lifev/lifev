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
/**
   \file ethiersteiman.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2011-03-08
 */

#ifndef __Ethiersteinman_H
#define __Ethiersteinman_H 1

#include <life/lifesolver/OseenSolver.hpp>
#include <life/lifemesh/RegionMesh3D.hpp>
#include <life/lifemesh/ElementShapes.hpp>

#include "life/lifefunctions/RossEthierSteinmanDec.hpp"
#include "life/lifefunctions/RossEthierSteinmanInc.hpp"

enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

/*!
 * \class Ethiersteinman
 * \brief 2D/3D Ethiersteinman Simulation class
 *
 *  @author Christophe Prud'homme
 *  @author Gwenol Grandperrin
 *  @see
 */

class Ethiersteinman
{
public:
    typedef LifeV::RegionMesh3D<LifeV::LinearTetra>       mesh_Type;
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

    Ethiersteinman( int argc,
                    char** argv );

    ~Ethiersteinman()
    {}

    //@}

    /** @name  Methods
     */
    //@{
    //! Computes L2 errors
    void run();

    //@}


private:
    enum TestType{None, Accuracy,SpaceConvergence};
    enum InitializationType{Projection,Interpolation};
    enum MeshSourceType{File,RegularMesh};

    struct RESULT_CHANGED_EXCEPTION
    {
    public:
        RESULT_CHANGED_EXCEPTION()
        {
            std::cout << "Some modifications led to changes in the l2 norm of the solution" << std::endl;
        }
    };

    void computeErrors(const vector_Type& velocityAndPressureSolution,
                       LifeV::Real& uL2Error, LifeV::Real& uRelError, feSpacePtr_Type& uFESpace,
                       LifeV::Real& pL2Error, LifeV::Real& pRelError, feSpacePtr_Type& pFESpace,
                       LifeV::Real time);

    bool checkConvergenceRate(const std::vector<std::string>& uFELabel,
                              const std::vector<std::vector<LifeV::Real> >& uL2Error,
                              const std::vector<LifeV::UInt>& uConvergenceOrder,
                              const std::vector<std::string>& pFELabel,
                              const std::vector<std::vector<LifeV::Real> > pL2Error,
                              const std::vector<LifeV::UInt>& pConvergenceOrder,
                              const std::vector<LifeV::UInt>& meshDiscretizations,
                              LifeV::Real convTolerance);


    struct Private;
    boost::shared_ptr<Private> d;

    std::vector<LifeV::UInt> meshDiscretization;
    std::vector<std::string> uFE;
    std::vector<std::string> pFE;
    std::vector<LifeV::UInt> uConvergenceOrder;
    std::vector<LifeV::UInt> pConvergenceOrder;

    // Test to be performed (accuracy or convergence in space)
    TestType    M_test;
    LifeV::Real M_convTol; // Tolerance of the test (should be <1)
                           // Actually for convTol=1, the test failed
                           // if the improvement of accuracy is less
                           // than predicted by the theory.
                           // convTol lower down the theoretical bounds
    LifeV::Real        M_accuracyTol;

    // Data related to norm export
    bool               M_exportNorms;
    std::ofstream      out_norm;

    // Data related to solution export
    bool               M_exportExactSolutions;

    // Initialization method
    InitializationType M_initMethod;

    // Mesh source
    MeshSourceType     M_meshSource;
};

#endif /* __Ethiersteinman_H */
