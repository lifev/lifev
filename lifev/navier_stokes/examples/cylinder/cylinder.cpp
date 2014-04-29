/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2005-04-19

  Copyright (C) 2005 EPFL

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
   \file cylinder.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-04-19
 */


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


//#include "life/lifesolver/NavierStokesSolver.hpp"
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/navier_stokes/solver/OseenData.hpp>
#include <lifev/navier_stokes/fem/TimeAdvanceBDFNavierStokes.hpp>
#include <lifev/core/fem/FESpace.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEnsight.hpp>

#include <lifev/navier_stokes/solver/OseenSolver.hpp>

#include "cylinder.hpp"
#include <iostream>



using namespace LifeV;

typedef RegionMesh<LinearTetra>                             mesh_Type;

const int INLET       = 2;
const int WALL        = 1;
const int OUTLET      = 3;
const int RINGIN      = 20;
const int RINGOUT     = 30;


Real zero_scalar ( const Real& /* t */,
                   const Real& /* x */,
                   const Real& /* y */,
                   const Real& /* z */,
                   const ID& /* i */ )
{
    return 0.;
}

Real u2 (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 2:
            if ( t <= 0.003 )
            {
                return 1.3332e4;
            }
            //      return 0.01;
            return 0.0;
            break;
        case 1:
            return 0.0;
            //      return 1.3332e4;
            //    else
            //      return 0.0;
            break;
    }
    return 0;
}

void
postProcessFluxesPressures ( OseenSolver< mesh_Type >& nssolver,
                             BCHandler& bcHandler,
                             const LifeV::Real& t, bool _verbose )
{
    LifeV::Real Q, P;
    UInt flag;

    for ( BCHandler::bcBaseIterator_Type it = bcHandler.begin();
            it != bcHandler.end(); ++it )
    {
        flag = it->flag();

        Q = nssolver.flux (flag);
        P = nssolver.pressure (flag);

        if ( _verbose )
        {
            std::ofstream outfile;
            std::stringstream filenamess;
            std::string filename;

            // file name contains the label
            filenamess << flag;
            // writing down fluxes
            filename = "flux_label" + filenamess.str() + ".m";
            outfile.open (filename.c_str(), std::ios::app);
            outfile << Q << " " << t << "\n";
            outfile.close();
            // writing down pressures
            filename = "pressure_label" + filenamess.str() + ".m";
            outfile.open (filename.c_str(), std::ios::app);
            outfile << P << " " << t << "\n";
            outfile.close();
            // reset ostringstream
            filenamess.str ("");
        }
    }

}


struct Cylinder::Private
{
    Private() :
        //check(false),
        nu (1),
        //rho(1),
        H (1), D (1)
        //H(20), D(1)
        //H(0.41), D(0.1)
    {}
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_Type;

    double Re;

    std::string data_file_name;

    double      nu;  /**< viscosity (in m^2/s) */
    //const double rho; /**< density is constant (in kg/m^3) */
    double      H;   /**< height and width of the domain (in m) */
    double      D;   /**< diameter of the cylinder (in m) */
    bool        centered; /**< true if the cylinder is at the origin */

    std::string initial_sol;

    boost::shared_ptr<Epetra_Comm>   comm;
    /**
     * get the characteristic velocity
     *
     * @return the characteristic velocity
     */
    double Ubar() const
    {
        return nu * Re / D;
    }

    /**
     * get the magnitude of the profile velocity
     *
     *
     * @return the magnitude of the profile velocity
     */
    double Um_3d() const
    {
        return 9 * Ubar() / 4;
    }

    double Um_2d() const
    {
        return 3 * Ubar() / 2;
    }


    /**
     * u3d 3D velocity profile.
     *
     * Define the velocity profile at the inlet for the 3D cylinder
     */
    Real u3d ( const Real& /* t */,
               const Real& /* x */,
               const Real& y,
               const Real& z,
               const ID&   id ) const
    {
        if ( id == 0 )
        {
            if ( centered )
            {
                return Um_3d() * (H + y) * (H - y) * (H + z) * (H - z) / std::pow (H, 4);
            }
            else
            {
                return 16 * Um_3d() * y * z * (H - y) * (H - z) / std::pow (H, 4);
            }
        }
        else
        {
            return 0;
        }
    }

    fct_Type getU_3d()
    {
        fct_Type f;
        f = boost::bind (&Cylinder::Private::u3d, this, _1, _2, _3, _4, _5);
        return f;
    }

    /**
     * u2d flat 2D velocity profile.
     *
     * Define the velocity profile at the inlet for the 2D cylinder
     */
    Real u2d ( const Real& t,
               const Real& /*x*/,
               const Real& /*y*/,
               const Real& /*z*/,
               const ID&   id ) const
    {

        switch (id)
        {
            case 0: // x component
                return 0.0;
                break;
            case 2: // z component
                if ( t <= 0.003 )
                {
                    return 1.3332e4;
                }
                //      return 0.01;
                return 0.0;
                break;
            case 1: // y component
                return 0.0;
                //      return 1.3332e4;
                //    else
                //      return 0.0;
                break;
        }
        return 0;
    }

    fct_Type getU_2d()
    {
        fct_Type f;
        f = boost::bind (&Cylinder::Private::u2d, this, _1, _2, _3, _4, _5);
        return f;
    }

    /**
     * one flat (1,1,1)
     *
     * Define the velocity profile at the inlet for the 2D cylinder
     */
    Real poiseuille ( const Real& /*t*/,
                      const Real& x,
                      const Real& y,
                      const Real& /*z*/,
                      const ID&   id ) const
    {
        double r = std::sqrt (x * x + y * y);

        if (id == 2)
        {
            return Um_2d() * 2 * ( (D / 2.) * (D / 2.) - r * r);
        }

        return 0.;
    }

    fct_Type getU_pois()
    {
        fct_Type f;
        f = boost::bind (&Cylinder::Private::poiseuille, this, _1, _2, _3, _4, _5);
        return f;
    }


    Real oneU ( const Real& /*t*/,
                const Real& /*x*/,
                const Real& /*y*/,
                const Real& /*z*/,
                const ID&   /*id*/ ) const
    {
        //            if (id == 3)
        return 10.;

        return -1.;
    }

    fct_Type getU_one()
    {
        fct_Type f;
        f = boost::bind (&Cylinder::Private::oneU, this, _1, _2, _3, _4, _5);
        return f;
    }


};

Cylinder::Cylinder ( int argc,
                     char** argv )
    :
    d ( new Private )
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );
    d->data_file_name = data_file_name;

    d->Re          = dataFile ( "fluid/problem/Re", 1. );
    d->nu          = dataFile ( "fluid/physics/viscosity", 1. ) /
                     dataFile ( "fluid/physics/density", 1. );
    d->H           = 20.;//dataFile( "fluid/problem/H", 20. );
    d->D           =               dataFile ( "fluid/problem/D", 1. );
    d->centered    = (bool)        dataFile ( "fluid/problem/centered", 0 );
    d->initial_sol = (std::string) dataFile ( "fluid/problem/initial_sol", "stokes");
    std::cout << d->initial_sol << std::endl;


#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    //    MPI_Init(&argc,&argv);

    int ntasks = 0;
    d->comm.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    if (!d->comm->MyPID() )
    {
        std::cout << "My PID = " << d->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
        std::cout << "Re = " << d->Re << std::endl
                  << "nu = " << d->nu << std::endl
                  << "H  = " << d->H  << std::endl
                  << "D  = " << d->D  << std::endl;
    }
    //    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    d->comm.reset ( new Epetra_SerialComm() );
#endif

}

void
Cylinder::run()

{
    typedef FESpace< mesh_Type, MapEpetra >       feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>       feSpacePtr_Type;
    typedef OseenSolver< mesh_Type >::vector_Type vector_Type;
    typedef boost::shared_ptr<vector_Type>        vectorPtr_Type;
    // Reading from data file
    //
    GetPot dataFile ( d->data_file_name );

    //    int save = dataFile("fluid/miscellaneous/save", 1);

    bool verbose = (d->comm->MyPID() == 0);

    // Boundary conditions
    BCHandler bcH;
    BCFunctionBase uZero ( zero_scalar );
    std::vector<ID> zComp (1);
    zComp[0] = 3;

    BCFunctionBase uIn  (  d->getU_2d() );
    BCFunctionBase uOne (  d->getU_one() );
    BCFunctionBase uPois (  d->getU_pois() );


    //BCFunctionBase unormal(  d->get_normal() );

    //cylinder

    bcH.addBC ( "Inlet",    INLET,    Essential,     Full,     uPois  , 3 );
    bcH.addBC ( "Ringin",   RINGIN,   Essential,     Full,     uZero  , 3 );
    bcH.addBC ( "Ringout",  RINGOUT,  Essential,     Full,     uZero  , 3 );
    bcH.addBC ( "Outlet",   OUTLET,   Natural,     Full,     uZero, 3 );
    bcH.addBC ( "Wall",     WALL,     Essential,   Full,     uZero, 3 );

    int numLM = 0;

    boost::shared_ptr<OseenData> oseenData (new OseenData() );
    oseenData->setup ( dataFile );

    MeshData meshData;
    meshData.setup (dataFile, "fluid/space_discretization");

    boost::shared_ptr<mesh_Type> fullMeshPtr ( new mesh_Type ( d->comm ) );
    readMesh (*fullMeshPtr, meshData);

    boost::shared_ptr<mesh_Type> meshPtr;
    {
        MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, d->comm);
        meshPtr = meshPart.meshPartition();
    }
    if (verbose)
    {
        std::cout << std::endl;
    }
    if (verbose)
    {
        std::cout << "Time discretization order " << oseenData->dataTimeAdvance()->orderBDF() << std::endl;
    }

    //oseenData.meshData()->setMesh(meshPtr);

    std::string uOrder =  dataFile ( "fluid/space_discretization/vel_order", "P1");
    if (verbose)
    {
        std::cout << "Building the velocity FE space ... " << std::flush;
    }

    feSpacePtr_Type uFESpacePtr ( new feSpace_Type (meshPtr, uOrder, 3, d->comm) );

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }


    std::string pOrder =  dataFile ( "fluid/space_discretization/press_order", "P1");

    if (verbose)
    {
        std::cout << "Building the pressure FE space ... " << std::flush;
    }

    feSpacePtr_Type pFESpacePtr ( new feSpace_Type (meshPtr, pOrder, 1, d->comm) );

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    UInt totalVelDof   = uFESpacePtr->map().map (Unique)->NumGlobalElements();
    UInt totalPressDof = pFESpacePtr->map().map (Unique)->NumGlobalElements();



    if (verbose)
    {
        std::cout << "Total Velocity DOF = " << totalVelDof << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total Pressure DOF = " << totalPressDof << std::endl;
    }

    if (verbose)
    {
        std::cout << "Calling the fluid constructor ... ";
    }

    bcH.setOffset ( "Inlet", totalVelDof + totalPressDof - 1 );

    OseenSolver< mesh_Type > fluid (oseenData,
                                    *uFESpacePtr,
                                    *pFESpacePtr,
                                    d->comm, numLM);
    MapEpetra fullMap (fluid.getMap() );

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    fluid.setUp (dataFile);
    fluid.buildSystem();

    MPI_Barrier (MPI_COMM_WORLD);

    // Initialization

    Real dt     = oseenData->dataTime()->timeStep();
    Real t0     = oseenData->dataTime()->initialTime();
    Real tFinal = oseenData->dataTime()->endTime();


    // bdf object to store the previous solutions

    TimeAdvanceBDFNavierStokes<vector_Type> bdf;
    bdf.setup (oseenData->dataTimeAdvance()->orderBDF() );

    vector_Type beta ( fullMap );
    vector_Type rhs ( fullMap );

#ifdef HAVE_HDF5
    ExporterHDF5<mesh_Type > ensight ( dataFile, meshPtr, "cylinder", d->comm->MyPID() );
#else
    ExporterEnsight<mesh_Type > ensight ( dataFile, meshPtr, "cylinder", d->comm->MyPID() );
#endif

    vectorPtr_Type velAndPressure ( new vector_Type (*fluid.solution(), ensight.mapType() ) );

    ensight.addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", uFESpacePtr,
                          velAndPressure, UInt (0) );

    ensight.addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpacePtr,
                          velAndPressure, UInt (3 * uFESpacePtr->dof().numTotalDof() ) );

    // initialization with stokes solution

    if (d->initial_sol == "stokes")
    {
        if (verbose)
        {
            std::cout << std::endl;
        }
        if (verbose)
        {
            std::cout << "Computing the stokes solution ... " << std::endl << std::endl;
        }

        oseenData->dataTime()->setTime (t0);

        MPI_Barrier (MPI_COMM_WORLD);

        beta *= 0.;
        rhs  *= 0.;

        fluid.updateSystem (0, beta, rhs );
        fluid.iterate ( bcH );

        //    fluid.postProcess();

        *velAndPressure = *fluid.solution();
        ensight.postProcess ( 0 );
        fluid.resetPreconditioner();
    }

    bdf.bdfVelocity().setInitialCondition ( *fluid.solution() );

    // Temporal loop

    LifeChrono chrono;
    int iter = 1;

    for ( Real time = t0 + dt ; time <= tFinal + dt / 2.; time += dt, iter++)
    {

        oseenData->dataTime()->setTime (time);

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "We are now at time " << oseenData->dataTime()->time() << " s. " << std::endl;
            std::cout << std::endl;
        }

        chrono.start();

        double alpha = bdf.bdfVelocity().coefficientFirstDerivative ( 0 ) / oseenData->dataTime()->timeStep();

        bdf.bdfVelocity().extrapolation (beta); // Extrapolation for the convective term

        bdf.bdfVelocity().updateRHSContribution ( oseenData->dataTime()->timeStep() );

        fluid.setVelocityRhs(bdf.bdfVelocity().rhsContributionFirstDerivative());

        if (oseenData->conservativeFormulation() )
            rhs  = fluid.matrixMass() * bdf.bdfVelocity().rhsContributionFirstDerivative();

        fluid.updateSystem ( alpha, beta, rhs );

        if (!oseenData->conservativeFormulation() )
            rhs  = fluid.matrixMass() * bdf.bdfVelocity().rhsContributionFirstDerivative();

        fluid.iterate ( bcH );

        bdf.bdfVelocity().shiftRight ( *fluid.solution() );

        //         if (((iter % save == 0) || (iter == 1 )))
        //         {
        *velAndPressure = *fluid.solution();
        ensight.postProcess ( time );
        //         }
        //         postProcessFluxesPressures(fluid, bcH, time, verbose);


        MPI_Barrier (MPI_COMM_WORLD);

        chrono.stop();
        if (verbose)
        {
            std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
        }
    }

}


//////////////////////


