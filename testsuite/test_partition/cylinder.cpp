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
    @file
    @brief

    @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
            Radu Popescu <radu.popescu@epfl.ch>
    @date 19-04-2005

 */




// #include <life/lifecore/life.hpp>

// #include <life/lifecore/life.hpp>
// #include <life/lifecore/GetPot.hpp>
// #include <life/lifecore/debug.hpp>

// #include <life/lifefilters/importer.hpp>

//#include "NavierStokesSolverBlockIP.hpp"

//#include "Epetra_SerialComm.h"
#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#include <mpi.h>
#else
#include <Epetra_SerialComm.h>
#endif
//#include "life/lifesolver/NavierStokesSolver.hpp"
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>
#include <life/lifefilters/HDF5Filter3DMesh.hpp>
#include <life/lifefilters/ensight.hpp>

#include <life/lifesolver/Oseen.hpp>

#include "cylinder.hpp"
#include <iostream>



using namespace LifeV;



const int INLET       = 2;
const int WALL        = 1;
const int OUTLET      = 3;

Real zero_scalar( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.;
}

Real u2(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
    case 1:
        return 0.0;
        break;
    case 3:
        if ( t <= 0.003 )
            return 1.3332e4;
//      return 0.01;
        return 0.0;
        break;
    case 2:
        return 0.0;
//      return 1.3332e4;
//    else
//      return 0.0;
        break;
    }
    return 0;
}

void
postProcessFluxesPressures( Oseen< RegionMesh3D<LinearTetra> >& nssolver,
                            BCHandler& bcHandler,
                            const LifeV::Real& t, bool _verbose )
{
    LifeV::Real Q, P;
    UInt flag;

    for ( BCHandler::bcBaseIterator_Type it = bcHandler.begin();
            it != bcHandler.end(); ++it )
    {
        flag = it->flag();

        Q = nssolver.flux(flag);
        P = nssolver.pressure(flag);

        if ( _verbose )
        {
            std::ofstream outfile;
            std::stringstream filenamess;
            std::string filename;

            // file name contains the label
            filenamess << flag;
            // writing down fluxes
            filename = "flux_label" + filenamess.str() + ".m";
            outfile.open(filename.c_str(),std::ios::app);
            outfile << Q << " " << t << "\n";
            outfile.close();
            // writing down pressures
            filename = "pressure_label" + filenamess.str() + ".m";
            outfile.open(filename.c_str(),std::ios::app);
            outfile << P << " " << t << "\n";
            outfile.close();
            // reset ostringstream
            filenamess.str("");
        }
    }

}


struct Cylinder::Private
{
    Private() :
            //check(false),
            nu(1),
            //rho(1),
            H(1), D(1)
            //H(20), D(1)
            //H(0.41), D(0.1)
    {}
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;

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
    double Ubar() const { return nu*Re/D; }

    /**
     * get the magnitude of the profile velocity
     *
     *
     * @return the magnitude of the profile velocity
     */
    double Um_3d() const { return 9*Ubar()/4; }

    double Um_2d() const { return 3*Ubar()/2; }


    /**
     * u3d 3D velocity profile.
     *
     * Define the velocity profile at the inlet for the 3D cylinder
     */
    Real u3d( const Real& /* t */,
              const Real& /* x */,
              const Real& y,
              const Real& z,
              const ID&   id ) const
    {
        if ( id == 1 )
        {
            if ( centered )
            {
                return Um_3d() * (H+y)*(H-y) * (H+z)*(H-z) / pow(H,4);
            }
            else
            {
                return 16 * Um_3d() * y * z * (H-y) * (H-z) / pow(H,4);
            }
        }
        else
        {
            return 0;
        }
    }

    fct_type getU_3d()
    {
        fct_type f;
        f = boost::bind(&Cylinder::Private::u3d, this, _1, _2, _3, _4, _5);
        return f;
    }

    /**
     * u2d flat 2D velocity profile.
     *
     * Define the velocity profile at the inlet for the 2D cylinder
     */
    Real u2d( const Real& t,
              const Real& x,
              const Real& y,
              const Real& z,
              const ID&   id ) const
    {

        switch (id)
        {
        case 1: // x component
            return 0.0;
            break;
        case 3: // z component
            if ( t <= 0.003 )
                return 1.3332e4;
            //      return 0.01;
            return 0.0;
            break;
        case 2: // y component
            return 0.0;
            //      return 1.3332e4;
            //    else
            //      return 0.0;
            break;
        }
        return 0;
    }

    fct_type getU_2d()
    {
        fct_type f;
        f = boost::bind(&Cylinder::Private::u2d, this, _1, _2, _3, _4, _5);
        return f;
    }

    /**
     * one flat (1,1,1)
     *
     * Define the velocity profile at the inlet for the 2D cylinder
     */
    Real poiseuille( const Real& t,
                     const Real& x,
                     const Real& y,
                     const Real& z,
                     const ID&   id ) const
    {
        double result(0.);
        double dist2 = D*D/4. - (x*x + y*y);

        if (id == 3)
        {
            // make sure it is exactely 0 on the cylinder surface.
            result = (fabs(dist2)<1e-5) ? 0. : Um_2d()*2*dist2;
        }

        return result;
    }

    fct_type getU_pois()
    {
        fct_type f;
        f = boost::bind(&Cylinder::Private::poiseuille, this, _1, _2, _3, _4, _5);
        return f;
    }


    Real oneU( const Real& /*t*/,
               const Real& /*x*/,
               const Real& /*y*/,
               const Real& /*z*/,
               const ID&   id ) const
    {
        //            if (id == 3)
        return 10.;

        return -1.;
    }

    fct_type getU_one()
    {
        fct_type f;
        f = boost::bind(&Cylinder::Private::oneU, this, _1, _2, _3, _4, _5);
        return f;
    }


};

Cylinder::Cylinder( int argc,
                    char** argv )
        :
        d( new Private )
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );
    d->data_file_name = data_file_name;

    d->Re          = dataFile( "fluid/problem/Re", 1. );
    d->nu          = dataFile( "fluid/physics/viscosity", 1. ) /
                     dataFile( "fluid/physics/density", 1. );
    d->H           = 20.;//dataFile( "fluid/problem/H", 20. );
    d->D           =               dataFile( "fluid/problem/D", 1. );
    d->centered    = (bool)        dataFile( "fluid/problem/centered", 0 );
    d->initial_sol = (std::string) dataFile( "fluid/problem/initial_sol", "stokes");
    std::cout << d->initial_sol << std::endl;


#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    //    MPI_Init(&argc,&argv);

    int ntasks = 0;
    d->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    if (!d->comm->MyPID())
    {
        std::cout << "My PID = " << d->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
        std::cout << "Re = " << d->Re << std::endl
                  << "nu = " << d->nu << std::endl
                  << "H  = " << d->H  << std::endl
                  << "D  = " << d->D  << std::endl;
    }
//    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    d->comm.reset( new Epetra_SerialComm() );
#endif

}

void
Cylinder::run()

{

    typedef Oseen< RegionMesh3D<LinearTetra> >::vector_Type  vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;
    // Reading from data file
    //
    GetPot dataFile( d->data_file_name );

    //    int save = dataFile("fluid/miscellaneous/save", 1);

    bool verbose = (d->comm->MyPID() == 0);

    // Boundary conditions
    BCHandler bcH;
    BCFunctionBase uZero( zero_scalar );
    std::vector<ID> zComp(1);
    zComp[0] = 3;

    BCFunctionBase uIn  (  d->getU_2d() );
    BCFunctionBase uOne (  d->getU_one() );
    BCFunctionBase uPois(  d->getU_pois() );


    //BCFunctionBase unormal(  d->get_normal() );

    //cylinder

    bcH.addBC( "Inlet",    INLET,    Essential,     Full,     uPois  , 3 );
    bcH.addBC( "Outlet",   OUTLET,   Natural,     Full,     uZero, 3 );
    bcH.addBC( "Wall",     WALL,     Essential,   Full,     uZero, 3 );

    int numLM = 0;

    boost::shared_ptr<DataNavierStokes> dataNavierStokes(new DataNavierStokes());
    dataNavierStokes->setup( dataFile );

    partitionMesh< RegionMesh3D<LinearTetra> >   meshPart;
    meshPart.setup(1, (d->comm));

    HDF5Filter3DMesh<RegionMesh3D<LinearTetra> > HDF5Input(dataFile, "cylinderPart");
    HDF5Input.setComm(d->comm);
//    HDF5Input.loadGraph(meshPart.elementDomains(), d->comm);
    meshPart.elementDomains() = HDF5Input.getGraph();
    meshPart.update();
    meshPart.meshPartition() = HDF5Input.getMeshPartition();
//    HDF5Input.loadMyPartition(meshPart.meshPartition(), d->comm);
    HDF5Input.closeFile();

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "Time discretization order " << dataNavierStokes->dataTime()->orderBDF() << std::endl;

    std::string uOrder =  dataFile( "fluid/space_discretization/vel_order", "P1");
    if (verbose)
        std::cout << "Building the velocity FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > uFESpace(meshPart,uOrder,3,d->comm);

    if (verbose)
        std::cout << "ok." << std::endl;


    std::string pOrder =  dataFile( "fluid/space_discretization/press_order", "P1");

    if (verbose)
        std::cout << "Building the pressure FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > pFESpace(meshPart,pOrder,1,d->comm);

    if (verbose)
        std::cout << "ok." << std::endl;

    UInt totalVelDof   = uFESpace.map().map(Unique)->NumGlobalElements();
    UInt totalPressDof = pFESpace.map().map(Unique)->NumGlobalElements();



    if (verbose) std::cout << "Total Velocity Dof = " << totalVelDof << std::endl;
    if (verbose) std::cout << "Total Pressure Dof = " << totalPressDof << std::endl;

    if (verbose) std::cout << "Calling the fluid constructor ... ";

    bcH.setOffset("Inlet", totalVelDof + totalPressDof);

    Oseen< RegionMesh3D<LinearTetra> > fluid (dataNavierStokes,
                                              uFESpace,
                                              pFESpace,
                                              d->comm, numLM);
    EpetraMap fullMap(fluid.getMap());

    if (verbose) std::cout << "ok." << std::endl;

    fluid.setUp(dataFile);
    fluid.buildSystem();

    MPI_Barrier(MPI_COMM_WORLD);

    // Initialization

    Real dt     = dataNavierStokes->dataTime()->timeStep();
    Real t0     = dataNavierStokes->dataTime()->initialTime();
    Real tFinal = dataNavierStokes->dataTime()->endTime();


    // bdf object to store the previous solutions

    BdfTNS<vector_type> bdf;
    bdf.setup(dataNavierStokes->dataTime()->orderBDF());

    vector_type beta( fullMap );
    vector_type rhs ( fullMap );

#ifdef HAVE_HDF5
    HDF5Filter3DMesh<RegionMesh3D<LinearTetra> > ensight( dataFile, meshPart.meshPartition(), "cylinder", d->comm->MyPID());
#else
    Ensight<RegionMesh3D<LinearTetra> > ensight( dataFile, meshPart.meshPartition(), "cylinder", d->comm->MyPID());
#endif

    vector_ptrtype velAndPressure ( new vector_type(*fluid.solution(), ensight.mapType() ) );

    ensight.addVariable( ExporterData::Vector, "velocity", velAndPressure,
                         UInt(0), uFESpace.dof().numTotalDof() );

    ensight.addVariable( ExporterData::Scalar, "pressure", velAndPressure,
                         UInt(3*uFESpace.dof().numTotalDof() ),
                         UInt(  pFESpace.dof().numTotalDof() ) );

    // initialization with stokes solution

    if (d->initial_sol == "stokes")
    {
        if (verbose) std::cout << std::endl;
        if (verbose) std::cout << "Computing the stokes solution ... " << std::endl << std::endl;

        dataNavierStokes->dataTime()->setTime(t0);

        MPI_Barrier(MPI_COMM_WORLD);

        beta *= 0.;
        rhs  *= 0.;

        fluid.updateSystem(0, beta, rhs );
        fluid.iterate( bcH );

//    fluid.postProcess();

        *velAndPressure = *fluid.solution();
        ensight.postProcess( 0 );
        fluid.resetPreconditioner();
    }

    bdf.bdfVelocity().setInitialCondition( *fluid.solution() );

    // Temporal loop

    Chrono chrono;
    int iter = 1;

    for ( Real time = t0 + dt ; time <= tFinal + dt/2.; time += dt, iter++)
    {

        dataNavierStokes->dataTime()->setTime(time);

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "We are now at time "<< dataNavierStokes->dataTime()->time() << " s. " << std::endl;
            std::cout << std::endl;
        }

        chrono.start();

        double alpha = bdf.bdfVelocity().coefficientFirstDerivative( 0 ) / dataNavierStokes->dataTime()->timeStep();

        beta = bdf.bdfVelocity().extrapolation();
        bdf.bdfVelocity().updateRHSContribution( dataNavierStokes->dataTime()->timeStep());
        rhs  = fluid.matrMass()*bdf.bdfVelocity().rhsContributionFirstDerivative();

        fluid.updateSystem( alpha, beta, rhs );
        fluid.iterate( bcH );

        bdf.bdfVelocity().shiftRight( *fluid.solution() );

//         if (((iter % save == 0) || (iter == 1 )))
//         {
        *velAndPressure = *fluid.solution();
        ensight.postProcess( time );
//         }
//         postProcessFluxesPressures(fluid, bcH, time, verbose);


        MPI_Barrier(MPI_COMM_WORLD);

        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    }

}


//////////////////////


