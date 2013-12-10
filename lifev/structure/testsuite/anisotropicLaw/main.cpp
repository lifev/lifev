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
/**
   \file main.cpp

   This test is the case of traction of a cube. It does not use the symmetry BCs
   This test uses the FESpace which is the standard in LifeV and the ETFESpace
   The FESpace is used for the BCs of Neumann type since in the ET branch there
   is not the integration over the boundary faces.

   \author Paolo Tricerri <paolo.tricerri@epfl.ch>
   \date 1861-03-17
 */

#ifdef TWODIM
#error test_structure cannot be compiled in 2D
#endif

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/fem/TimeAdvance.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialLinear.hpp>
#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialNonLinear.hpp>
#include <lifev/structure/solver/isotropic/ExponentialMaterialNonLinear.hpp>
#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialNonLinearPenalized.hpp>
#include <lifev/structure/solver/isotropic/SecondOrderExponentialMaterialNonLinear.hpp>
#include <lifev/structure/solver/isotropic/NeoHookeanMaterialNonLinear.hpp>


#include <lifev/structure/solver/anisotropic/StructuralAnisotropicConstitutiveLaw.hpp>
#include <lifev/structure/solver/anisotropic/HolzapfelMaterialNonLinear.hpp>
#include <lifev/structure/solver/anisotropic/DistributedHolzapfelMaterialNonLinear.hpp>


#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

//Includes for the Expression Template
#include <lifev/eta/fem/ETFESpace.hpp>

#include <iostream>
#include <fstream>

#include "ud_functions.hpp"


using namespace LifeV;

int returnValue = EXIT_FAILURE;
enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}


std::set<UInt> parseList ( const std::string& list )
{
    std::string stringList = list;
    std::set<UInt> setList;
    if ( list == "" )
    {
        return setList;
    }
    size_t commaPos = 0;
    while ( commaPos != std::string::npos )
    {
        commaPos = stringList.find ( "," );
        setList.insert ( atoi ( stringList.substr ( 0, commaPos ).c_str() ) );
        stringList = stringList.substr ( commaPos + 1 );
    }
    setList.insert ( atoi ( stringList.c_str() ) );
    return setList;
}


class Structure
{
public:

    // Public typedefs
    typedef StructuralOperator< RegionMesh<LinearTetra> >::vector_Type  vector_Type;
    typedef boost::shared_ptr<vector_Type>                              vectorPtr_Type;
    typedef boost::shared_ptr< TimeAdvance< vector_Type > >             timeAdvance_Type;
    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;

    // typedefs for fibers
    // Boost function for fiber direction
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fiberFunction_Type;
    typedef boost::shared_ptr<fiberFunction_Type> fiberFunctionPtr_Type;

    typedef std::vector<fiberFunctionPtr_Type>                          vectorFiberFunction_Type;
    typedef boost::shared_ptr<vectorFiberFunction_Type>                 vectorFiberFunctionPtr_Type;

    typedef std::vector<vectorPtr_Type>                                 listOfFiberDirections_Type;

    /** @name Constructors, destructor
     */
    //@{
    Structure ( int                                   argc,
                char**                                argv,
                boost::shared_ptr<Epetra_Comm>        structComm );

    ~Structure()
    {}
    //@}

    //@{
    void run()
    {
        run3d();
    }
    void CheckResultHolzapfelModel (const Real& dispNorm, const Real& time);
    void CheckResultDistributedModel (const Real& dispNorm, const Real& time);
    void resultChanged (Real time);
    //@}

protected:

private:

    /**
     * run the 2D driven cylinder simulation
     */
    void run2d();

    /**
     * run the 3D driven cylinder simulation
     */
    void run3d();

private:
    struct Private;
    boost::shared_ptr<Private> parameters;
};



struct Structure::Private
{
    Private() :
        rho (1), poisson (1), young (1), bulk (1), alpha (1), gamma (1)
    {}
    double rho, poisson, young, bulk, alpha, gamma;

    std::string data_file_name;

    boost::shared_ptr<Epetra_Comm>     comm;
};



Structure::Structure ( int                                   argc,
                       char**                                argv,
                       boost::shared_ptr<Epetra_Comm>        structComm) :
    parameters ( new Private() )
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );
    parameters->data_file_name = data_file_name;

    parameters->comm = structComm;
    int ntasks = parameters->comm->NumProc();

    if (!parameters->comm->MyPID() )
    {
        std::cout << "My PID = " << parameters->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
    }
}



void
Structure::run2d()
{
    std::cout << "2D cylinder test case is not available yet\n";
}



void
Structure::run3d()
{
    bool verbose = (parameters->comm->MyPID() == 0);

    //! Number of boundary conditions for the velocity and mesh motion
    boost::shared_ptr<BCHandler> BCh ( new BCHandler() );

    //! dataElasticStructure
    GetPot dataFile ( parameters->data_file_name.c_str() );

    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    dataStructure->showMe();

    MeshData             meshData;
    meshData.setup (dataFile, "solid/space_discretization");

    boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr (new RegionMesh<LinearTetra> (  parameters->comm  ) );
    readMesh (*fullMeshPtr, meshData);

    MeshPartitioner< RegionMesh<LinearTetra> > meshPart ( fullMeshPtr, parameters->comm );

    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");

    //Mainly used for BCs assembling (Neumann type)
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (meshPart, dOrder, 3, parameters->comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (meshPart, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );

    // Setting the fibers
    vectorFiberFunctionPtr_Type pointerToVectorOfFamilies( new vectorFiberFunction_Type( ) );
    (*pointerToVectorOfFamilies).resize( dataStructure->numberFibersFamilies( ) );

    listOfFiberDirections_Type fiberDirections;
    fiberDirections.resize( dataStructure->numberFibersFamilies( ) );

    std::cout << "Size of the number of families: " << (*pointerToVectorOfFamilies).size() << std::endl;

    fibersDirectionList setOfFiberFunctions;
    setOfFiberFunctions.setupFiberDefinitions( dataStructure->numberFibersFamilies( ) );

    if (verbose)
    {
        std::cout << std::endl;
    }

    std::string timeAdvanceMethod =  dataFile ( "solid/time_discretization/method", "BDF");

    timeAdvance_Type  timeAdvance ( TimeAdvanceFactory::instance().createObject ( timeAdvanceMethod ) );

    UInt OrderDev = 2;

    timeAdvance->setup (dataStructure->dataTimeAdvance()->orderBDF() , OrderDev);

    timeAdvance->setTimeStep (dataStructure->dataTime()->timeStep() );
    //timeAdvance->showMe();

    //! #################################################################################
    //! BOUNDARY CONDITIONS
    //! #################################################################################
    vector <ID> compx (1), compy (1), compz (1), compxy (2), compxz (2), compyz (2);
    compx[0] = 0;
    compy[0] = 1, compz[0] = 2;
    compxy[0] = 0;
    compxy[1] = 1;
    compxz[0] = 0;
    compxz[1] = 2;
    compyz[0] = 1;
    compyz[1] = 2;

    BCFunctionBase zero (bcZero);
    BCFunctionBase nonZero;

    nonZero.setFunction (bcNonZero);

    //! =================================================================================
    //! BC for StructuredCube4_test_structuralsolver.mesh
    //! =================================================================================
    BCh->addBC ("EdgesIn",      20,  Natural,   Component, nonZero, compy);
    BCh->addBC ("EdgesIn",      40,  Essential, Component, zero,    compy);

    //! Symmetry BC
    BCh->addBC ("EdgesIn",      50,   EssentialVertices, Component, zero, compxy);
    BCh->addBC ("EdgesIn",      30,   EssentialVertices, Component, zero, compyz);
    BCh->addBC ("EdgesIn",      80,   EssentialVertices, Component, zero, compxz);
    BCh->addBC ("EdgesIn",      100,  EssentialVertices,  Full, zero, 3);

    BCh->addBC ("EdgesIn",      7, Essential, Component , zero, compx);
    BCh->addBC ("EdgesIn",      3, Essential, Component , zero, compz);
    //! =================================================================================

    //! 1. Constructor of the structuralSolver
    StructuralOperator< RegionMesh<LinearTetra> > solid;

    //! 2. Setup of the structuralSolver
    solid.setup (dataStructure,
                 dFESpace,
                 dETFESpace,
                 BCh,
                 parameters->comm);

    //! 3. Setting data from getPot
    solid.setDataFromGetPot (dataFile);

    // Setting the vector of fibers functions
    for( UInt k(1); k <= pointerToVectorOfFamilies->size( ); k++ )
    {
        // Setting up the name of the function to define the family
        std::string family="Family";
        // adding the number of the family
        std::string familyNumber;
        std::ostringstream number;
        number << ( k );
        familyNumber = number.str();

        // Name of the function to create
        std::string creationString = family + familyNumber;
        (*pointerToVectorOfFamilies)[ k-1 ].reset( new fiberFunction_Type() );
        (*pointerToVectorOfFamilies)[ k-1 ] = setOfFiberFunctions.fiberDefinition( creationString );

        fiberDirections[ k-1 ].reset( new vector_Type(solid.displacement(), Unique) );
    }

    //! 3.b Setting the fibers in the abstract class of Anisotropic materials
    solid.material()->anisotropicLaw()->setupFiberDirections( pointerToVectorOfFamilies );

    //! 4. Building system using TimeAdvance class
    double timeAdvanceCoefficient = timeAdvance->coefficientSecondDerivative ( 0 ) / (dataStructure->dataTime()->timeStep() * dataStructure->dataTime()->timeStep() );
    solid.buildSystem (timeAdvanceCoefficient);


    //dataStructure->showMe();
    //! =================================================================================
    //! Temporal data and initial conditions
    //! =================================================================================

    //! 5. Initial data
    Real dt = dataStructure->dataTime()->timeStep();

    vectorPtr_Type rhs (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type disp (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type vel (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type acc (new vector_Type (solid.displacement(), Unique) );

    if (verbose)
    {
        std::cout << "S- initialization ... ";
    }

    std::vector<vectorPtr_Type> uv0;

    if (timeAdvanceMethod == "Newmark")
    {
        uv0.push_back (disp);
        uv0.push_back (vel);
        uv0.push_back (acc);
    }

    vectorPtr_Type initialDisplacement (new vector_Type (solid.displacement(), Unique) );

    if ( !dataStructure->solidTypeIsotropic().compare ("secondOrderExponential") )
    {
        dFESpace->interpolate ( static_cast<solidFESpace_Type::function_Type> ( d0 ), *initialDisplacement, 0.0 );
    }

    if (timeAdvanceMethod == "BDF")
    {
        Real tZero = dataStructure->dataTime()->initialTime();

        for ( UInt previousPass = 0; previousPass < timeAdvance->size() ; previousPass++)
        {
            Real previousTimeStep = tZero - previousPass * dt;
            std::cout << "BDF " << previousTimeStep << "\n";
            if ( !dataStructure->solidTypeIsotropic().compare ("secondOrderExponential") )
            {
                uv0.push_back (initialDisplacement);
            }
            else
            {
                uv0.push_back (disp);
            }
        }
    }

    timeAdvance->setInitialCondition (uv0);

    timeAdvance->setTimeStep ( dt );

    timeAdvance->updateRHSContribution ( dt );

    if ( !dataStructure->solidTypeIsotropic().compare ("secondOrderExponential") )
    {
        solid.initialize ( initialDisplacement );
    }
    else
    {
        solid.initialize ( disp );
    }

    MPI_Barrier (MPI_COMM_WORLD);

    if (verbose )
    {
        std::cout << "ok." << std::endl;
    }

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile ( "exporter/type", "ensight");
#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "structure" ) );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporter.reset ( new ExporterEmpty<RegionMesh<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID() ) );
        }

        else
        {
            exporter.reset ( new ExporterEnsight<RegionMesh<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID() ) );
        }
    }

    exporter->setPostDir ( "./" );
    exporter->setMeshProcId ( meshPart.meshPartition(), parameters->comm->MyPID() );

    vectorPtr_Type solidDisp ( new vector_Type (solid.displacement(),  exporter->mapType() ) );
    vectorPtr_Type solidVel  ( new vector_Type (solid.displacement(),  exporter->mapType() ) );
    vectorPtr_Type solidAcc  ( new vector_Type (solid.displacement(),  exporter->mapType() ) );
    //vectorPtr_Type rhsVector ( new vector_Type (solid.displacement(),  exporter->mapType() ) );

    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "velocity",     dFESpace, solidVel,  UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "acceleration", dFESpace, solidAcc,  UInt (0) );
    //exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "forcingTerm", dFESpace, rhsVector,  UInt (0) );

    // Adding the fibers vectors
    // Setting the vector of fibers functions
    for( UInt k(1); k <= pointerToVectorOfFamilies->size( ); k++ )
    {
        // Setting up the name of the function to define the family
        std::string family="Family-";
        // adding the number of the family
        std::string familyNumber;
        std::ostringstream number;
        number << ( k );
        familyNumber = number.str();

        // Name of the function to create
        std::string creationString = family + familyNumber;
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationString, dFESpace, fiberDirections[ k-1 ], UInt (0) );

        // Extracting the fibers vectors
        *(fiberDirections[ k-1 ]) = solid.material()->anisotropicLaw()->ithFiberVector( k );
    }


    exporter->postProcess ( 0 );
    cout.precision(16);

    // int IDPointX = 618;
    // int IDPointY = 331;

    /*
    //!--------------------------------------------------------------------------------------------
    //! MATLAB FILE WITH DISPLACEMENT OF A CHOSEN POINT
    //!--------------------------------------------------------------------------------------------

    ofstream file_comp( "Displacement_components_NL.m" );
    if ( !file_comp )
    {
      std::cout <<" Unable to open file! You need to specify the output folder in the data file " << std::endl;
    }


    //int IDPoint = 401; // StructuredCube8
    //int IDPoint = 2593; // StructuredCube16

    //int IDPoint = 74;// cube4x4.mesh
    //int IDPoint = 315;// cube8x8.mesh
    //int IDPoint = 1526;// cube16x16.mesh

    file_comp << " % TEST NONLINEAR ELASTICITY" << endl;
    file_comp << " % Displacement components of ID node  " << IDPoint << " :" << endl;
    file_comp << " % Each row is a time step" << endl;
    file_comp << " % First column = comp x, second = comp y and third = comp z. " << endl;
    file_comp <<  endl;
    file_comp << " SolidDisp_NL = [ " ;

    for ( UInt k = IDPoint - 1; k <= solid.displacement().size() - 1; k = k + solid.displacement().size()/nDimensions )
    {
    file_comp<< solid.displacement()[ k ] << " ";
    }

    file_comp<< endl;
    */
    //!--------------------------------------------------------------------------------------------
    //!The update of the RHS is done by the TimeAdvance class
    //solid.updateSystem();
    //! =================================================================================

    Real normVect;
    normVect =  solid.displacement().norm2();

    if( verbose )
    {
        std::cout << "The norm 2 of the displacement field is: " << normVect << std::endl;
    }

    //! =============================================================================
    //! Temporal loop
    //! =============================================================================
    //    for (Real time = dt; time <= T; time += dt)
    for (dataStructure->dataTime()->setTime ( dt ) ; dataStructure->dataTime()->canAdvance( ); dataStructure->dataTime()->updateTime( ) )
    {
        returnValue = EXIT_FAILURE;

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "S- Now we are at time " << dataStructure->dataTime()->time() << " s." << std::endl;
        }

        //! 6. Updating right-hand side
        *rhs *= 0;
        timeAdvance->updateRHSContribution ( dt );
        *rhs += *solid.massMatrix() * timeAdvance->rhsContributionSecondDerivative() / timeAdvanceCoefficient;
        solid.setRightHandSide ( *rhs );

        //! 7. Iterate --> Calling Newton
        solid.iterate ( BCh );

        timeAdvance->shiftRight ( solid.displacement() );

        *solidDisp = solid.displacement();
        *solidVel  = timeAdvance->firstDerivative();
        *solidAcc  = timeAdvance->secondDerivative();

        // This vector is to export the forcing term
        //*rhsVector = solid.rhsCopy();

        exporter->postProcess ( dataStructure->dataTime()->time() );

        /* This part lets to save the displacement at one point of the mesh and to check the result
           w.r.t. manufactured solution.
           //!--------------------------------------------------------------------------------------------------
           //! MATLAB FILE WITH DISPLACEMENT OF A CHOOSEN POINT
           //!--------------------------------------------------------------------------------------------------
           cout <<"*******  DISPLACEMENT COMPONENTS of ID node "<< IDPoint << " *******"<< std::endl;
           for ( UInt k = IDPoint - 1 ; k <= solid.displacement().size() - 1; k = k + solid.displacement().size()/nDimensions )
           {
           file_comp<< solid.displacement()[ k ] << " ";
           cout.precision(16);
           cout <<"*********************************************************"<< std::endl;
           cout <<" solid.disp()[ "<< k <<" ] = "<<  solid.displacement()[ k ]  << std::endl;
           cout <<"*********************************************************"<< std::endl;
           }
           file_comp<< endl;

           cout <<"*********************************************************"<< std::endl;
           cout <<" solid.disp()[ "<< IDPointX - 1  <<" ] = "<<  solid.displacement()[ IDPointX - 1 ]  << std::endl;
           cout <<" solid.disp()[ "<< IDPointX - 1 + dFESpace->dof().numTotalDof() <<" ] = "<<  solid.displacement()[ IDPointX - 1 + dFESpace->dof().numTotalDof() ]  << std::endl;
           cout <<"*********************************************************"<< std::endl;
        */

        Real normVect;
        normVect =  solid.displacement().norm2();
        std::cout << "The norm 2 of the displacement field is: " << normVect << std::endl;

        // Check results
        if ( !dataStructure->solidTypeAnisotropic().compare ("holzapfel") )
            CheckResultHolzapfelModel (normVect, dataStructure->dataTime()->time() );
        else
            CheckResultDistributedModel (normVect, dataStructure->dataTime()->time() );

        //!--------------------------------------------------------------------------------------------------

        MPI_Barrier (MPI_COMM_WORLD);
    }

}

void Structure::CheckResultHolzapfelModel (const Real& dispNorm, const Real& time)
{
    if ( time == 0.05  && std::fabs (dispNorm - 0.84960668) <= 1e-7 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.1   && std::fabs (dispNorm - 0.84981715) <= 1e-7 )
    {
        this->resultChanged (time);
    }
}

void Structure::CheckResultDistributedModel (const Real& dispNorm, const Real& time)
{
    if ( time == 0.05  && std::fabs (dispNorm - 1.7747561302) <= 1e-7 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.1   && std::fabs (dispNorm - 1.7757263996) <= 1e-7 )
    {
        this->resultChanged (time);
    }
}

void Structure::resultChanged (Real time)
{
    std::cout << "Correct value at time: " << time << std::endl;
    returnValue = EXIT_SUCCESS;
}



int
main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_MpiComm> Comm (new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    if ( Comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }
#else
    boost::shared_ptr<Epetra_SerialComm> Comm ( new Epetra_SerialComm() );
    cout << "% using serial Version" << endl;
#endif

    Structure structure ( argc, argv, Comm );
    structure.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return returnValue ;
}
