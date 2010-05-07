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
   \file heart.cpp
   \author Lucia Mirabella <lucia.mirabella@mail.polimi.it> and Mauro Perego <mauro.perego@polimi.it>
   \date 2007-11
   \mantainer R. Ruiz
   \date 2010-04
 */
#include <Epetra_ConfigDefs.h>
#include <Epetra_MpiComm.h>

#include <heart.hpp>

using namespace LifeV;



//! Identifiers for heart boundaries
const int EPICARDIUM    = 40;
const int ENDOCARDIUM   = 60;
const int TRUNC_SEC 	= 50;

Real zero_scalar( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.;
}

Real minus84_scalar( const Real& /* t */,
                     const Real& /* x */,
                     const Real& /* y */,
                     const Real& /* z */,
                     const ID& /* i */ )
{
    return  -84.0;
}


Heart::Heart( int argc,
              char** argv,
              LifeV::AboutData const& /*ad*/,
              LifeV::po::options_description const& /*od*/ )
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile(data_file_name);

    //! Pointer to access functors
    M_heart_fct.reset(new HeartFunctors( dataFile));
    ion_model=dataFile("electric/physics/ion_model",1);
    std::cout << "mpi initialization ... " << std::endl;
    M_heart_fct->comm = new Epetra_MpiComm( MPI_COMM_WORLD );

    if (!M_heart_fct->comm->MyPID()) {
        std::cout << "My PID = " << M_heart_fct->comm->MyPID() << std::endl;
    }
}

void
Heart::run()
{
    Chrono chronoinitialsettings;
    Chrono chronototaliterations;
    chronoinitialsettings.start();
    Real normu;
    Real meanu;
    Real minu;
    //! Construction of data classes

#ifdef MONODOMAIN
    DataMonodomain<RegionMesh3D<LinearTetra> > _data(M_heart_fct);
#else
    DataBidomain<RegionMesh3D<LinearTetra> > _data(M_heart_fct);
#endif
    DataIonic<RegionMesh3D<LinearTetra> > _dataIonic(M_heart_fct->_dataFile);

    bool verbose = (M_heart_fct->comm->MyPID() == 0);

    //! Boundary conditions handler and function
    BCFunctionBase uZero( zero_scalar );
    BCHandler bcH( 3, BCHandler::HINT_BC_NONE );
    bcH.addBC( "Endo",   	ENDOCARDIUM,	Natural,	Full,	uZero,  1 );
    bcH.addBC( "Epi",   	EPICARDIUM, 	Natural,   	Full,   uZero, 	1 );
    bcH.addBC( "Trunc",    	TRUNC_SEC,  	Natural, 	Full,   uZero, 	1 );

    const RefFE*    refFE_w;
    const QuadRule* qR_w;
    const QuadRule* bdQr_w;

    const RefFE*    refFE_u;
    const QuadRule* qR_u;
    const QuadRule* bdQr_u;


    //! Construction of the partitioned mesh
    partitionMesh< RegionMesh3D<LinearTetra> >   meshPart(*_data.mesh(), *M_heart_fct->comm);
    std::string uOrder =  M_heart_fct->_dataFile( "electric/space_discretization/u_order", "P1");

    //! Initialization of the FE type and quadrature rules for both the variables
    if ( uOrder.compare("P1") == 0 )
    {
        if (verbose) std::cout << "P1 potential " << std::flush;
        refFE_u = &feTetraP1;
        qR_u    = &quadRuleTetra15pt;
        bdQr_u  = &quadRuleTria3pt;
    }
    else
    {
        cout<<"\n "<<uOrder<<" finite element not implemented yet \n";
        exit(1);
    }

    Dof uDof(*_data.mesh(), *refFE_u);
    std::string wOrder =  M_heart_fct->_dataFile( "electric/space_discretization/w_order", "P1");
    if ( wOrder.compare("P1") == 0 )
    {
        if (verbose) std::cout << "P1 recovery variable " << std::flush;
        refFE_w = &feTetraP1;
        qR_w    = &quadRuleTetra4pt;
        bdQr_w  = &quadRuleTria3pt;
    }
    else
    {
        cout<<"\n "<<wOrder<<" finite element not implemented yet \n";
        exit(1);
    }

    _data.setMesh(meshPart.mesh());

    //! Construction of the FE spaces
    if (verbose)
        std::cout << "Building the potential FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > uFESpace(meshPart,
                                                             *refFE_u,
                                                             *qR_u,
                                                             *bdQr_u,
                                                             1,
                                                             *M_heart_fct->comm);

#ifdef BIDOMAIN
    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > _FESpace(meshPart,
                                                             *refFE_u,
                                                             *qR_u,
                                                             *bdQr_u,
                                                             2,
                                                             *M_heart_fct->comm);
#endif
    if (verbose)
        std::cout << "ok." << std::endl;
    if (verbose)
        std::cout << "Building the recovery variable FE space ... " << std::flush;
    if (verbose)
        std::cout << "ok." << std::endl;

    UInt totalUDof  = uFESpace.map().getMap(Unique)->NumGlobalElements();
    if (verbose) std::cout << "Total Potential Dof = " << totalUDof << std::endl;
    if (verbose) std::cout << "Calling the electric model constructor ... ";

#ifdef MONODOMAIN
    MonodomainSolver< RegionMesh3D<LinearTetra> > electricModel (_data, uFESpace, bcH, *M_heart_fct->comm);
#else
    BidomainSolver< RegionMesh3D<LinearTetra> > electricModel (_data, _FESpace, uFESpace, bcH, *M_heart_fct->comm);
#endif

    if (verbose) std::cout << "ok." << std::endl;
    EpetraMap fullMap(electricModel.getMap());
    vector_type rhs ( fullMap);
    electricModel.setUp( M_heart_fct->_dataFile );std::cout<<"setup ok"<<std::endl;

    if (verbose) std::cout << "Calling the ionic model constructor ... ";
    boost::shared_ptr< IonicSolver< RegionMesh3D<LinearTetra> > > ionicModel;
    if(ion_model==1)
    {
    	std::cout<<"Ion Model = Rogers-McCulloch"<<std::endl<<std::flush;
    	ionicModel.reset(new Rogers_McCulloch< RegionMesh3D<LinearTetra> >(_dataIonic, uFESpace, *M_heart_fct->comm));
    }else if(ion_model==2)
    {
    	std::cout<<"Ion Model = Luo-Rudy"<<std::endl<<std::flush;
    	ionicModel.reset(new Luo_Rudy< RegionMesh3D<LinearTetra> >(_dataIonic, uFESpace, *M_heart_fct->comm));
    }
#ifdef MONODOMAIN
    electricModel.initialize( M_heart_fct->get_initial_scalar());
#else
    electricModel.initialize( M_heart_fct->get_initial_scalar(), M_heart_fct->get_zero_scalar() );
#endif

    if (verbose) std::cout << "ok." << std::endl;

    ionicModel->initialize( );

    //! Building time-independent part of the system
    electricModel.buildSystem( );std::cout<<"buildsystem ok"<<std::endl;
    //! Initialization
    Real dt     = _data.getTimeStep();
    Real t0     = 0;
    Real tFinal = _data.getEndTime ();
    MPI_Barrier(MPI_COMM_WORLD);

    if (verbose) std::cout << "Setting the initial solution ... " << std::endl << std::endl;
    _data.setTime(t0);
    electricModel.resetPrec();
    if (verbose) std::cout << " ok "<< std::endl;

    //! Setting Ensight postprocessing
    Ensight<RegionMesh3D<LinearTetra> > ensight( M_heart_fct->_dataFile, meshPart.mesh(), "heart", M_heart_fct->comm->MyPID());
    vector_ptrtype Uptr( new vector_type(electricModel.solution_u(), Repeated ) );

    ensight.addVariable( ExporterData::Scalar,  "potential", Uptr,
                         UInt(0), uFESpace.dof().numTotalDof() );

#ifdef BIDOMAIN
    vector_ptrtype Ueptr( new vector_type(electricModel.solution_ue(), Repeated ) );
    ensight.addVariable( ExporterData::Scalar,  "potential_e", Ueptr,
                         UInt(0), uFESpace.dof().numTotalDof() );
#endif

    vector_ptrtype Fptr( new vector_type(electricModel.fiber_vector(), Repeated ) );

    if(_data.has_fibers() )
        ensight.addVariable( ExporterData::Vector, "fibers", Fptr, UInt(0), uFESpace.dof().numTotalDof(), 1 );
    ensight.postProcess( 0 );

    MPI_Barrier(MPI_COMM_WORLD);
    chronoinitialsettings.stop();

    //! Temporal loop
    Chrono chrono;
    int iter = 1;
    chronototaliterations.start();
    for ( Real time = t0 + dt ; time <= tFinal + dt/2.; time += dt, iter++)
    {
        _data.setTime(time);
        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "We are now at time "<< _data.getTime() << " s. " << std::endl;
            std::cout << std::endl;
        }
        chrono.start();
        MPI_Barrier(MPI_COMM_WORLD);
        ionicModel->ionModelSolve( electricModel.solution_u(), _data.getTimeStep() );
        rhs*=0;
        computeRhs( rhs, electricModel, ionicModel, _data );
        //! Updating the PDE system
        electricModel.updatePDESystem( rhs );  //if matrix is time independent (otherwise calling an overload of updatePDESystem)

        //! Solving the system
        electricModel.PDEiterate( bcH );

        normu=electricModel.solution_u().Norm2();
        //       	Real* meanu;
        electricModel.solution_u().getEpetraVector().MeanValue(&meanu);
        //       	Real* maxu;
        electricModel.solution_u().getEpetraVector().MaxValue(&minu);
        if (verbose)
        {
            std::cout << "norm u " << normu << std::endl;
            std::cout << "mean u " << meanu << std::endl;
            std::cout << "max u " << minu << std::endl<<std::flush;
        }

        //! Ensight postprocess
        *Uptr = electricModel.solution_u();
#ifdef BIDOMAIN
        *Ueptr = electricModel.solution_ue();
#endif

        ensight.postProcess( time );
        MPI_Barrier(MPI_COMM_WORLD);
        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
        chronototaliterations.stop();
    }

    if (verbose) std::cout << "Total iterations time " << chronototaliterations.diff() << " s." << std::endl;
    if (verbose) std::cout << "Total initial settings time " << chronoinitialsettings.diff() << " s." << std::endl;
    if (verbose) std::cout << "Total execution time " << chronoinitialsettings.diff()+chronototaliterations.diff() << " s." << std::endl;
}

#ifdef MONODOMAIN
void Heart::computeRhs( vector_type& rhs, MonodomainSolver< RegionMesh3D<LinearTetra> >& electricModel,
                        boost::shared_ptr< IonicSolver< RegionMesh3D<LinearTetra> > > ionicModel, DataMonodomain< RegionMesh3D<LinearTetra> >& data )
{
    bool verbose = (M_heart_fct->comm->MyPID() == 0);
    Real lambda = data.lambda();
    if (verbose) std::cout << "  f-  Computing Rhs ...        "<<"\n"<<std::flush;
    Chrono chrono;
    chrono.start();

    //! u, w with repeated map
    vector_type uVecRep(electricModel.solution_u(), Repeated);
    ionicModel->updateRepeated();
    ElemVec elvec_Iapp( electricModel.potentialFESpace().fe().nbNode, 2 ),
    elvec_u( electricModel.potentialFESpace().fe().nbNode, 1 ),
    elvec_Iion( electricModel.potentialFESpace().fe().nbNode, 1 );

    for(UInt iVol=1; iVol<=electricModel.potentialFESpace().mesh()->numVolumes(); ++iVol)
    {
        electricModel.potentialFESpace().fe().updateJacQuadPt( electricModel.potentialFESpace().mesh()->volumeList( iVol ) );
        elvec_Iapp.zero();
        elvec_u.zero();
        elvec_Iion.zero();
        int ig;
        UInt eleIDu = electricModel.potentialFESpace().fe().currentLocalId();
        UInt nbNode = ( UInt ) electricModel.potentialFESpace().fe().nbNode;

        //! Filling local elvec_u with potential values in the nodes
        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            ig = electricModel.potentialFESpace().dof().localToGlobal( eleIDu, iNode + 1 );
            elvec_u.vec()[ iNode ] = uVecRep[ig];
        }

        ionicModel->updateElvec(eleIDu);
        ionicModel->computeIion(data.Cm(), elvec_Iion, elvec_u, electricModel.potentialFESpace());

        //! Computing the current source of the righthand side
        source(M_heart_fct->get_stim(), elvec_Iapp, electricModel.potentialFESpace().fe(), data.getTime(), 0);
        source(M_heart_fct->get_stim(), elvec_Iapp, electricModel.potentialFESpace().fe(), data.getTime(), 1);

        //! Assembling the righthand side
        for ( int i = 0 ; i < electricModel.potentialFESpace().fe().nbNode ; i++ )
        {
            ig = electricModel.potentialFESpace().dof().localToGlobal( eleIDu, i + 1 );
            rhs.sumIntoGlobalValues (ig, (lambda*elvec_Iapp.vec()[i]+elvec_Iapp.vec()[i+nbNode])/(1+lambda)+data.Chi()*elvec_Iion.vec()[i] );
        }
    }
    rhs.GlobalAssemble();
    Real coeff= data.Chi()*data.Cm()/ data.getTimeStep();
    vector_type tmpvec(electricModel.solution_u());
    tmpvec*=coeff;
    rhs+=electricModel.matrMass()*tmpvec;
    MPI_Barrier(MPI_COMM_WORLD);
    chrono.stop();
    if (verbose) std::cout << "done in " << chrono.diff() << " s." << std::endl;
}
#else
void Heart::computeRhs( vector_type& rhs, BidomainSolver< RegionMesh3D<LinearTetra> >& electricModel,
                        boost::shared_ptr< IonicSolver< RegionMesh3D<LinearTetra> > > ionicModel, DataBidomain< RegionMesh3D<LinearTetra> >& data )
{
    bool verbose = (M_heart_fct->comm->MyPID() == 0);
    if (verbose) std::cout << "  f-  Computing Rhs ...        "<<"\n"<<std::flush;
    Chrono chrono;
    chrono.start();

    //! u, w with repeated map
    vector_type uVecRep(electricModel.solution_u(), Repeated);
    ionicModel->updateRepeated();

    ElemVec elvec_Iapp( electricModel.potentialFESpace().fe().nbNode, 2 ),
        elvec_u( electricModel.potentialFESpace().fe().nbNode, 1 ),
        elvec_Iion( electricModel.potentialFESpace().fe().nbNode, 1 );
    for(UInt iVol=1; iVol<=electricModel.potentialFESpace().mesh()->numVolumes(); ++iVol)
    {
        electricModel.potentialFESpace().fe().updateJacQuadPt( electricModel.potentialFESpace().mesh()->volumeList( iVol ) );
        elvec_u.zero();
        elvec_Iion.zero();
        elvec_Iapp.zero();

        UInt eleIDu = electricModel.potentialFESpace().fe().currentLocalId();
        UInt nbNode = ( UInt ) electricModel.potentialFESpace().fe().nbNode;
        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            int ig = electricModel.potentialFESpace().dof().localToGlobal( eleIDu, iNode + 1 );
            elvec_u.vec()[ iNode ] = uVecRep[ig];
        }

        UInt eleID = electricModel.potentialFESpace().fe().currentLocalId();
        ionicModel->updateElvec(eleID);
        ionicModel->computeIion(data.Cm(), elvec_Iion, elvec_u, electricModel.potentialFESpace());

        //! Computing Iapp
        source(M_heart_fct->get_stim(), elvec_Iapp, electricModel.potentialFESpace().fe(), data.getTime(), 0);
        source(M_heart_fct->get_stim(), elvec_Iapp, electricModel.potentialFESpace().fe(), data.getTime(), 1);
        UInt totalUDof  = electricModel.potentialFESpace().map().getMap(Unique)->NumGlobalElements();

        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            int ig = electricModel.potentialFESpace().dof().localToGlobal( eleIDu, iNode + 1 );
            rhs.sumIntoGlobalValues (ig, elvec_Iapp.vec()[iNode]+data.Chi()*elvec_Iion.vec()[iNode] ); // or
            rhs.sumIntoGlobalValues (ig+totalUDof, -elvec_Iapp.vec()[iNode+nbNode]-data.Chi()*elvec_Iion.vec()[iNode] );
        }
    }
    rhs.GlobalAssemble();

    rhs+=electricModel.matrMass()*data.Chi()*data.Cm()*electricModel.bdf_uiue().time_der(data.getTimeStep());

    MPI_Barrier(MPI_COMM_WORLD);

    chrono.stop();
    if (verbose) std::cout << "done in " << chrono.diff() << " s." << std::endl;
}
#endif
