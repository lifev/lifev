/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <lifemc/lifesolver/Monolithic.hpp>

namespace LifeV
{



Monolithic::Monolithic():
    super(),
    M_monolithicMap(),
    M_interfaceMap(),
    M_beta(),
    M_monolithicMatrix(),
    M_precPtr(),
    M_rhsFull(),
    M_BCh_flux(),
    M_BCh_Robin(),
    M_fluxes(0),
    M_BChWSS(),
    M_offset(0),
    M_solidAndFluidDim(0),
    M_fluidBlock(),
    M_solidBlock(),
    M_solidBlockPrec(),
    M_linearSolver(),
    M_numerationInterface(),
    M_BChs(),
    M_FESpaces(),
    //end of protected attributes
    M_PAAP(),
#ifdef OBSOLETE
    M_rhsShapeDerivatives(),
#endif
    M_diagonalScale(false),
    M_reusePrec(true),
    M_resetPrec(true),
    M_maxIterSolver(-1),
    M_restarts(false)
{}

// Destructor

Monolithic::~Monolithic()
{
}


void
Monolithic::setupFEspace()
{
	super::setupFEspace();

	// Monolitic: In the beginning I need a non-partitioned mesh. later we will do the partitioning
    M_dFESpace.reset( new FESpace<mesh_type, EpetraMap>( M_solidMesh,
                                                         M_data->dataSolid()->order(),
                                                         nDimensions,
                                                         M_epetraComm));
}


void
Monolithic::setupDOF( void )
{

    M_dofStructureToHarmonicExtension    .reset( new DofInterface3Dto3D );

	M_dofStructureToHarmonicExtension->setup(   M_uFESpace->refFE(), M_uFESpace->dof(),
											    M_dFESpace->refFE(), M_dFESpace->dof() );
	M_dofStructureToHarmonicExtension->update( *M_uFESpace->mesh(),  M_data->structureInterfaceFlag(),
											   *M_dFESpace->mesh(),  M_data->harmonicInterfaceFlag(),
											    M_data->interfaceTolerance() );

                                               createInterfaceMaps(M_dofStructureToHarmonicExtension);
}


void
Monolithic::setupFluidSolid()
{
    super::setupFluidSolid();

    // Note: up to now it works only with matching grids (and poly order) on the interface
    assert(M_fluidInterfaceMap->getMap(Unique)->NumGlobalElements() == M_solidInterfaceMap->getMap(Unique)->NumGlobalElements());

    M_interfaceMap = M_solidInterfaceMap;

    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;

    M_monolithicMap.reset(new EpetraMap(M_uFESpace->map()));
    M_BCh_flux = M_BCh_u; // For the moment M_BCh_u contains only the fluxes.
    M_fluxes=M_BCh_flux->size( );

    std::string opertype = M_dataFile("problem/blockOper", "AdditiveSchwarz");

    createOperator( opertype);

    *M_monolithicMap+= M_pFESpace->map();
    *M_monolithicMap+= M_fluxes;
    *M_monolithicMap+= M_dFESpace->map();

    M_monolithicMatrix->createInterfaceMap( *M_interfaceMap, M_dofStructureToHarmonicExtension->locDofMap(), M_dFESpace->map().getMap(Unique)->NumGlobalElements()/nDimensions, M_epetraWorldComm );
    *M_monolithicMap += *M_monolithicMatrix->getInterfaceMap();

    //the map for the interface coupling matrices should be done with respect to the coarser mesh.
    M_monolithicMatrix->getNumerationInterface(M_numerationInterface);
    M_beta.reset  (new vector_type(/*M_monolithicMap*/M_uFESpace->map()));

    M_offset = M_uFESpace->dof().numTotalDof()*nDimensions + M_fluxes +  M_pFESpace->dof().numTotalDof();
    M_solidAndFluidDim= M_offset + M_dFESpace->dof().numTotalDof()*nDimensions;
    M_BCh_d->setOffset(M_offset);
    M_BCh_flux->setOffset(M_offset-M_fluxes);
    std::vector<BCBase>::iterator fluxIt = M_BCh_flux->begin( );
    for ( UInt i = 0; i < M_fluxes; ++i, ++fluxIt )
        fluxIt->setOffset( i );

    M_BCh_Robin = M_BCh_d;// For the moment M_BCh_d contains only the Robin b.c..
    M_BCh_Robin->setOffset(M_offset);

}//end setup



void
Monolithic::setDataFile( const GetPot& dataFile )
{
    super::setDataFile( dataFile );

    M_diagonalScale    = dataFile( "linear_system/prec/diagonalScaling",  false );
    M_restarts         = dataFile( "exporter/start"  ,  0   );
}


void
Monolithic::buildSystem()
{
    M_solidBlock.reset(new matrix_type(*M_monolithicMap, 1));//since it is constant, we keep this throughout the simulation
    M_solid->buildSystem(M_solidBlock);
    M_solidBlock->GlobalAssemble();
    *M_solidBlock *= (M_data->dataSolid()->dataTime()->getTimeStep()*M_solid->rescaleFactor());
    M_solid->rescaleMatrices();
}



void
Monolithic::couplingRhs(vector_ptrtype rhs, vector_ptrtype un) // not working with non-matching grids
{
    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;
    //    UInt solidDim=M_dFESpace->map().getMap(Unique)->NumGlobalElements()/nDimensions;

    vector_type lambda(*M_interfaceMap, Unique);
    this->monolithicToInterface(lambda, *un);
    UInt interface(M_monolithicMatrix->getInterface());
    Real rescale(M_solid->rescaleFactor());
    UInt totalDofs(M_dFESpace->dof().numTotalDof());


    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap->getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                if(rhs.get())
                    (*rhs)[  (int)(*M_numerationInterface)[ITrow->second ] + dim*interface +M_solidAndFluidDim ] = -lambda( ITrow->second + dim*totalDofs )*rescale;
            }
        }
    }
}


void
Monolithic::updateSystem()
{

    vector_type solution(*this->M_monolithicMap);
    monolithicToX(*this->M_un, solution, M_uFESpace->map(), UInt(0));
    this->M_bdf->shift_right(solution);

    M_meshMotion->updateSystem();

    this->fluid().updateUn(*this->M_un);
    *M_rhs*=0;
    *M_rhsFull*=0;
    this->M_fluid->resetStab();
}


void
Monolithic::monolithicToX(const vector_type& disp, vector_type& dispFluid, EpetraMap& map, UInt offset)
{
    if(disp.getMaptype()== Repeated)
    {
        vector_type dispUnique(disp, Unique);
        monolithicToX(dispUnique, dispFluid, map, offset);
        dispFluid = dispUnique;
        return;
    }
    dispFluid.subset(disp, map, offset, offset);
}


void
Monolithic::monolithicToInterface(vector_type& lambdaSolid, const vector_type& disp)
{
    if (disp.getMaptype() == Repeated)
    {
        vector_type const  dispUnique(disp, Unique);
        monolithicToInterface(lambdaSolid, dispUnique);
        return;
    }
    if (lambdaSolid.getMaptype() == Repeated)
    {
        vector_type  lambdaSolidUn(lambdaSolid.getMap(), Unique);
        monolithicToInterface( lambdaSolidUn, disp);
        lambdaSolid = lambdaSolidUn;
        return;
    }
    /* UInt MyOffset(M_uFESpace->map().getMap(Unique)->NumMyElements()+M_pFESpace->map().getMap(Unique)->NumMyElements());
       vector_type subDisp(this->M_dFESpace->map(), Unique);
       subDisp.mySubset(disp, MyOffset);
       lambdaSolid=subDisp;*/

    EpetraMap subMap(*disp.getMap().getMap(Unique), M_offset,disp.getMap().getMap(Unique)->NumGlobalElements() );
    vector_type subDisp(subMap, Unique);
    subDisp.subset(disp, M_offset);
    lambdaSolid=subDisp;
}



void
Monolithic::setDispSolid(const vector_type &sol)
{
    vector_type disp(*M_monolithicMap);
    monolithicToX(sol, disp, M_dFESpace->map(), M_offset);
    this->M_solid->setDisp(disp);
}



void
Monolithic::evalResidual( vector_type&       res,
                          const vector_type& disp,
                          const UInt          iter )
{

    if((iter==0)|| !this->M_data->dataFluid()->isSemiImplicit())
    {
        setDispSolid(disp);
        vector_type lambdaFluid(*M_interfaceMap, Unique);

        monolithicToInterface(lambdaFluid, disp);

        lambdaFluid *= (M_data->dataFluid()->dataTime()->getTimeStep()*(M_solid->rescaleFactor()));//because of the matrix scaling
        this->setLambdaFluid(lambdaFluid); // it must be _disp restricted to the interface
        M_meshMotion->iterate(*M_BCh_mesh);
        M_meshMotion->updateDispDiff();

        M_beta.reset(new vector_type(M_uFESpace->map()));
        vector_type meshDispDiff( M_meshMotion->disp(), Repeated );

        this->moveMesh(meshDispDiff);//initialize the mesh position with the total displacement

        meshDispDiff=M_meshMotion->dispDiff();//repeating the mesh dispDiff
        this->interpolateVelocity(meshDispDiff, *this->M_beta);

        double alpha = 1./M_data->dataFluid()->dataTime()->getTimeStep();//mesh velocity w

        *this->M_beta *= -alpha;

        vector_ptrtype fluid(new vector_type(this->M_uFESpace->map()));
        fluid->subset(*M_un, (UInt)0);
        *this->M_beta += *fluid/*M_un*/;//relative velocity beta=un-w

        //M_monolithicMatrix.reset(new matrix_type(*M_monolithicMap));

        assembleFluidBlock(iter);

        M_solidBlockPrec.reset(new matrix_type(*M_monolithicMap, 1));
        *M_solidBlockPrec += *M_solidBlock;

        if ( !M_BCh_u->bdUpdateDone() )
            M_BCh_u->bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
        M_BCh_d->setOffset(M_offset);
        if ( !M_BCh_d->bdUpdateDone() )
            M_BCh_d->bdUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );

        M_monolithicMatrix->setRobin( M_robinCoupling, M_rhsFull );
        M_precPtr->setRobin(M_robinCoupling, M_rhsFull);

        if(!this->M_monolithicMatrix->set())
        {
            M_BChs.push_back(M_BCh_d);
            M_BChs.push_back(M_BCh_u);
            M_FESpaces.push_back(M_dFESpace);
            M_FESpaces.push_back(M_uFESpace);

            M_monolithicMatrix->push_back_matrix(M_solidBlockPrec, false);
            M_monolithicMatrix->push_back_matrix(M_fluidBlock, true);
            M_monolithicMatrix->setConditions(M_BChs);
            M_monolithicMatrix->setSpaces(M_FESpaces);
            M_monolithicMatrix->setOffsets(2, M_offset, 0);
            M_monolithicMatrix->coupler(M_monolithicMap, M_dofStructureToHarmonicExtension->locDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->getTimeStep());
        }
        else
        {
            M_monolithicMatrix->replace_matrix(M_fluidBlock, 1);
            M_monolithicMatrix->replace_matrix(M_solidBlockPrec, 0);
        }

        M_monolithicMatrix->blockAssembling();
        M_monolithicMatrix->applyBoundaryConditions(dataFluid().dataTime()->getTime(), M_rhsFull);

        M_monolithicMatrix->GlobalAssemble();
        //M_monolithicMatrix->getMatrix()->spy("M");

        //NOTE: M_monolithic->GlobalAssemble has to be called before M_precPtr->blockAssembling(), because they hold
        //shared pointers to the same blocks

        M_nbEval++ ;
        //M_precPtr->reset();
        //M_precPtr->resetBlocks();
        setupBlockPrec( *M_rhsFull );

        M_precPtr->blockAssembling();
        M_precPtr->applyBoundaryConditions(dataFluid().dataTime()->getTime());
        M_precPtr->GlobalAssemble();

    }
    evalResidual( disp,  M_rhsFull, res, M_diagonalScale);
}

int Monolithic::setupBlockPrec(vector_type& /*rhs*/)
{
     if(!(M_precPtr->set()))
     {
        M_precPtr->push_back_matrix(M_solidBlockPrec, false);
        M_precPtr->push_back_matrix(M_fluidBlock, true);
        M_precPtr->setConditions(M_BChs);
        M_precPtr->setSpaces(M_FESpaces);
        M_precPtr->setOffsets(2, M_offset, 0);
        M_precPtr->coupler(M_monolithicMap, M_dofStructureToHarmonicExtension->locDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->getTimeStep());
     }
     else
     {
         M_precPtr->replace_matrix(M_fluidBlock, 1);
         M_precPtr->replace_matrix(M_solidBlockPrec, 0);
     }
}

void
Monolithic::
evalResidual( const vector_type& sol,const vector_ptrtype& rhs, vector_type& res, bool diagonalScaling)
{
    //M_monolithicMatrix->getMatrix()->GlobalAssemble();
    if(diagonalScaling)
        diagonalScale(*rhs, M_monolithicMatrix->getMatrix());
    res = *(M_monolithicMatrix->getMatrix())*sol;
    res -= *rhs;
    // Ax-b
}




void
Monolithic::solveJac(vector_type         &_step,
                     const vector_type   &_res,
                     const Real         /*_linearRelTol*/)
{
    M_solid->getDisplayer().leaderPrint("  M-  Jacobian NormInf res:                    ", _res.NormInf(), "\n");
    M_solid->getDisplayer().leaderPrint("  M-  Solving Jacobian system ...              \n" );

    this->iterateMonolithic(_res, _step);

    M_solid->getDisplayer().leaderPrint("  M-  Jacobian NormInf res:                    ", _step.NormInf(), "\n");
}


void
Monolithic::iterateMesh(const vector_type& disp)
{
    vector_type lambdaFluid(*M_interfaceMap, Unique);

    monolithicToInterface(lambdaFluid, disp);

    lambdaFluid *= (M_data->dataFluid()->dataTime()->getTimeStep()*(M_solid->rescaleFactor()));

    this->setLambdaFluid(lambdaFluid); // it must be _disp restricted to the interface

    M_meshMotion->iterate(*M_BCh_mesh);

}


void
Monolithic::variablesInit(const std::string& dOrder)
{
    //    EpetraMap interfaceMap(*M_solidInterfaceMap);
    //M_solidMeshPart.reset( new  partitionMesh< FSIOperator::mesh_type > (*M_solidMesh, *M_epetraComm/*, M_solidInterfaceMap->getMap(Unique).get(), M_solidInterfaceMap->getMap(Repeated).get()*/));

    M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(*M_solidMeshPart,
                                                       dOrder,
                                                       3,
                                                       M_epetraComm));

    // INITIALIZATION OF THE VARIABLES
    M_lambdaFluid.reset(new vector_type(*M_fluidInterfaceMap, Unique) );
    M_lambdaFluidRepeated.reset(new vector_type(*M_fluidInterfaceMap, Repeated) );
}



void
Monolithic::
updateSolidSystem( vector_ptrtype & rhsFluidCoupling )
{
    M_solid->updateSystem();
    *rhsFluidCoupling += *M_solid->rhsWithoutBC();
}


void
Monolithic::setupSystem( )
{
    M_fluid->setUp( M_dataFile );
    setUp( M_dataFile );
}


void
Monolithic::setUp( const GetPot& dataFile )
{

    M_linearSolver.reset(new solver_type(M_epetraComm));

    M_linearSolver->setDataFromGetPot( dataFile, "linear_system/solver" );
    M_linearSolver->setUpPrec(dataFile, "linear_system/prec");//to avoid if we have already a prec.
    std::string prectype = dataFile("problem/DDBlockPrec", "NULL");
    std::string opertype = dataFile("problem/blockOper", "NULL");

    M_precPtr.reset(BlockPrecFactory::instance().createObject( prectype ));

    std::string section("linear_system/prec");
    M_precPtr->setComm(M_epetraComm);
    M_precPtr->setDataFromGetPot(dataFile, section);//to avoid if we build the prec from a matrix.
    M_monolithicMatrix->setComm(M_epetraComm);
    M_monolithicMatrix->setDataFromGetPot(dataFile, section);//to avoid if we build the prec from a matrix.
    M_reusePrec     = dataFile( "linear_system/prec/reuse", true);
    M_maxIterSolver = dataFile( "linear_system/solver/max_iter", -1);
}


void
Monolithic::
diagonalScale(vector_type& rhs, matrix_ptrtype matrFull)
{
    Epetra_Vector diagonal(*rhs.getMap().getMap(Unique));
    //M_matrFull->getMatrixPtr()->InvRowSums(diagonal);
    //M_matrFull->getMatrixPtr()->InvRowMaxs(diagonal);
    //M_matrFull->getMatrixPtr()->InvColSums(diagonal);
    matrFull->getMatrixPtr()->InvColMaxs(diagonal);
    matrFull->getMatrixPtr()->LeftScale(diagonal);
    rhs.getEpetraVector().Multiply(1, rhs.getEpetraVector(), diagonal,0);
}



void
Monolithic::solidInit(std::string const& dOrder)
{   // Monolitic: In the beginning I need a non-partitioned mesh. later we will do the partitioning
    M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(M_solidMesh,
                                                       dOrder,
                                                       nDimensions,
                                                       M_epetraComm));
}






void
Monolithic::initialize( FSIOperator::fluid_type::value_type::Function const& u0,
                             FSIOperator::solid_type::value_type::Function const& p0,
                             FSIOperator::solid_type::value_type::Function const& d0,
                             FSIOperator::solid_type::value_type::Function const& w0,
                             FSIOperator::solid_type::value_type::Function const& /*w0*/ )
{
    vector_type u(M_uFESpace->map());
    M_uFESpace->interpolate(u0, u, M_data->dataFluid()->dataTime()->getTime());

    vector_type p(M_pFESpace->map());
    M_pFESpace->interpolate(p0, p, M_data->dataFluid()->dataTime()->getTime());

    vector_type d(M_dFESpace->map());
    M_dFESpace->interpolate(d0, d, M_data->dataSolid()->dataTime()->getTime());

    initialize(u, p, d);
}


void
Monolithic::initialize( const vector_type& u0, const vector_type& p0, const vector_type& d0)
{
    *M_un=u0;
    M_un->add(p0, nDimensions*M_uFESpace->dof().numTotalDof());
    M_un->add(d0, M_offset);
}

#ifdef HAVE_TRILINOS_ANASAZI

Real&
Monolithic::computeMaxSingularValue( )
{
    typedef Epetra_Operator                                                operator_type;

    M_PAAP.reset(new ComposedPreconditioner<Epetra_Operator>(M_epetraComm));

    boost::shared_ptr<operator_type>  ComposedPrecPtr(M_linearSolver->getPrec()->getPrec());

    //M_monolithicMatrix->getMatrixPtr()->OptimizeStorage();
    boost::shared_ptr<Epetra_FECrsMatrix> matrCrsPtr(new Epetra_FECrsMatrix(*M_monolithicMatrix->getMatrix()->getMatrixPtr()));

    M_PAAP->push_back(boost::dynamic_pointer_cast<operator_type>(/*ComposedPrecPtr*/matrCrsPtr));
    M_PAAP->push_back(boost::dynamic_pointer_cast<operator_type>(ComposedPrecPtr/*matrCrsPtr*/),  true);
    M_PAAP->push_back(boost::dynamic_pointer_cast<operator_type>(ComposedPrecPtr/*matrCrsPtr*/), true, true);
    M_PAAP->push_back(boost::dynamic_pointer_cast<operator_type>(/*ComposedPrecPtr*/matrCrsPtr), false, true);

    std::vector<LifeV::Real> real;
    std::vector<LifeV::Real> imaginary;

    boost::shared_ptr<EigenSolver> eig;

    UInt nev = M_dataFile("eigensolver/nevec", 10);//number of eigenvectors
    if(nev)
    {
        eig.reset(new EigenSolver(M_PAAP, M_PAAP->OperatorDomainMap(), nev));
        eig->setDataFromGetPot(M_dataFile, "eigensolver/");
        eig->solve();
        eig->eigenvalues(real, imaginary);
    }
    else
    {
        throw UNDEF_EIGENSOLVER_EXCEPTION();
    }
    for (int i=0; i<real.size(); ++i)
    {
        displayer().leaderPrint("\n real part ", real[i]);
        displayer().leaderPrint("\n imaginary part ", imaginary[i]);
    }
    return real[0];
}
#endif



void
Monolithic::computeFNormals( vector_type& normals)
{
    BCNormalManager<mesh_type, matrix_type> normalManager(*M_uFESpace->mesh());
    if ( !M_BChWSS->bdUpdateDone() )//possibly to avoid
        M_BChWSS->bdUpdate(*M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
    normalManager.init((*M_BChWSS)[0], 0.);
    normalManager.computeIntegratedNormals(M_uFESpace->dof(), M_uFESpace->feBd(), normals, *M_uFESpace->mesh());
}


void
Monolithic::initializeMesh(vector_ptrtype fluid_dispOld)
{
    meshMotion().setDisplacement(*fluid_dispOld);
}


void
Monolithic::
iterateMonolithic(const vector_type& rhs, vector_type& step)
{
    Chrono chrono;

    displayer().leaderPrint("    Solving the system ... \n" );

    displayer().leaderPrint("    Updating the boundary conditions ... ");

    //M_monolithicMatrix->GlobalAssemble();
    //necessary if we did not imposed Dirichlet b.c.
    M_linearSolver->setOperator(*M_monolithicMatrix->getMatrix()->getMatrixPtr());

    M_linearSolver->setReusePreconditioner( (M_reusePrec) && (!M_resetPrec) );

    int numIter = M_precPtr->solveSystem( rhs, step, M_linearSolver );

    if (numIter < 0)
        {
            chrono.start();

            M_solid->getDisplayer().leaderPrint("   Iterative solver failed, numiter = ", -numIter );

            if (numIter <= -M_maxIterSolver)
                M_solid->getDisplayer().leaderPrint("   ERROR: Iterative solver failed.\n");
        }

    M_solid->getDisplayer().leaderPrint("   system solved.\n ");
}



void
Monolithic::assembleFluidBlock(UInt iter)
{
    double alpha = 1./M_data->dataFluid()->dataTime()->getTimeStep();//mesh velocity w
    matrix_ptrtype newMatrix(new matrix_type(*M_monolithicMap));
    M_fluid->updateSystem(alpha,*this->M_beta, *this->M_rhs, newMatrix );
    this->M_fluid->updateStab(*newMatrix);
    newMatrix->GlobalAssemble();

    M_fluidBlock.reset(new matrix_type(*M_monolithicMap));
    *M_fluidBlock += *newMatrix;


        if(iter==0)
        {
            M_nbEval = 0; // new time step
            M_resetPrec=true;
            *this->M_rhs               += M_fluid->matrMass()*M_bdf->time_der( M_data->dataFluid()->dataTime()->getTimeStep() );
            couplingRhs(this->M_rhs, M_un);

            if (!M_restarts)
            {
                this->M_solid->updateVel();
                M_restarts = false;
            }
            updateSolidSystem(this->M_rhs);
        }

        *M_rhsFull = *M_rhs;

}

}
