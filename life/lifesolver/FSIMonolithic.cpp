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

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <EpetraExt_Reindex_MultiVector.h>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/LifeV.hpp>
#include <life/lifesolver/FSIMonolithic.hpp>

namespace LifeV
{

// ===================================================
// Constructors, Destructor
// ===================================================
FSIMonolithic::FSIMonolithic():
        super_Type(),
        M_monolithicMap(),
        M_interfaceMap(),
        M_beta(),
        M_monolithicMatrix(),
        M_precPtr(),
        M_rhsFull(),
        M_BCh_flux(),
        M_BCh_Robin(),
        M_BChWS(),
        M_offset(0),
        M_solidAndFluidDim(0),
        M_fluidBlock(),
        M_solidBlock(),
        M_solidBlockPrec(),
        M_linearSolver(),
        M_numerationInterface(),
        M_BChs(),
        M_FESpaces(),
        M_diagonalScale(false),
        M_reusePrec(true),
        M_resetPrec(true),
        M_maxIterSolver(-1),
        M_restarts(false),
        //end of protected attributes
        M_preconditionedSymmetrizedMatrix(),
        M_stress(),
        M_fluxes(0)
#ifdef OBSOLETE
        ,M_rhsShapeDerivatives()
#endif
{}

FSIMonolithic::~FSIMonolithic()
{
}

// ===================================================
// Setup Methods
// ===================================================
void
FSIMonolithic::setupFEspace()
{
    super_Type::setupFEspace();

    // Monolitic: In the beginning I need a non-partitioned mesh. later we will do the partitioning
    M_dFESpace.reset( new FESpace<mesh_Type, MapEpetra>( M_solidMesh,
                                                         M_data->dataSolid()->order(),
                                                         nDimensions,
                                                         M_epetraComm));
}

void
FSIMonolithic::setupDOF( void )
{

    M_dofStructureToHarmonicExtension    .reset( new DOFInterface3Dto3D );
    M_dofStructureToFluid    .reset( new DOFInterface3Dto3D );

    M_dofStructureToHarmonicExtension->setup(   M_mmFESpace->refFE(), M_mmFESpace->dof(),
                                                M_dFESpace->refFE(), M_dFESpace->dof() );
    M_dofStructureToHarmonicExtension->update( *M_mmFESpace->mesh(),  M_data->fluidInterfaceFlag(),
                                               *M_dFESpace->mesh(),  M_data->structureInterfaceFlag(),
                                               M_data->interfaceTolerance(),
                                               M_data->fluidInterfaceVertexFlag() );

    M_dofStructureToFluid->setup(   M_uFESpace->refFE(), M_uFESpace->dof(),
                                    M_dFESpace->refFE(), M_dFESpace->dof() );
    M_dofStructureToFluid->update( *M_uFESpace->mesh(),  M_data->fluidInterfaceFlag(),
                                   *M_dFESpace->mesh(),  M_data->structureInterfaceFlag(),
                                   M_data->interfaceTolerance(),
                                   M_data->fluidInterfaceVertexFlag() );


    createInterfaceMaps(M_dofStructureToHarmonicExtension->localDofMap());
}

void
FSIMonolithic::setupDOF( meshFilter_Type& filterMesh )
{
    createInterfaceMaps(*filterMesh.getStoredInterface(0));
}


void
FSIMonolithic::setupSystem( )
{
    M_fluid->setUp( M_dataFile );
    setUp( M_dataFile );
}

void
FSIMonolithic::setUp( const GetPot& dataFile )
{

    M_linearSolver.reset(new solver_Type(M_epetraComm));

    M_linearSolver->setDataFromGetPot( dataFile, "linear_system/solver" );
    M_linearSolver->setupPreconditioner(dataFile, "linear_system/prec");//to avoid if we have already a prec.
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
    M_diagonalScale    = dataFile( "linear_system/prec/diagonalScaling",  false );
    M_restarts         = dataFile( "exporter/start"  ,  0   );
}

void
FSIMonolithic::setupFluidSolid( )
{
    M_BCh_flux = M_BCh_u; // For the moment M_BCh_u contains only the fluxes.
    M_fluxes=M_BCh_u->size( );

    setupFluidSolid( M_fluxes );

    M_BCh_flux->setOffset(M_offset-M_fluxes);
    std::vector<BCBase>::iterator fluxIt = M_BCh_flux->begin( );
    for ( UInt i = 0; i < M_fluxes; ++i, ++fluxIt )
        fluxIt->setOffset( i );

    M_BCh_Robin = M_BCh_d;// For the moment M_BCh_d contains only the Robin b.c..
    M_BCh_Robin->setOffset(M_offset);
}

void
FSIMonolithic::setupFluidSolid( UInt const fluxes )
{

    // Note: up to now it works only with matching grids (and poly order) on the interface
    assert(M_fluidInterfaceMap->map(Unique)->NumGlobalElements() == M_solidInterfaceMap->map(Unique)->NumGlobalElements());

    M_interfaceMap = M_solidInterfaceMap;

    //std::map<ID, ID> const& localDofMap = M_dofStructureToHarmonicExtension->localDofMap();
    std::map<ID, ID>::const_iterator ITrow;

    M_monolithicMap.reset(new MapEpetra(M_uFESpace->map()));

    std::string opertype = M_dataFile("problem/blockOper", "AdditiveSchwarz");

    createOperator( opertype );

    *M_monolithicMap+= M_pFESpace->map();
    *M_monolithicMap+= fluxes;
    *M_monolithicMap+= M_dFESpace->map();

    M_monolithicMatrix->createInterfaceMap( *M_interfaceMap, M_dofStructureToHarmonicExtension->localDofMap(), M_dFESpace->map().map(Unique)->NumGlobalElements()/nDimensions, M_epetraWorldComm );
    *M_monolithicMap += *M_monolithicMatrix->interfaceMap();

    //the map for the interface coupling matrices should be done with respect to the coarser mesh.
    M_monolithicMatrix->numerationInterface(M_numerationInterface);
    M_beta.reset  (new vector_Type(/*M_monolithicMap*/M_uFESpace->map()));

    M_offset = M_uFESpace->dof().numTotalDof()*nDimensions + fluxes +  M_pFESpace->dof().numTotalDof();
    M_solidAndFluidDim= M_offset + M_dFESpace->dof().numTotalDof()*nDimensions;
    M_BCh_d->setOffset(M_offset);
}

// ===================================================
// Public Methods
// ===================================================
void
FSIMonolithic::monolithicToInterface(vector_Type& lambdaSolid, const vector_Type& disp)
{
    if (disp.mapType() == Repeated)
    {
        vector_Type const  dispUnique(disp, Unique);
        monolithicToInterface(lambdaSolid, dispUnique);
        return;
    }
    if (lambdaSolid.mapType() == Repeated)
    {
        vector_Type  lambdaSolidUn(lambdaSolid.map(), Unique);
        monolithicToInterface( lambdaSolidUn, disp);
        lambdaSolid = lambdaSolidUn;
        return;
    }
    /* UInt MyOffset(M_uFESpace->map().getMap(Unique)->NumMyElements()+M_pFESpace->map().getMap(Unique)->NumMyElements());
       vector_Type subDisp(this->M_dFESpace->map(), Unique);
       subDisp.mySubset(disp, MyOffset);
       lambdaSolid=subDisp;*/

    MapEpetra subMap(*disp.map().map(Unique), M_offset,disp.map().map(Unique)->NumGlobalElements() );
    vector_Type subDisp(subMap, Unique);
    subDisp.subset(disp, M_offset);
    lambdaSolid=subDisp;
}

void
FSIMonolithic::monolithicToX(const vector_Type& disp, vector_Type& dispFluid, MapEpetra& map, UInt offset)
{
    if(disp.mapType()== Repeated)
    {
        vector_Type dispUnique(disp, Unique);
        monolithicToX(dispUnique, dispFluid, map, offset);
        dispFluid = dispUnique;
        return;
    }
    dispFluid.subset(disp, map, offset, offset);
}

void
FSIMonolithic::setDispSolid( const vector_Type& solution )
{
    vector_Type disp(*M_monolithicMap);
    monolithicToX(solution, disp, M_dFESpace->map(), M_offset);
    this->M_solid->setDisp(disp);
}

void
FSIMonolithic::buildSystem()
{
    M_solidBlock.reset(new matrix_Type(*M_monolithicMap, 1));//since it is constant, we keep this throughout the simulation
    M_solid->buildSystem(M_solidBlock, M_data->dataSolid()->dataTime()->timeStep()*M_solid->rescaleFactor());//M_data->dataSolid()->rescaleFactor());
    M_solidBlock->globalAssemble();
    M_solid->rescaleMatrices();
}

#ifdef HAVE_TRILINOS_ANASAZI
Real&
FSIMonolithic::computeMaxSingularValue( )
{
    typedef Epetra_Operator                                                operatorPtr_Type;

    M_preconditionedSymmetrizedMatrix.reset(new ComposedOperator<Epetra_Operator>(M_epetraComm));

    boost::shared_ptr<operatorPtr_Type>  ComposedPrecPtr(M_linearSolver->preconditioner()->preconditioner());

    //M_monolithicMatrix->getMatrixPtr()->OptimizeStorage();
    boost::shared_ptr<Epetra_FECrsMatrix> matrCrsPtr(new Epetra_FECrsMatrix(*M_monolithicMatrix->matrix()->matrixPtr()));
    M_preconditionedSymmetrizedMatrix->push_back(boost::dynamic_pointer_cast<operatorPtr_Type>(/*ComposedPrecPtr*/matrCrsPtr));
    M_preconditionedSymmetrizedMatrix->push_back(boost::dynamic_pointer_cast<operatorPtr_Type>(ComposedPrecPtr/*matrCrsPtr*/), true);
    M_preconditionedSymmetrizedMatrix->push_back(boost::dynamic_pointer_cast<operatorPtr_Type>(ComposedPrecPtr/*matrCrsPtr*/), true, true);
    M_preconditionedSymmetrizedMatrix->push_back(boost::dynamic_pointer_cast<operatorPtr_Type>(/*ComposedPrecPtr*/matrCrsPtr), false, true);

    std::vector<LifeV::Real> real;
    std::vector<LifeV::Real> imaginary;

    boost::shared_ptr<EigenSolver> eig;

    UInt nev = M_dataFile("eigensolver/nevec", 10);//number of eigenvectors
    if(nev)
    {
        eig.reset(new EigenSolver(M_preconditionedSymmetrizedMatrix, M_preconditionedSymmetrizedMatrix->OperatorDomainMap(), nev));
        eig->setDataFromGetPot(M_dataFile, "eigensolver/");
        eig->solve();
        eig->eigenvalues(real, imaginary);
    }
    else
    {
        throw UNDEF_EIGENSOLVER_EXCEPTION();
    }
    for (UInt i=0; i<real.size(); ++i)
    {
        displayer().leaderPrint("\n real part ", real[i]);
        displayer().leaderPrint("\n imaginary part ", imaginary[i]);
    }
    return real[0];
}
#endif

void
FSIMonolithic::computeFluidNormals( vector_Type& normals)
{
    BCManageNormal<matrix_Type> normalManager;
    if ( !M_BChWS->bcUpdateDone() )//possibly to avoid
        M_BChWS->bcUpdate(*M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
    normalManager.init((*M_BChWS)[0], 0.);
    normalManager.computeIntegratedNormals(M_uFESpace->dof(), M_uFESpace->feBd(), normals, *M_uFESpace->mesh());
}

void
FSIMonolithic::solveJac( vector_Type& step, const vector_Type& res, const Real /*linearRelTol*/ )
{
    setupBlockPrec( );

    checkIfChangedFluxBC( M_precPtr );

    M_precPtr->blockAssembling();
    M_precPtr->applyBoundaryConditions( dataFluid()->dataTime()->time() );
    M_precPtr->GlobalAssemble();

#ifdef HAVE_LIFEV_DEBUG
    M_solid->displayer().leaderPrint("  M-  Residual NormInf:                        ", res.normInf(), "\n");
#endif

    iterateMonolithic(res, step);

#ifdef HAVE_LIFEV_DEBUG
    M_solid->displayer().leaderPrint("  M-  Solution NormInf:                        ", step.normInf(), "\n");
#endif
}

void
FSIMonolithic::updateSystem()
{
    vector_Type solution(*this->M_monolithicMap);
    monolithicToX(*this->M_un, solution, M_uFESpace->map(), UInt(0));
    this->M_bdf->shiftRight(solution);

    M_meshMotion->updateSystem();

    this->fluid().updateUn(*this->M_un);
    *M_rhs*=0;
    *M_rhsFull*=0;
    this->M_fluid->resetStabilization();
}

void
FSIMonolithic::initialize( fluidPtr_Type::value_type::function_Type const& u0,
                           fluidPtr_Type::value_type::function_Type const& p0,
                           solidPtr_Type::value_type::Function const& d0,
                           solidPtr_Type::value_type::Function const& /*w0*/,
                           fluidPtr_Type::value_type::function_Type const& /*df0*/ )
{
    vector_Type u(M_uFESpace->map());
    M_uFESpace->interpolate(u0, u, M_data->dataFluid()->dataTime()->time());

    vector_Type p(M_pFESpace->map());
    M_pFESpace->interpolate(p0, p, M_data->dataFluid()->dataTime()->time());

    vector_Type d(M_dFESpace->map());
    M_dFESpace->interpolate(d0, d, M_data->dataSolid()->dataTime()->time());

    initialize(u, p, d);
}

void
FSIMonolithic::initialize( const vectorPtr_Type& fluidVelocityAndPressure,
                           const vectorPtr_Type& fluidDisplacement,
                           const vectorPtr_Type& solidVelocity,
                           const vectorPtr_Type& solidDisplacement )
{
    // Solution
    setSolution( *fluidVelocityAndPressure );

    // Fluid
    M_fluid->initialize( *fluidVelocityAndPressure );
    M_meshMotion->initialize( *fluidDisplacement );
    initializeBDF( *fluidVelocityAndPressure );

    // Solid
    // Extend the external solid vectors to have the monolithic map
    vectorPtr_Type extendedSolidDisplacement( new vector_Type( *M_monolithicMap ) );
    vectorPtr_Type extendedSolidVelocity( new vector_Type( *M_monolithicMap ) );

    extendedSolidDisplacement->subset( *solidDisplacement, solidDisplacement->map(), static_cast <UInt> ( 0 ), M_offset );
    extendedSolidVelocity->subset( *solidVelocity, solidVelocity->map(), static_cast <UInt> ( 0 ), M_offset );

    // Rescale the quantities
    *extendedSolidDisplacement /= M_data->dataFluid()->dataTime()->timeStep() * M_solid->rescaleFactor();
    *extendedSolidVelocity /= M_data->dataFluid()->dataTime()->timeStep() * M_solid->rescaleFactor();

    M_solid->initialize( extendedSolidDisplacement, extendedSolidVelocity );
}

void
FSIMonolithic::initializeMesh(vectorPtr_Type fluid_dispOld)
{
    meshMotion().setDisplacement(*fluid_dispOld);
}

void
FSIMonolithic::initialize( const vector_Type& u0, const vector_Type& p0, const vector_Type& d0)
{
    *M_un=u0;
    M_un->add(p0, nDimensions*M_uFESpace->dof().numTotalDof());
    M_un->add(d0, M_offset);
}



// ===================================================
// Protected Methods
// ===================================================
void
FSIMonolithic::iterateMonolithic(const vector_Type& rhs, vector_Type& step)
{
    LifeChrono chrono;

    displayer().leaderPrint("  M-  Solving the system ... \n" );

    //M_monolithicMatrix->GlobalAssemble();
    //necessary if we did not imposed Dirichlet b.c.
    M_linearSolver->setOperator(*M_monolithicMatrix->matrix()->matrixPtr());

    M_linearSolver->setReusePreconditioner( (M_reusePrec) && (!M_resetPrec) );

    int numIter = M_precPtr->solveSystem( rhs, step, M_linearSolver );

    if (numIter < 0)
    {
        chrono.start();

        M_solid->displayer().leaderPrint("   Iterative solver failed, numiter = ", -numIter );

        if (numIter <= -M_maxIterSolver)
            M_solid->displayer().leaderPrint("   ERROR: Iterative solver failed.\n");
    }

    //M_solid->getDisplayer().leaderPrint("  M-  System solved.\n" );
}

void
FSIMonolithic::couplingRhs(vectorPtr_Type rhs, vectorPtr_Type un) // not working with non-matching grids
{
    std::map<ID, ID> const& localDofMap = M_dofStructureToHarmonicExtension->localDofMap();
    std::map<ID, ID>::const_iterator ITrow;
    //    UInt solidDim=M_dFESpace->map().getMap(Unique)->NumGlobalElements()/nDimensions;

    vector_Type lambda(*M_interfaceMap, Unique);
    this->monolithicToInterface(lambda, *un);
    UInt interface(M_monolithicMatrix->interface());
    //Real rescale(M_solid->rescaleFactor());
    UInt totalDofs(M_dFESpace->dof().numTotalDof());


    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = localDofMap.begin(); ITrow != localDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap->map(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                if(rhs.get())
                    (*rhs)[  (int)(*M_numerationInterface)[ITrow->second ] + dim*interface +M_solidAndFluidDim ] = -lambda( ITrow->second + dim*totalDofs )/**rescale*/;
            }
        }
    }
}

void
FSIMonolithic::
evalResidual( const vector_Type& sol,const vectorPtr_Type& rhs, vector_Type& res, bool diagonalScaling)
{
    if( diagonalScaling )
        diagonalScale(*rhs, M_monolithicMatrix->matrix());

    res = *(M_monolithicMatrix->matrix())*sol;
    res -= *rhs; // Ax-b
}

void
FSIMonolithic::
updateSolidSystem( vectorPtr_Type & rhsFluidCoupling )
{
    M_solid->updateSystem();
    *rhsFluidCoupling += *M_solid->rhsWithoutBC();
}

void
FSIMonolithic::
diagonalScale(vector_Type& rhs, matrixPtr_Type matrFull)
{
    Epetra_Vector diagonal(*rhs.map().map(Unique));
    //M_matrFull->matrixPtr()->InvRowSums(diagonal);
    //M_matrFull->matrixPtr()->InvRowMaxs(diagonal);
    //M_matrFull->matrixPtr()->InvColSums(diagonal);
    matrFull->matrixPtr()->InvColMaxs(diagonal);
    matrFull->matrixPtr()->LeftScale(diagonal);
    rhs.epetraVector().Multiply(1, rhs.epetraVector(), diagonal,0);
}

void
FSIMonolithic::solidInit(std::string const& dOrder)
{   // Monolitic: In the beginning I need a non-partitioned mesh. later we will do the partitioning
    M_dFESpace.reset(new FESpace<mesh_Type, MapEpetra>(M_solidMesh,
                                                       dOrder,
                                                       nDimensions,
                                                       M_epetraComm));
}

void
FSIMonolithic::variablesInit(const std::string& dOrder)
{
    M_dFESpace.reset(new FESpace<mesh_Type, MapEpetra>(*M_solidMeshPart,
                                                       dOrder,
                                                       3,
                                                       M_epetraComm));
    // INITIALIZATION OF THE VARIABLES
    M_lambdaFluid.reset(new vector_Type(*M_fluidInterfaceMap, Unique) );
    M_lambdaFluidRepeated.reset(new vector_Type(*M_fluidInterfaceMap, Repeated) );
}

void FSIMonolithic::setupBlockPrec( )
{
    if(!(M_precPtr->set()))
     {
         M_precPtr->push_back_matrix(M_solidBlockPrec, false);
         M_precPtr->push_back_matrix(M_fluidBlock, true);
         M_precPtr->setConditions(M_BChs);
         M_precPtr->setSpaces(M_FESpaces);
         M_precPtr->setOffsets(2, M_offset, 0);
         M_precPtr->coupler(M_monolithicMap, M_dofStructureToHarmonicExtension->localDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->timeStep());
     }
    else
    {
        M_precPtr->replace_matrix(M_fluidBlock, 1);
        M_precPtr->replace_matrix(M_solidBlockPrec, 0);
    }
}

void
FSIMonolithic::assembleSolidBlock( UInt iter, vectorPtr_Type& solution )
{
    if (iter == 0)
    {
        updateSolidSystem(this->M_rhs);
    }
    else
    {
        M_solid->computeMatrix( *solution, 1.);
    }

    M_solid->getSolidMatrix( M_solidBlock );
    M_solidBlockPrec.reset( new matrix_Type( *M_monolithicMap, 1 ) );
    *M_solidBlockPrec += *M_solidBlock;
}

void
FSIMonolithic::assembleFluidBlock(UInt iter, vectorPtr_Type& solution)
{
    M_fluidBlock.reset(new matrix_Type(*M_monolithicMap));

    Real alpha = 1./M_data->dataFluid()->dataTime()->timeStep();//mesh velocity w
    M_fluid->updateSystem(alpha,*this->M_beta, *this->M_rhs, M_fluidBlock, solution );
    this->M_fluid->updateStabilization(*M_fluidBlock);

    if (iter==0)
    {
        M_resetPrec=true;
        M_bdf->updateRHSContribution( M_data->dataFluid()->dataTime()->timeStep() );
        *M_rhs += M_fluid->matrixMass()*(*M_un)/M_data->dataFluid()->dataTime()->timeStep();//(M_bdf->rhsContributionFirstDerivative()) ;
        couplingRhs(M_rhs, M_un);
    }
}

void FSIMonolithic::updateRHS()
{
    M_bdf->updateRHSContribution( M_data->dataFluid()->dataTime()->timeStep() );
    *M_rhs += M_fluid->matrixMass()*(*M_un)*1/M_data->dataFluid()->dataTime()->timeStep();//M_bdf->rhsContributionFirstDerivative() ;
    couplingRhs(M_rhs, M_un);
    //M_solid->updateVel();
    updateSolidSystem(M_rhs);
    *M_rhsFull = *M_rhs;
}

namespace
{
static Real fZero(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{return 0.;}
static Real fOne(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{return 1.;}
}

void FSIMonolithic::enableStressComputation(UInt flag)
{
    M_BChWS.reset(new solidBchandler_Type());
    BCFunctionBase bcfZero(fZero);
    BCFunctionBase bcfOne(fOne);
    M_bcfWs.setFunctions_Robin(bcfOne,bcfOne);

    M_BChWS->addBC("WS", (UInt) flag, Robin, Full, M_bcfWs, 3);
}

FSIMonolithic::vectorPtr_Type FSIMonolithic::computeStress()
{
    vector_Type lambda(M_monolithicMatrix->interfaceMap());
    lambda.subset(*M_un, M_solidAndFluidDim);

    M_boundaryMass.reset(new matrix_Type(*M_interfaceMap));
    if ( !M_BChWS->bcUpdateDone() )
        M_BChWS->bcUpdate(*M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
    bcManageMatrix(*M_boundaryMass, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BChWS, M_dFESpace->feBd(), 1., dataSolid()->getdataTime()->time() );
    M_boundaryMass->globalAssemble();

    solver_Type solverMass(M_epetraComm);
    solverMass.setDataFromGetPot( M_dataFile, "solid/solver" );
    solverMass.setupPreconditioner(M_dataFile, "solid/prec");//to avoid if we have already a prec.

    boost::shared_ptr<PreconditionerIfpack> P(new PreconditionerIfpack());

    vectorPtr_Type sol(new vector_Type(M_monolithicMatrix->interfaceMap()));
    solverMass.setMatrix(*M_boundaryMass);
    solverMass.setReusePreconditioner(false);
    solverMass.solveSystem( lambda, *sol, M_boundaryMass);

    EpetraExt::MultiVector_Reindex reindexMV(*M_interfaceMap->map(Unique));
    boost::shared_ptr<MapEpetra> newMap(new MapEpetra( *M_interfaceMap ));
    M_stress.reset(new vector_Type(reindexMV(lambda.epetraVector()), newMap, Unique));
    return M_stress;
}

void
FSIMonolithic::checkIfChangedFluxBC( precPtr_Type oper )
{
   UInt nfluxes(M_BChs[1]->numberOfBCWithType(Flux));
    if(M_fluxes != nfluxes)
    {
        //std::vector<bcName_Type> names = M_BChs[1]->findAllBCWithType(Flux);
        for (UInt i=0; i<M_fluxes; ++i)
        {
            const BCBase* bc (M_BChs[1]->findBCWithName(M_BCFluxNames[i]));
            if(bc->type() != Flux)
            {
                oper->addToCoupling(1., M_fluxOffset[i], M_fluxOffset[i],1 );
            }
        }
    }
}


}
