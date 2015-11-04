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


#include <EpetraExt_Reindex_MultiVector.h>


#include <lifev/core/LifeV.hpp>
#include <lifev/fsi/solver/FSIMonolithic.hpp>

namespace LifeV
{

// ===================================================
// Constructors, Destructor
// ===================================================
FSIMonolithic::FSIMonolithic() :
    super_Type(),
    M_monolithicMap(),
    M_interfaceMap(),
    M_beta(),
    M_monolithicMatrix(),
    M_precPtr(),
    M_rhsFull(),
    M_BCh_flux(),
    M_BChWS(),
    M_offset (0),
    M_solidAndFluidDim (0),
    M_fluidBlock(),
    M_solidBlockPrec(),
    M_linearSolver(),
    M_numerationInterface(),
    M_BChs(),
    M_FESpaces(),
    M_diagonalScale (false),
    M_reusePrec (true),
    M_resetPrec (true),
    M_maxIterSolver (-1),
    M_restarts (false),
    //end of protected attributes
    M_preconditionedSymmetrizedMatrix(),
    M_stress(),
    M_fluxes (0)
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
    M_dFESpace.reset ( new FESpace<mesh_Type, MapEpetra> ( M_solidMesh,
                                                           M_data->dataSolid()->order(),
                                                           nDimensions,
                                                           M_epetraComm) );

}

void
FSIMonolithic::setupDOF ( void )
{
    M_dofStructureToHarmonicExtension    .reset ( new DOFInterface3Dto3D );
    M_dofStructureToFluid    .reset ( new DOFInterface3Dto3D );

    M_dofStructureToHarmonicExtension->setup (   M_mmFESpace->refFE(), M_mmFESpace->dof(),
                                                 M_dFESpace->refFE(), M_dFESpace->dof() );
    M_dofStructureToHarmonicExtension->update ( *M_mmFESpace->mesh(),  M_data->fluidInterfaceFlag(),
                                                *M_dFESpace->mesh(),  M_data->structureInterfaceFlag(),
                                                M_data->interfaceTolerance(),
                                                M_data->fluidInterfaceVertexFlag() );

    M_dofStructureToFluid->setup (   M_uFESpace->refFE(), M_uFESpace->dof(),
                                     M_dFESpace->refFE(), M_dFESpace->dof() );
    M_dofStructureToFluid->update ( *M_uFESpace->mesh(),  M_data->fluidInterfaceFlag(),
                                    *M_dFESpace->mesh(),  M_data->structureInterfaceFlag(),
                                    M_data->interfaceTolerance(),
                                    M_data->fluidInterfaceVertexFlag() );

    createInterfaceMaps (M_dofStructureToFluid/*HarmonicExtension*/->localDofMap() );

    M_fluidMesh.reset();
    M_solidMesh.reset();
}

#ifdef HAVE_HDF5
void
FSIMonolithic::setupDOF ( meshFilter_Type& filterMesh )
{
    createInterfaceMaps (*filterMesh.getStoredInterface (0) );
}
#endif

void
FSIMonolithic::setupSystem( )
{
    M_fluid->setUp ( M_dataFile );
    setup ( M_dataFile );
}

void
FSIMonolithic::setup ( const GetPot& dataFile )
{

    M_linearSolver.reset (new solver_Type (M_epetraComm) );

    M_linearSolver->setDataFromGetPot ( dataFile, "linear_system/solver" );
    std::string prectype = dataFile ("problem/DDBlockPrec", "PRECTYPE_UNDEFINED");
    std::string opertype = dataFile ("problem/blockOper", "PRECTYPE_UNDEFINED");

    M_precPtr.reset (BlockPrecFactory::instance().createObject ( prectype ) );

    M_precPtr->setupSolver (*M_linearSolver, dataFile);
    std::string section ("linear_system/prec");
    M_precPtr->setComm (M_epetraComm);
    M_precPtr->setDataFromGetPot (dataFile, section); //to avoid if we build the prec from a matrix.
    M_monolithicMatrix->setComm (M_epetraComm);
    M_monolithicMatrix->setDataFromGetPot (dataFile, section); //to avoid if we build the prec from a matrix.
    M_reusePrec     = dataFile ( "linear_system/prec/reuse", true);
    M_maxIterSolver = dataFile ( "linear_system/solver/max_iter", -1);
    M_diagonalScale    = dataFile ( "linear_system/prec/diagonalScaling",  false );
    M_restarts         = dataFile ( "exporter/start"  ,  0   );
}

void
FSIMonolithic::setupFluidSolid( )
{
    M_BCh_flux = M_BCh_u; // For the moment M_BCh_u contains only the fluxes.
    M_fluxes = M_BCh_u->size( );

    setupFluidSolid ( M_fluxes );

    M_BCh_flux->setOffset (M_offset - M_fluxes);
    std::vector<BCBase>::iterator fluxIt = M_BCh_flux->begin( );
    for ( UInt i = 0; i < M_fluxes; ++i, ++fluxIt )
    {
        fluxIt->setOffset ( i );
    }

}

void
FSIMonolithic::setupFluidSolid ( UInt const fluxes )
{

    // Note: up to now it works only with matching grids (and poly order) on the interface
    assert (M_fluidInterfaceMap->map (Unique)->NumGlobalElements() == M_solidInterfaceMap->map (Unique)->NumGlobalElements() );

    M_interfaceMap = M_solidInterfaceMap;

    //std::map<ID, ID> const& localDofMap = M_dofStructureToHarmonicExtension->localDofMap();
    std::map<ID, ID>::const_iterator ITrow;

    M_monolithicMap.reset (new MapEpetra (M_uFESpace->map() ) );

    std::string opertype = M_dataFile ("problem/blockOper", "AdditiveSchwarz");

    createOperator ( opertype );

    *M_monolithicMap += M_pFESpace->map();
    *M_monolithicMap += fluxes;
    *M_monolithicMap += M_dFESpace->map();

    M_monolithicMatrix->createInterfaceMap ( *M_interfaceMap, M_dofStructureToFluid->localDofMap(), M_dFESpace->map().map (Unique)->NumGlobalElements() / nDimensions, M_epetraWorldComm );
    *M_monolithicMap += *M_monolithicMatrix->interfaceMap();

    //the map for the interface coupling matrices should be done with respect to the coarser mesh.
    M_monolithicMatrix->numerationInterface (M_numerationInterface);
    M_beta.reset  (new vector_Type (M_uFESpace->map() ) );

    M_offset = M_uFESpace->dof().numTotalDof() * nDimensions + fluxes +  M_pFESpace->dof().numTotalDof();
    M_solidAndFluidDim = M_offset + M_dFESpace->dof().numTotalDof() * nDimensions;
    M_BCh_d->setOffset (M_offset);

    M_fluxes = fluxes;
}

// ===================================================
// Public Methods
// ===================================================
void
FSIMonolithic::monolithicToInterface (vector_Type& lambdaSolid, const vector_Type& disp)
{
    if (disp.mapType() == Repeated)
    {
        vector_Type const  dispUnique (disp, Unique);
        monolithicToInterface (lambdaSolid, dispUnique);
        return;
    }
    if (lambdaSolid.mapType() == Repeated)
    {
        vector_Type  lambdaSolidUn (lambdaSolid.map(), Unique);
        monolithicToInterface ( lambdaSolidUn, disp);
        lambdaSolid = lambdaSolidUn;
        return;
    }

    MapEpetra subMap (*disp.map().map (Unique), M_offset, disp.map().map (Unique)->NumGlobalElements() );
    vector_Type subDisp (subMap, Unique);
    subDisp.subset (disp, M_offset);
    lambdaSolid = subDisp;
}

void
FSIMonolithic::monolithicToX (const vector_Type& disp, vector_Type& dispFluid, MapEpetra& map, UInt offset)
{
    if (disp.mapType() == Repeated)
    {
        vector_Type dispUnique (disp, Unique);
        monolithicToX (dispUnique, dispFluid, map, offset);
        dispFluid = dispUnique;
        return;
    }
    dispFluid.subset (disp, map, offset, offset);
}


void
FSIMonolithic::buildSystem()
{
    M_solid->buildSystem ( M_solidTimeAdvance->coefficientSecondDerivative ( 0 ) / (M_data->dataSolid()->dataTime()->timeStep() *M_data->dataSolid()->dataTime()->timeStep() ) );
}

#ifdef HAVE_TRILINOS_ANASAZI
Real&
FSIMonolithic::computeMaxSingularValue( )
{
    typedef Epetra_Operator                                                operatorPtr_Type;

    M_preconditionedSymmetrizedMatrix.reset (new ComposedOperator<Epetra_Operator> (M_epetraComm) );

    std::shared_ptr<operatorPtr_Type>  ComposedPrecPtr (M_linearSolver->preconditioner()->preconditioner() );

    //M_monolithicMatrix->getMatrixPtr()->OptimizeStorage();
    std::shared_ptr<Epetra_FECrsMatrix> matrCrsPtr (new Epetra_FECrsMatrix (*M_monolithicMatrix->matrix()->matrixPtr() ) );
    M_preconditionedSymmetrizedMatrix->push_back (std::static_pointer_cast<operatorPtr_Type> (/*ComposedPrecPtr*/matrCrsPtr) );
    M_preconditionedSymmetrizedMatrix->push_back ( (ComposedPrecPtr/*matrCrsPtr*/), true);
    M_preconditionedSymmetrizedMatrix->push_back ( (ComposedPrecPtr/*matrCrsPtr*/), true, true);
    M_preconditionedSymmetrizedMatrix->push_back (std::static_pointer_cast<operatorPtr_Type> (/*ComposedPrecPtr*/matrCrsPtr), false, true);

    std::vector<LifeV::Real> real;
    std::vector<LifeV::Real> imaginary;

    std::shared_ptr<EigenSolver> eig;

    UInt nev = M_dataFile ("eigensolver/nevec", 10); //number of eigenvectors
    if (nev)
    {
        eig.reset (new EigenSolver (M_preconditionedSymmetrizedMatrix, M_preconditionedSymmetrizedMatrix->OperatorDomainMap(), nev) );
        eig->setDataFromGetPot (M_dataFile, "eigensolver/");
        eig->solve();
        eig->eigenvalues (real, imaginary);
    }
    else
    {
        throw UNDEF_EIGENSOLVER_EXCEPTION();
    }
    for (UInt i = 0; i < real.size(); ++i)
    {
        displayer().leaderPrint ("\n real part ", real[i]);
        displayer().leaderPrint ("\n imaginary part ", imaginary[i]);
    }
    return real[0];
}
#endif

void
FSIMonolithic::computeFluidNormals ( vector_Type& normals)
{
    BCManageNormal<matrix_Type> normalManager;
    if ( !M_BChWS->bcUpdateDone() ) //possibly to avoid
    {
        M_BChWS->bcUpdate (*M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
    }
    normalManager.init ( (*M_BChWS) [0], 0.);
    normalManager.computeIntegratedNormals (M_uFESpace->dof(), M_uFESpace->feBd(), normals, *M_uFESpace->mesh() );
}

void
FSIMonolithic::solveJac ( vector_Type& step, const vector_Type& res, const Real /*linearRelTol*/ )
{
    setupBlockPrec( );

    checkIfChangedFluxBC ( M_precPtr );

    M_precPtr->blockAssembling();
    M_precPtr->applyBoundaryConditions ( dataFluid()->dataTime()->time() );
    M_precPtr->GlobalAssemble();

#ifdef HAVE_LIFEV_DEBUG
    M_solid->displayer().leaderPrint ("  M-  Residual NormInf:                        ", res.normInf(), "\n");
#endif
    iterateMonolithic (res, step);
#ifdef HAVE_LIFEV_DEBUG
    M_solid->displayer().leaderPrint ("  M-  Solution NormInf:                        ", step.normInf(), "\n");
#endif
}

void
FSIMonolithic::updateSystem()
{
    *M_rhs *= 0;
    *M_rhsFull *= 0;
    this->M_fluid->resetStabilization();
}


// ===================================================
// Protected Methods
// ===================================================
void
FSIMonolithic::iterateMonolithic (const vector_Type& rhs, vector_Type& step)
{
    LifeChrono chrono;

    displayer().leaderPrint ("  M-  Solving the system ... \n" );

    //M_monolithicMatrix->GlobalAssemble();
    //necessary if we did not imposed Dirichlet b.c.

    M_linearSolver->setOperator (*M_monolithicMatrix->matrix()->matrixPtr() );

    M_linearSolver->setReusePreconditioner ( (M_reusePrec) && (!M_resetPrec) );

    int numIter = M_precPtr->solveSystem ( rhs, step, M_linearSolver );

    if (numIter < 0)
    {
        chrono.start();

        M_solid->displayer().leaderPrint ("   Iterative solver failed, numiter = ", -numIter );

        if (numIter <= -M_maxIterSolver)
        {
            M_solid->displayer().leaderPrint ("   ERROR: Iterative solver failed.\n");
        }
    }

    //M_solid->getDisplayer().leaderPrint("  M-  System solved.\n" );
}

void
FSIMonolithic::couplingRhs (vectorPtr_Type rhs) // not working with non-matching grids
{
    std::map<ID, ID> const& localDofMap = M_dofStructureToFluid->localDofMap();
    std::map<ID, ID>::const_iterator ITrow;

    vector_Type rhsStructureVelocity (M_solidTimeAdvance->rhsContributionFirstDerivative() *M_solid->rescaleFactor(), Unique, Add);
    vector_Type lambda (*M_interfaceMap, Unique);

    this->monolithicToInterface (lambda, rhsStructureVelocity);

    UInt interface (M_monolithicMatrix->interface() );
    UInt totalDofs (M_dFESpace->dof().numTotalDof() );


    for (UInt dim = 0; dim < nDimensions; ++dim)
    {
        for ( ITrow = localDofMap.begin(); ITrow != localDofMap.end() ; ++ITrow)
        {
            if (M_interfaceMap->map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->second /*+ dim*solidDim*/) ) >= 0 ) //to avoid repeated stuff
            {
                if (rhs.get() )
                {
                    (*rhs) [  (int) (*M_numerationInterface) [ITrow->second ] + dim * interface + M_solidAndFluidDim ] = -lambda ( ITrow->second + dim * totalDofs );
                }
            }
        }
    }
}

void
FSIMonolithic::
evalResidual ( const vector_Type& sol, const vectorPtr_Type& rhs, vector_Type& res, bool diagonalScaling)
{
    if ( diagonalScaling )
    {
        diagonalScale (*rhs, M_monolithicMatrix->matrix() );
    }

    // case of nonlinear laws
    if ( ! ( M_data->dataSolid()->lawType().compare ("nonlinear") ) )
    {
        // need to set correctly the vectors to remove offset
        // Extract the right proportion
        vectorPtr_Type solidPart ( new vector_Type ( M_dFESpace->map() ) );
        solidPart->subset (sol, M_offset );

        vectorPtr_Type resSolidPart ( new vector_Type ( M_dFESpace->map() ) );
        resSolidPart->subset (res, M_offset );

        // Multiplying it by the rescale factor
        *solidPart *= M_solid->rescaleFactor();

        // Computing residual
        M_solid->apply (*solidPart, *resSolidPart);

        resSolidPart->globalAssemble();

        // reassembling them in the right places
        // Only the residual is needed since the sol is not touched inside the solid part
        // sol.subset( *solidPart, solidPart->map(), UInt(0), M_offset);
        res.subset ( *resSolidPart, resSolidPart->map(), UInt (0), M_offset);

        res.globalAssemble();

        M_fluidBlock->globalAssemble();

        res += ( (*M_fluidBlock) * sol);

        M_monolithicMatrix->coupling()->globalAssemble();

        res += *M_monolithicMatrix->coupling() * sol;
    }
    else
    {
        // this works for the linear elastic case where the matrix is not touched
        res = * (M_monolithicMatrix->matrix() ) * sol;

        res -= *rhs; // Ax-b
    }
}

void
FSIMonolithic::
updateSolidSystem ( vectorPtr_Type& rhsFluidCoupling )
{
    Real coefficient ( M_data->dataSolid()->dataTime()->timeStep() * M_data->dataSolid()->dataTime()->timeStep() * M_solid->rescaleFactor() /  M_solidTimeAdvance->coefficientSecondDerivative ( 0 ) );

    M_solidTimeAdvance->updateRHSContribution ( M_data->dataSolid()->dataTime()->timeStep() );

    // Extract the right portion from th right hand side contribution
    vectorPtr_Type solidPortionRHSSecondDerivative ( new vector_Type ( M_dFESpace->map() ) );
    solidPortionRHSSecondDerivative->subset (M_solidTimeAdvance->rhsContributionSecondDerivative(), M_offset );

    // Create a matrix of the size of the structure
    std::shared_ptr<MatrixEpetra<Real> >  solidMassMatrix (new MatrixEpetra<Real> ( M_solid->map(), 1 ) );
    *solidMassMatrix *= 0.0;
    *solidMassMatrix += *M_solid->massMatrix();
    solidMassMatrix->globalAssemble();

    vectorPtr_Type solidPortionRHSFluidCoupling ( new vector_Type ( M_dFESpace->map() ) );
    *solidPortionRHSFluidCoupling *= 0.0;

    // Computing the vector
    *solidPortionRHSFluidCoupling = *solidMassMatrix * ( *solidPortionRHSSecondDerivative * coefficient );

    // Putting it in the right place
    rhsFluidCoupling->subset ( *solidPortionRHSFluidCoupling, solidPortionRHSFluidCoupling->map(), UInt (0), M_offset);
}


void FSIMonolithic::setVectorInStencils ( const vectorPtr_Type& vel,
                                          const vectorPtr_Type& pressure,
                                          const vectorPtr_Type& solidDisp,
                                          //const vectorPtr_Type& fluidDisp,
                                          const UInt iter)
{
    setFluidVectorInStencil (vel, pressure, iter);
    setSolidVectorInStencil (solidDisp, iter);
    //  setALEVectorInStencil(fluidDisp, iter);

}

void FSIMonolithic::setFluidVectorInStencil ( const vectorPtr_Type& vel,
                                              const vectorPtr_Type& pressure,
                                              const UInt /*iter*/ )
{

    //The fluid and solid TimeAdvance classes have a stencil of dimension
    //as big as the coupled problem.

    //Fluid Problem
    vectorPtr_Type vectorMonolithicFluidVelocity (new vector_Type (*M_monolithicMap, Unique, Zero) );
    vectorPtr_Type vectorMonolithicFluidPressure (new vector_Type (*M_monolithicMap, Unique, Zero) );

    *vectorMonolithicFluidVelocity *= 0.0;
    *vectorMonolithicFluidPressure *= 0.0;

    vectorMonolithicFluidVelocity->subset (*vel, vel->map(), UInt (0), UInt (0) ) ;
    vectorMonolithicFluidPressure->subset ( *pressure, pressure->map(), UInt (0), (UInt) 3 * M_uFESpace->dof().numTotalDof() );

    *vectorMonolithicFluidVelocity += *vectorMonolithicFluidPressure;

    vector_Type* normalPointerToFluidVector ( new vector_Type (*vectorMonolithicFluidVelocity) );
    (M_fluidTimeAdvance->stencil() ).push_back ( normalPointerToFluidVector );
}


void FSIMonolithic::setSolidVectorInStencil ( const vectorPtr_Type& solidDisp,
                                              const UInt iter)
{
    //Solid problem
    vectorPtr_Type vectorMonolithicSolidDisplacement (new vector_Type (*M_monolithicMap, Unique, Zero) );
    *vectorMonolithicSolidDisplacement *= 0.0;
    vectorMonolithicSolidDisplacement->subset ( *solidDisp, solidDisp->map(), (UInt) 0, M_offset);
    *vectorMonolithicSolidDisplacement *= 1.0 / M_solid->rescaleFactor();

    //The fluid timeAdvance has size = orderBDF because it is seen as an equation frist order in time.
    //We need to add the solidVector to the fluidVector in the fluid TimeAdvance because we have the
    //extrapolation on it.
    if ( iter <= M_fluidTimeAdvance->size() - 1 )
    {
        * ( M_fluidTimeAdvance->stencil() [ iter ] ) += *vectorMonolithicSolidDisplacement;
    }

    vector_Type* normalPointerToSolidVector ( new vector_Type (*vectorMonolithicSolidDisplacement) );
    (M_solidTimeAdvance->stencil() ).push_back ( normalPointerToSolidVector );

}

void FSIMonolithic::finalizeRestart( )
{
    //Set the initialRHS for the TimeAdvance classes
    vector_Type zeroFluidSolid (*M_monolithicMap, Unique, Zero);
    vector_Type zeroALE (M_mmFESpace->map(), Unique, Zero);

    zeroFluidSolid *= 0.0;
    zeroALE *= 0.0;

    M_fluidTimeAdvance->setInitialRHS (zeroFluidSolid);
    M_solidTimeAdvance->setInitialRHS (zeroFluidSolid);
    M_ALETimeAdvance->setInitialRHS (zeroALE);

    //This updates at the current value (the one when the simulation was stopped) the RHScontribution
    //of the first derivative which is use to compute the velocity in TimeAdvance::velocity().
    //Please note that, even if it is ugly, at this stage, the fluidTimeAdvance is leading the Time Discretization
    //and this is why there  is the dataFluid class to get the dt.
    M_ALETimeAdvance->updateRHSFirstDerivative ( M_data->dataFluid()->dataTime()->timeStep() );
}

void FSIMonolithic::initializeMonolithicOperator ( std::vector< vectorPtr_Type> u0,
                                                   std::vector< vectorPtr_Type> ds0,
                                                   std::vector< vectorPtr_Type> df0)
{
    UInt i;
    if (!u0.size() || !ds0.size() || !df0.size() )
    {
        if ( this->isFluid() )
        {
            for (i = 0; i < M_fluidTimeAdvance->size(); ++i)
            {
                vectorPtr_Type vec (new vector_Type ( *M_monolithicMap ) );
                u0.push_back (vec); // couplingVariableMap()
            }
            for (i = 0; i < M_ALETimeAdvance->size(); ++i)
            {
                vectorPtr_Type vec (new vector_Type ( M_mmFESpace->map() ) );
                df0.push_back (vec); // couplingVariableMap()
            }
        }
        if ( this->isSolid() )
        {
            for (i = 0; i < M_solidTimeAdvance->size(); ++i)
            {
                vectorPtr_Type vec (new vector_Type ( *M_monolithicMap ) );
                ds0.push_back (vec); // couplingVariableMap()
            }
        }
        initializeTimeAdvance (u0, ds0, df0);
        //  M_oper->initializeBDF(*u0);
    }
    else
    {
        initializeTimeAdvance (u0, ds0, df0); // couplingVariableMap()//copy
    }

}

void
FSIMonolithic::
diagonalScale (vector_Type& rhs, matrixPtr_Type matrFull)
{
    Epetra_Vector diagonal (*rhs.map().map (Unique) );
    //M_matrFull->matrixPtr()->InvRowSums(diagonal);
    //M_matrFull->matrixPtr()->InvRowMaxs(diagonal);
    //M_matrFull->matrixPtr()->InvColSums(diagonal);
    matrFull->matrixPtr()->InvColMaxs (diagonal);
    matrFull->matrixPtr()->LeftScale (diagonal);
    rhs.epetraVector().Multiply (1, rhs.epetraVector(), diagonal, 0);
}

void
FSIMonolithic::variablesInit (const std::string& dOrder)
{
    M_dFESpace.reset (new FESpace<mesh_Type, MapEpetra> (M_solidLocalMesh,
                                                         dOrder,
                                                         3,
                                                         M_epetraComm) );

    M_dETFESpace.reset (new ETFESpace<mesh_Type, MapEpetra, 3, 3> (M_solidLocalMesh,
                                                                   & (M_dFESpace->refFE() ),
                                                                   & (M_dFESpace->fe().geoMap() ),
                                                                   M_epetraComm) );


    // INITIALIZATION OF THE VARIABLES
    M_lambdaFluid.reset (new vector_Type (*M_fluidInterfaceMap, Unique) );
    M_lambdaFluidRepeated.reset (new vector_Type (*M_fluidInterfaceMap, Repeated) );
}

void FSIMonolithic::setupBlockPrec( )
{
    if (! (M_precPtr->set() ) )
    {
        M_precPtr->push_back_matrix (M_solidBlockPrec, M_structureNonLinear);
        M_precPtr->push_back_matrix (M_fluidBlock, true);
        M_precPtr->setConditions (M_BChs);
        M_precPtr->setSpaces (M_FESpaces);
        M_precPtr->setOffsets (2, M_offset, 0);
        M_precPtr->coupler (M_monolithicMap, M_dofStructureToFluid->localDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->timeStep(), M_solidTimeAdvance->coefficientFirstDerivative ( 0 ), M_solid->rescaleFactor() );
    }
    else
    {
        M_precPtr->replace_matrix (M_fluidBlock, 1);
        M_precPtr->replace_matrix (M_solidBlockPrec, 0);
    }

#ifdef HAVE_NS_PREC
    /* This shall be enabled once the branch about NS_PREC is going to be merged to master*/
    std::shared_ptr<MonolithicBlockComposed> blockPrecPointer ( std::dynamic_pointer_cast<MonolithicBlockComposed> M_precPtr);

    if (blockPrecPointer.get() != 0)
    {
        UInt fluidPosition = blockPrecPointer->whereIsBlock (MonolithicBlockComposed::fluid);
        ASSERT (blockPrecPointer->blockPrecs().size() < fluidPosition, "The preconditioner corresponding to the fluid block is probably not PCD. Check in the data file");
        std::shared_ptr<PreconditionerPCD> prec_PCD ( std::dynamic_pointer_cast<PreconditionerPCD> blockPrecPointer->blockPrecs() [fluidPosition] );


        if (prec_PCD.get() != 0)
        {
            prec_PCD->setFESpace (M_uFESpace, M_pFESpace);
            prec_PCD->setBCHandler (M_BCh_u);
            prec_PCD->setTimestep (M_data->dataFluid()->dataTime()->timeStep() );
            prec_PCD->setViscosity (M_data->dataFluid()->viscosity() );
            prec_PCD->setDensity (M_data->dataFluid()->density() );
            prec_PCD->setCouplingMatrixView (M_precPtr->couplingVector() [MonolithicBlockComposed::fluid]);
            prec_PCD->setMapStructure (&M_dFESpace->map() );
            prec_PCD->updateBeta (M_fluidTimeAdvance->singleElement (0) );
        }
    }
#endif
}

void
FSIMonolithic::assembleSolidBlock ( UInt iter, const vector_Type& solution )
{
    if (iter == 0)
    {
        updateSolidSystem (this->M_rhs);
    }

    // Resetting the solidBlockPrec term
    M_solidBlockPrec.reset (new matrix_Type (*M_monolithicMap, 1) );
    *M_solidBlockPrec *= 0.0;

    // When ET for structures is used, there is not offset parameter. This is why
    // we need to extract portions of vector and mount them in the proper parts.
    // Extractig the right portion from the total solution of the solid part
    vectorPtr_Type solidPortion ( new vector_Type ( M_dFESpace->map() ) );
    solidPortion->subset (solution, M_offset );

    //Multiplying it by the rescale factor
    *solidPortion *= M_solid->rescaleFactor();

    // Start building a block matrix to handle easily the assembled part
    // in the structure module.
    //1. Create the correct global map for fluxes
    MapEpetra mapEpetraFluxes ( M_fluxes, M_epetraComm );

    //2. Construct a global MapVector with all the maps
    // The globalMatrix will be adjusted according to the monolithic solver
    // that is used.
    matrixBlockPtr_Type globalMatrixWithOnlyStructure;

    //The map is composed of ( u , p , fluxes, solid, interface, d_f )
    if ( !M_data->method().compare ("monolithicGI") )
    {
        globalMatrixWithOnlyStructure.reset (new matrixBlock_Type ( M_uFESpace->map() | M_pFESpace->map() | mapEpetraFluxes | M_dFESpace->map()
                                                                    | * (M_monolithicMatrix->interfaceMap() ) |  M_mmFESpace->map() ) );
    }
    //The map is composed of ( u , p , fluxes, solid, interface )
    else
    {
        globalMatrixWithOnlyStructure.reset (new matrixBlock_Type ( M_uFESpace->map() | M_pFESpace->map() | mapEpetraFluxes | M_dFESpace->map()
                                                                    | * (M_monolithicMatrix->interfaceMap() ) ) );
    }

    // In the case of LE it does not do anything.
    M_solid->material()->updateJacobianMatrix ( *solidPortion, dataSolid(), M_solid->mapMarkersVolumes(), M_solid->mapMarkersIndexes(), M_solid->displayerPtr() ); // computing the derivatives if nonlinear (comment this for inexact Newton);

    //Need to inglobe it into a boost::shared to avoid memeory leak
    std::shared_ptr<MatrixEpetra<Real> >  solidMatrix (new MatrixEpetra<Real> ( M_solid->map(), 1 ) );
    *solidMatrix *= 0.0;

    *solidMatrix += *M_solid->massMatrix();
    *solidMatrix += * (M_solid->material()->jacobian() );

    solidMatrix->globalAssemble();

    MatrixEpetra<Real>* rawPointerToMatrix = new MatrixEpetra<Real> ( *solidMatrix );

    matrixBlockView_Type structurePart (* (globalMatrixWithOnlyStructure->block (3, 3) ) );

    structurePart.setup (UInt (0), UInt (0), UInt (3 * M_dFESpace->dof().numTotalDof() ), UInt (3 * M_dFESpace->dof().numTotalDof() ), rawPointerToMatrix);

    using namespace MatrixEpetraStructuredUtility;

    //3. Put the matrix assembled in the solid in the proper vector
    copyBlock ( structurePart, * (globalMatrixWithOnlyStructure->block (3, 3) ) );

    globalMatrixWithOnlyStructure->globalAssemble();

    //Summing the local matrix into the global
    *M_solidBlockPrec += *globalMatrixWithOnlyStructure;

    M_solidBlockPrec->globalAssemble();

    *M_solidBlockPrec *= M_solid->rescaleFactor();

    delete rawPointerToMatrix;
}

void
FSIMonolithic::assembleFluidBlock (UInt iter, const vector_Type& solution)
{
    M_fluidBlock.reset (new  FSIOperator::fluid_Type::matrix_Type (*M_monolithicMap) );

    Real alpha = M_fluidTimeAdvance->coefficientFirstDerivative ( 0 ) / M_data->dataFluid()->dataTime()->timeStep(); //mesh velocity w

    //This line is based on the hypothesis that the conservativeFormulation flag is set on FALSE
    M_fluid->updateSystem (alpha, *this->M_beta, *this->M_rhs, M_fluidBlock, solution );

    //This in the case of conservativeFormulation == true
    // else
    //   if (! M_fluid->matrixMassPtr().get() )
    //  M_fluid->buildSystem( );

    if (iter == 0)
    {
        M_resetPrec = true;
        M_fluidTimeAdvance->updateRHSContribution ( M_data->dataFluid()->dataTime()->timeStep() );
        if (!M_data->dataFluid()->conservativeFormulation() )
        {
            *M_rhs += M_fluid->matrixMass() * (M_fluidTimeAdvance->rhsContributionFirstDerivative() );
        }
        else
        {
            *M_rhs += (M_fluidMassTimeAdvance->rhsContributionFirstDerivative() );
        }
        couplingRhs (M_rhs/*, M_fluidTimeAdvance->stencil()[0]*/);
    }
    //the conservative formulation as it is now is of order 1. To have higher order (TODO) we need to store the mass matrices computed at the previous time steps.
    //At the moment, the flag conservativeFormulation should be always kept on FALSE
    if (M_data->dataFluid()->conservativeFormulation() )
    {
        M_fluid->updateSystem (alpha, *this->M_beta, *this->M_rhs, M_fluidBlock, solution );
    }
    this->M_fluid->updateStabilization (*M_fluidBlock);
}

void FSIMonolithic::updateRHS()
{
    // Update fluid (iter == 0)
    M_fluidTimeAdvance->updateRHSContribution ( M_data->dataFluid()->dataTime()->timeStep() );
    *M_rhs += M_fluid->matrixMass() * (M_fluidTimeAdvance->rhsContributionFirstDerivative() );
    couplingRhs (M_rhs);

    // Update solid (iter == 0)
    updateSolidSystem (M_rhs);

    // Update RHS
    *M_rhsFull = *M_rhs;
}

namespace
{
static Real fZero (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.;
}
static Real fOne (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 1.;
}
}

void FSIMonolithic::enableStressComputation (UInt flag)
{
    M_BChWS.reset (new solidBchandler_Type() );
    BCFunctionBase bcfZero (fZero);
    BCFunctionBase bcfOne (fOne);
    M_bcfWs.setFunctions_Robin (bcfOne, bcfOne);

    M_BChWS->addBC ("WS", (UInt) flag, Robin, Full, M_bcfWs, 3);
}

FSIMonolithic::vectorPtr_Type FSIMonolithic::computeStress()
{
    vector_Type lambda (M_monolithicMatrix->interfaceMap() );
    lambda.subset (M_fluidTimeAdvance->singleElement (0), M_solidAndFluidDim);

    M_boundaryMass.reset (new matrix_Type (*M_interfaceMap) );
    if ( !M_BChWS->bcUpdateDone() )
    {
        M_BChWS->bcUpdate (*M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
    }
    bcManageMatrix (*M_boundaryMass, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BChWS, M_dFESpace->feBd(), 1., dataSolid()->dataTime()->time() );
    M_boundaryMass->globalAssemble();

    solver_Type solverMass (M_epetraComm);
    solverMass.setDataFromGetPot ( M_dataFile, "solid/solver" );
    solverMass.setupPreconditioner (M_dataFile, "solid/prec"); //to avoid if we have already a prec.

    std::shared_ptr<PreconditionerIfpack> P (new PreconditionerIfpack() );

    vectorPtr_Type sol (new vector_Type (M_monolithicMatrix->interfaceMap() ) );
    solverMass.setMatrix (*M_boundaryMass);
    solverMass.setReusePreconditioner (false);
    solverMass.solveSystem ( lambda, *sol, M_boundaryMass);

    EpetraExt::MultiVector_Reindex reindexMV (*M_interfaceMap->map (Unique) );
    std::shared_ptr<MapEpetra> newMap (new MapEpetra ( *M_interfaceMap ) );
    M_stress.reset (new vector_Type (reindexMV (lambda.epetraVector() ), newMap, Unique) );
    return M_stress;
}

void
FSIMonolithic::checkIfChangedFluxBC ( precPtr_Type oper )
{
    UInt nfluxes (M_BChs[1]->numberOfBCWithType (Flux) );
    if (M_fluxes != nfluxes)
    {
        for (UInt i = 0; i < M_fluxes; ++i)
        {
            const BCBase* bc (M_BChs[1]->findBCWithName (M_BCFluxNames[i]) );
            if (bc->type() != Flux)
            {
                oper->addToCoupling (1., M_fluxOffset[i], M_fluxOffset[i], 1 );
            }
        }
    }
}


}
