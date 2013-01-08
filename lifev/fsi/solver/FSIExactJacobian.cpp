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

#include <lifev/fsi/solver/FSIExactJacobian.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
FSIExactJacobian::FSIExactJacobian():
        super(),
        M_rhsNew(),
        M_beta(),
        M_aitkFS(),
        M_linearSolver(),
        M_epetraOper(),
        M_matrShapeDer(),
        M_recomputeShapeDer(true)
{

}



FSIExactJacobian::~FSIExactJacobian()
{}

// ===================================================
// Methods
// ===================================================
void  FSIExactJacobian::solveJac(vector_Type         &_muk,
                              const vector_Type   &_res,
                              const Real         _linearRelTol)
{
    if (this->isFluid() && this->isLeader()) std::cout << "  f- ";
    if (this->isSolid() && this->isLeader()) std::cout << "  s- ";

    this->displayer().leaderPrint( "FSI --- solveJac: NormInf res " , _res.normInf(), "\n" );
    _muk *= 0.;

    M_linearSolver.setTolerance( _linearRelTol );
    M_linearSolver.setMaxNumIterations( 100 );

    vector_Type res(_res);

    M_linearSolver.setOperator(M_epetraOper);

    this->displayer().leaderPrint( "FSI --- Solving Jacobian FSI system... " );

    M_recomputeShapeDer=true;
    M_linearSolver.solve(_muk, res);

    this->displayer().leaderPrint( "FSI --- Solving the Jacobian FSI system done.\n" );
}

void FSIExactJacobian::evalResidual(vector_Type&       res,
                                 const vector_Type& disp,
                                 const UInt          iter)
{
    if (this->isSolid())
    {
        std::cout << "FSI ---      Residual computation g(x_" << iter <<" )\n";
    }

    setLambdaSolidOld(disp);


    eval(disp, iter);

    res  = this->lambdaSolid();
    res -=  disp;

    if (this->isSolid())
      M_solid->updateJacobian(M_solid->displacement());

    this->displayer().leaderPrint("FSI ---      NormInf res        =                     " , res.normInf(), "\n" );
    if (this->isSolid())
        this->displayer().leaderPrint("FSI ---      NormInf res_d      =                     " , this->solid().residual().normInf(), "\n" );

}


void  FSIExactJacobian::solveLinearFluid()
{
    //vector_Type dispFluidDomainRep( M_fluid->matrixNoBC().getMap(), Repeated);
    vector_Type dispFluidDomain( M_fluid->matrixNoBC().map(), Unique, Zero);
    dispFluidDomain.setCombineMode(Zero);
    vector_Type dispFluidMesh(this->derVeloFluidMesh().map(), Repeated);
//if statement: in order not to iterate the mesh for each linear residual calculation, needed just for exact Jac case.
    if (false && this->M_data->dataFluid()->isSemiImplicit()==true)
    {//to be corrected: up to now also in the semi implicit case the harmonic extension eq.
        //is solved at each GMRES iteration

        //vector_Type solidDisplacementRepeated(this->M_solid->disp(), Repeated);

        //            this->transferSolidOnInterface(sldsp, lambdaSolid());

        //            vector_Type repLambdaSolid(lambdaSolid(), Repeated);

        //this->transferInterfaceOnFluid(repLambdaSolid, dispFluidMesh);
        this->transferSolidOnFluid(this->M_solid->displacement(), dispFluidMesh);
    }
    else
    {
        this->transferMeshMotionOnFluid(M_meshMotion->disp(),
                                        dispFluidMesh);
    }
    //dispFluidDomainRep = dispFluidMesh;//import
    dispFluidDomain=dispFluidMesh;//import


    this->derVeloFluidMesh() = dispFluidMesh;
    this->derVeloFluidMesh() *= M_fluidTimeAdvance->coefficientFirstDerivative(0)/(M_data->dataFluid()->dataTime()->timeStep());
    this->displayer().leaderPrint( " norm inf dw = " , derVeloFluidMesh().normInf(), "\n" );
    *M_rhsNew *= 0.;

    double alpha = this->M_fluidTimeAdvance->coefficientFirstDerivative( 0 ) / M_data->dataFluid()->dataTime()->timeStep();

    if (!this->M_fluid->stabilization())//if using P1Bubble
    {

        this->M_fluid->updateLinearSystem( M_fluid->matrixNoBC(),
                                           alpha,
                                           M_fluidTimeAdvance->singleElement(0),
                                           *M_fluid->solution(),
                                           dispFluidMesh,
                                           veloFluidMesh(),
                                           derVeloFluidMesh(),
                                           *M_rhsNew );
    }
    else
    {
        if (M_recomputeShapeDer)
        {
            M_recomputeShapeDer=false;
            M_matrShapeDer.reset(new matrix_Type(M_fluid->matrixNoBC().map()/*, M_mmFESpace->map()*/));
            this->M_fluid->updateShapeDerivatives(
                                                  *M_matrShapeDer,
                                                  alpha,
                                                  M_fluidTimeAdvance->singleElement(0),
                                                  *M_fluid->solution(),
                                                  //dispFluidMesh,
                                                  veloFluidMesh(),
                                                  (UInt)0,
                                                  *M_mmFESpace,
                                                  true,
                                                  false
                                                  //this->derVeloFluidMesh(),
                                                  //*M_rhsNew
                                                  );
            M_matrShapeDer->globalAssemble();
            *M_matrShapeDer*=-1;
            //M_matrShapeDer->spy("matrsd");
        }
        *M_rhsNew=(*M_matrShapeDer)*dispFluidDomain;
        M_fluid->updateLinearRightHandSideNoBC(*M_rhsNew);
    }

    this->M_fluid->solveLinearSystem( *this->M_BCh_du );
}


void  FSIExactJacobian::solveLinearSolid()
{
    M_solid->iterateLin( M_BCh_dz );
}

void
FSIExactJacobian::setupFEspace()
{
    setLinearFluid(true);
    setLinearSolid(true);

    super::setupFEspace();
}

void
FSIExactJacobian::setupFluidSolid()
{
    super::setupFluidSolid();

    M_epetraOper.setOperator(this);

    if ( this->isFluid() )
    {
        M_rhsNew.reset(new vector_Type(this->M_fluid->getMap()));
        M_beta.reset  (new vector_Type(this->M_fluid->getMap()));
    }

}

void
FSIExactJacobian::setDataFile( const GetPot& dataFile )
{
    super::setDataFile( dataFile );

    M_aitkFS.setDefaultOmega( M_data->defaultOmega(), 0.001 );
    M_aitkFS.setOmegaRange( M_data->OmegaRange() );

    M_linearSolver.setCommunicator( M_epetraComm );
    M_linearSolver.setDataFromGetPot(dataFile, "jacobian");

}


void FSIExactJacobian::registerMyProducts( )
{
    FSIFactory_Type::instance().registerProduct( "exactJacobian", &createEJ );

    solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "linearVenantKirchhoff", &super::createVenantKirchhoffLinear );
    solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "nonLinearVenantKirchhoff", &super::createVenantKirchhoffNonLinear );
}

// ===================================================
// Private Methods
// ===================================================


UInt
FSIExactJacobian::imposeFlux( void )
{
    if ( this->isFluid() )
    {
        UInt numLM = super::imposedFluxes();

        std::vector<bcName_Type> fluxVector = M_BCh_du->findAllBCWithType( Flux );
        if ( numLM != ( static_cast<UInt> ( fluxVector.size() ) ) )
        {
            ERROR_MSG("Different number of fluxes imposed on Fluid and on LinearFluid");
        }

        UInt offset = M_uFESpace->map().map(Unique)->NumGlobalElements()
                      + M_pFESpace->map().map(Unique)->NumGlobalElements();

        for ( UInt i = 0; i < numLM; ++i )
            M_BCh_du->setOffset( fluxVector[i], offset + i );

        return numLM;
    }
    else
        return 0;
}

//
// Residual computation
//

void FSIExactJacobian::eval(const vector_Type& _disp,
                         const UInt          iter)
{
    LifeChrono chronoFluid, chronoSolid, chronoInterface;

    bool recomputeMatrices ( iter == 0 || ( !this->M_data->dataFluid()->isSemiImplicit() &&
                                            ( M_data->updateEvery() > 0 &&
                                              (iter % M_data->updateEvery() == 0) ) ) );

    if (iter == 0)
    {
        if (isFluid())
            this->M_fluid->resetPreconditioner();
        //this->M_solid->resetPrec();
    }


    this->setLambdaFluid(_disp);


    vector_Type sigmaFluidUnique (this->sigmaFluid(), Unique);

    M_epetraWorldComm->Barrier();
    chronoFluid.start();

    if (isFluid())
    {
        M_meshMotion->iterate(*M_BCh_mesh);
        // M_meshMotion->updateDispDiff();


            // copying displacement to a repeated indeces displacement, otherwise the mesh wont know
            // the value of the displacement for some points


        if ( iter==0 || !this->M_data->dataFluid()->isSemiImplicit() )
        {
            *M_beta *= 0.;

            if( iter==0 )
            {
                M_ALETimeAdvance->updateRHSFirstDerivative(M_data->dataFluid()->dataTime()->timeStep());
                M_ALETimeAdvance->shiftRight(M_meshMotion->disp());
                M_ALETimeAdvance->extrapolation(M_meshMotion->disp());//closer initial solution
            }
            else
            {
                M_ALETimeAdvance->setSolution(M_meshMotion->disp());
            }

            vector_Type meshDispRepeated( M_meshMotion->disp(), Repeated );
            this->moveMesh(meshDispRepeated);
            vector_Type vel ( this->M_ALETimeAdvance->firstDerivative( ) );
            transferMeshMotionOnFluid( vel, veloFluidMesh() );
            M_fluidTimeAdvance->extrapolation(*M_beta);//explicit
            *M_beta -= veloFluidMesh();//implicit


	    //the conservative formulation as it is now is of order 1. To have higher order (TODO) we need to store the mass matrices computed at the previous time steps.
	    if(M_data->dataFluid()->conservativeFormulation())
	      //*M_rhs = M_fluid->matrixMass()*M_fluidTimeAdvance->rhsContributionFirstDerivative();
	      *M_rhs = M_fluidMassTimeAdvance->rhsContributionFirstDerivative();


            if (recomputeMatrices)
            {
                double alpha = M_fluidTimeAdvance->coefficientFirstDerivative(0)/M_data->dataFluid()->dataTime()->timeStep();
                this->M_fluid->updateSystem( alpha, *M_beta, *M_rhs );
            }
            else
            {
                this->M_fluid->updateRightHandSide( *M_rhs );
            }
	    if(!M_data->dataFluid()->conservativeFormulation())
	      {
		*M_rhs = M_fluid->matrixMass()*M_fluidTimeAdvance->rhsContributionFirstDerivative();
		this->M_fluid->updateRightHandSide( *M_rhs );
	      }
        }

        this->M_fluid->iterate( *M_BCh_u );

        this->transferFluidOnInterface(this->M_fluid->residual(), sigmaFluidUnique);
    }

    M_epetraWorldComm->Barrier();
    chronoFluid.stop();
    this->displayer().leaderPrintMax("      Fluid solution total time:               ", chronoFluid.diff() );

    if ( false && this->isFluid() )
    {
        vector_Type vel  (this->fluid().velocityFESpace().map());
        vector_Type press(this->fluid().pressureFESpace().map());

        vel.subset(*this->M_fluid->solution());
        press.subset(*this->M_fluid->solution(), this->fluid().velocityFESpace().dim()*this->fluid().pressureFESpace().fieldDim());

        std::cout << "norm_inf( vel ) " << vel.normInf() << std::endl;
        std::cout << "norm_inf( press ) " << press.normInf() << std::endl;

        //this->M_fluid->postProcess();
    }

    chronoInterface.start();

    this->setSigmaFluid( sigmaFluidUnique );
    this->setSigmaSolid( sigmaFluidUnique );



    vector_Type lambdaSolidUnique   (this->lambdaSolid(),    Unique);
    vector_Type lambdaDotSolidUnique(this->lambdaDotSolid(), Unique);
    vector_Type sigmaSolidUnique    (this->sigmaSolid(),     Unique);

    chronoInterface.stop();
    M_epetraWorldComm->Barrier();
    chronoSolid.start();

    if (this->isSolid())
    {
        this->M_solid->iterate( M_BCh_d );
    }

    M_epetraWorldComm->Barrier();
    chronoSolid.stop();
    this->displayer().leaderPrintMax("      Solid solution total time:               ", chronoSolid.diff() );

    chronoInterface.start();

    if (this->isSolid())
    {
        this->transferSolidOnInterface(this->M_solid->displacement(),     lambdaSolidUnique);
        this->transferSolidOnInterface( this->M_solidTimeAdvance->firstDerivative( this->M_solid->displacement()), lambdaDotSolidUnique );
        this->transferSolidOnInterface(this->M_solid->residual(), sigmaSolidUnique);
    }

    this->setLambdaSolid(    lambdaSolidUnique);
    this->setLambdaDotSolid( lambdaDotSolidUnique);
    this->setSigmaSolid(     sigmaSolidUnique);

    chronoInterface.stop();
    this->displayer().leaderPrintMax("      Interface transfer total time:           ", chronoInterface.diffCumul() );


// possibly unsafe when using more cpus, since both has repeated maps

    this->displayer().leaderPrint("      Norm(disp     )    =                     ", _disp.normInf(), "\n");
    this->displayer().leaderPrint("      Norm(dispNew  )    =                     " , this->lambdaSolid().normInf(), "\n" );
    this->displayer().leaderPrint("      Norm(velo     )    =                     " , this->lambdaDotSolid().normInf(), "\n" );
    this->displayer().leaderPrint("      Max Residual Fluid =                     " , this->sigmaFluid().normInf(), "\n" );
    this->displayer().leaderPrint("      Max Residual Solid =                     " , this->sigmaSolid().normInf(), "\n" );

    if ( this->isFluid() )
        this->displayer().leaderPrint("      Max ResidualF      =                     " , M_fluid->residual().normInf(), "\n" );
    if ( this->isSolid() )
    {
        this->displayer().leaderPrint("      NL2 Diplacement S. =                     " , M_solid->displacement().norm2(), "\n" );
        this->displayer().leaderPrint("      Max Residual Solid =                     " , M_solid->residual().normInf(), "\n" );
    }
}






// ===================================================
// Epetra_ExactJacobian
// ===================================================

// ===================================================
// Constructors & Destructorxs
// ===================================================

FSIExactJacobian::Epetra_ExactJacobian::Epetra_ExactJacobian() :
    M_ej(),
    M_operatorDomainMap(),
    M_operatorRangeMap(),
    M_comm()
{
}

// ===================================================
// Methods
// ===================================================

void FSIExactJacobian::Epetra_ExactJacobian::setOperator(FSIExactJacobian* ej)
{
    ASSERT(ej != 0, "passing NULL pointer to se operator");
    M_ej                = ej;
    M_operatorDomainMap = M_ej->solidInterfaceMap()->map(Repeated);
    M_operatorRangeMap  = M_ej->solidInterfaceMap()->map(Repeated);
    M_comm              = M_ej->worldComm();
}

int FSIExactJacobian::Epetra_ExactJacobian::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{

    LifeChrono chronoFluid, chronoSolid, chronoInterface;

    M_comm->Barrier();

    double xnorm = 0.;
    X.NormInf(&xnorm);

    Epetra_FEVector  dz(Y.Map());

    if (M_ej->isSolid())
        std::cout << "\n ***** norm (z)= " << xnorm << std::endl << std::endl;


    if ( xnorm == 0.0 )
    {
        Y.Scale(0.);
        dz.Scale(0.);
    }
    else
    {
        vector_Type const z(X,  M_ej->solidInterfaceMap(), Unique);

        M_ej->displayer().leaderPrint( "NormInf res   " , z.normInf(), "\n" );

        //M_ej->solid().residual() *= 0.;
        //M_ej->transferInterfaceOnSolid(z, M_ej->solid().residual());

        //std::cout << "NormInf res_d " << M_ej->solid().residual().NormInf() << std::endl;

        //if (M_ej->isSolid())
        //    M_ej->solid().postProcess();

        M_ej->setLambdaFluid(z);

        //M_ej->transferInterfaceOnSolid(z, M_ej->solid().disp());

        chronoInterface.start();
        vector_Type sigmaFluidUnique (M_ej->sigmaFluid(), Unique);
        chronoInterface.stop();

        M_comm->Barrier();
        chronoFluid.start();

        if (M_ej->isFluid())
        {

            //to be used when we correct the other lines
            if (true || ( !this->M_ej->dataFluid()->isSemiImplicit() /*|| this->M_ej->dataFluid().semiImplicit()==-1*/))
            {
                M_ej->meshMotion().iterate(*M_ej->BCh_harmonicExtension());
                //std::cout<<" mesh motion iterated!!!"<<std::endl;
            }

            M_ej->displayer().leaderPrint( " norm inf dx = " , M_ej->meshMotion().disp().normInf(), "\n" );

            M_ej->solveLinearFluid();

            M_ej->transferFluidOnInterface(M_ej->fluid().residual(), sigmaFluidUnique);

            //M_ej->fluidPostProcess();
        }

        M_comm->Barrier();
        chronoFluid.stop();
        M_ej->displayer().leaderPrintMax( "Fluid linear solution: total time : ", chronoFluid.diff() );


        chronoInterface.start();
        // M_ej->setSigmaFluid(sigmaFluidUnique);
        M_ej->setSigmaSolid(sigmaFluidUnique);


        vector_Type lambdaSolidUnique (M_ej->lambdaSolid(), Unique);
        chronoInterface.stop();

        M_comm->Barrier();
        chronoFluid.start();

        if (M_ej->isSolid())
        {
            M_ej->solveLinearSolid();
            M_ej->transferSolidOnInterface(M_ej->solid().displacement(), lambdaSolidUnique);
        }

        M_comm->Barrier();
        chronoSolid.stop();
        M_ej->displayer().leaderPrintMax( "Solid linear solution: total time : " , chronoSolid.diff() );

        chronoInterface.start();
        M_ej->setLambdaSolid(lambdaSolidUnique);

        chronoInterface.stop();
        M_ej->displayer().leaderPrintMax( "Interface linear transfer: total time : " , chronoInterface.diffCumul() );

        dz = lambdaSolidUnique.epetraVector();
    }


    Y = X;
    Y.Update(1., dz, -1.);

    double ynorm;
    Y.NormInf(&ynorm);

    if (M_ej->isSolid())
        std::cout << "\n\n ***** norm (Jz)= " << ynorm
                  << std::endl << std::endl;

    return 0;
}

}
