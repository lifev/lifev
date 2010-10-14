/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

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


#ifndef TWODIM
#include <life/lifesolver/exactJacobianBase.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>
//#include <life/lifesolver/reducedLinFluid.hpp>


namespace LifeV
{

Real fzeroEJ(const Real& /*t*/,
             const Real& /*x*/,
             const Real& /*y*/,
             const Real& /*z*/,
             const ID& /*i*/)
{return 0.0;}


exactJacobian::exactJacobian():
    super       (),
    M_matrShapeDer()
//    M_epetraOper(this)
{

//     M_epetraOper = new Epetra_ExactJacobian();
//    M_epetraOper->setOperator(this);
}



exactJacobian::~exactJacobian()
{}

void
exactJacobian::setDataFile( const GetPot& dataFile )
{
    super::setDataFile( dataFile );

    M_aitkFS.setDefaultOmega( M_data->defaultOmega(), 0.001 );
    M_aitkFS.setOmegaRange( M_data->OmegaRange() );

    M_linearSolver.setDataFromGetPot(dataFile, "jacobian");

}

void
exactJacobian::setupFEspace()
{
    setLinearFluid(true);
    setLinearSolid(true);

    super::setupFEspace();
}

void
exactJacobian::setupFluidSolid()
{
    super::setupFluidSolid();

    M_epetraOper.reset( new Epetra_ExactJacobian(this));

    if ( this->isFluid() )
    {
        M_rhsNew.reset(new vector_type(this->M_fluid->getMap()));
        M_beta.reset  (new vector_type(this->M_fluid->getMap()));
    }


//    M_reducedLinFluid.reset(new reducedLinFluid(this, M_fluid, M_solid));
}


UInt
exactJacobian::imposeFlux( void )
{
    if ( this->isFluid() )
    {
        UInt numLM = super::imposeFlux();

        std::vector<BCName> fluxVector = M_BCh_du->getBCWithType( Flux );
        if( numLM != ( static_cast<UInt> ( fluxVector.size() ) ) )
        {
            ERROR_MSG("Different number of fluxes imposed on Fluid and on LinearFluid");
        }

        UInt offset = M_uFESpace->map().getMap(Unique)->NumGlobalElements()
                    + M_pFESpace->map().getMap(Unique)->NumGlobalElements();

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

void exactJacobian::eval(const vector_type& _disp,
                         const UInt          iter)
{
    Chrono chronoFluid, chronoSolid, chronoInterface;

    bool recomputeMatrices ( iter == 0 || ( !this->M_data->dataFluid()->isSemiImplicit() &&
                                            ( M_data->updateEvery() > 0 &&
                                              (iter % M_data->updateEvery() == 0) ) ) );

    if(iter == 0)
    {
        if (isFluid())
            this->M_fluid->resetPrec();
        //this->M_solid->resetPrec();
    }


    this->setLambdaFluid(_disp);


    vector_type sigmaFluidUnique (this->sigmaFluid(), Unique);

    M_epetraWorldComm->Barrier();
    chronoFluid.start();

    if (this->isFluid())
    {
        this->M_meshMotion->iterate(*M_BCh_mesh);
        this->M_meshMotion->updateDispDiff();

        this->transferMeshMotionOnFluid(M_meshMotion->disp(),
                                        this->veloFluidMesh());

        this->veloFluidMesh()    -= dispFluidMeshOld();
        this->veloFluidMesh()    *= 1./(M_data->dataFluid()->dataTime()->getTimeStep());

        if( iter==0 || !this->M_data->dataFluid()->isSemiImplicit() )
            {
                // copying displacement to a repeated indeces displacement, otherwise the mesh wont know
                // the value of the displacement for some points


                *M_beta *= 0.;

                vector_type meshDisp( M_meshMotion->disp(), Repeated );

                this->moveMesh(meshDisp);

                vector_type meshDispDiff( M_meshMotion->dispDiff(), Repeated );
                this->interpolateVelocity(meshDispDiff, *M_beta);

                *M_beta *= -1./M_data->dataFluid()->dataTime()->getTimeStep();

                *M_beta  += *this->M_un;

                if(recomputeMatrices)
                    {
                        double alpha = 1./M_data->dataFluid()->dataTime()->getTimeStep();
                        this->M_fluid->updateSystem( alpha, *M_beta, *M_rhs );
                    }
                else
                    {
                        this->M_fluid->updateRHS( *M_rhs );
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
        vector_type vel  (this->fluid().velFESpace().map());
        vector_type press(this->fluid().pressFESpace().map());

        vel.subset(*this->M_fluid->solution());
        press.subset(*this->M_fluid->solution(), this->fluid().velFESpace().dim()*this->fluid().pressFESpace().fieldDim());

        std::cout << "norm_inf( vel ) " << vel.NormInf() << std::endl;
        std::cout << "norm_inf( press ) " << press.NormInf() << std::endl;

        //this->M_fluid->postProcess();
    }

    chronoInterface.start();

    this->setSigmaFluid( sigmaFluidUnique );
    this->setSigmaSolid( sigmaFluidUnique );



    vector_type lambdaSolidUnique   (this->lambdaSolid(),    Unique);
    vector_type lambdaDotSolidUnique(this->lambdaDotSolid(), Unique);
    vector_type sigmaSolidUnique    (this->sigmaSolid(),     Unique);

    chronoInterface.stop();
    M_epetraWorldComm->Barrier();
    chronoSolid.start();

    if (this->isSolid())
    {
        this->M_solid->iterate( M_BCh_d );
//         this->transferSolidOnInterface(this->M_solid->disp(),     lambdaSolidUnique);
//         this->transferSolidOnInterface(this->M_solid->vel(),      lambdaDotSolidUnique);
//         this->transferSolidOnInterface(this->M_solid->residual(), sigmaSolidUnique);
    }

    M_epetraWorldComm->Barrier();
    chronoSolid.stop();
    this->displayer().leaderPrintMax("      Solid solution total time:               ", chronoSolid.diff() );

    chronoInterface.start();

    if (this->isSolid())
        {
            this->transferSolidOnInterface(this->M_solid->disp(),     lambdaSolidUnique);
            this->transferSolidOnInterface(this->M_solid->vel(),      lambdaDotSolidUnique);
            this->transferSolidOnInterface(this->M_solid->residual(), sigmaSolidUnique);
        }

    this->setLambdaSolid(    lambdaSolidUnique);
    this->setLambdaDotSolid( lambdaDotSolidUnique);
    this->setSigmaSolid(     sigmaSolidUnique);

    chronoInterface.stop();
    this->displayer().leaderPrintMax("      Interface transfer total time:           ", chronoInterface.diff_cumul() );

    if ( false && this->isSolid() )
    {
        //this->solid().postProcess();
    }

// possibly unsafe when using more cpus, since both has repeated maps

    this->displayer().leaderPrint("      Norm(disp     )    =                     ", _disp.NormInf(), "\n");
    this->displayer().leaderPrint("      Norm(dispNew  )    =                     " , this->lambdaSolid().NormInf(), "\n" );
    this->displayer().leaderPrint("      Norm(velo     )    =                     " , this->lambdaDotSolid().NormInf(), "\n" );
    this->displayer().leaderPrint("      Max Residual Fluid =                     " , this->sigmaFluid().NormInf(), "\n" );
    this->displayer().leaderPrint("      Max Residual Solid =                     " , this->sigmaSolid().NormInf(), "\n" );

    if ( this->isFluid() )
        this->displayer().leaderPrint("      Max ResidualF      =                     " , M_fluid->residual().NormInf(), "\n" );
    if ( this->isSolid() )
    {
        this->displayer().leaderPrint("      NL2 Diplacement S. =                     " , M_solid->disp().Norm2(), "\n" );
        this->displayer().leaderPrint("      Max Residual Solid =                     " , M_solid->residual().NormInf(), "\n" );
    }
}

void exactJacobian::evalResidual(vector_type&       res,
                                 const vector_type& disp,
                                 const UInt          iter)
{
    if (this->isSolid())
    {
        std::cout << "      Residual computation g(x_" << iter <<" )\n";
    }

    this->setLambdaSolidOld(disp);


    eval(disp, iter);

    res  = this->lambdaSolid();
    res -=  disp;

    this->displayer().leaderPrint("      NormInf res        =                     " , res.NormInf(), "\n" );
    if (this->isSolid())
        this->displayer().leaderPrint("      NormInf res_d      =                     " , this->solid().residual().NormInf(), "\n" );

}


//
// new step computation resolution
//


void  exactJacobian::solveJac(vector_type         &_muk,
                              const vector_type   &_res,
                              const Real         _linearRelTol)
{
    if (this->isFluid() && this->isLeader()) std::cout << "  f- ";
    if (this->isSolid() && this->isLeader()) std::cout << "  s- ";

    this->displayer().leaderPrint( "solveJac: NormInf res " , _res.NormInf(), "\n" );
    _muk *= 0.;

    M_linearSolver.setTolMaxiter(_linearRelTol, 100);

    vector_type res(_res);

    M_linearSolver.setOperator(*M_epetraOper);

    this->displayer().leaderPrint( "Solving Jacobian system... " );

    M_recomputeShapeDer=true;
    M_linearSolver.solve(_muk, res);

    this->displayer().leaderPrint( "Solving the Jacobian system done.\n" );
}


void  exactJacobian::solveLinearFluid()
{
    //vector_type dispFluidDomainRep( M_fluid->matrNoBC().getMap(), Repeated);
    vector_type dispFluidDomain( M_fluid->matrNoBC().getMap(), Unique, Zero);
    dispFluidDomain.setCombineMode(Zero);
    vector_type dispFluidMesh(this->derVeloFluidMesh().getMap(), Repeated);
//if statement: in order not to iterate the mesh for each linear residual calculation, needed just for exact Jac case.
    if(false && this->M_data->dataFluid()->isSemiImplicit()==true)// not working in parallel
        {//to be corrected: up to now also in the semi implicit case the harmonic extension eq.
            //is solved at each GMRES iteration

            //vector_type solidDisplacementRepeated(this->M_solid->disp(), Repeated);

            //            this->transferSolidOnInterface(sldsp, lambdaSolid());

            //            vector_type repLambdaSolid(lambdaSolid(), Repeated);

            //this->transferInterfaceOnFluid(repLambdaSolid, dispFluidMesh);
            this->transferSolidOnFluid(this->M_solid->disp(), dispFluidMesh);
        }
    else
        {

            this->transferMeshMotionOnFluid(M_meshMotion->disp(),
                                            dispFluidMesh);
        }
    //dispFluidDomainRep = dispFluidMesh;//import
    dispFluidDomain=dispFluidMesh;//import
    this->derVeloFluidMesh() = dispFluidMesh;
    this->derVeloFluidMesh() *= 1./(M_data->dataFluid()->dataTime()->getTimeStep());
    this->displayer().leaderPrint( " norm inf dw = " , this->derVeloFluidMesh().NormInf(), "\n" );
    *M_rhsNew *= 0.;

    double alpha = this->M_bdf->coeff_der( 0 ) / M_data->dataFluid()->dataTime()->getTimeStep();

    if(!this->M_fluid->stab())//if using P1Bubble
    {

        this->M_fluid->updateLinearSystem( M_fluid->matrNoBC(),
                                           alpha,
                                           *M_un,
                                           *M_fluid->solution(),
                                           dispFluidMesh,
                                           this->veloFluidMesh(),
                                           this->derVeloFluidMesh(),
                                           *M_rhsNew );

    }
    else
     {
        if(M_recomputeShapeDer)
        {
            M_recomputeShapeDer=false;
            M_matrShapeDer.reset(new matrix_type(M_fluid->matrNoBC().getMap()/*, M_mmFESpace->map()*/));
            this->M_fluid->updateShapeDerivatives(
                                                  *M_matrShapeDer,
                                                  alpha,
                                                  *M_un,
                                                  *M_fluid->solution(),
                                                  //dispFluidMesh,
                                                  this->veloFluidMesh(),
                                                  (UInt)0,
                                                  *M_mmFESpace,
                                                  true,
                                                  false
                                                  //this->derVeloFluidMesh(),
                                                  //*M_rhsNew
                                                  );
            M_matrShapeDer->GlobalAssemble();
            *M_matrShapeDer*=-1;
            //M_matrShapeDer->spy("matrsd");
        }
        *M_rhsNew=(*M_matrShapeDer)*dispFluidDomain;
        M_fluid->updateRhsLinNoBC(*M_rhsNew);
    }

//    //DEBUG:
//     vector_type rhs_debug(M_fluid->matrNoBC().getMap());
//     rhs_debug=*M_matrShapeDer*dispFluidDomain;
//     std::cout<<"normInf1 "<<rhs_debug.NormInf()<<std::endl;
//     //    std::cout<<"normInf2 "<<M_rhsLin->NormInf()<<std::endl;
//     std::cout<<"normInf displ "<<dispFluidMesh.NormInf()<<std::endl;
//     std::cout<<"normInf solution "<<M_fluid->solution().NormInf()<<std::endl;
//     std::cout<<"normInf sould be 0 "<<rhs_debug.NormInf()<<std::endl;
//     //END DEBUG

    this->M_fluid->iterateLin( *this->M_BCh_du );
}


//


void  exactJacobian::solveLinearSolid()
{
    this->M_solid->iterateLin( M_BCh_dz );
}



int Epetra_ExactJacobian::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{

    Chrono chronoFluid, chronoSolid, chronoInterface;

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
        FSIOperator::vector_type const z(X,  M_ej->solidInterfaceMap(), Unique);

        M_ej->displayer().leaderPrint( "NormInf res   " , z.NormInf(), "\n" );

        //M_ej->solid().residual() *= 0.;
        //M_ej->transferInterfaceOnSolid(z, M_ej->solid().residual());

        //std::cout << "NormInf res_d " << M_ej->solid().residual().NormInf() << std::endl;

        //if (M_ej->isSolid())
        //    M_ej->solid().postProcess();

        M_ej->setLambdaFluid(z);

        //M_ej->transferInterfaceOnSolid(z, M_ej->solid().disp());

        chronoInterface.start();
        vector_type sigmaFluidUnique (M_ej->sigmaFluid(), Unique);
        chronoInterface.stop();

        M_comm->Barrier();
        chronoFluid.start();

        if (M_ej->isFluid())
            {

                //to be used when we correct the other lines
                if(true || ( !this->M_ej->dataFluid()->isSemiImplicit() /*|| this->M_ej->dataFluid().semiImplicit()==-1*/))
                    {
                        M_ej->meshMotion().iterate(*M_ej->BCh_harmonicExtension());
                        //std::cout<<" mesh motion iterated!!!"<<std::endl;
                    }

                M_ej->displayer().leaderPrint( " norm inf dx = " , M_ej->meshMotion().disp().NormInf(), "\n" );

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


        vector_type lambdaSolidUnique (M_ej->lambdaSolid(), Unique);
        chronoInterface.stop();

        M_comm->Barrier();
        chronoFluid.start();

        if (M_ej->isSolid())
            {
                M_ej->solveLinearSolid();
                M_ej->transferSolidOnInterface(M_ej->solid().disp(), lambdaSolidUnique);
            }

        M_comm->Barrier();
        chronoSolid.stop();
        M_ej->displayer().leaderPrintMax( "Solid linear solution: total time : " , chronoSolid.diff() );

        chronoInterface.start();
        M_ej->setLambdaSolid(lambdaSolidUnique);

        chronoInterface.stop();
        M_ej->displayer().leaderPrintMax( "Interface linear transfer: total time : " , chronoInterface.diff_cumul() );

        dz = lambdaSolidUnique.getEpetraVector();

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

void  exactJacobian::bcManageVec( fluid_bchandler_type& /*bch*/, vector_type& /*rhs*/ )
{

}

void exactJacobian::registerMyProducts( )
{
FSIFactory::instance().registerProduct( "exactJacobian", &createEJ );
solid_raw_type::StructureSolverFactory::instance().registerProduct( "LinearVenantKirchhof", &createLinearStructure );
solid_raw_type::StructureSolverFactory::instance().registerProduct( "NonLinearVenantKirchhof", &createNonLinearStructure );
}


}
#endif
