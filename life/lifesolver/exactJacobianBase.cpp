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
    M_updateEvery(0),
    M_nbEvalAux(0)
//    M_epetraOper(this)
{

//     M_epetraOper = new Epetra_ExactJacobian();
//    M_epetraOper->setOperator(this);
}



exactJacobian::~exactJacobian()
{}

void
exactJacobian::setDataFromGetPot( GetPot const& data )
{
    // 0-order transpiration bnd conditions
    // If M_updateEvery == 1, normal Newton algorithm
    // If M_updateEvery  > 1, recompute computational domain every M_updateEvery iterations (transpiration)
    // If M_updateEvery <= 0, recompute computational domain and matrices only at first subiteration (semi-implicit)
      M_updateEvery = data("problem/updateEvery", 0);

     // call the super class to setup the data from getpot file if needed
    super::setDataFromGetPot( data );

    M_linearSolver.setDataFromGetPot(data, "jacobian");

    if (false && this->isFluid())
        M_ensightFluid.reset( new  Ensight<RegionMesh3D<LinearTetra> > ( data, "fixedPtFluidInnerIterations") );



    M_defOmega = data("problem/defOmega", 0.001);

    this->M_dataFluid->setSemiImplicit( data("problem/ifSemiImplicit", 0) );

    Debug( 6205 ) << "fixedPoint::setDataFromGetPot(GetPot) OmegaS = " << M_defOmega << "\n";

    M_aitkFS.setDefault(M_defOmega, 0.001);

}


void
exactJacobian::setup()
{
    setLinearFluid(true);
    setLinearSolid(true);

    super::setup();

    M_epetraOper.reset( new Epetra_ExactJacobian(this));

    if ( this->isFluid() )
    {
        M_rhsNew.reset(new vector_type(this->M_fluid->getMap()));
        M_beta.reset  (new vector_type(this->M_fluid->getMap()));
    }

    if ( false && this->isFluid())
        {
            if (M_ensightFluid.get() == 0)
                M_ensightFluid.reset( new  Ensight<RegionMesh3D<LinearTetra> > ( "", "fixedPtFluidInnerIterations") );

            assert( M_uFESpace.get() );
            assert( M_uFESpace->mesh().get() );

            M_ensightFluid->setMeshProcId(M_uFESpace->mesh(), M_uFESpace->map().Comm().MyPID());

            assert( this->M_fluid.get() );
            M_velAndPressure.reset( new vector_type( this->M_fluid->getMap(), Repeated ));

            M_ensightFluid->addVariable( ExporterData::Vector, "velocityInner", M_velAndPressure,
                                         UInt(0), M_uFESpace->dof().numTotalDof() );

            M_ensightFluid->addVariable( ExporterData::Scalar, "pressureInner", M_velAndPressure,
                                         UInt(3*M_uFESpace->dof().numTotalDof()),
                                         UInt(3*M_uFESpace->dof().numTotalDof()+M_pFESpace->dof().numTotalDof()) );

        }


//    M_reducedLinFluid.reset(new reducedLinFluid(this, M_fluid, M_solid));
}
//
// Residual computation
//

void exactJacobian::eval(const vector_type& _disp,
                         const int          iter)
{
    Chrono chronoFluid, chronoSolid, chronoInterface;

    bool recomputeMatrices ( M_updateEvery > 0 && iter % M_updateEvery == 0 );

    if(iter == 0)
        {
            M_nbEval = 0; // new time step
            this->M_fluid->resetPrec();
            //this->M_solid->resetPrec();
        }

    M_nbEval++ ;

    // possibly unsafe when using more cpus, since both has repeated maps
    this->setLambdaFluid(_disp);


    vector_type sigmaFluidUnique (this->sigmaFluid(), Unique);

    M_epetraWorldComm->Barrier();
    chronoFluid.start();

    if (this->isFluid())
    {
        this->M_meshMotion->iterate();

        this->transferMeshMotionOnFluid(M_meshMotion->disp(),
                                        this->veloFluidMesh());

        this->veloFluidMesh()    -= dispFluidMeshOld();
        this->veloFluidMesh()    *= 1./(M_dataFluid->timestep());

        if((iter==0 && this->M_dataFluid->semiImplicit())|| !this->M_dataFluid->semiImplicit() || this->M_dataFluid->semiImplicit() == -1)
            {
                // copying displacement to a repeated indeces displacement, otherwise the mesh wont know
                // the value of the displacement for some points


                *M_beta *= 0.;

                vector_type const meshDispDiff( M_meshMotion->dispDiff(), Repeated );

                this->moveMesh(meshDispDiff);

        this->interpolateVelocity(meshDispDiff, *M_beta);

        *M_beta *= -1./M_dataFluid->timestep();

        *M_beta  += *this->M_un;

                double alpha = 1./M_dataFluid->timestep();
                if(recomputeMatrices)
                    {
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
    this->leaderPrintMax( "Fluid solution: total time : " , chronoFluid.diff() );

    if ( false && this->isFluid() )
    {
        vector_type vel  (this->fluid().velFESpace().map());
        vector_type press(this->fluid().pressFESpace().map());

        vel.subset(this->M_fluid->solution());
        press.subset(this->M_fluid->solution(), this->fluid().velFESpace().dim()*this->fluid().pressFESpace().fieldDim());

        std::cout << "norm_inf( vel ) " << vel.NormInf() << std::endl;
        std::cout << "norm_inf( press ) " << press.NormInf() << std::endl;

        this->M_fluid->postProcess();
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
        this->M_solid->iterate( *M_BCh_d );
        this->transferSolidOnInterface(this->M_solid->disp(),     lambdaSolidUnique);
        this->transferSolidOnInterface(this->M_solid->vel(),      lambdaDotSolidUnique);
        this->transferSolidOnInterface(this->M_solid->residual(), sigmaSolidUnique);
    }

    M_epetraWorldComm->Barrier();
    chronoSolid.stop();
    this->leaderPrintMax( "Solid solution: total time : " , chronoSolid.diff() );
    chronoInterface.start();

    this->setLambdaSolid(    lambdaSolidUnique);
    this->setLambdaDotSolid( lambdaDotSolidUnique);
    this->setSigmaSolid(     sigmaSolidUnique);

    chronoInterface.stop();
    this->leaderPrintMax( "Interface transfer: total time : " , chronoInterface.diff_cumul() );

    if ( false && this->isSolid() )
    {
        this->solid().postProcess();
    }

// possibly unsafe when using more cpus, since both has repeated maps



    this->leaderPrint( " ::: norm(disp     )    = ", _disp.NormInf());
    this->leaderPrint( " ::: norm(dispNew  )    = " , this->lambdaSolid().NormInf() );
    this->leaderPrint( " ::: norm(velo     )    = " , this->lambdaDotSolid().NormInf() );
    this->leaderPrint( " ::: max Residual Fluid = " , this->sigmaFluid().NormInf() );
    this->leaderPrint( " ::: max Residual Solid = " , this->sigmaSolid().NormInf() );

    if (this->isFluid())
        this->leaderPrint( "Max ResidualF        = " , M_fluid->residual().NormInf() );
    if (this->isSolid())
        {
            this->leaderPrint( "NL2 DiplacementS     = " , M_solid->disp().Norm2() );
            this->leaderPrint( "Max ResidualS        = " , M_solid->residual().NormInf() );
        }
}

void exactJacobian::evalResidual(vector_type&       res,
                                 const vector_type& disp,
                                 const int          iter)
{
    if (this->isSolid())
    {
        std::cout << "*** Residual computation g(x_" << iter <<" )";
        if (iter == 0) std::cout << " [NEW TIME STEP] ";
        std::cout << std::endl;
    }

    this->setLambdaSolidOld(disp);


    eval(disp, iter);

    res  = this->lambdaSolid();
    res -=  disp;

    this->leaderPrint( "NormInf res   " , res.NormInf() );
    this->leaderPrint( "NormInf res_d " , this->solid().residual().NormInf() );

}


//
// new step computation resolution
//


void  exactJacobian::solveJac(vector_type         &_muk,
                              const vector_type   &_res,
                              const double         _linearRelTol)
{
    if (this->isFluid() && this->isLeader()) std::cout << "  f- ";
    if (this->isSolid() && this->isLeader()) std::cout << "  s- ";

    this->leaderPrint( "solveJac: NormInf res " , _res.NormInf() );
    _muk *= 0.;

    M_linearSolver.setTolMaxiter(_linearRelTol, 100);

    vector_type res(_res);

    M_linearSolver.setOperator(*M_epetraOper);

    this->leaderPrint( "Solving Jacobian system... " );
    M_linearSolver.solve(_muk, res);

    this->leaderPrint( "done.\n" );
}


void  exactJacobian::solveLinearFluid()
{

    double alpha = this->M_bdf->coeff_der( 0 ) / M_dataFluid->timestep();

    vector_type dispFluidMesh(this->derVeloFluidMesh().getMap(), Repeated);
//if statement: in order not to iterate the mesh for each linear residual calculation, needed just for exact Jac case.
    if(false && this->M_dataFluid->semiImplicit()==1)// not working in parallel
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

    this->derVeloFluidMesh() = dispFluidMesh;

    this->derVeloFluidMesh() *= 1./(M_dataFluid->timestep());

    this->leaderPrint( " norm inf dw = " , this->derVeloFluidMesh().NormInf() );

    *M_rhsNew *= 0.;

    this->M_fluid->updateLinearSystem( M_fluid->matrNoBC(),
                                       alpha,
                                       *M_un,
                                       M_fluid->solution(),
                                       dispFluidMesh,
                                       this->veloFluidMesh(),
                                       this->derVeloFluidMesh(),
                                       *M_rhsNew );

    this->M_fluid->iterateLin( *this->M_BCh_du );
}


//


void  exactJacobian::solveLinearSolid()
{
    this->M_solid->iterateLin( *M_BCh_dz );
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

        M_ej->leaderPrint( "NormInf res   " , z.NormInf() );

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
                if(true || ( !this->M_ej->dataFluid().semiImplicit() || this->M_ej->dataFluid().semiImplicit()==-1))
                    M_ej->meshMotion().iterate();

                M_ej->leaderPrint( " norm inf dx = " , M_ej->meshMotion().disp().NormInf() );

                M_ej->solveLinearFluid();

                M_ej->transferFluidOnInterface(M_ej->fluid().residual(), sigmaFluidUnique);

                M_ej->fluidPostProcess();
            }

        M_comm->Barrier();
        chronoFluid.stop();
        M_ej->leaderPrintMax( "Fluid linear solution: total time : ", chronoFluid.diff() );


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
        M_ej->leaderPrintMax( "Solid linear solution: total time : " , chronoSolid.diff() );

        chronoInterface.start();
        M_ej->setLambdaSolid(lambdaSolidUnique);

        chronoInterface.stop();
        M_ej->leaderPrintMax( "Interface linear transfer: total time : " , chronoInterface.diff_cumul() );

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


namespace
{
FSIOperator* createEJ(){ return new exactJacobian(); }
static bool reg = FSIFactory::instance().registerProduct( "exactJacobian", &createEJ );
}

}
