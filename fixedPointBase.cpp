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


#include <life/lifesolver/fixedPointBase.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>


namespace LifeV
{
fixedPoint::fixedPoint():
    super(),
    M_defOmega( 0.001 ),
    M_updateEvery(1),
    M_aitkFS(),
    //    M_displacement(),
    //    M_stress(),
    //    M_velocity(),
    //    M_residualFSI(),
    M_rhsNew(),
    M_beta()
{
    M_aitkFS.setDefault( M_defOmega );
}


fixedPoint::~fixedPoint()
{}


void
fixedPoint::setDataFromGetPot( GetPot const& data )
{
    // call the super class to setup the data from getpot file if needed
    super::setDataFromGetPot( data );

    M_defOmega = data("problem/defOmega", 0.001);

    Debug( 6205 ) << "fixedPoint::setDataFromGetPot(GetPot) OmegaS = " << M_defOmega << "\n";

    M_aitkFS.setDefault(M_defOmega, 0.001);

    // 0-order transpiration bnd conditions
    // If M_updateEvery == 1, normal fixedPoint algorithm
    // If M_updateEvery  > 1, recompute computational domain every M_updateEvery iterations (transpiration)
    // If M_updateEvery <= 0, recompute computational domain and matrices only at first subiteration (semi-implicit)
    M_updateEvery = data("problem/updateEvery", 1);

    if(this->algorithm()=="RobinNeumann")  M_aitkFS.setDefault(-1, 1); //RN

    if (false && this->isFluid())
        M_ensightFluid.reset( new  Ensight<RegionMesh3D<LinearTetra> > ( data, "fixedPtFluidInnerIterations") );

}

void
fixedPoint::setup()
{
    // call FSIOperator setup()

    Debug( 6205 ) << "Setting up the FSI problem \n";

    super::setLinearFluid(false);
    super::setLinearSolid(false);

    super::setup();

    if ( this->isFluid() )
        {
            M_rhsNew.reset(new vector_type(this->M_fluid->getMap()));
            M_beta.reset(new vector_type(this->M_fluid->getMap()));
        }

    if ( false && this->isFluid())
        {
            if (M_ensightFluid.get() == 0)
                M_ensightFluid.reset( new  Ensight<RegionMesh3D<LinearTetra> > ( "", "fixedPtFluidInnerIterations") );

            assert( M_uFESpace.get() );
            assert( M_uFESpace->mesh().get() );

            M_ensightFluid->setMeshProcId(M_uFESpace->mesh(), M_uFESpace->map().Comm().MyPID());

            assert( this->M_fluid.get() );
            M_velAndPressure.reset( new vector_type( this->M_fluid->getMap() , Repeated ));

            M_ensightFluid->addVariable( ExporterData::Vector, "velocityInner", M_velAndPressure,
                                         UInt(0), M_uFESpace->dof().numTotalDof() );

            M_ensightFluid->addVariable( ExporterData::Scalar, "pressureInner", M_velAndPressure,
                                         UInt(3*M_uFESpace->dof().numTotalDof()),
                                         UInt(3*M_uFESpace->dof().numTotalDof()+M_pFESpace->dof().numTotalDof()) );

        }

//@    setUpBC();
}

void fixedPoint::eval( const vector_type& _disp,
                       int                iter)
{
    // If M_updateEvery == 1, normal fixedPoint algorithm
    // If M_updateEvery  > 1, recompute computational domain every M_updateEvery iterations (transpiration)
    // If M_updateEvery <= 0, recompute computational domain and matrices only at first subiteration (semi-implicit)
    bool recomputeMatrices ( iter == 0 || ( M_updateEvery > 0 && iter % M_updateEvery == 0 ) );

    std::cout << "recomputeMatrices = " << recomputeMatrices << " == " << true
              << "; iter = " << iter
              << "; M_updateEvery = " << M_updateEvery << std::endl;

    if (iter == 0)
        {
            M_nbEval = 0; // new time step
            this->M_fluid->resetPrec();
            //this->M_solid->resetPrec();
        }
    else
        {
            if(!recomputeMatrices)
                this->M_fluid->resetPrec(false);
        }

    M_nbEval++ ;


    MPI_Barrier(MPI_COMM_WORLD);

    // possibly unsafe when using more cpus,
    this->setLambdaFluid(_disp);

    //Change in sign in the residual for RN:
    // if(this->algorithm()=="RobinNeumann")   this->setMinusSigmaFluid( sigmaSolidUnique );
    if(this->algorithm()=="RobinNeumann")   this->setMinusSigmaFluid( this->sigmaSolid() );


    vector_type sigmaFluidUnique (this->sigmaFluid().getMap(), Unique);

    if (this->isFluid())
    {
        this->M_meshMotion->iterate();

        this->transferMeshMotionOnFluid(M_meshMotion->disp(),
                                        this->veloFluidMesh());

        this->veloFluidMesh()    -= this->dispFluidMeshOld();
        this->veloFluidMesh()    *= 1./(M_dataFluid->timestep());

        // copying displacement to a repeated indeces displacement, otherwise the mesh wont know
        // the value of the displacement for some points

        vector_type const meshDisplacement( M_meshMotion->disp(), Repeated );
        this->moveMesh(meshDisplacement);
        /*
        vector_type const meshDisplacement( M_meshMotion->dispDiff(), Repeated );
        this->moveMesh(meshDispDiff);
        */

        *M_beta *= 0.;

        this->transferMeshMotionOnFluid(M_meshMotion->dispDiff(),  *M_beta);

        *M_beta *= -1./M_dataFluid->timestep();

        *M_beta += this->M_bdf->extrap();

        double alpha = this->M_bdf->coeff_der( 0 ) / M_dataFluid->timestep();

        //*M_rhsNew   = *this->M_rhs;
        //*M_rhsNew  *= alpha;


        if (recomputeMatrices)
            {

                this->M_fluid->updateSystem( alpha, *M_beta, *M_rhs );
            }
        else
            {
                this->M_fluid->updateRHS( *M_rhs );
            }

        //	if(this->algorithm()=="RobinNeumann") this->updatealphaf(this->veloFluidMesh());// this->setAlphaf();

        this->M_fluid->iterate( *M_BCh_u );

        std::cout << "Finished iterate, transfering to interface \n" << std::flush;

        this->transferFluidOnInterface(this->M_fluid->residual(), sigmaFluidUnique);

    }

    MPI_Barrier(MPI_COMM_WORLD);

    if ( true && this->isFluid() )
        {
            vector_type vel  (this->fluid().velFESpace().map());
            vector_type press(this->fluid().pressFESpace().map());

            vel.subset(this->M_fluid->solution());
            press.subset(this->M_fluid->solution(), this->fluid().velFESpace().dim()*this->fluid().pressFESpace().fieldDim());

            std::cout << "norm_inf( vel ) " << vel.NormInf() << std::endl;
            std::cout << "norm_inf( press ) " << press.NormInf() << std::endl;

            //              this->M_fluid->postProcess();
            //              *M_velAndPressure = this->M_fluid->solution();
            //              M_ensightFluid->postProcess( M_nbEval );


        }



     this->setSigmaFluid( sigmaFluidUnique );
     this->setSigmaSolid( sigmaFluidUnique );
    //*this->M_sigmaSolid = -1.0 * *this->M_sigmaFluid;

//    this->M_sigmaSolid->spy("sigmasolid");
    MPI_Barrier(MPI_COMM_WORLD);


    vector_type lambdaSolidUnique   (this->lambdaSolid().getMap(),    Unique);
    vector_type lambdaDotSolidUnique(this->lambdaDotSolid().getMap(), Unique);
    vector_type sigmaSolidUnique    (this->sigmaSolid().getMap(),     Unique);

    if (this->isSolid())
    {
        this->M_solid->iterate( *M_BCh_d );
        this->transferSolidOnInterface(this->M_solid->disp(),     lambdaSolidUnique);
        this->transferSolidOnInterface(this->M_solid->vel(),      lambdaDotSolidUnique);
        this->transferSolidOnInterface(this->M_solid->residual(), sigmaSolidUnique);
    }

    this->setLambdaSolid( lambdaSolidUnique );
    this->setLambdaDotSolid( lambdaDotSolidUnique );
    this->setSigmaSolid( sigmaSolidUnique );

//     dispNew = *this->M_lambdaSolid;
//     velo    = *this->M_lambdaDotSolid;

    MPI_Barrier(MPI_COMM_WORLD);


    // Some displays:
    Real norm;

    norm = _disp.NormInf();
    if (this->isSolid() && this->isLeader())
        std::cout << " ::: norm(disp     )    = " << norm << std::endl;

    norm = this->lambdaSolid().NormInf();
    if (this->isSolid() && this->isLeader())
        std::cout << " ::: norm(dispNew  )    = " << norm << std::endl;

    norm = this->lambdaDotSolid().NormInf();
    if (this->isSolid() && this->isLeader())
        std::cout << " ::: norm(velo     )    = " << norm << std::endl;

    norm = this->sigmaFluid().NormInf();
    if (this->isSolid() && this->isLeader())
        std::cout << " ::: max Residual Fluid = " << norm << std::endl;

    norm = this->sigmaSolid().NormInf();
    if (this->isSolid() && this->isLeader())
        std::cout << " ::: max Residual Solid = " << norm  << std::endl;

    if (this->isFluid())
    {
        norm = M_fluid->residual().NormInf();
        if (this->isLeader())
            std::cout << "Max ResidualF        = " << norm << std::endl;
    }
    if (this->isSolid())
    {
        norm = M_solid->disp().Norm2();
        if (this->isLeader())
            std::cout << "NL2 DiplacementS     = " << norm << std::endl;

        norm = M_solid->residual().NormInf();
        if (this->isLeader())
            std::cout << "Max ResidualS        = " << norm << std::endl;
    }

}

// Residual evaluation
//
void fixedPoint::evalResidual(vector_type &res, const vector_type& disp, int iter)
{


    if (this->isSolid())
    {
        std::cout << "*** Residual computation g(x_" << iter <<" )";
        if (iter == 0) std::cout << " [NEW TIME STEP] ";
        std::cout << std::endl;
    }

    this->setLambdaSolidOld(disp);

    eval(disp, iter);

    res  = disp;
    res -= this->lambdaSolid();
//     res  = this->lambdaSolid();
//     res -= disp;
}


//
// new step computation resolution
//


void  fixedPoint::solveJac(vector_type        &muk,
                           const vector_type  &res,
                           const double   /*_linearRelTol*/)
{
    if (M_nbEval == 1) M_aitkFS.restart();

    if(this->algorithm()=="RobinNeumann")
    {
        muk = M_aitkFS.computeDeltaLambda(this->lambdaSolidOld(), res);
    }
    else
    {
        muk = M_aitkFS.computeDeltaLambda(this->lambdaSolidOld(), -1.*res);
    }
}


//
// add fixedPoint to factory
//


// namespace
// {
// FSIOperator* createFP(){ return new fixedPoint(); }
// bool reg = FSIFactory::instance().registerProduct( "fixedPoint", &createFP );
// }

}
