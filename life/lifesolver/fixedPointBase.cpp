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

#include <life/lifesolver/fixedPointBase.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

fixedPoint::fixedPoint():
        super(),
        M_aitkFS(),
        M_rhsNew(),
        M_beta()
{
}


fixedPoint::~fixedPoint()
{}


// ===================================================
// Methods
// ===================================================

void  fixedPoint::solveJac(vector_Type        &muk,
                           const vector_Type  &res,
                           const Real   /*_linearRelTol*/)
{
    M_aitkFS.restart();

    if (M_data->algorithm()=="RobinNeumann")
    {
        muk = M_aitkFS.computeDeltaLambdaScalar(this->lambdaSolidOld(), res);
    }
    else
    {
        muk = M_aitkFS.computeDeltaLambdaScalar(this->lambdaSolidOld(), -1.*res);
    }
}

void fixedPoint::evalResidual(vector_Type &res, const vector_Type& disp, UInt iter)
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
}


void
fixedPoint::setupFEspace()
{
    super::setLinearFluid(false);
    super::setLinearSolid(false);

    super::setupFEspace();
}



void
fixedPoint::setupFluidSolid()
{
    // call FSIOperator setup()

    Debug( 6205 ) << "Setting up the FSI problem \n";

    super::setLinearFluid(false);
    super::setLinearSolid(false);

    super::setupFluidSolid();

    if ( this->isFluid() )
    {
        M_rhsNew.reset(new vector_Type(this->M_fluid->getMap()));
        M_beta.reset(new vector_Type(this->M_fluid->getMap()));
    }

}

void
fixedPoint::setDataFile( GetPot const& dataFile )
{
    super::setDataFile( dataFile );

    M_aitkFS.setDefaultOmega(M_data->defaultOmega(), 0.001);
    M_aitkFS.setOmegaRange( M_data->OmegaRange() );

    if ( M_data->algorithm() == "RobinNeumann" )
        M_aitkFS.setDefaultOmega(-1, 1);

}

// ===================================================
// Private Methods
// ===================================================

void fixedPoint::eval( const vector_Type& _disp,
                       UInt                iter)
{
    // If M_data->updateEvery() == 1, normal fixedPoint algorithm
    // If M_data->updateEvery()  > 1, recompute computational domain every updateEvery iterations (transpiration)
    // If M_data->updateEvery() <= 0, recompute computational domain and matrices only at first subiteration (semi-implicit)
    bool recomputeMatrices ( iter == 0 || ( M_data->updateEvery() > 0 && iter % M_data->updateEvery() == 0 ) );

    std::cout << "recomputeMatrices = " << recomputeMatrices << " == " << true
              << "; iter = " << iter
              << "; updateEvery = " << M_data->updateEvery() << std::endl;

    if (iter == 0 && this->isFluid())
    {
        this->M_fluid->resetPreconditioner();
        //this->M_solid->resetPrec();
    }



    MPI_Barrier(MPI_COMM_WORLD);

    // possibly unsafe when using more cpus,
    this->setLambdaFluid(_disp);

    //Change in sign in the residual for RN:
    // if(this->algorithm()=="RobinNeumann")   this->setMinusSigmaFluid( sigmaSolidUnique );
    if (M_data->algorithm()=="RobinNeumann")   this->setMinusSigmaFluid( this->sigmaSolid() );


    vector_Type sigmaFluidUnique (this->sigmaFluid().map(), Unique);

    if (this->isFluid())
    {
        this->M_meshMotion->iterate(*M_BCh_mesh);

        this->transferMeshMotionOnFluid(M_meshMotion->disp(),
                                        this->veloFluidMesh());

        this->veloFluidMesh()    -= this->dispFluidMeshOld();
        this->veloFluidMesh()    *= 1./(M_data->dataFluid()->dataTime()->getTimeStep());

        // copying displacement to a repeated indeces displacement, otherwise the mesh wont know
        // the value of the displacement for some points

        vector_Type const meshDisplacement( M_meshMotion->disp(), Repeated );
        this->moveMesh(meshDisplacement);
        /*
        vector_Type const meshDisplacement( M_meshMotion->dispDiff(), Repeated );
        this->moveMesh(meshDispDiff);
        */

        *M_beta *= 0.;

        this->transferMeshMotionOnFluid(M_meshMotion->dispDiff(),  *M_beta);

        *M_beta *= -1./M_data->dataFluid()->dataTime()->getTimeStep();

        *M_beta += this->M_bdf->extrap();

        double alpha = this->M_bdf->coeff_der( 0 ) / M_data->dataFluid()->dataTime()->getTimeStep();

        //*M_rhsNew   = *this->M_rhs;
        //*M_rhsNew  *= alpha;


        if (recomputeMatrices)
        {

            this->M_fluid->updateSystem( alpha, *M_beta, *M_rhs );
        }
        else
        {
            this->M_fluid->updateRightHandSide( *M_rhs );
        }

        //	if(this->algorithm()=="RobinNeumann") this->updatealphaf(this->veloFluidMesh());// this->setAlphaf();

        this->M_fluid->iterate( *M_BCh_u );

        std::cout << "Finished iterate, transfering to interface \n" << std::flush;

        this->transferFluidOnInterface(this->M_fluid->residual(), sigmaFluidUnique);

    }

    MPI_Barrier(MPI_COMM_WORLD);

    if ( true && this->isFluid() )
    {
        vector_Type vel  (this->fluid().velocityFESpace().map());
        vector_Type press(this->fluid().pressureFESpace().map());

        vel.subset(*this->M_fluid->solution());
        press.subset(*this->M_fluid->solution(), this->fluid().velocityFESpace().dim()*this->fluid().pressureFESpace().fieldDim());

        std::cout << "norm_inf( vel ) " << vel.normInf() << std::endl;
        std::cout << "norm_inf( press ) " << press.normInf() << std::endl;

    }



    this->setSigmaFluid( sigmaFluidUnique );
    this->setSigmaSolid( sigmaFluidUnique );

    MPI_Barrier(MPI_COMM_WORLD);


    vector_Type lambdaSolidUnique   (this->lambdaSolid().map(),    Unique);
    vector_Type lambdaDotSolidUnique(this->lambdaDotSolid().map(), Unique);
    vector_Type sigmaSolidUnique    (this->sigmaSolid().map(),     Unique);

    if (this->isSolid())
    {
        this->M_solid->iterate( M_BCh_d );
        this->transferSolidOnInterface(this->M_solid->disp(),     lambdaSolidUnique);
        this->transferSolidOnInterface(this->M_solid->vel(),      lambdaDotSolidUnique);
        this->transferSolidOnInterface(this->M_solid->residual(), sigmaSolidUnique);
    }

    this->setLambdaSolid( lambdaSolidUnique );
    this->setLambdaDotSolid( lambdaDotSolidUnique );
    this->setSigmaSolid( sigmaSolidUnique );

    MPI_Barrier(MPI_COMM_WORLD);


    // Some displays:
    Real norm;

    norm = _disp.normInf();
    if (this->isSolid() && this->isLeader())
        std::cout << " ::: norm(disp     )    = " << norm << std::endl;

    norm = this->lambdaSolid().normInf();
    if (this->isSolid() && this->isLeader())
        std::cout << " ::: norm(dispNew  )    = " << norm << std::endl;

    norm = this->lambdaDotSolid().normInf();
    if (this->isSolid() && this->isLeader())
        std::cout << " ::: norm(velo     )    = " << norm << std::endl;

    norm = this->sigmaFluid().normInf();
    if (this->isSolid() && this->isLeader())
        std::cout << " ::: max Residual Fluid = " << norm << std::endl;

    norm = this->sigmaSolid().normInf();
    if (this->isSolid() && this->isLeader())
        std::cout << " ::: max Residual Solid = " << norm  << std::endl;

    if (this->isFluid())
    {
        norm = M_fluid->residual().normInf();
        if (this->isLeader())
            std::cout << "Max ResidualF        = " << norm << std::endl;
    }
    if (this->isSolid())
    {
        norm = M_solid->disp().norm2();
        if (this->isLeader())
            std::cout << "NL2 DiplacementS     = " << norm << std::endl;

        norm = M_solid->residual().normInf();
        if (this->isLeader())
            std::cout << "Max ResidualS        = " << norm << std::endl;
    }

}




void fixedPoint::registerMyProducts( )
{
    FSIFactory_Type::instance().registerProduct( "fixedPoint", &createFP );
    solid_Type::StructureSolverFactory::instance().registerProduct( "LinearVenantKirchhof", &FSIOperator::createLinearStructure );
    //solid_Type::StructureSolverFactory::instance().registerProduct( "NonLinearVenantKirchhof", &FSIOperator::createNonLinearStructure );
}

}   // Namespace LifeV
