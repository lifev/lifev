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
            M_velAndPressure.reset( new vector_type( this->M_fluid->getRepeatedEpetraMap() ));

            M_ensightFluid->addVariable( ExporterData::Vector, "velocityInner", M_velAndPressure,
                                         UInt(0), M_uFESpace->dof().numTotalDof() );

            M_ensightFluid->addVariable( ExporterData::Scalar, "pressureInner", M_velAndPressure,
                                         UInt(3*M_uFESpace->dof().numTotalDof()),
                                         UInt(3*M_uFESpace->dof().numTotalDof()+M_pFESpace->dof().numTotalDof()) );

        }

//@    setUpBC();
}

void fixedPoint::eval(vector_type&       dispNew,
                      vector_type&       velo,
                      const vector_type& disp,
                      int                status)
{
    if(status)
        {
            M_nbEval = 0; // new time step
            this->M_fluid->resetPrec();
            this->M_solid->resetPrec();
        }

    M_nbEval++ ;


    MPI_Barrier(MPI_COMM_WORLD);

    // possibly unsafe when using more cpus, since both has repeated maps
    *M_lambdaFluid       = disp;


    *this->M_sigmaFluid *= 0.;

    if (this->isFluid())
    {
        this->M_meshMotion->iterate();

        vector_type const meshDisplacement( M_meshMotion->displacement(), this->M_meshMotion->getRepeatedEpetraMap() );

        this->transferMeshMotionOnFluid(meshDisplacement,
                                        *this->M_veloFluidMesh);

        *this->M_veloFluidMesh    -= *M_dispFluidMeshOld;
        *this->M_veloFluidMesh    *= 1./(M_dataFluid->timestep());

        // copying displacement to a repeated indeces displacement, otherwise the mesh wont know
        // the value of the displacement for some points

        this->moveMesh(meshDisplacement);

        *M_beta *= 0.;

        vector_type const meshDispDiff( M_meshMotion->dispDiff(), this->M_meshMotion->getRepeatedEpetraMap() );

        this->interpolateVelocity(meshDispDiff, *M_beta);

        *M_beta *= -1./M_dataFluid->timestep();

        *M_beta += this->M_bdf->extrap();

        double alpha = this->M_bdf->coeff_der( 0 ) / M_dataFluid->timestep();

        *M_rhsNew   = *this->M_rhs;
        *M_rhsNew  *= alpha;

        this->M_fluid->updateSystem( alpha, *M_beta, *M_rhsNew );

        this->M_fluid->iterate( *M_BCh_u );

        this->transferFluidOnInterface(this->M_fluid->residual(), *this->M_sigmaFluid);

//        this->shiftSolution();
    }

    MPI_Barrier(MPI_COMM_WORLD);

     if ( false && this->isFluid() )
         {
             *M_velAndPressure = this->M_fluid->solution();
             M_ensightFluid->postProcess( M_nbEval );
         }



    // possibly unsafe when using more cpus, since both has repeated maps
    *this->M_sigmaSolid = *this->M_sigmaFluid;
    //*this->M_sigmaSolid = -1.0 * *this->M_sigmaFluid;

//    this->M_sigmaSolid->spy("sigmasolid");
    MPI_Barrier(MPI_COMM_WORLD);


    if (this->isSolid())
    {
        this->M_solid->iterate( *M_BCh_d );
        this->transferSolidOnInterface(this->M_solid->disp(), *this->M_lambdaSolid);
        this->transferSolidOnInterface(this->M_solid->vel() , *this->M_lambdaDotSolid);
        this->transferSolidOnInterface(this->M_solid->residual() , *this->M_sigmaSolid);
    }

    // possibly unsafe when using more cpus, since both has repeated maps
//    *this->M_lambdaFluid    = *this->M_lambdaSolid;

//    this->M_lambdaFluid->getEpetraVector().Print(std::cout);

//    *this->M_lambdaDotFluid = *this->M_lambdaDotSolid;


    dispNew = *this->M_lambdaSolid;
    velo    = *this->M_lambdaDotSolid;

    MPI_Barrier(MPI_COMM_WORLD);


//      if (this->isSolid())
//      {
    std::cout << " ::: norm(disp     )    = " << disp.NormInf() << std::endl;
    std::cout << " ::: norm(dispNew  )    = " << dispNew.NormInf() << std::endl;
    std::cout << " ::: norm(velo     )    = " << velo.NormInf() << std::endl;
//      }
    std::cout << " ::: max Residual Fluid = " << M_sigmaFluid->NormInf() << std::endl;
    std::cout << " ::: max Residual Solid = " << M_sigmaSolid->NormInf() << std::endl;

    if (this->isFluid())
        std::cout << "Max ResidualF        = " << M_fluid->residual().NormInf() << std::endl;
    if (this->isSolid())
        {
            std::cout << "NL2 DiplacementS     = " << M_solid->disp().Norm2() << std::endl;
            std::cout << "Max ResidualS        = " << M_solid->residual().NormInf() << std::endl;
        }

}

// Residual evaluation
//
void fixedPoint::evalResidual(vector_type &res, const vector_type& disp, int iter)
{
    int status = 0;
    if(iter == 0) status = 1;


    if (this->isSolid())
    {
        std::cout << "*** Residual computation g(x_" << iter <<" )";
        if (status) std::cout << " [NEW TIME STEP] ";
        std::cout << std::endl;
    }

    M_lambdaSolidOld.reset(new vector_type(disp));

    eval(*M_lambdaSolid, *M_lambdaDotSolid, disp, status);

    res  = *M_lambdaSolid;
    res -= disp;

//    res.getEpetraVector().Print(std::cout);
//     transferOnInterface(res,
//                         M_solid->BC_solid(),
//                         "Interface",
//                         res);
}


//
// new step computation resolution
//


void  fixedPoint::solveJac(vector_type        &_muk,
                           const vector_type  &_res,
                           const double   /*_linearRelTol*/)
{
    if (M_nbEval == 1) M_aitkFS.restart();
    //    _muk = 0.001*_res;
    _muk = M_aitkFS.computeDeltaLambda(*M_lambdaSolidOld, -1.*_res);
}


//
// add fixedPoint to factory
//


namespace
{
FSIOperator* createFP(){ return new fixedPoint(); }
}
bool LifeV::fixedPoint::reg = FSIFactory::instance().registerProduct( "fixedPoint", &createFP );

}
