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

/*!
 *  @file
 *  @brief This file contains solvers for St. Venant-Kirchhof materials (linear for the moment)
 *
 *  @version 1.0
 *  @date 01-06-2003
 *  @author Miguel Angel Fernandez
 *
 *  @version 1.1
 *  @date 01-09-2010
 *  @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#incude<life/lifesolver/VenantKirchhofSolver.hpp>

namespace LifeV
{


//====================================
// Constructor
//=====================================

template <typename Mesh, typename SolverType>
VenantKirchhofSolver<Mesh, SolverType>::VenantKirchhofSolver( ):
        M_data                       ( ),
        M_FESpace                    ( ),
        M_Displayer                  ( ),
        M_me                         ( 0 ),
        M_linearSolver               ( ),
        M_elmatK                     ( ),
        M_elmatM                     ( ),
        M_elmatC                     ( ),
        M_disp                       ( ),
        M_vel                        ( ),
        M_rhs                        ( /*new vector_type(M_localMap)*/),//useful
        M_rhsW                       ( ),
        M_rhsNoBC                    ( ),
        M_f                          ( ),//useless
        M_residual_d                 ( ),//useless
        M_sxx                        (/*M_localMap*/),//useless
        M_syy                        (/*M_localMap*/),//useless
        M_szz                        (/*M_localMap*/),//useless
        M_out_iter                   ( "out_iter_solid" ),
        M_out_res                    ( "out_res_solid" ),
        M_BCh                        (),
        M_localMap                   ( ),
        M_mass                       ( ),
        M_massStiff                  ( /*new matrix_type(monolithicMap) */),//constructed outside
        M_stiff                      ( /*new matrix_type(M_localMap)*/ ),
        M_linearStiff                ( ),
        M_jacobian                   (/*M_localMap*/),
        _recur                       (),
        M_source                     (),
        M_count                      ( 0 ),//useless
        M_offset                     ( 0 ),
        M_rescaleFactor              ( 1. ),
        M_matrFull                   ( )//useless
{
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::setup(boost::shared_ptr<data_type>          data,
					      const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
					      bchandler_type&                BCh,
					      boost::shared_ptr<Epetra_Comm>&              comm)
{
    setup(data, dFESpace, comm);
    M_BCh = BCh;
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::setup(boost::shared_ptr<data_type>        data,
					      const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
					      boost::shared_ptr<Epetra_Comm>&     comm)
{
    setup( data, dFESpace, comm, dFESpace->mapPtr(), (UInt)0 );
    M_stiff.reset                      ( new matrix_type(*M_localMap) );
    M_massStiff.reset                  ( new matrix_type(*M_localMap) );
    M_jacobian.reset                   ( new matrix_type(*M_localMap) );
    M_rhs.reset                        ( new vector_type(*M_localMap));
    M_f.reset                          ( new vector_type(*M_localMap));
    M_residual_d.reset                 ( new vector_type(*M_localMap));
    M_sxx.reset                        ( new vector_type(*M_localMap) );
    M_syy.reset                        ( new vector_type(*M_localMap) );
    M_szz.reset                        ( new vector_type(*M_localMap) );
    M_linearSolver.reset               ( new SolverType( comm ) );

}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::setup(boost::shared_ptr<data_type>        data,
					      const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
					      boost::shared_ptr<Epetra_Comm>&     comm,
					      const boost::shared_ptr<const EpetraMap>&  monolithicMap,
					      UInt                                offset)
{
    M_data                            = data;
    M_FESpace                         = dFESpace;
    M_Displayer.reset                 (new Displayer(comm));
    M_me                              = comm->MyPID();
    M_elmatK.reset                    ( new ElemMat( M_FESpace->fe().nbNode, nDimensions, nDimensions ) );
    M_elmatM.reset                    ( new ElemMat( M_FESpace->fe().nbNode, nDimensions, nDimensions ) );
    M_elmatC.reset                    ( new ElemMat( M_FESpace->fe().nbNode, nDimensions, nDimensions ) );
    M_localMap                        = monolithicMap;
    M_disp.reset                      (new vector_type(*M_localMap));
    M_vel.reset                       (new vector_type(*M_localMap));
    M_rhsW.reset                      ( new vector_type(*M_localMap) );
    M_rhsNoBC.reset                   ( new vector_type(*M_localMap) );
    M_mass.reset                      (new matrix_type(*M_localMap));
    M_linearStiff.reset               (new matrix_type(*M_localMap));
    M_offset                          = offset;
}

template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::updateSystem(  )
{
    updateSystem(M_stiff);
}

template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::updateSystem( matrix_ptrtype& /*stiff*/ )
{
    M_Displayer->leaderPrint("  S-  Updating mass term on right hand side... ");

    Chrono chrono;
    chrono.start();

    //M_stiff = M_linearStiff;

    // Number of displacement components
//    UInt nc = nDimensions;

    Real coef;

    coef = (Real) M_data->dataTime()->getTimeStep();

    vector_type _z = *M_disp;
    _z            += coef*(*M_vel);

    std::cout<< "M_disp in solid" << M_disp->Norm2()<<std::endl;

    *M_rhsNoBC  = *M_mass*_z;

    std::cout<< "rhsNoBC in solid 1" << M_rhsNoBC->Norm2()<<std::endl;

    *M_rhsNoBC -= *M_linearStiff*(*M_disp);

    std::cout<< "rhsNoBC in solid 2" << M_rhsNoBC->Norm2()<<std::endl;

    coef = 2.0/M_data->dataTime()->getTimeStep();

    *M_rhsW  = coef*(*M_disp);
    *M_rhsW += *M_vel;

    chrono.stop();

    M_Displayer->leaderPrintMax("done in ", chrono.diff());
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::buildSystem()
{
    M_Displayer->leaderPrint("  S-  Computing constant matrices ...          ");
    Chrono chrono;
    chrono.start();

    M_massStiff.reset(new matrix_type(*M_localMap));
    buildSystem(M_massStiff);
    M_massStiff->GlobalAssemble();

    chrono.stop();
    M_Displayer->leaderPrintMax( "done in ", chrono.diff() );
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::buildSystem(matrix_ptrtype massStiff, const Real& factor)
{
    UInt totalDof = M_FESpace->dof().numTotalDof();

    // Number of displacement components
    UInt nc = nDimensions;

    //inverse of dt:
    Real dti2 = 2.0 / ( M_data->dataTime()->getTimeStep() * M_data->dataTime()->getTimeStep() );

    // Elementary computation and matrix assembling
    // Loop on elements
    for ( UInt i = 1; i <= M_FESpace->mesh()->numVolumes(); i++ )
    {

        M_FESpace->fe().updateFirstDerivQuadPt( M_FESpace->mesh()->volumeList( i ) );

        M_elmatK->zero();
        M_elmatM->zero();


        // stiffness
        UInt marker = M_FESpace->mesh()->volumeList( i ).marker();
        stiff_strain(     M_data->mu(marker), *M_elmatK, M_FESpace->fe() );
        stiff_div   ( 0.5*M_data->lambda(marker), *M_elmatK, M_FESpace->fe() );

        M_elmatC->mat() = M_elmatK->mat();

        // mass
        mass( dti2 * M_data->rho(), *M_elmatM, M_FESpace->fe(), 0, 0, nDimensions );

        M_elmatC->mat() += M_elmatM->mat();

        // assembling
        for ( UInt ic = 0; ic < nc; ic++ )
        {
            for ( UInt jc = 0; jc < nc; jc++ )
            {
                assembleMatrix( *M_linearStiff, *M_elmatK, M_FESpace->fe(), M_FESpace->fe(), M_FESpace->dof(), M_FESpace->dof(),  ic,  jc,  M_offset +ic*totalDof, M_offset + jc*totalDof );

                assembleMatrix( *massStiff  , *M_elmatC, M_FESpace->fe(), M_FESpace->fe(), M_FESpace->dof(), M_FESpace->dof(),  ic,jc, M_offset +ic*totalDof,M_offset + jc*totalDof );
            }

            //mass
            assembleMatrix( *M_mass, *M_elmatM, M_FESpace->fe(), M_FESpace->dof(),ic,ic, M_offset +  ic*totalDof, M_offset +  ic*totalDof);
        }
    }

    comm()->Barrier();

    M_linearStiff->GlobalAssemble();
    massStiff->GlobalAssemble();
    M_mass->GlobalAssemble();
    *massStiff *= factor; //M_data.dataTime()->getTimeStep() * M_rescaleFactor;
}

template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::iterate(vector_type &_sol)
{
    int status;

    vector_type sol(*M_localMap);

    *M_vel  = ( 2.0 / M_data->dataTime()->getTimeStep() ) * (*M_disp);
    *M_vel -= *M_rhsW;

    //M_Displayer->leaderPrint("sol norm = ", norm(this->sol));

    *M_residual_d  = M_massStiff*sol;
    *M_residual_d -= *M_rhsNoBC;
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::iterate( bchandler_type& bch )
{
    // matrix and vector assembling communication
    matrix_ptrtype matrFull( new matrix_type( *M_localMap, M_massStiff->getMeanNumEntries()));
    *matrFull += *M_massStiff;

    M_rhsNoBC->GlobalAssemble();
    M_rhsW->GlobalAssemble();

    vector_type rhsFull (*M_rhsNoBC);

    // boundary conditions update
    //M_comm->Barrier();

    M_Displayer->leaderPrint("  S-  Applying boundary conditions ...         ");
    Chrono chrono;
    chrono.start();

    applyBoundaryConditions( *matrFull, rhsFull, bch);

    chrono.stop();
    M_Displayer->leaderPrintMax("done in " , chrono.diff());


    // solving the system
    M_linearSolver->setMatrix(*matrFull);
    //matrFull->spy( "matrix.data" );
    //rhsFull.spy( "rhs.dat" );

    M_linearSolver->solveSystem( rhsFull, *M_disp, matrFull);

    //M_disp.spy( "sol.dat" );

    *M_vel  = ( 2.0 / M_data->dataTime()->getTimeStep() ) * (*M_disp);
    *M_vel -= *M_rhsW;

    *M_residual_d =  *M_massStiff * (*M_disp);
    *M_residual_d -= *M_rhsNoBC;

}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::iterateLin( bchandler_type& bch )
{

    matrix_ptrtype matrFull( new matrix_type( *M_localMap, M_massStiff->getMeanNumEntries()));
    *matrFull += *M_massStiff;

    M_rhsNoBC->GlobalAssemble();
    M_rhsW->GlobalAssemble();

    vector_type rhsFull (M_rhsNoBC->map());

    M_Displayer->leaderPrint(" LS-  Applying boundary conditions ...         ");
    Chrono chrono;
    chrono.start();

    // boundary conditions update
    applyBoundaryConditions( *matrFull, rhsFull, bch);

    chrono.stop();
    M_Displayer->leaderPrintMax("done in " , chrono.diff());

    //M_Displayer->leaderPrint("rhs_dz norm = ", rhsFull.NormInf() );

    M_linearSolver->setMatrix(*matrFull);
    int numIter = M_linearSolver->solveSystem(  rhsFull, *M_disp, matrFull );


    //M_Displayer->leaderPrintMax("dz norm     = " , M_disp.NormInf() );

    numIter = abs(numIter);

    *M_residual_d =  *M_massStiff*(*M_disp);
//    M_residual_d -= M_rhsNoBC;

}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::evalResidual( vector_type &res, const vector_type& sol, int /*iter*/)
{
    M_Displayer->leaderPrint("  S-  Computing residual ...                   ");
    Chrono chrono;
    chrono.start();

    // Matrices initialization
    M_stiff = M_massStiff;


    M_Displayer->leaderPrint("updating the boundary conditions ... ");
    if ( !M_BCh->bdUpdateDone() )
        M_BCh->bdUpdate( M_FESpace->mesh(), M_FESpace->feBd(), M_FESpace->dof() );

    bcManageMatrix( M_stiff, *M_FESpace->mesh(), M_FESpace->dof(), *M_BCh, M_FESpace->feBd(), 1.0 );

    *M_rhs = *M_rhsNoBC;

    bcManageVector( *M_rhs, *M_FESpace->mesh(), M_FESpace->dof(), *M_BCh, M_FESpace->feBd(), M_data->dataTime()->getTime(), 1.0 );

    res  = M_stiff * sol;
//    res -= M_rhs;

    chrono.stop();
    M_Displayer->leaderPrintMax("done in ", chrono.diff() );
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::evalConstraintTensor()
{
    vector_type count(*M_localMap);

    *M_sxx *= 0.;
    *M_syy *= 0.;
    *M_szz *= 0.;

    for ( UInt ielem = 1; ielem <= M_FESpace->mesh()->numVolumes(); ielem++ )
    {
        //UInt elem = M_FESpace->mesh()->volumeList( ielem ).id();
        M_FESpace->fe().updateFirstDerivQuadPt( M_FESpace->mesh()->volumeList( ielem ) );

        //int    marker = M_FESpace->mesh()->volumeList( ielem ).marker();
        Real s      = 0;
        Real volume = M_FESpace->fe().detJac(0);

        for ( int ig = 0; ig < M_FESpace->fe().nbQuadPt; ++ig )
        {
            for ( int k = 0; k < M_FESpace->fe().nbNode; ++k )
            {
                int i    = M_FESpace->fe().patternFirst(k);
                int idof = M_FESpace->dof().localToGlobal(M_FESpace->fe().currentLocalId(), i + 1);

                s+= (2*M_data->mu() + M_data->lambda())*
                    M_FESpace->fe().weightDet( ig )*
                    M_FESpace->fe().phiDer( k, 0 , ig )*
                    (*M_disp)[idof + 0*M_FESpace->dim()];

                s+= M_data->lambda()*
                    M_FESpace->fe().weightDet( ig )*
                    M_FESpace->fe().phiDer( k, 1 , ig )*
                    (*M_disp)[idof + 1*M_FESpace->dim()];

                s+= M_data->lambda()*
                    M_FESpace->fe().weightDet( ig )*
                    M_FESpace->fe().phiDer( k, 2 , ig )*
                    (*M_disp)[idof + 2*M_FESpace->dim()];

                count[idof]++;
            }
        }

        for ( int k = 0; k < M_FESpace->fe().nbNode; ++k )
        {
            int i    = M_FESpace->fe().patternFirst(k);
            int idof = M_FESpace->dof().localToGlobal(M_FESpace->fe().currentLocalId(), i + 1);

            (*M_sxx)[idof] += s/M_FESpace->fe().detJac(0);
        }

        s = 0;

        for ( int ig = 0; ig < M_FESpace->fe().nbQuadPt; ++ig )
        {
            for ( int k = 0; k < M_FESpace->fe().nbNode; ++k )
            {
                int i    = M_FESpace->fe().patternFirst(k);
                int idof = M_FESpace->dof().localToGlobal(M_FESpace->fe().currentLocalId(), i + 1);

                s += M_data->lambda()*
                     M_FESpace->fe().weightDet( ig )*
                     M_FESpace->fe().phiDer( k, 0 , ig )*
                     (*M_disp)[idof + 0*M_FESpace->dim()];

                s += (2*M_data->mu() + M_data->lambda())*
                     M_FESpace->fe().weightDet( ig )*
                     M_FESpace->fe().phiDer( k, 1 , ig )*
                     (*M_disp)[idof + 1*M_FESpace->dim()];

                s += M_data->lambda()*
                     M_FESpace->fe().weightDet( ig )*
                     M_FESpace->fe().phiDer( k, 2 , ig )*
                     (*M_disp)[idof + 2*M_FESpace->dim()];
//                         M_sxx[idof] += s;

//                M_syy[idof] += s/volume;
            }
        }

        for ( int k = 0; k < M_FESpace->fe().nbNode; ++k )
        {
            int i    = M_FESpace->fe().patternFirst(k);
            int idof = M_FESpace->dof().localToGlobal(M_FESpace->fe().currentLocalId(), i + 1);

            (*M_syy)[idof] += s/volume;
        }


        s = 0;

        for ( int ig = 0; ig < M_FESpace->fe().nbQuadPt; ++ig )
        {
            for ( int k = 0; k < M_FESpace->fe().nbNode; ++k )
            {
                int i    = M_FESpace->fe().patternFirst(k);
                int idof = M_FESpace->dof().localToGlobal(M_FESpace->fe().currentLocalId(), i + 1);

                s += M_data->lambda()*
                     M_FESpace->fe().weightDet( ig )*
                     M_FESpace->fe().phiDer( k, 0 , ig )*
                     (*M_disp)[idof + 0*M_FESpace->dim()];

                s += M_data->lambda()*
                     M_FESpace->fe().weightDet( ig )*
                     M_FESpace->fe().phiDer( k, 1 , ig )*
                     (*M_disp)[idof + 1*M_FESpace->dim()];

                s += (2*M_data->mu() + M_data->lambda())*
                     M_FESpace->fe().weightDet( ig )*
                     M_FESpace->fe().phiDer( k, 2 , ig )*
                     (*M_disp)[idof + 2*M_FESpace->dim()];


//                         M_sxx[idof] += s;
            }
        }

        for ( int k = 0; k < M_FESpace->fe().nbNode; ++k )
        {
            int i    = M_FESpace->fe().patternFirst(k);
            int idof = M_FESpace->dof().localToGlobal(M_FESpace->fe().currentLocalId(), i + 1);

            (*M_szz)[idof] += s/M_FESpace->fe().detJac(0);
        }

    }

    for (UInt ii = 1; ii <= M_FESpace->dim(); ++ii)
    {
        (*M_sxx)[ii] /= count[ii];
        (*M_syy)[ii] /= count[ii];
        (*M_szz)[ii] /= count[ii];
    }
}


template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::initialize( vector_ptrtype disp, vector_ptrtype vel, vector_ptrtype /*acc*/)
{
    *M_disp = *disp;
    if (vel.get())
        initializeVel(*vel);
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::initializeVel( const vector_type& vel)
{
    *M_vel = vel;
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::initialize( const Function& d0, const Function& w0, const Function& /*a0*/ )
{
    M_FESpace->interpolate(d0, *M_disp, 0.0);
    M_FESpace->interpolate(w0, *M_vel , 0.0);
}


template <typename Mesh, typename SolverType> // for monolithic
void
VenantKirchhofSolver<Mesh, SolverType>::updateVel()
{
    *M_vel  = ( 2.0 / M_data->dataTime()->getTimeStep() ) * (*M_disp);
    *M_vel -= *M_rhsW;
}

template<typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::reduceSolution( Vector& disp, Vector& vel )
{
    vector_type displacement(*M_disp, 0);
    vector_type velocity    (*M_vel , 0);

    if ( comm()->MyPID() == 0 )
    {
        for ( UInt iDof = 0; iDof < nDimensions*dim(); ++iDof )
        {
            disp[ iDof ] = displacement[ iDof + 1 ];
            vel [ iDof ] = velocity    [ iDof + 1 ];
        }
    }
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::rescaleMatrices()
{
    *M_mass *=(M_data->dataTime()->getTimeStep()*M_rescaleFactor);
    *M_linearStiff *= (M_data->dataTime()->getTimeStep()*M_rescaleFactor);
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::setDataFromGetPot( const GetPot& dataFile )
{
    M_linearSolver->setDataFromGetPot( dataFile, "solid/solver" );
    M_linearSolver->setUpPrec(dataFile, "solid/prec");

    UInt marker = M_FESpace->mesh()->volumeList( 1 ).marker();
    if (!M_data->young(marker))
        M_data->setYoung(dataFile( "solid/physics/young", 0. ), marker);
    if (!M_data->poisson(marker))
        M_data->setPoisson(dataFile( "solid/physics/poisson", 0. ), marker);
}


template<typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::applyBoundaryConditions( matrix_type&        matrix,
                                                                 vector_type&        rhs,
                                                                 bchandler_type&     BCh,
                                                                 UInt                offset)
{
    // BC manage for the velocity
    if (offset)
        BCh->setOffset(offset);
    if ( !BCh->bdUpdateDone() )
        BCh->bdUpdate( *M_FESpace->mesh(), M_FESpace->feBd(), M_FESpace->dof() );

    // vector_type rhsFull(rhs, Repeated, Zero); // ignoring non-local entries, Otherwise they are summed up lately
    vector_type rhsFull(rhs, Unique);  // bcManages now manages the also repeated parts

    bcManage( matrix, rhsFull, *M_FESpace->mesh(), M_FESpace->dof(), *BCh, M_FESpace->feBd(), 1.,
              M_data->dataTime()->getTime() );

    // matrix should be GlobalAssembled by  bcManage

    rhs = rhsFull;

}

} /*namespace LifeV*/
