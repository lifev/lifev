//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing a solver class for the 1D model.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @author Tiziano Passerini
 *  @author Lucia Mirabella
 *  @date 01-10-2006
 *
 *  @version 2.0
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *  @date 01-08-2009
 *
 *  @version 2.1
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 21-04-2010
 *
 */

#include <lifemc/lifesolver/OneDimensionalModel_Solver.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_Solver::OneDimensionalModel_Solver():
        M_Physics               (),
        M_Flux                  (),
        M_Source                (),
        M_FESpace               (),
        M_Comm                  (),
        M_leftNodeId            (),
        M_leftInternalNodeId    (),
        M_rightNodeId           (),
        M_rightInternalNodeId   (),
        M_coeffMass             (),
        M_coeffStiff            (),
        M_coeffGrad             (),
        M_coeffDiv              (),
        M_elmatMass             (),
        M_elmatStiff            (),
        M_elmatGrad             (),
        M_elmatDiv              (),
        M_U_thistime            ( new Solution_Type() ),
        M_U_prevtime            (),
        M_U_2prevtime           (),
        M_rhs                   (),
        M_FluxVector            (),
        M_diffFlux              ( 4 ),
        M_SourceVector          (),
        M_diffSrc               ( 4 ),
        M_massMatrix            (),
        M_massMatrixDiffSrc     ( 4  ),
        M_stiffMatrixDiffFlux   ( 4  ),
        M_gradMatrix            (),
        M_gradMatrixDiffFlux    ( 4  ),
        M_divMatrixDiffSrc      ( 4 ),
        M_linearSolver          (),
        M_post_process_buffer   (),
        M_post_process_buffer_offset (),
        M_bcDirLeft             (),
        M_bcDirRight            (),
        M_variable_index_map    (),
        M_variable_filter_map   ()
{
}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalModel_Solver::setup()
{
    // These maps allow to use the switch... case construct with string variables
    M_oneDstring2initializeVarMap["A"]      = OneDInitArea;
    M_oneDstring2initializeVarMap["Q"]      = OneDInitFlux;
    M_oneDstring2initializeVarMap["W1"]     = OneDInitRiemann1;
    M_oneDstring2initializeVarMap["W2"]     = OneDInitRiemann2;
    M_oneDstring2initializeVarMap["P"]      = OneDInitPressure;

    Debug( 6310 ) << "[OneDimensionalModel_Solver] O-  Nb of unknowns: " << M_FESpace->dim()     << "\n";
    Debug( 6310 ) << "[OneDimensionalModel_Solver] O-  Computing the constant matrices... \n";
    Debug( 6310 ) << "[OneDimensionalModel_Solver] O-  Adopting a"
                  << std::string( M_Physics->Data()->viscoelasticWall() ? " viscoelastic " : "n elastic " )
                  << "model for vessel wall... \n";

    Chrono chrono;
    chrono.start();

    // Vector and Matrices initialization
    // insert variables!
    UInt nvar(0);
    M_variable_index_map.insert( make_pair("A",  nvar++ ) );
    M_variable_index_map.insert( make_pair("Q",  nvar++ ) );
    M_variable_index_map.insert( make_pair("W1", nvar++ ) );
    M_variable_index_map.insert( make_pair("W2", nvar++ ) );

    M_variable_index_map.insert( make_pair("P",  nvar++ ) );

    // correction flux with viscoelastic term
    if( M_Physics->Data()->viscoelasticWall() )
    {
        // correction on the flux due to the viscoelastic term
        M_variable_index_map.insert( make_pair("Q_visc",  nvar++ ) );
        // viscoelastic contribution to the pressure
        M_variable_index_map.insert( make_pair("P_visc",  nvar++ ) );
        // elastic contribution to the pressure
        M_variable_index_map.insert( make_pair("P_elast", nvar++ ) );
        // time derivative of the section area
        M_variable_index_map.insert( make_pair("dA_dt",   nvar++ ) );
    }

    // flux second derivative
    if( M_Physics->Data()->fluxSecondDer() )
        M_variable_index_map.insert( make_pair("d2Q_dx2", nvar++ ) );

    // correction flux with inertial term
    if( M_Physics->Data()->inertialWall() )
        M_variable_index_map.insert( make_pair("Q_inert", nvar++ ) );

    // correction flux with longitudinal term
    if( M_Physics->Data()->longitudinalWall() )
        M_variable_index_map.insert( make_pair("Q_long", nvar++ ) );

    // activate the export filters
    // matlab postprocessing

    M_variable_filter_map.insert( make_pair(".m",
                                            &LifeV::OneDimensionalModel_Solver::output_to_matlab) );

    M_U_thistime->resize(nvar, Vector_Type(M_FESpace->map()));

    for (UInt ii = 0; ii < nvar; ++ii)
        (*M_U_thistime)[ii] = 0.;

    // initialize matrices
    std::fill( M_diffFlux.begin(),    M_diffFlux.end(),    ublas::zero_vector<Real>(M_FESpace->dim()) );
    std::fill( M_diffSrc.begin() ,    M_diffSrc.end() ,    ublas::zero_vector<Real>(M_FESpace->dim()) );

    for( UInt i = 0; i < 4; ++i )
    {
        M_massMatrixDiffSrc[i].reset(new Matrix_Type( M_FESpace->map() ));
        M_stiffMatrixDiffFlux[i].reset( new Matrix_Type( M_FESpace->map() ));
        M_gradMatrixDiffFlux[i].reset( new Matrix_Type( M_FESpace->map() ));
        M_divMatrixDiffSrc[i].reset( new Matrix_Type( M_FESpace->map() ));
    }

    // set the coeff to 1.
    M_coeffMass = 1.;
    M_coeffGrad = 1.;

    // Elementary computation and matrix assembling
    // Loop on elements
    for(UInt iedge = 1; iedge <= M_FESpace->mesh()->numEdges(); iedge++)
    {
        // set the elementary matrices to 0.
        M_elmatMass->zero();
        M_elmatGrad->zero();

        // update the current element

        M_FESpace->fe().updateFirstDerivQuadPt(M_FESpace->mesh()->edgeList(iedge));

        // update the mass and grad matrices
        mass( M_coeffMass, *M_elmatMass, M_FESpace->fe(),0, 0 );
        grad( 0 , - M_coeffGrad, *M_elmatGrad, M_FESpace->fe(), M_FESpace->fe(), 0, 0 );
        // assemble the mass and grad matrices
        assemb_mat( *M_massMatrix, *M_elmatMass, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );
        //assemb_mat( M_factorMassMatrix, M_elmatMass, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );
        assemb_mat( *M_gradMatrix, *M_elmatGrad, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );
    }

    // Dirichlet boundary conditions set in the mass matrix
    M_massMatrix->GlobalAssemble();
    M_gradMatrix->GlobalAssemble();

    chrono.stop();

    Debug( 6310 ) << "[OneDimensionalModel_Solver] \tdone in " << chrono.diff() << " s.\n";
}

void
OneDimensionalModel_Solver::initialize()
{
    ASSERT_PRE( (M_leftNodeId <= M_rightNodeId) && (M_rightNodeId <= M_FESpace->dim()),
                "[initialize] outside tube boundaries" );

    Real value1, value2;

    Vector exponent( M_FESpace->dim() );
    for ( UInt inode=M_leftNodeId; inode <= M_rightNodeId ; ++inode )
    {
        exponent[inode - 1] *= 0.;

        // first half of a gaussian signal, centered in firstnode;
        // second half of a gaussian signal, centered in lastnode;
        // width represents the total duration of the gaussian signal
        // (rise + decay)
        if( (inode < M_leftNodeId) || (inode > M_rightNodeId) )
        {
            exponent[inode - 1]  = - std::pow( Real(int( inode - M_leftNodeId )), 2 );
            exponent[inode - 1] /= 2*std::pow( M_Physics->Data()->Length(), 2 );
        }
    }

    switch( M_oneDstring2initializeVarMap[M_Physics->Data()->initialVariable()] )
    {
        // case 1, 2: initialize physical variables to desired value
        case OneDInitPressure:
            //std::cout << "OneDInitPressure" << std::endl;
            // this is a pressure value! has to be converted in area value
            Debug( 6310 ) << "[initialize] 0- OneDInitPressure\n";

            value1 = M_Physics->Data()->restValue();
            // HYPOTHESIS: when initializing pressure, flux is imposed constant = 0
            value2 = 0;

            Debug( 6310 ) << "[initialize] pressure " << value1 << "\n";
            (*M_U_thistime)[1] = Vector_Type(M_FESpace->map());
            (*M_U_thistime)[1] = value2;

            Debug( 6310 ) << "[initialize] Q done\n";

            for ( UInt inode = M_leftNodeId; inode <= M_rightNodeId ; ++inode )
            {
                // reusing value2 as help variable
                value2 = value1*( 1 + M_Physics->Data()->multiplier()*std::exp( exponent[inode - 1] ) );
                (*M_U_thistime)[0][inode] = M_Physics->A_from_P( value2 );

                Debug( 6310 ) << "[initialize] A(" << inode <<") done\n";
                M_Physics->W_from_U( (*M_U_thistime)[2][inode], (*M_U_thistime)[3][inode],
                                     (*M_U_thistime)[0][inode], (*M_U_thistime)[1][inode],
                                     inode - 1);
                Debug( 6310 ) << "[initialize] W_i(" << inode <<") done\n";
            }
            break;

        case OneDInitArea:
            Debug( 6310 ) << "[initialize] 0- OneDInitArea\n";

            //std::cout << "OneDInitArea" << std::endl;
            value1 = M_Physics->Data()->initialValue();
            value2 = 0;

            //ScalarVector( (*M_U_thistime)[1].size(), value2 );
            //(*M_U_thistime)[0] = Vector_Type(M_FESpace->map());
            (*M_U_thistime)[1] = value2;


            for ( UInt inode=M_leftNodeId; inode <= M_rightNodeId ; ++inode )
            {
                (*M_U_thistime)[0][inode] = value1 *
                  ( 1 + M_Physics->Data()->multiplier() * std::exp( exponent[inode - 1] ) );

                M_Physics->W_from_U( (*M_U_thistime)[2][inode], (*M_U_thistime)[3][inode],
                                     (*M_U_thistime)[0][inode], (*M_U_thistime)[1][inode],
                                     inode - 1 );
            }
            break;

        case OneDInitFlux:
            // HYPOTHESIS: when initializing flux, area is equal to Area0
            Debug( 6310 ) << "[initialize] 0- OneDInitFlux\n";

            value1 = M_Physics->Data()->Area0(0); // this if Area0 is constant
            value2 = M_Physics->Data()->initialValue();

            for ( UInt inode = M_leftNodeId; inode <= M_rightNodeId ; ++inode )
            {
                (*M_U_thistime)[0][inode] = M_Physics->Data()->Area0(inode - 1);
                (*M_U_thistime)[1][inode] = value2*( 1 + M_Physics->Data()->multiplier()*std::exp( exponent[inode - 1] ) );

                M_Physics->W_from_U( (*M_U_thistime)[2][inode], (*M_U_thistime)[3][inode],
                                     (*M_U_thistime)[0][inode], (*M_U_thistime)[1][inode],
                                     inode - 1 );
            }
            break;

        case OneDInitRiemann1:
            Debug( 6310 ) << "[initialize] 0- OneDInitRiemann1\n";
            value1 = M_Physics->Data()->initialValue();
            value2 = -value1;
            std::cout << "[initialize] WARNING! Initializing W2 = - W1"
                      << " (assuming Q = 0)" << std::endl;

            for ( UInt inode=M_leftNodeId; inode <= M_rightNodeId ; ++inode )
            {
                (*M_U_thistime)[2][inode] = value1 *
                  ( 1 + M_Physics->Data()->multiplier() * std::exp( exponent[inode - 1] ) );
                (*M_U_thistime)[3][inode] = value2 *
                  ( 1 + M_Physics->Data()->multiplier() * std::exp( exponent[inode - 1] ) );
                M_Physics->U_from_W( (*M_U_thistime)[0][inode],
                                     (*M_U_thistime)[1][inode],
                                     (*M_U_thistime)[2][inode],
                                     (*M_U_thistime)[3][inode],
                                     inode - 1);
            }

            break;

        case OneDInitRiemann2:
            Debug( 6310 ) << "[initialize] 0- OneDInitRiemann2\n";
            value1 = M_Physics->Data()->initialValue();
            value2 = - value1;
            std::cout << "[initialize] WARNING! Initializing W1 = - W2"
                      << " (assuming Q = 0)" << std::endl;

            for (UInt inode = M_leftNodeId; inode <= M_rightNodeId; ++inode )
            {
                (*M_U_thistime)[3][inode] = value1 *
                  ( 1 + M_Physics->Data()->multiplier() * std::exp( exponent[inode - 1] ) );
                (*M_U_thistime)[2][inode] = value2 *
                  ( 1 + M_Physics->Data()->multiplier() * std::exp( exponent[inode - 1] ) );
                M_Physics->U_from_W( (*M_U_thistime)[0][inode],
                                     (*M_U_thistime)[1][inode],
                                     (*M_U_thistime)[2][inode],
                                     (*M_U_thistime)[3][inode],
                                     inode );
            }

            break;

        default:
            ERROR_MSG("No such initializing option.");
    }

    Debug( 6310 ) << "[initialize]\t\tvalue1         = " << value1 << "\n";
    Debug( 6310 ) << "[initialize]\t\tvalue1_step    = " << value1 * M_Physics->Data()->multiplier() << "\n";
    Debug( 6310 ) << "[initialize]\t\tvalue2         = " << value2 << "\n";

    for( UInt i = 0; i < 2; ++i )
    {
        M_U_prevtime [i] = (*M_U_thistime)[i];
        M_U_2prevtime[i] = M_U_prevtime[i];
    }

    // Compute initial pressure (taking into account the viscoelastic wall)
    _updatePressure( M_Physics->Data()->dataTime()->getTimeStep() );

    // Prepare the buffers
    openFileBuffers();
}



void
OneDimensionalModel_Solver::initialize( const Real& u10, const Real& u20, const std::string& var )
{
    Debug( 6310 ) << "[OneDimensionalModel_Solver::initialize] O- Initialize: var " << var << "\n";

    if( var == "physical")
    {
        Debug( 6310 ) << "[OneDimensionalModel_Solver::initialize] O- Imposing real values ... \n";
        Debug( 6310 ) << "[OneDimensionalModel_Solver::initialize] O- A0 = " << u10 << "\n";
        Debug( 6310 ) << "[OneDimensionalModel_Solver::initialize] O- Q0 = " << u20 << "\n";

        (*M_U_thistime)[0] = Vector_Type( M_FESpace->map() );
        //(*M_U_thistime)[0][LeftNodeId()] = u10;

        (*M_U_thistime)[1] = Vector_Type( M_FESpace->map() );
        //(*M_U_thistime)[1][LeftNodeId()] = u20;

        for (UInt ielem = 0; ielem < M_FESpace->dim() ; ielem++ )
        {
            (*M_U_thistime)[0][ielem + 1] = u10;
            (*M_U_thistime)[1][ielem + 1] = u20;

            //              std::cout << ielem << " " << (*M_U_thistime)[0][ielem + 1] << std::endl;
            M_Physics->W_from_U( (*M_U_thistime)[2][ielem + 1], (*M_U_thistime)[3][ielem + 1],
                                 (*M_U_thistime)[0][ielem + 1], (*M_U_thistime)[1][ielem + 1],
                                 ielem );
        }
        Debug( 6310 ) << "[OneDimensionalModel_Solver::initialize] O- ok\n";
    }
    else if( var == "Riemann" )
    {
        //        for( UInt i=0; i<2; ++i )
        (*M_U_thistime)[2] = Vector_Type( M_FESpace->map() );
        (*M_U_thistime)[2] = u10;
        (*M_U_thistime)[3] = Vector_Type( M_FESpace->map() );
        (*M_U_thistime)[3] = u20;

        for (UInt ielem = 0; ielem < M_FESpace->dim() ; ielem++ )
            M_Physics->U_from_W( (*M_U_thistime)[0][ielem + 1], (*M_U_thistime)[1][ielem + 1],
                                   (*M_U_thistime)[2][ielem + 1], (*M_U_thistime)[3][ielem + 1],
                                   ielem + 1 ); // WARNING the +1 is not debugged yet (GF 12/2009)
    }
    else
    {
        std::cout << "[initialize] trying to initialize " << var << " variables!" << std::endl;
        abort();
    }

    for( UInt i = 0; i < 2; ++i )
    {
        M_U_prevtime [i] = (*M_U_thistime)[i];
        M_U_2prevtime[i] = M_U_prevtime[i];
    }

    // Compute initial pressure (taking into account the viscoelastic wall)
    _updatePressure( M_Physics->Data()->dataTime()->getTimeStep() );

    // Prepare the buffers
    openFileBuffers();
}

void
OneDimensionalModel_Solver::initialize( const Vector_Type& u10, const Vector_Type& u20 )
{
    (*M_U_thistime)[0] = u10;
    (*M_U_thistime)[1] = u20;

    for (UInt ielem=0; ielem <= M_FESpace->dim() ; ielem++ )
    {
        M_Physics->W_from_U( (*M_U_thistime)[2][ielem], (*M_U_thistime)[3][ielem],
                             (*M_U_thistime)[0][ielem], (*M_U_thistime)[1][ielem], ielem );
    }

    for( UInt i=0; i<2; ++i )
    {
        M_U_prevtime[i] = (*M_U_thistime)[i];
        M_U_2prevtime[i] = M_U_prevtime[i];
    }

    // Compute initial pressure (taking into account the viscoelastic wall)
    _updatePressure( M_Physics->Data()->dataTime()->getTimeStep() );

    // Prepare the buffers
    openFileBuffers();
}

void
OneDimensionalModel_Solver::initialize( const Real& u20 )
{
    (*M_U_thistime)[1] = Vector_Type( M_FESpace->map() );
    (*M_U_thistime)[1] = u20;

    for (UInt ielem = 0; ielem <= M_FESpace->dim() ; ielem++ )
    {
        (*M_U_thistime)[0][ielem]=M_Physics->Data()->Area0(ielem);
        M_Physics->W_from_U( (*M_U_thistime)[2][ielem], (*M_U_thistime)[3][ielem],
                             (*M_U_thistime)[0][ielem], (*M_U_thistime)[1][ielem], ielem );
    }

    for( UInt i=0; i<2; ++i )
    {
        M_U_prevtime[i] = (*M_U_thistime)[i];
        M_U_2prevtime[i] = M_U_prevtime[i];
    }

    // Compute initial pressure (taking into account the viscoelastic wall)
    _updatePressure( M_Physics->Data()->dataTime()->getTimeStep() );

    // Prepare the buffers
    openFileBuffers();
}

void
OneDimensionalModel_Solver::iterate( OneDimensionalModel_BCHandler& bcH, const Real& Time, const Real& TimeStep )
{
    // Update the RHS and the matrices
    _updateSystem( TimeStep );

    //---------------------------------------------------
    // 3/ take into account the BOUNDARY CONDITIONS
    //---------------------------------------------------
    // compute the values for the boundary conditions
    Debug( 6310 ) << "[timeAdvance] \tcompute BC\n";
    bcH.applyBC( Time, TimeStep, M_U_thistime, M_Flux, M_bcDirLeft, M_bcDirRight );

    // take into account the bc
    Debug( 6310 ) << "[timeAdvance] \tcompute BC dirichlet vector\n";
    _updateBCDirichletVector();

    Debug( 6310 ) << "[timeAdvance] \trhs0 norm2 = " << M_rhs[0].Norm2() << "\n";
    Debug( 6310 ) << "[timeAdvance] \trhs1 norm2 = " << M_rhs[1].Norm2() << "\n";

    Debug( 6310 ) << "[iterate] o-  Solving the system... t = " << Time << " ... \n";
    Chrono chrono;

    Chrono chrono1;
    Chrono chrono2;
    Chrono chrono3;
    Chrono chrono4;
    Chrono chrono5;

    if( M_Physics->Data()->UW() )
    {
        Real Ainode, Qinode;

        Real lambda1_plus  = 0.;
        Real lambda2_plus  = 0.;

        Real lambda1_minus = 0.;
        Real lambda2_minus = 0.;

        Real eigval1, eigval2;
        Real tmp11, tmp12, tmp21, tmp22;

        //Real deltaX = M_FESpace->mesh()->edgeList( 1 ).point(1).x() - M_FESpace->mesh()->edgeList( 1 ).point(2).x();
        Real delta  = - std::sqrt( ( ( M_FESpace->mesh()->edgeList[ 1 ].point( 2 ) ).x() - ( M_FESpace->mesh()->edgeList[ 1 ].point( 1 ) ).x() ) *
                                   ( ( M_FESpace->mesh()->edgeList[ 1 ].point( 2 ) ).x() - ( M_FESpace->mesh()->edgeList[ 1 ].point( 1 ) ).x() ) +
                                   ( ( M_FESpace->mesh()->edgeList[ 1 ].point( 2 ) ).y() - ( M_FESpace->mesh()->edgeList[ 1 ].point( 1 ) ).y() ) *
                                   ( ( M_FESpace->mesh()->edgeList[ 1 ].point( 2 ) ).y() - ( M_FESpace->mesh()->edgeList[ 1 ].point( 1 ) ).y() ) +
                                   ( ( M_FESpace->mesh()->edgeList[ 1 ].point( 2 ) ).z() - ( M_FESpace->mesh()->edgeList[ 1 ].point( 1 ) ).z() ) *
                                   ( ( M_FESpace->mesh()->edgeList[ 1 ].point( 2 ) ).z() - ( M_FESpace->mesh()->edgeList[ 1 ].point( 1 ) ).z() )
                                 );

        // working on riemann invariants
        Vector_Type W1_UW(M_FESpace->map());
        Vector_Type W2_UW(M_FESpace->map());

        // converting boundary conditions on physical variables
        // in boundary conditions on characteristic variables
        M_Physics->W_from_U( W1_UW(0), W2_UW(0),
                               M_rhs[0][0], M_rhs[1][0], 0 );
        M_Physics->W_from_U( W1_UW(M_FESpace->dim()-1), W2_UW(M_FESpace->dim()-1),
                               M_rhs[0][M_FESpace->dim()-1], M_rhs[1][M_FESpace->dim()-1],
                               M_FESpace->dim() - 1 );

        for ( UInt ii=1; ii < (M_FESpace->dim()-1) ; ii++ )
        {
            // compute the eigenvalues at node
            Ainode = (*M_U_thistime)[0]( ii );
            Qinode = (*M_U_thistime)[1]( ii );
            M_Flux->jacobian_EigenValues_Vectors( Ainode, Qinode,
                                                  eigval1, eigval2,
                                                  tmp11, tmp12,
                                                  tmp21, tmp22,
                                                  ii );

            lambda1_plus = std::max<Real>( eigval1, 0. );
            lambda1_minus = std::min<Real>( eigval1, 0. );
            lambda2_plus = std::max<Real>( eigval2, 0. );
            lambda2_minus = std::min<Real>( eigval2, 0. );
            // update the solution for the next time step
            W1_UW[ii] = (*M_U_thistime)[2][ii]
                      - (TimeStep / delta)  * lambda1_plus  * ( (*M_U_thistime)[2][ii]   - (*M_U_thistime)[2][ii-1])
                      - (TimeStep / delta)  * lambda1_minus * ( (*M_U_thistime)[2][ii+1] - (*M_U_thistime)[2][ii])
                      -  TimeStep * ( tmp11 * M_SourceVector[0][ii] + tmp12 * M_SourceVector[1][ii] );
            W2_UW[ii] = (*M_U_thistime)[3][ii]
                      - (TimeStep / delta)  * lambda2_plus  * ( (*M_U_thistime)[3][ii]   - (*M_U_thistime)[3][ii-1])
                      - (TimeStep / delta)  * lambda2_minus * ( (*M_U_thistime)[3][ii+1] - (*M_U_thistime)[3][ii])
                      -  TimeStep * ( tmp21 * M_SourceVector[0][ii] + tmp22 * M_SourceVector[1][ii] );
        }

        (*M_U_thistime)[2] = W1_UW;
        (*M_U_thistime)[3] = W2_UW;

        for (UInt ielem=0; ielem <= M_FESpace->dim() ; ielem++ )
            M_Physics->U_from_W( (*M_U_thistime)[0][ielem], (*M_U_thistime)[1][ielem],
                                 (*M_U_thistime)[2][ielem], (*M_U_thistime)[3][ielem],
                                 ielem );
    }
    else
    {
        chrono.start();
        // cholesky or lapack lu solve
        // solve the system: rhs1 = massFactor^{-1} * rhs1

        chrono1.start();
        Vector_Type sol0(M_rhs[0]);

        //Matrix_PtrType matrFull( new Matrix_Type( M_FESpace->map(), M_factorMassMatrix.getMeanNumEntries()));
        Matrix_PtrType matrFull( new Matrix_Type( *M_massMatrix ));
        //M_massMatrix.spy("mass");
        _updateBCDirichletMatrix( *matrFull );
        chrono1.stop();
        //M_factorMassMatrix.GlobalAssemble();
        //*matrFull = M_massMatrix;
        //matrFull->GlobalAssemble();

        chrono2.start();
        //matrFull->spy("matr");
        //M_rhs[0].ShowMe();
        M_linearSolver->setMatrix(*matrFull);
        M_linearSolver->setReusePreconditioner(false);
        M_linearSolver->solveSystem( M_rhs[0], sol0, matrFull );
        //std::cout << "sol0 norm2 = " << sol0.Norm2() << std::endl;
        chrono2.stop();

        chrono3.start();
        M_rhs[0] = sol0;
        Vector_Type sol1(M_rhs[1]);

        //std::cout << "rhs1 norm2 = " << M_rhs[1].Norm2() << std::endl;
//             // solve the system: rhs2 = massFactor^{-1} * rhs2
//            M_linearSolver->setMatrix
        M_linearSolver->setReusePreconditioner(false);
        M_linearSolver->solveSystem( M_rhs[1], sol1, matrFull );
        //std::cout << "sol1 norm2 = " << sol1.Norm2() << std::endl;

        M_rhs[1] = sol1;

        chrono3.stop();
        // correct flux with inertial term
        chrono4.start();
        if( M_Physics->Data()->inertialWall() )
        {
            (*M_U_thistime)[M_variable_index_map.find("Q_inert")->second] =
                _correct_flux_inertial( M_rhs[1] );
            M_rhs[1] += (*M_U_thistime)[M_variable_index_map.find("Q_inert")->second];
        }
        // correct flux with viscoelastic term
        if( M_Physics->Data()->viscoelasticWall() )
        {
            (*M_U_thistime)[M_variable_index_map.find("Q_visc")->second] =
                _correct_flux_viscoelastic( M_rhs[1], TimeStep );
            M_rhs[1] += (*M_U_thistime)[M_variable_index_map.find("Q_visc")->second];
        }
        // compute L2 projection of d2Q_dx2
        //    if( M_Physics->Data()->fluxSecondDer() )
        //      M_d2_U2_dx2 = _compute_d2Q_dx2( M_rhs[1] );

        // correct flux with longitudinal term
        if( M_Physics->Data()->longitudinalWall() )
        {
//                     (*M_U_thistime)[M_variable_index_map.find("Q_long")->second] =
//                         _correct_flux_longitudinal(  );
            M_rhs[1] += (*M_U_thistime)[M_variable_index_map.find("Q_long")->second];
        }

        // store solution at previous timesteps & update the solution for the next time step
        for( UInt i=0; i<2; ++i )
        {
            M_U_2prevtime[i] = M_U_prevtime[i];
            M_U_prevtime [i] = (*M_U_thistime)[i];
            (*M_U_thistime)[i] = M_rhs[i];
        }

        chrono4.stop();
        chrono5.start();
        for (UInt ielem = 0; ielem < M_FESpace->dim() ; ielem++ )
        {
            M_Physics->W_from_U( (*M_U_thistime)[2][ielem + 1], (*M_U_thistime)[3][ielem + 1],
                                 (*M_U_thistime)[0][ielem + 1], (*M_U_thistime)[1][ielem + 1],
                                  ielem);
        }

        // Update the pressure (taking into account the viscoelastic wall)
        _updatePressure( TimeStep );

        chrono5.stop();
        chrono.stop();
    }

    Debug( 6310 ) << "chrono1 " << chrono1.diff() << "s.\n";
    Debug( 6310 ) << "chrono2 " << chrono2.diff() << "s.\n";
    Debug( 6310 ) << "chrono3 " << chrono3.diff() << "s.\n";
    Debug( 6310 ) << "chrono4 " << chrono4.diff() << "s.\n";
    Debug( 6310 ) << "chrono5 " << chrono5.diff() << "s.\n";
    Debug( 6310 ) << "chrono  " << chrono.diff() << " s.\n";

    Debug( 6310 ) << "[iterate] \tdone in " << chrono.diff() << " s.\n";
}

Real
OneDimensionalModel_Solver::ComputeCFL( const Real& timeStep ) const
{
    Real CFL = 0.;

    Real lambdaMax = 0.;

    Real eigval1, eigval2;
    Real tmp11, tmp12, tmp21, tmp22;

    for ( UInt inode(1); inode <= M_FESpace->dim() ; ++inode )
    {
        // compute the eigenvalues at node
        M_Flux->jacobian_EigenValues_Vectors( (*M_U_thistime)[0]( inode ), //A
                                              (*M_U_thistime)[1]( inode ), //Q
                                              eigval1, eigval2,
                                              tmp11, tmp12, tmp21, tmp22, inode );

        lambdaMax = std::max<Real>( std::max<Real>( std::fabs(eigval1), std::fabs(eigval2) ), lambdaMax );
    }

    return CFL = lambdaMax * timeStep / M_FESpace->mesh()->minH();
}

void
OneDimensionalModel_Solver::savesol()
{
    M_U_prevtime[0] = (*M_U_thistime)[0];
    M_U_prevtime[1] = (*M_U_thistime)[1];
}

void
OneDimensionalModel_Solver::loadsol()
{
    Container2D_Type U_leftbd, U_rightbd;
    for( UInt i=0; i<2; ++i )
    {
        U_leftbd[i] = (*M_U_thistime)[i]( M_leftNodeId );
        U_rightbd[i] = (*M_U_thistime)[i]( M_rightNodeId );
    }

    (*M_U_thistime)[0] = M_U_prevtime[0];
    (*M_U_thistime)[1] = M_U_prevtime[1];

    for( UInt i=0; i<2; ++i )
    {
        (*M_U_thistime)[i]( M_leftNodeId ) = U_leftbd[i];
        (*M_U_thistime)[i]( M_rightNodeId ) = U_rightbd[i];
    }

    for ( UInt inode=M_leftNodeId; inode <= M_rightNodeId ; ++inode )
        M_Physics->W_from_U( (*M_U_thistime)[2][inode], (*M_U_thistime)[3][inode],
                             (*M_U_thistime)[0][inode], (*M_U_thistime)[1][inode],
                             inode - 1 );
}

void
OneDimensionalModel_Solver::postProcess( const Real& /*time_val*/ )
{
    std::string str;

    std::ofstream outfile;

    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator it;
    std::map< std::string, UInt>::iterator iter;

    Debug( 6310 ) << "[postProcess] o- Dumping solutions on files (1d)!" << "/n";

    it  = M_post_process_buffer.begin();
    for( iter = M_variable_index_map.begin(); iter != M_variable_index_map.end(); ++iter, ++it)
    {
        std::string varname  = iter->first;

        int offset           = iter->second;

        std::string filename = it->first;
        Debug( 6310 ) << "Writing file " << filename << " with var " << offset << "/n";

        outfile.open( filename.c_str(), std::ios::app );
        //outfile << "# time = " << time_val << std::endl;

        for ( UInt ii(0); ii < M_FESpace->dim(); ++ii )
        {
            outfile << (*M_U_thistime)[offset](ii + 1) << " ";
            //outfile << (*(*iter).second).str();
        }
        outfile << std::endl;
        outfile.close();
    }
    resetFileBuffers();
}

void
OneDimensionalModel_Solver::output_to_matlab( std::string  fname,
                                              Real         /*time_val*/,
                                              const        Vector_Type& U,
                                              std::string  /*vname*/ )
{
#if 0

    fstream filestr;

    filestr.open (fname.c_str(), fstream::in | fstream::out | fstream::app);

    //filestr << "# time " << time_val << std::endl;
    for ( int ii = LeftNodeId(); ii <= RightNodeId() ; ++ii )
        {
            //        (*buf) << U[ii] << "; ";
            //std::cout << U[ii] << " " << std::endl;
            filestr << U[ii] << " ";
            //test << U[ii] << std::endl;
        }
    filestr << "\n";
    filestr.close();

#else
    boost::shared_ptr<std::ostringstream> buf;

    //std::cout << "buffer = " << fname << std::endl;
    buf = M_post_process_buffer[fname];
    //buf = &std::cout;
//     std::ostringstream test;
//     test.open("test.txt");

    //    (*buf) << vname << M_Physics->Data()->PostFile()
    //           << "( " << (static_cast<int>( std::floor( time_val/M_Physics->Data()->dataTime()->getTimeStep() + 0.5 ) )
    //              /  M_Physics->Data()->verbose() )+1
    //           << ", : ) = [ " << std::flush;

    //std::cout << "Norm2 U = " << U.Norm2() << std::endl;

    for ( UInt ii(LeftNodeId()); ii <= RightNodeId() ; ++ii )
    {
        //        (*buf) << U[ii] << "; ";
        //std::cout << U[ii] << std::endl;
        (*buf) << U[ii] << " ";
        //std::cout << buf->str();
        //test << U[ii] << std::endl;
    }
    //    (*buf) << "]';\n" << std::endl;
    (*buf) << "\n";
    //std::cout << (*buf) << std::endl;
#endif
}

void
OneDimensionalModel_Solver::create_movie_file()
{
    std::ofstream outfile;
    std::string file_output;
    file_output = M_Physics->Data()->PostDirectory() + "/" + "areamovie"+M_Physics->Data()->PostFile()+".m";
    outfile.open(file_output.c_str(), std::ios::out);
    outfile <<"Area"<<M_Physics->Data()->PostFile()<<";\n"
            <<"[n,m]=size(A" << M_Physics->Data()->PostFile() << ");\n"
            <<"Max=max(max(A" << M_Physics->Data()->PostFile() << "));\n"
            <<"Min=min(min(A" << M_Physics->Data()->PostFile() << "));\n"
            <<"for i=1:n\n"
            <<"  plot(A" << M_Physics->Data()->PostFile() << "(i,:));\n"
            <<"  title((i-1)*"<<M_Physics->Data()->dataTime()->getTimeStep()<<"*"<<M_Physics->Data()->verbose()<<");\n"
            <<"  axis([0 m Min-(Min/100) Max+(Min/100)]);\n"
            <<"  %pause;\n"
            <<"  F(i) = getframe;\n"
            <<"end\n"
            <<"movie(F)";

    outfile.close();
    file_output = M_Physics->Data()->PostDirectory() + "/" + "portatamovie"+M_Physics->Data()->PostFile()+".m";
    outfile.open(file_output.c_str(), std::ios::out);
    outfile <<"Portata"<<M_Physics->Data()->PostFile()<<";\n"
            <<"[n,m]=size(Q" << M_Physics->Data()->PostFile() << ");\n"
            <<"Max=max(max(Q" << M_Physics->Data()->PostFile() << "));\n"
            <<"Min=min(min(Q" << M_Physics->Data()->PostFile() << "));\n"
            <<"for i=1:n\n"
            <<"  plot(Q" << M_Physics->Data()->PostFile() << "(i,:));\n"
            <<"  title((i-1)*"<< M_Physics->Data()->dataTime()->getTimeStep() <<"*"<<M_Physics->Data()->verbose()<<");\n"
            <<"  axis([0 m Min-(Min/100) Max+(Min/100)]);\n"
            <<"  %pause;\n"
            <<"  F(i) = getframe;\n"
            <<"end\n"
            <<"movie(F)";
    outfile.close();
}

void
OneDimensionalModel_Solver::openFileBuffers()
{
    std::string file_output;

    boost::shared_ptr<std::ostringstream> buf;

    M_variable_index_iter iter_variable;
    M_variable_filter_iter iter_suffix;

    for( iter_variable = M_variable_index_map.begin(); iter_variable != M_variable_index_map.end(); ++iter_variable )
    {

        for( iter_suffix = M_variable_filter_map.begin(); iter_suffix != M_variable_filter_map.end(); ++iter_suffix )
        {

            file_output = M_Physics->Data()->PostDirectory() + "/" + M_Physics->Data()->PostFile() + iter_variable->first + iter_suffix->first;
            Debug( 6310 ) << "[openFileBuffers] setting output for file " << file_output << "\n";

            buf.reset( new std::ostringstream() );
            buf->setf(std::ios_base::scientific);
            buf->precision(5);
            buf->width(13);

            M_post_process_buffer.insert( std::map<std::string,
                                          boost::shared_ptr<std::ostringstream> >::value_type( file_output, buf ) );

            M_post_process_buffer_offset.insert( std::map<std::string,
                                                 long>::value_type( file_output, buf->tellp() ) );
        }
    }
}

void
OneDimensionalModel_Solver::output2FileBuffers( const Real& time_val )
{
    Debug( 6310 ) << "[output2FileBuffers] begin \n";

    if( !( static_cast<int>( std::floor( (time_val/M_Physics->Data()->dataTime()->getTimeStep()) + .5 ) ) % M_Physics->Data()->verbose() ) )
    {

        Debug( 6310 ) << "[output2FileBuffers] writting output \n";

        std::string file_output;

        M_variable_index_iter  iter_variable;
        M_variable_filter_iter iter_suffix;

        for( iter_variable = M_variable_index_map.begin(); iter_variable != M_variable_index_map.end(); ++iter_variable )
        {
            for( iter_suffix = M_variable_filter_map.begin(); iter_suffix != M_variable_filter_map.end(); ++iter_suffix )
            {
                file_output = M_Physics->Data()->PostDirectory() + "/" + M_Physics->Data()->PostFile() +
                    iter_variable->first + iter_suffix->first;

                Debug( 6310 ) << "[output2FileBuffers] setting output for file "
                              << file_output << ", writing variable "
                              << iter_variable->first << "\n";
                Debug( 6310 ) << "[output2FileBuffers] variable size = "
                              << (*M_U_thistime)[iter_variable->second].size()
                              << " variable norm2 = " << (*M_U_thistime)[iter_variable->second].Norm2() << "\n";

                (this->*(iter_suffix->second))
                    ( file_output, time_val, (*M_U_thistime)[iter_variable->second], iter_variable->first );
            }
        }
    }
}

void
OneDimensionalModel_Solver::resetFileBuffers()
{
    closeFileBuffers();
    openFileBuffers();
}

void
OneDimensionalModel_Solver::seekpFileBuffers()
{
    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator iter;

    for( iter = M_post_process_buffer.begin(); iter != M_post_process_buffer.end(); iter++ )
        (*iter).second->seekp( M_post_process_buffer_offset[(*iter).first] );
}

void
OneDimensionalModel_Solver::tellpFileBuffers()
{
    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator iter;

    for( iter = M_post_process_buffer.begin(); iter != M_post_process_buffer.end(); iter++ )
        M_post_process_buffer_offset[(*iter).first] = (*iter).second->tellp();
}

void
OneDimensionalModel_Solver::closeFileBuffers()
{
    // as I have a boost::shared_ptr, I expect the objects to be deallocated
    // now that the pointers are destroyed
    M_post_process_buffer.erase( M_post_process_buffer.begin(),
        M_post_process_buffer.end() );
    M_post_process_buffer_offset.erase( M_post_process_buffer_offset.begin(),
        M_post_process_buffer_offset.end() );
}

void
OneDimensionalModel_Solver::showMe( std::ostream& output ) const
{
    M_Physics->Data()->showMe( output );
}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_Solver::setProblem( const Physics_PtrType Physics,
                                        const Flux_PtrType    Flux,
                                        const Source_PtrType  Source )
{
    M_Physics = Physics;
    M_Flux    = Flux;
    M_Source  = Source;
}

void
OneDimensionalModel_Solver::setCommunicator( const Comm_PtrType Comm )
{
    M_Comm = Comm;
    M_Displayer.SetCommunicator( *Comm );
}

void
OneDimensionalModel_Solver::setFESpace( const FESpace_PtrType FESpace )
{
    M_FESpace = FESpace;

    //Id of left and right bc nodes
    M_leftNodeId          = 1;
    M_leftInternalNodeId  = M_leftNodeId  + 1;
    M_rightNodeId         = M_FESpace->dim();
    M_rightInternalNodeId = M_rightNodeId - 1;

    //Elementary Matrices
    M_elmatMass.reset(  new ElemMat ( M_FESpace->fe().nbNode, 1, 1 ) );
    M_elmatStiff.reset( new ElemMat ( M_FESpace->fe().nbNode, 1, 1 ) );
    M_elmatGrad.reset(  new ElemMat ( M_FESpace->fe().nbNode, 1, 1 ) );
    M_elmatDiv.reset(   new ElemMat ( M_FESpace->fe().nbNode, 1, 1 ) );

    //Vectors
    M_U_thistime->resize(   5, Vector_Type( M_FESpace->map() ) );
    M_U_prevtime.resize(   2, Vector_Type( M_FESpace->map() ) );
    M_U_2prevtime.resize(  2, Vector_Type( M_FESpace->map() ) );
    M_rhs.resize(          2, Vector_Type( M_FESpace->map() ) );
    M_FluxVector.resize(   2, Vector_Type( M_FESpace->map() ) );
    M_SourceVector.resize( 2, Vector_Type( M_FESpace->map() ) );

    //Matrix
    M_massMatrix.reset( new Matrix_Type( M_FESpace->map() ) );
    M_gradMatrix.reset( new Matrix_Type( M_FESpace->map() ) );
}

void
OneDimensionalModel_Solver::setLinearSolver( const LinearSolver_PtrType linearSolver )
{
    M_linearSolver = linearSolver;
}

void
OneDimensionalModel_Solver::setBCValuesLeft( const Real& bcL1, const Real& bcL2 )
{
    M_bcDirLeft[0] = bcL1;
    M_bcDirLeft[1] = bcL2;
}

void
OneDimensionalModel_Solver::setBCValuesRight( const Real& bcR1, const Real& bcR2 )
{
    M_bcDirRight[0] = bcR1;
    M_bcDirRight[1] = bcR2;
}

// ===================================================
// Get Methods
// ===================================================
const OneDimensionalModel_Solver::Solution_PtrType
OneDimensionalModel_Solver::U_thistime() const
{
    return M_U_thistime;
}

const OneDimensionalModel_Solver::Vector_Type&
OneDimensionalModel_Solver::U1_thistime() const
{
    return (*M_U_thistime)[0];
}

const OneDimensionalModel_Solver::Vector_Type&
OneDimensionalModel_Solver::U2_thistime() const
{
    return (*M_U_thistime)[1];
}

const OneDimensionalModel_Solver::Vector_Type&
OneDimensionalModel_Solver::W1_thistime() const
{
    return (*M_U_thistime)[2];
}

const OneDimensionalModel_Solver::Vector_Type&
OneDimensionalModel_Solver::W2_thistime() const
{
    return (*M_U_thistime)[3];
}

const OneDimensionalModel_Solver::Vector_Type&
OneDimensionalModel_Solver::P_thistime() const
{
    return (*M_U_thistime)[4];
}

const OneDimensionalModel_Solver::Physics_PtrType&
OneDimensionalModel_Solver::Physics() const
{
    return M_Physics;
}

const OneDimensionalModel_Solver::Flux_PtrType&
OneDimensionalModel_Solver::Flux() const
{
    return M_Flux;
}

const OneDimensionalModel_Solver::Source_PtrType&
OneDimensionalModel_Solver::Source() const
{
    return M_Source;
}

const UInt&
OneDimensionalModel_Solver::LeftNodeId() const
{
    return M_leftNodeId;
}

const UInt&
OneDimensionalModel_Solver::LeftInternalNodeId() const
{
    return M_leftInternalNodeId;
}

const UInt&
OneDimensionalModel_Solver::RightNodeId() const
{
    return M_rightNodeId;
}

const UInt&
OneDimensionalModel_Solver::RightInternalNodeId() const
{
    return M_rightInternalNodeId;
}

Container2D_Type
OneDimensionalModel_Solver::BCValuesLeft() const
{
    Container2D_Type temp;

    temp[0] = (*M_U_thistime)[0]( LeftNodeId() );
    temp[1] = (*M_U_thistime)[1]( LeftNodeId() );

    return temp;
}

Container2D_Type
OneDimensionalModel_Solver::BCValuesInternalLeft() const
{
    Container2D_Type temp;

    temp[0] = (*M_U_thistime)[0]( LeftInternalNodeId() );
    temp[1] = (*M_U_thistime)[1]( LeftInternalNodeId() );

    return temp;
}

Container2D_Type
OneDimensionalModel_Solver::BCValuesRight() const
{
    Container2D_Type temp;

    temp[0] = (*M_U_thistime)[0]( RightNodeId() );
    temp[1] = (*M_U_thistime)[1]( RightNodeId() );

    return temp;
}

Container2D_Type
OneDimensionalModel_Solver::BCValuesInternalRight() const
{
    Container2D_Type temp;

    temp[0] = (*M_U_thistime)[0]( RightInternalNodeId() );
    temp[1] = (*M_U_thistime)[1]( RightInternalNodeId() );

    return temp;
}

Real
OneDimensionalModel_Solver::value( std::string var, UInt pos ) const
{
    std::map< std::string, UInt>::const_iterator it = M_variable_index_map.find( var );

    return (*M_U_thistime)[it->second]( pos );
}

Real
OneDimensionalModel_Solver::BoundaryValue( const OneD_BC& bcType, const OneD_BCSide& bcSide ) const
{
    UInt boundaryDof;

    switch( static_cast<OneD_BCSide> (bcSide) )
    {
        case OneD_left:
            boundaryDof = 1;
        break;
        case OneD_right:
            boundaryDof = M_Flux->Physics()->Data()->nbElem() + 1;
        break;
        default:
            std::cout << "Warning: bcSide \"" << bcSide << "\" not available!" << std::endl;
            return 0;
        break;
    }

    switch( bcType )
    {
        case OneD_A:
            return (*M_U_thistime)[0]( boundaryDof );
        case OneD_Q:
            // Flow rate is positive with respect to the outgoing normal
            return (*M_U_thistime)[1]( boundaryDof ) * ( ( bcSide == OneD_left ) ? -1. : 1. );
        case OneD_W1:
            return (*M_U_thistime)[2]( boundaryDof );
        case OneD_W2:
            return (*M_U_thistime)[3]( boundaryDof );
        case OneD_P:
            // Pressure is positive with respect to the outgoing normal
            return (*M_U_thistime)[4]( boundaryDof ) * ( ( bcSide == OneD_left ) ? -1. : 1. );
        default:
            std::cout << "Warning: bcType \"" << bcType << "\"not available!" << std::endl;
            return 0.;
    }
}

// ===================================================
// Private Methods
// ===================================================
void
OneDimensionalModel_Solver::_updatePressure( const Real& TimeStep )
{
    Vector pressures(4*M_FESpace->dim());

    for ( UInt i(0); i < M_FESpace->dim() ; ++i )
    {
        ublas::subrange(pressures, 4*i, 4 + 4*i) = M_Physics->pressure( (*M_U_thistime) [0][i+1],
                                                                          M_U_prevtime  [0][i+1],
                                                                          M_U_2prevtime [0][i+1],
                                                                          i,
                                                                          TimeStep );
        (*M_U_thistime)[4][i + 1] = pressures(4*i);
    }

    if( M_Physics->Data()->viscoelasticWall() )
    {
        for ( UInt i(0); i < M_FESpace->dim() ; ++i )
        {
            (*M_U_thistime)[M_variable_index_map.find("P_elast")->second][i + 1] = pressures(1 + 4*i);
            (*M_U_thistime)[M_variable_index_map.find("P_visc")->second][i + 1]  = pressures(2 + 4*i);
            (*M_U_thistime)[M_variable_index_map.find("dA_dt")->second][i + 1]   = pressures(3 + 4*i);
        }
    }
}

void
OneDimensionalModel_Solver::_updateSystem( const Real& TimeStep )
{
    Chrono chrono; chrono.start();

    // Update the vector containing the values of the flux at the nodes and its jacobian
    _updateFluxDer();

    // Update the vector containing the values of the source term at the nodes and its jacobian
    _updateSourceDer();

    // Update the matrices for the non-linear terms
    _updateMatrices();

    // Taylor-Galerkin scheme: (explicit, U = [U1,U2]^T )
    // (Un+1, phi) =                                // massFactor^{-1} * Un+1
    // (Un, phi)                                    //            mass * U
    // + dt     * (       Fh(Un), dphi/dz )         //            grad * F(U)
    // - dt^2/2 * (diffFh(Un) Sh(Un), dphi/dz )     // gradDiffFlux(U) * S(U)
    // + dt^2/2 * (diffSh(Un) dFh/dz(Un), phi )     //   divDiffSrc(U) * F(U)
    // - dt^2/2 * (diffFh(Un) dFh/dz(Un), dphi/dz ) //stiffDiffFlux(U) * F(U)
    // - dt     * (       Sh(Un), phi )             //            mass * S(U)
    // + dt^2/2 * (diffSh(Un) Sh(Un), phi )         //  massDiffSrc(U) * S(U)

    // Compute M_rhs[0], M_rhs[1] (systems in U1=A, U2=Q)
    Real dt2over2 = TimeStep * TimeStep * 0.5;
    for( UInt i(0); i < 2; ++i )
    {
        // rhs = mass * Un
        M_rhs[i]  = (*M_massMatrix) * (*M_U_thistime)[i];
        M_rhs[i] += (*M_gradMatrix) * (TimeStep*M_FluxVector[i]);

        // rhs = rhs - dt * mass * S(Un)
        M_rhs[i] += (*M_massMatrix) * (- TimeStep*M_SourceVector[i]);

        for( UInt j(0); j<2; ++j )
        {
            // rhs = rhs - dt^2/2 * gradDiffFlux * S(Un)
            M_gradMatrixDiffFlux[2*i + j]->GlobalAssemble();
            M_rhs[i] += *M_gradMatrixDiffFlux[2*i + j]*(- dt2over2*M_SourceVector[j]);

            // rhs = rhs + dt^2/2 * divDiffSrc * F(Un)
            M_divMatrixDiffSrc[2*i + j]->GlobalAssemble();
            M_rhs[i] += *M_divMatrixDiffSrc[2*i + j]*(dt2over2*M_FluxVector[j]);

            // rhs = rhs - dt^2/2 * stiffDiffFlux * F(Un)
            M_stiffMatrixDiffFlux[2*i + j]->GlobalAssemble();
            M_rhs[i] += *M_stiffMatrixDiffFlux[2*i + j]*(- dt2over2*M_FluxVector[j]);

            // rhs = rhs + dt^2/2 * massDiffSrc * S(Un)
            M_massMatrixDiffSrc[2*i + j]->GlobalAssemble();
            M_rhs[i] += *M_massMatrixDiffSrc[2*i + j]*(dt2over2*M_SourceVector[j]);
        }
    }

    chrono.stop(); Debug( 6310 ) << "[timeAdvance] overall computation time " << chrono.diff() << " s.\n";
}

void
OneDimensionalModel_Solver::_updateFlux()
{
    Real Aii, Qii;

    for ( UInt ii = 1; ii <= M_FESpace->dim() ; ii++ )
    {
        Aii = (*M_U_thistime)[0]( ii );
        Qii = (*M_U_thistime)[1]( ii );
        M_FluxVector[0]( ii ) = ( *M_Flux )( Aii, Qii, 1, ii - 1 );
        M_FluxVector[1]( ii ) = ( *M_Flux )( Aii, Qii, 2, ii - 1 );
    }
}

void
OneDimensionalModel_Solver::_updateFluxDer()
{
    // first update the Flux vector
    _updateFlux();

    // then update the derivative of the Flux vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;
    UInt ii, iip1;

    for ( UInt ielem = 0; ielem < M_FESpace->dim() - 1; ielem++ )
    {
        // for P1Seg and appropriate mesh only!
        ii    = ielem;      // left node of current element
        iip1  = ielem + 1;  // right node of current element

        Aii   = (*M_U_thistime)[0]( ielem + 1 );
        Qii   = (*M_U_thistime)[1]( ielem + 1 );

        Aiip1 = (*M_U_thistime)[0]( ielem + 2 );
        Qiip1 = (*M_U_thistime)[1]( ielem + 2 );

        for( UInt ii=1; ii<3; ++ii )
        {
            for( UInt jj=1; jj<3; ++jj )
            {
                tmp  = M_Flux->diff(   Aii,   Qii, ii, jj, ii );
                tmp += M_Flux->diff( Aiip1, Qiip1, ii, jj, iip1 );

                M_diffFlux[ 2*(ii - 1) + jj - 1 ]( ielem ) = 0.5*tmp;
            }
        }
    }
}

void
OneDimensionalModel_Solver::_updateSource( )
{
    Real Aii, Qii;

    for ( UInt ii = 1; ii <= M_FESpace->dim() ; ii++ )
    {
        Aii = (*M_U_thistime)[0]( ii );
        Qii = (*M_U_thistime)[1]( ii );
        for(UInt k=0; k<2; ++k)
            M_SourceVector[k]( ii ) = ( *M_Source )( Aii, Qii, k+1, ii - 1);
    }
}

void
OneDimensionalModel_Solver::_updateSourceDer( )
{
    // first update the Source vector
    _updateSource();

    // then update the derivative of the Source vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;
    UInt ii, iip1;

    for ( UInt ielem=0; ielem < M_FESpace->dim() - 1 ; ielem++ )
    {
        // for P1Seg and appropriate mesh only!
        ii = ielem;        // left node of current element
        iip1 = ielem + 1;  // right node of current element

        Aii   = (*M_U_thistime)[0]( ielem + 1);
        Qii   = (*M_U_thistime)[1]( ielem + 1);
        Aiip1 = (*M_U_thistime)[0]( ielem + 2 );
        Qiip1 = (*M_U_thistime)[1]( ielem + 2 );

        for( UInt ii=1; ii<3; ++ii )
        {
            for( UInt jj=1; jj<3; ++jj )
            {
                tmp =  M_Source->diff(   Aii,   Qii, ii, jj, ii );
                tmp += M_Source->diff( Aiip1, Qiip1, ii, jj, iip1 );
                M_diffSrc[ 2*(ii - 1) + jj - 1 ]( ielem ) = 0.5 * tmp;
            }
        }
    }
}

void
OneDimensionalModel_Solver::_updateMatrices()
{
    //--------------------------------------------------------
    // Chrono chrono;
    // chrono.start();
    // std::cout << "o-loop over the matrices INIT... ";

    // Matrices initialization
    for( UInt i = 0; i < 4; ++i )
    {
//             *M_massMatrixDiffSrc[i]   *= 0.;
//             *M_stiffMatrixDiffFlux[i] *= 0.;
//             *M_gradMatrixDiffFlux[i]  *= 0.;
//             *M_divMatrixDiffSrc[i]    *= 0.;
        M_massMatrixDiffSrc[i].reset(new Matrix_Type( M_FESpace->map() ));
        //M_stiffMatrixDiffFlux.push_back( new Matrix_Type( M_FESpace->map() ));
        M_stiffMatrixDiffFlux[i].reset( new Matrix_Type( M_FESpace->map() ));
        //M_gradMatrixDiffFlux. push_back( new Matrix_Type( M_FESpace->map() ));
        M_gradMatrixDiffFlux[i].reset( new Matrix_Type( M_FESpace->map() ));
        //M_divMatrixDiffSrc.   push_back( new Matrix_Type( M_FESpace->map() ));
        M_divMatrixDiffSrc[i].reset( new Matrix_Type( M_FESpace->map() ));
    }

    /*
      chrono.stop();
      std::cout << "done in " << chrono.diff() << " s." << std::endl;
      chrono.start();
      std::cout << "o-loop over the matrices... ";
    */

    // Elementary computation and matrix assembling
    // Loop on elements
    for(UInt iedge = 1; iedge <= M_FESpace->mesh()->numEdges(); iedge++)
    {
        // update the current element
        M_FESpace->fe().updateFirstDerivQuadPt(M_FESpace->mesh()->edgeList(iedge));
        //  std::cout << M_FESpace->fe().currentId() << std::endl;

        for(UInt ii = 1; ii <= 2; ii ++)
        {
            for(UInt jj = 1; jj <= 2; jj ++)
            {
                // update the M_coeff*
                _updateMatrixCoefficients( ii , jj, iedge);

                // update the M_elmat*
                _updateElemMatrices();

                // assemble the global matrices
                _assemble_matrices( ii, jj );
            }
        }
    } // end loop on elements

    /*
      chrono.stop();
      std::cout << "done in " << chrono.diff() << " s." << std::endl;
      chrono.start();
      std::cout << "o-loop over the matrices BC DIR... ";
    */
}

void
OneDimensionalModel_Solver::_updateMatrixCoefficients( const UInt& ii, const UInt& jj , const UInt& iedge )
{
    Real dFluxdUelem(0), dSrcdUelem(0);

    ASSERT_BD( 0 < ii && ii < 3 && 0 < jj && jj < 3 );

    dFluxdUelem = M_diffFlux[ 2*(ii-1) + jj-1 ]( iedge - 1 ); // iedge starts from 1...
    dSrcdUelem  = M_diffSrc [ 2*(ii-1) + jj-1 ]( iedge - 1 );

    M_coeffGrad  = dFluxdUelem; // term gradDiffFlux(U) [* S(U)]
    M_coeffDiv   = dSrcdUelem;  // term  divDiffSrc(U) [* F(U)]
    M_coeffStiff = dFluxdUelem; // term stiffDiffFlux(U) [* F(U)]
    M_coeffMass  = dSrcdUelem;  // term  massDiffSrc(U) [* S(U)]
}

void
OneDimensionalModel_Solver::_updateElemMatrices()
{
    // set the elementary matrices to 0.
    M_elmatMass->zero();
    M_elmatStiff->zero();
    M_elmatGrad->zero();
    M_elmatDiv->zero();

    // update the mass matrix
    mass( M_coeffMass, *M_elmatMass, M_FESpace->fe(),0, 0 );
    //  M_elmatMass->showMe( std::cout );

    // update the stiffness matrix
    stiff( M_coeffStiff, *M_elmatStiff, M_FESpace->fe(),0 ,0 );
    // std::cout << "Elem Stiff matrix :" << std::endl;
    // M_elmatStiff->showMe( std::cout );

    /*! update the gradient matrix
      gradient operator:
      grad_{ij} = \int_{fe} coeff \phi_j \frac{d \phi_i}{d x}

      BEWARE :
      \param 0: the first argument "0" corresponds to the first
      and only coordinate (1D!), and HERE it starts from 0... (Damm'!)

      \param - M_coeffGrad: the sign "-" in the second argument
      is added to correspond to the described operator.
      (There is a minus in the elemOper implementation).
    */
    grad( 0 , - M_coeffGrad, *M_elmatGrad, M_FESpace->fe(), M_FESpace->fe(), 0, 0 );
    //  std::cout << "Elem Grad matrix :" << std::endl;
    //  M_elmatGrad->showMe( std::cout );

    /*! update the divergence matrix
      divergence operator: (transpose of the gradient)
      div_{ij} = \int_{fe} coeff \frac{d \phi_j}{d x} \phi_i

      \note formally this M_elmatDiv is not necessary
      as it is the transpose of the M_elmatGrad.
      But for the sake of clarity, I prefer to keep it. (low cost!)

      BEWARE : same remarks as grad (see above).
    */
    div( 0 , - M_coeffDiv, *M_elmatDiv, M_FESpace->fe(), M_FESpace->fe(), 0, 0 );
    //  std::cout << "Elem Div matrix :" << std::endl;
    //  M_elmatDiv->showMe( std::cout );
}

int
OneDimensionalModel_Solver::_assemble_matrices(const UInt& ii, const UInt& jj )
{
    ASSERT_BD( 0 < ii && ii < 3 && 0 < jj && jj < 3 );

    // assemble the mass matrix
    assemb_mat( *M_massMatrixDiffSrc[ 2*(ii-1) + jj-1 ]  , *M_elmatMass, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );

    // assemble the stiffness matrix
    assemb_mat( *M_stiffMatrixDiffFlux[ 2*(ii-1) + jj-1 ], *M_elmatStiff, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );

    // assemble the gradient matrix
    assemb_mat( *M_gradMatrixDiffFlux[ 2*(ii-1) + jj-1 ] , *M_elmatGrad, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );

    // assemble the divergence matrix
    assemb_mat( *M_divMatrixDiffSrc[ 2*(ii-1) + jj-1 ]   , *M_elmatDiv, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );

    return 0;

    ERROR_MSG("Invalid values for the _assemble_matrices method.");
    return 1;
}

void
OneDimensionalModel_Solver::_updateBCDirichletVector()
{
    UInt firstDof = M_leftNodeId;
    UInt lastDof  = M_rightNodeId;

    Debug( 6310 ) << "[_updateBCDirichletVector] \t firstDof = " << firstDof
                  << " lastDof = " << lastDof << ";\n";

    Debug( 6310 ) << "[_updateBCDirichletVector] \t bcDirLeft[0] = " << M_bcDirLeft[0]
                  << " bcDirLeft[1] = " << M_bcDirLeft[1] << ";\n";

    Debug( 6310 ) << "[_updateBCDirichletVector] \t bcDirRight[0] = " << M_bcDirRight[0]
                  << " bcDirRight[1] = " << M_bcDirRight[1] << ";\n";

    // unsymmetric treatment (LU must be used!)
    // first row modified
    M_rhs[0]( firstDof ) = M_bcDirLeft[0];
    M_rhs[1]( firstDof ) = M_bcDirLeft[1];

    // last row modified
    M_rhs[0]( lastDof ) = M_bcDirRight[0];
    M_rhs[1]( lastDof ) = M_bcDirRight[1];

    // symmetric treatment (cholesky can be used)
//     for(UInt i=0; i<2; ++i) {
//         // first row modified (Dirichlet)
//         M_rhs[i]( firstDof ) = M_bcDirLeft[i];
//         // second row modified (for symmetry)
//         M_rhs[i]( firstDof + 1 ) += - M_massMatrix.LowDiag()( firstDof ) * M_bcDirLeft[i];
//         // last row modified (Dirichlet)
//         M_rhs[i]( lastDof ) = M_bcDirRight[i];
//         // penultimate row modified (for symmetry)
//         M_rhs[i]( lastDof - 1 ) += - M_massMatrix.UpDiag()( lastDof - 1 ) * M_bcDirRight[i];
    // }
}

void
OneDimensionalModel_Solver::_updateBCDirichletMatrix( Matrix_Type& mat )
{
    UInt firstDof = M_leftNodeId;
    UInt lastDof  = M_rightNodeId;

    // unsymmetric treatment (LU must be used!)
    // modify the first row

    mat.diagonalize(firstDof - 1, 1, 0);

//     mat.set_mat_inc(firstDof, firstDof    , 1. );
//     mat.set_mat_inc(firstDof, firstDof + 1, 0. );

    // modify the last row

    mat.diagonalize(lastDof - 1, 1, 0);
//     mat.set_mat_inc(lastDof,  lastDof      , 1.);
//     mat.set_mat_inc(lastDof,  lastDof  - 1 , 0.);

    // symmetric treatment (cholesky can be used)
    // modify the first row
//     mat.Diag()( firstDof )    = 1.;
//     mat.UpDiag()( firstDof )  = 0.;
//     mat.LowDiag()( firstDof ) = 0.; //and second row

//     // modify the last row
//     mat.Diag()( lastDof )      = 1.;
//     mat.UpDiag()( lastDof-1 )  = 0.; //and penultimate row
//     mat.LowDiag()( lastDof-1 ) = 0.;
}

OneDimensionalModel_Solver::Vector_Type
OneDimensionalModel_Solver::_correct_flux_inertial( const Vector_Type& flux )
{
    Matrix_Type _matrixLHS(M_FESpace->map());
    Matrix_Type _stiffRHS (M_FESpace->map());

    ElemMat _elmatMassLHS  (M_FESpace->fe().nbNode, 1, 1);
    ElemMat _elmatStiffLHS (M_FESpace->fe().nbNode, 1, 1);
    ElemMat _elmatStiffRHS (M_FESpace->fe().nbNode, 1, 1);

    Vector_Type _rhs(M_FESpace->map());

    Real _coeffMass;
    Real _coeffStiff;

    Real m, meanA0;

    //    std::ostringstream output;

    _matrixLHS *= 0.;
     _stiffRHS  *= 0.;

    // Elementary computation and matrix assembling
    // Loop on elements
    for(UInt iedge = 1; iedge <= M_FESpace->mesh()->numEdges(); iedge++)
    {
        // set the elementary matrices to 0.
        _elmatMassLHS. zero();
        _elmatStiffLHS.zero();
        _elmatStiffRHS.zero();

        _coeffMass  = M_rhs[0]( iedge - 1 ) + M_rhs[0]( iedge );
        _coeffMass /= 2;
        _coeffMass  = 1./_coeffMass;

        meanA0  = M_Physics->Data()->Area0(iedge - 1) + M_Physics->Data()->Area0(iedge);
        meanA0 /= 2;

        m = M_Physics->Data()->DensityWall()*M_Physics->Data()->Thickness()/
            ( 2*std::sqrt(4*std::atan(1))*std::sqrt(meanA0) );

        _coeffStiff = m/M_Physics->Data()->DensityRho();

        // update the current element
        M_FESpace->fe().updateFirstDerivQuadPt(M_FESpace->mesh()->edgeList(iedge));

        mass (   _coeffMass,  _elmatMassLHS,  M_FESpace->fe(), 0, 0 );
        stiff(   _coeffStiff, _elmatStiffLHS, M_FESpace->fe(), 0, 0 );
        stiff( - _coeffStiff, _elmatStiffRHS, M_FESpace->fe(), 0, 0 );

        // assemble the mass and grad matrices
        assemb_mat( _matrixLHS, _elmatMassLHS, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );
        assemb_mat( _matrixLHS, _elmatStiffLHS, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );

        assemb_mat( _stiffRHS, _elmatStiffRHS, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );

        Debug( 6310 ) << "\n\tm = "           << m
                      << "\n\t_coeffMass = "  << _coeffMass
                      << "\n\t_coeffStiff = " << _coeffStiff << "\n";
    } // end loop on elements

    // update rhs
    //_stiffRHS.Axpy( 1., flux , 0., _rhs );

    _rhs = _stiffRHS*flux;

    UInt firstDof = 1;
    UInt lastDof  = _rhs.size();

    // symmetric treatment (cholesky can be used)
    // first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    // second row modified (for symmetry)
    //  _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    // last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    // penultimate row modified (for symmetry)
    //  _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    _updateBCDirichletMatrix(_matrixLHS);

    //_tridiagsolver.Factor( _matrixLHS );

    // cholesky or lapack lu solve
    // solve the system: rhs1 = massFactor^{-1} * rhs1
    //_tridiagsolver.Solve( _matrixLHS, _rhs );

    Vector_Type _sol(_rhs);

    M_linearSolver->setMatrix(_matrixLHS);
    //@int numIter = M_linearSolver->solveSystem( _rhs, _sol, _matrixLHS, true);

    //std::cout <<" iterations number :  " << numIter << std::endl;

    return _sol;
}

OneDimensionalModel_Solver::Vector_Type
OneDimensionalModel_Solver::_correct_flux_viscoelastic( const Vector_Type& flux, const Real& TimeStep )
{
    Matrix_Type _matrixLHS(M_FESpace->map());
    Matrix_Type _stiffRHS (M_FESpace->map());

    //TriDiagCholesky< Real, Matrix_Type, Vector > _tridiagsolver(M_FESpace->dim());

    ElemMat _elmatMassLHS  (M_FESpace->fe().nbNode,1,1);
    ElemMat _elmatStiffLHS (M_FESpace->fe().nbNode,1,1);
    ElemMat _elmatStiffRHS (M_FESpace->fe().nbNode,1,1);

    Vector_Type _rhs(M_FESpace->map());

    Real _coeffMass;
    Real _coeffStiff;

    Real gamma, meanA0(1.);

    _matrixLHS *= 0.;
    _stiffRHS  *= 0.;

    // Elementary computation and matrix assembling
    // Loop on elements
    for(UInt iedge = 1; iedge <= M_FESpace->mesh()->numEdges(); iedge++)
    {
        // set the elementary matrices to 0.
        _elmatMassLHS.zero();
        _elmatStiffLHS.zero();
        _elmatStiffRHS.zero();

        // this comes from the exact derivation of generalized string model
        // + Voigt viscoelasticity
        _coeffMass = M_rhs[0]( iedge - 1 ) + M_rhs[0]( iedge );
        _coeffMass *= 0.5;
        _coeffMass = 1./std::sqrt(_coeffMass);

        if(M_Physics->Data()->linearizeStringModel())
        {
            // this is to recover the linearized version (_coeffMass = 1/A)
            _coeffMass *= _coeffMass;

            meanA0  = M_Physics->Data()->Area0( iedge - 1 ) + M_Physics->Data()->Area0( iedge );
            meanA0 *= 0.5;

            if(M_Physics->Data()->linearizeEquations())
            {
                // when using linearized equations, A \simeq A0
                _coeffMass = 1./ meanA0;
            }
        }
        gamma = M_Physics->Data()->ViscoelasticModulus() / ( 2 * std::sqrt(4*std::atan(1)) );
        gamma *= 1 / std::sqrt( meanA0 );
        _coeffStiff = TimeStep * gamma / M_Physics->Data()->DensityRho();

        // update the current element
        M_FESpace->fe().updateFirstDerivQuadPt(M_FESpace->mesh()->edgeList(iedge));
        //std::cout << M_FESpace->fe().currentId() << std::endl;

        mass (   _coeffMass,  _elmatMassLHS,  M_FESpace->fe(),0, 0 );
        stiff(   _coeffStiff, _elmatStiffLHS, M_FESpace->fe(),0, 0 );
        stiff( - _coeffStiff, _elmatStiffRHS, M_FESpace->fe(),0, 0 );

        // assemble the mass and grad matrices
        assemb_mat( _matrixLHS, _elmatMassLHS, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );
        assemb_mat( _matrixLHS, _elmatStiffLHS, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );

        assemb_mat( _stiffRHS, _elmatStiffRHS, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );

        Debug( 6310 ) << "\n\tgamma = " << gamma
                      << "\n\t_coeffMass = " << _coeffMass
                      << "\n\t_coeffStiff = " << _coeffStiff << "\n";

    } // end loop on elements

    // update rhs
    // rhs = _stiffRHS * rhs
    _rhs = _stiffRHS*flux;

    // NOTE I should add to rhs the value of boundary integral ("neumann-like")
    // BUT it's useless since I'm going to impose dirichlet conditions on all boundaries

    UInt firstDof = 1;
    UInt lastDof  = _rhs.size();

    // symmetric treatment (cholesky can be used)
    // first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    // second row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    // _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    // last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    // penultimate row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    _updateBCDirichletMatrix(_matrixLHS);

    Vector_Type _sol(_rhs);

    M_linearSolver->setMatrix(_matrixLHS);
    //int numIter = M_linearSolver->solveSystem( _rhs, _sol, _matrixLHS, true);

    return _sol;
}

OneDimensionalModel_Solver::Vector_Type
OneDimensionalModel_Solver::_correct_flux_longitudinal( )
{
    Matrix_Type _massLHS(M_FESpace->map());
    Matrix_Type _massRHS(M_FESpace->map());

    //TriDiagCholesky< Real, Matrix_Type, Vector > _tridiagsolver(M_FESpace->dim());

    ElemMat _elmatMassLHS (M_FESpace->fe().nbNode,1,1);
    ElemMat _elmatMassRHS (M_FESpace->fe().nbNode,1,1);

    Vector_Type _rhs(M_FESpace->map());
    // let g = sqrt(A) - sqrt(A0)
    // now f = _d3g_dz3
    Vector_Type _g(M_FESpace->map());
    Vector_Type _f(M_FESpace->map());

    //          _g = M_rhs[0];
    for( UInt i=0; i<M_FESpace->dim(); ++i )
        _g(i) = std::sqrt(M_rhs[0](i)) - std::sqrt(M_Physics->Data()->Area0(i));

    UInt inode;

    Real _coeffMassLHS;
    Real _coeffMassRHS;

    Real _a;
    //    std::ostringstream output;
    _massLHS *= 0.;
    _massRHS *= 0.;

    Real _h( M_Physics->Data()->Length() / static_cast<Real>(M_Physics->Data()->nbElem() - 1) );

    // Elementary computation and matrix assembling
    // Loop on elements
    for(UInt iedge = 1; iedge <= M_FESpace->mesh()->numEdges(); iedge++)
    {
        inode = iedge - 1;

        // set the elementary matrices to 0.
        _elmatMassLHS.zero();
        _elmatMassRHS.zero();

        // coeff (1/A) (average values over the element)
        _coeffMassLHS = M_rhs[0]( inode ) + M_rhs[0]( inode+1 );
        _coeffMassLHS /= 2;
        _coeffMassLHS = 1./_coeffMassLHS;

        _a = M_Physics->Data()->InertialModulus() / std::sqrt(4*std::atan(1));
        _coeffMassRHS = M_Physics->Data()->dataTime()->getTimeStep() * _a / M_Physics->Data()->DensityRho();

        // backward differentiation when near to the left boundary
        // d3 A (xi)/ dz3 = (1/h^3) * ( -A(xi) + A(xi+3) - 3A(xi+2) + 3A(xi+1) )

        // forward differentiation when near to the right boundary
        // d3 A (xi)/ dz3 = (1/h^3) * ( A(xi) - A(xi-3) + 3A(xi-2) - 3A(xi-1) )

        // central differentiation otherwise
        // d3 A (xi)/ dz3 = (1/h^3) * ( -A(xi-2) + 2A(xi-1) - 2A(xi+1) + A(xi+2) )

        Debug( 6310 ) << "\ninode = " << inode << "\n";
        if(inode<2)
        { // backward differentiation
            _f( inode ) = -_g( inode ) + _g( inode+3 )
                - 3*_g( inode+2 ) + 3*_g( inode+1 );
            Debug( 6310 ) << "\n\tbackward differentiation = " << _coeffMassLHS << "\n";
        }
        else if(inode>M_FESpace->mesh()->numEdges()-2)
        { // forward differentiation
            _f( inode ) = _g( inode ) - _g( inode-3 )
                + 3*_g( inode-2 ) - 3*_g( inode-1 );
            Debug( 6310 ) << "\n\forward differentiation = " << _coeffMassLHS << "\n";
        }
        else
        { // central differentiation
            _f( inode ) = -_g( inode - 2 ) + 2*_g( inode-1 )
                - 2*_g( inode+1 ) + _g( inode+2 );
            Debug( 6310 ) << "\n\tcentral differentiation = " << _coeffMassLHS << "\n";
        }

        _f(inode) *= 1/(2*std::pow(_h,3));

        // update the current element
        M_FESpace->fe().updateFirstDerivQuadPt(M_FESpace->mesh()->edgeList(iedge));

        mass( _coeffMassLHS, _elmatMassLHS, M_FESpace->fe(),0, 0 );
        mass( _coeffMassRHS, _elmatMassRHS, M_FESpace->fe(),0, 0 );

        // assemble the mass and grad matrices
        assemb_mat( _massLHS, _elmatMassLHS, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );
        assemb_mat( _massRHS, _elmatMassRHS, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );

        Debug( 6310 ) << "\n\t_coeffMassLHS = " << _coeffMassLHS << "\n";
        Debug( 6310 ) << "\n\t_coeffMassRHS = " << _coeffMassRHS << "\n";

    } // end loop on elements

    // update rhs
    //    _massLHS.Axpy( 1., (*M_U_thistime)[4] , 0., _rhs );
    _rhs = _massRHS*_f;

    UInt firstDof = 1;
    UInt lastDof  = _rhs.size();

    // symmetric treatment (cholesky can be used)
    // first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    // second row modified (for symmetry)
    //    _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    // last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    // penultimate row modified (for symmetry)
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    _updateBCDirichletMatrix(_massLHS);

    //@_tridiagsolver.Factor( _massLHS );

    // cholesky or lapack lu solve
    // solve the system: rhs1 = massFactor^{-1} * rhs1
    //@_tridiagsolver.Solve( _massLHS, _rhs );

    return _rhs;
}

/*
template< class Params, class Flux, class Source >
ScalVec
OneDimensionalModel_Solver::_compute_d2Q_dx2( const ScalVec& flux )
{
    Matrix_Type _massLHS (M_FESpace->dim());
    Matrix_Type _stiffRHS(M_FESpace->dim());

    TriDiagCholesky< Real, Matrix_Type, Vector > _tridiagsolver(M_FESpace->dim());

    ElemMat _elmatMassLHS (M_FESpace->fe().nbNode,1,1);
    ElemMat _elmatStiffRHS (M_FESpace->fe().nbNode,1,1);

    ScalVec _rhs(M_FESpace->dim());

    _massLHS  *= 0.;
    _stiffRHS *= 0.;

    // Elementary computation and matrix assembling
    // Loop on elements
    for(UInt iedge = 1; iedge <= M_FESpace->mesh()->numEdges(); iedge++){

        // set the elementary matrices to 0.
        _elmatMassLHS *= 0.;
        _elmatStiffRHS *= 0.;

        // update the current element
        M_FESpace->fe().updateFirstDerivQuadPt(M_FESpace->mesh()->edgeList(iedge));
        //std::cout << M_FESpace->fe().currentId() << std::endl;

        mass( 1., _elmatMassLHS, M_FESpace->fe(),0, 0 );
        stiff( -1., _elmatStiffRHS, M_FESpace->fe(),0, 0 );

        // assemble the mass and grad matrices
        assemb_mat( _massLHS, _elmatMassLHS, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );
        assemb_mat( _stiffRHS, _elmatStiffRHS, M_FESpace->fe(), M_FESpace->dof() , 0, 0 );

    } // end loop on elements

    // update rhs
    _rhs = _stiffRHS*flux;

    UInt firstDof = 0;
    UInt lastDof  = _rhs.size()-1;

    // symmetric treatment (cholesky can be used)
    // first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    // second row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    // _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    // last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    // penultimate row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    _updateBCDirichletMatrix(_massLHS);

//     Debug( 6310 ) << "\n[_compute_d2Q_dx2] _massLHS.Diag():\n";
//     for(int i=0; i<_massLHS.OrderMatrix(); ++i)
//         Debug( 6310 ) << "\t" << _massLHS.Diag()(i);
//     Debug( 6310 ) << "\n[_compute_d2Q_dx2] _massLHS.UpDiag():\n";
//     for(int i=0; i<_massLHS.OrderMatrix()-1; ++i)
//         Debug( 6310 ) << "\t" << _massLHS.UpDiag()(i);
//     Debug( 6310 ) << "\n[_compute_d2Q_dx2] _massLHS.LowDiag():\n";
//     for(int i=0; i<_massLHS.OrderMatrix()-1; ++i)
//         Debug( 6310 ) << "\t" << _massLHS.LowDiag()(i);
//     Debug( 6310 ) << "\n";

//     _tridiagsolver.Factor( _massLHS );

//     Debug( 6310 ) << "\n[_compute_d2Q_dx2] solving with rhs:\n";
//     for(int i=0; i<_massLHS.OrderMatrix(); ++i)
//         Debug( 6310 ) << "\t" << _rhs(i);
//     Debug( 6310 ) << "\n";

    // cholesky or lapack lu solve
    // solve the system: rhs1 = massFactor^{-1} * rhs1
    _tridiagsolver.Solve( _massLHS, _rhs );

    return _rhs;
}
*/

}
