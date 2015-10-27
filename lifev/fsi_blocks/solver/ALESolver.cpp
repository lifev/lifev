
#include <lifev/fsi_blocks/solver/ALESolver.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================


ALESolver::ALESolver ( FESpace<mesh_Type, MapEpetra>& mmFESpace,
                          boost::shared_ptr<Epetra_Comm>    comm ) :
    M_FESpace               ( mmFESpace ),
    M_localMap              ( M_FESpace.map() ),
    M_matrHE                ( new matrix_Type (M_localMap ) ),
    M_displayer              ( comm ),
    M_me                    ( comm->MyPID() ),
    M_verbose               ( M_me == 0 ),
    M_elmat                 ( M_FESpace.fe().nbFEDof(), nDimensions, nDimensions ),
    M_disp                  ( ),
    M_secondRHS             ( ),
    M_diffusion             ( 1. ),
    M_offset                (0)
{
}

ALESolver::ALESolver ( FESpace<mesh_Type, MapEpetra>& mmFESpace,
                          boost::shared_ptr<Epetra_Comm>              comm ,
                          MapEpetra& localMap,
                          UInt offset) :
    M_FESpace               ( mmFESpace ),
    M_localMap              ( localMap),
    M_matrHE                ( new matrix_Type (M_localMap ) ),
    M_displayer              ( comm ),
    M_me                    ( comm->MyPID() ),
    M_verbose               ( M_me == 0 ),
    M_elmat                 ( M_FESpace.fe().nbFEDof(), nDimensions, nDimensions ),
    M_secondRHS             ( ),
    M_diffusion             ( 1. ),
    M_offset                (offset)
{
}

// ===================================================
// Methods
// ===================================================


void ALESolver::setUp ( const GetPot& dataFile )
{
    M_diffusion = dataFile ("mesh_motion/diffusion", 1.0);

    computeMatrix( );

    M_secondRHS.reset (new vector_Type (M_FESpace.map() ) );
    M_disp.reset (new vector_Type (M_FESpace.map() ) );
}


void ALESolver::iterate ( BCHandler& BCh )
{
    LifeChrono chrono;

    // matrix and vector assembling communication
    M_displayer.leaderPrint (" ALE-  Applying boundary conditions ...         ");
    chrono.start();

    *M_secondRHS *= 0.;
    applyBoundaryConditions (*M_secondRHS, BCh);

    chrono.stop();
    M_displayer.leaderPrintMax ("done in " , chrono.diff() );
}

void ALESolver::applyBoundaryConditions (vector_Type& rhs, BCHandler& BCh)
{
    // CHANGED BY S. QUINODOZ !
    // "if" exchanged

    if (  ! BCh.bcUpdateDone() )
    {
        // BC boundary information update
        BCh.bcUpdate ( *M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );
    }

    if (M_offset) //mans that this is the fullMonolithic case
    {
        BCh.setOffset (M_offset);
    }
    else
    {
        bcManageRhs (rhs, *M_FESpace.mesh(), M_FESpace.dof(), BCh, M_FESpace.feBd(), 1., 0.0);
    }

    bcManageMatrix ( *M_matrHE, *M_FESpace.mesh(), M_FESpace.dof(), BCh, M_FESpace.feBd(), 1.0, 0. );
}

void ALESolver::applyBoundaryConditions (BCHandler& BCh)
{

    if (  ! BCh.bcUpdateDone() )
    {
    	BCh.bcUpdate ( *M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );
    }

    bcManageMatrix ( *M_matrHE_BC, *M_FESpace.mesh(), M_FESpace.dof(), BCh, M_FESpace.feBd(), 1.0, 0. );
}

void ALESolver::computeMatrix( )
{
    LifeChrono chrono;
    chrono.start();
    M_displayer.leaderPrint (" ALE-  Computing constant matrices ...          ");

    M_matrHE.reset ( new matrix_Type (M_localMap ) );

    UInt totalDof   = M_FESpace.dof().numTotalDof();
    // Loop on elements
    for ( UInt i = 0; i < M_FESpace.mesh()->numVolumes(); ++i )
    {
        // Updating derivatives
        M_FESpace.fe().updateFirstDerivQuadPt ( M_FESpace.mesh()->volumeList ( i ) );
        M_elmat.zero();
        AssemblyElemental::stiff ( M_diffusion, M_elmat, M_FESpace.fe(), 0, 0, 3 );
        // Assembling
        for ( UInt j = 0; j < M_FESpace.fieldDim(); ++j )
        {
            assembleMatrix ( *M_matrHE, M_elmat, M_FESpace.fe(), M_FESpace.dof(), j, j, j * totalDof + M_offset, j * totalDof + M_offset );
        }
    }

    M_matrHE->globalAssemble();

    M_matrHE_BC.reset ( new matrix_Type (M_localMap ) );
    M_matrHE_BC->zero();
    *M_matrHE_BC += *M_matrHE;

    chrono.stop();
    M_displayer.leaderPrintMax ("done in " , chrono.diff() );

}

void ALESolver::updateShapeDerivatives ( Real&                          alpha,
                                         const Real&                    density,
                                         const Real&                    viscosity,
                                         const vector_Type&             un,
                                         const vector_Type&             uk,
                                         //const vector_Type&           disp,
                                         const vector_Type&             w,
                                         FESpace<mesh_Type, MapEpetra>& velocityFESpace,
                                         FESpace<mesh_Type, MapEpetra>& pressureFESpace,
                                         bool                           wImplicit,
                                         bool                           convectiveTermDerivative,
                                         BCHandler& 					BCh )

{
    LifeChrono chrono;
    chrono.start();
    M_displayer.leaderPrint (" ALE-  Computing shapeDerivatives ...          ");

    UInt numVelocityComponent = nDimensions;

    // Velocity + Pressure. Probably needs also FESpaces for u and p
    M_matrShapeDerVel.reset ( new matrix_Type ( velocityFESpace.map() ) );
    M_matrShapeDerPressure.reset ( new matrix_Type ( pressureFESpace.map() ) );

    VectorElemental elementConvectionVelocity( velocityFESpace.fe().nbFEDof(), nDimensions );
    VectorElemental elementMeshVelocity      ( velocityFESpace.fe().nbFEDof(), nDimensions );
    VectorElemental u_loc                    ( velocityFESpace.fe().nbFEDof(), nDimensions );
    VectorElemental elementPressure          ( pressureFESpace.fe().nbFEDof(), 1 );
    VectorElemental elementVelocity          ( velocityFESpace.fe().nbFEDof(), nDimensions );


        vector_Type unRepeated  ( un  , Repeated );
        vector_Type ukRepeated  ( uk  , Repeated );
        // vector_Type dispRepeated( disp, Repeated );
        vector_Type wRepeated   ( w   , Repeated );
        // vector_Type dwRepeated  ( dw  , Repeated );


        for ( UInt i = 0; i < this->M_FESpace.mesh()->numVolumes(); i++ )
        {

            pressureFESpace.fe().update ( pressureFESpace.mesh()->volumeList ( i ) );
            velocityFESpace.fe().updateFirstDerivQuadPt ( velocityFESpace.mesh()->volumeList ( i ) );
            pressureFESpace.fe().updateFirstDerivQuadPt ( velocityFESpace.mesh()->volumeList ( i ) );

            // just to provide the id number in the assem_mat_mixed
            //pressureFESpace.fe().updateFirstDeriv( velocityFESpace.mesh()->volumeList( i ) );
            //as updateFirstDer
            //velocityFESpace.fe().updateFirstDeriv( velocityFESpace.mesh()->volumeList( i ) );
            M_FESpace.fe().updateFirstDerivQuadPt ( M_FESpace.mesh()->volumeList ( i ) );

            // initialization of elementary vectors
            boost::shared_ptr<MatrixElemental> elementMatrixPressure ( new MatrixElemental ( pressureFESpace.fe().nbFEDof(),
                                                                       1,
                                                                       0,
                                                                       M_FESpace.fe().nbFEDof(),
                                                                       0,
                                                                       nDimensions ) );
            boost::shared_ptr<MatrixElemental> elementMatrixVelocity ( new MatrixElemental ( velocityFESpace.fe().nbFEDof(),
                                                                       nDimensions,
                                                                       0,
                                                                       velocityFESpace.fe().nbFEDof(),
                                                                       0,
                                                                       nDimensions ) );
            boost::shared_ptr<MatrixElemental> elementMatrixConvective;

            if ( convectiveTermDerivative )
            {
                elementMatrixConvective.reset ( new MatrixElemental ( velocityFESpace.fe().nbFEDof(),
                                                                      nDimensions,
                                                                      0,
                                                                      velocityFESpace.fe().nbFEDof(),
                                                                      0,
                                                                      nDimensions ) );
                elementMatrixConvective->zero();
            }

            elementMatrixPressure->zero();
            elementMatrixVelocity->zero();

            for ( UInt iNode = 0 ; iNode < velocityFESpace.fe().nbFEDof() ; iNode++ )
            {
                UInt iLocal = velocityFESpace.fe().patternFirst ( iNode ); // iLocal = iNode

                for ( UInt iComponent = 0; iComponent < numVelocityComponent; ++iComponent )
                {
                    UInt iGlobal = velocityFESpace.dof().localToGlobalMap ( i, iLocal ) + iComponent * velocityFESpace.dim();

                    // if(!wImplicit)
                    // u^n - w^iNode local
                    elementConvectionVelocity.vec() [ iLocal + iComponent * velocityFESpace.fe().nbFEDof() ] = unRepeated (iGlobal)
                            - wRepeated ( iGlobal );
                    // else
                    // u^n - w^iNode local
                    // elementConvectionVelocity.vec() [ iLocal + iComponent*velocityFESpace.fe().nbFEDof() ] = ukRepeated(iGlobal)
                    // - wRepeated(iGlobal);
                    // w^iNode local
                    elementMeshVelocity.vec( )  [ iLocal + iComponent * velocityFESpace.fe().nbFEDof() ] = wRepeated ( iGlobal );
                    // u^iNode local
                    elementVelocity.vec( ) [ iLocal + iComponent * velocityFESpace.fe().nbFEDof() ] = ukRepeated ( iGlobal );
                    // dw local
                    //M_elementDisplacement.vec( ) [ iLocal + iComponent*velocityFESpace.fe().nbFEDof() ] = dispRepeated( iGlobal );
                    // dw local
                    //elementVelocityRightHandSide.vec( ) [ iLocal + iComponent*velocityFESpace.fe().nbFEDof() ] = dwRepeated( iGlobal );
                    // un local
                    u_loc.vec()   [ iLocal + iComponent * velocityFESpace.fe().nbFEDof() ] = unRepeated ( iGlobal );
                }
            }
            /*
            std::cout << elementConvectionVelocity.vec() << std::endl;
            std::cout << elementMeshVelocity.vec() << std::endl;
            std::cout << elementVelocity.vec() << std::endl;
            std::cout << M_elementDisplacement.vec() << std::endl;
            std::cout << elementVelocityRightHandSide.vec() << std::endl;
            std::cout << u_loc.vec() << std::endl;
            */
            for ( UInt iNode = 0 ; iNode < pressureFESpace.fe().nbFEDof() ; iNode++ )
            {
                // iLocal = iNode
                UInt iLocal = pressureFESpace.fe().patternFirst ( iNode );
                UInt iGlobal = pressureFESpace.dof().localToGlobalMap ( i, iLocal ) + numVelocityComponent * velocityFESpace.dim();
                // p^iNode local
                elementPressure[ iLocal ] = ukRepeated[ iGlobal ];
            }

            AssemblyElemental::shape_terms ( //M_elementDisplacement,
                density,
                viscosity,
                u_loc,
                elementVelocity,
                elementMeshVelocity,
                elementConvectionVelocity,
                elementPressure,
                *elementMatrixVelocity,
                M_FESpace.fe(),
                pressureFESpace.fe(),
                (ID) M_FESpace.fe().nbFEDof(),
                *elementMatrixPressure,
                0,
                wImplicit,
                alpha//,
                //elementMatrixConvective
            );

            //elementMatrixVelocity->showMe(std::cout);

            /*
            source_mass2( density,
                          elementVelocity,
                          *M_elementMatrixConvective,
                          velocityFESpace.fe(),
                          alpha );
            */

            AssemblyElemental::source_press ( 1.0,
                                              elementVelocity,
                                              *elementMatrixPressure,
                                              M_FESpace.fe(),
                                              pressureFESpace.fe(),
                                              (ID) M_FESpace.fe().nbFEDof() );

            //derivative of the convective term
            if ( convectiveTermDerivative )
                AssemblyElemental::mass_gradu ( density,
                                                elementVelocity,
                                                *elementMatrixConvective,
                                                velocityFESpace.fe() );
            /*
              std::cout << "source_press -> norm_inf( M_elementVectorVelocity )"  << std::endl;
            M_elementVectorPressure.showMe( std::cout );
            */
            //
            // Assembling
            //
            /*
            std::cout << "debut ====================" << std::endl;
            M_elementVectorPressure.showMe( std::cout );
            M_elementVectorVelocity.showMe( std::cout );
            std::cout << "fin   ====================" << std::endl;
            */
            UInt const velocityTotalDof ( velocityFESpace.dof().numTotalDof() );
            UInt const meshTotalDof ( M_FESpace.dof().numTotalDof() );
            for ( UInt iComponent = 0; iComponent < numVelocityComponent; ++iComponent )
            {
                for ( UInt jComponent = 0; jComponent < numVelocityComponent; ++jComponent )
                {
                    assembleMatrix ( *M_matrShapeDerVel,
                                     *elementMatrixVelocity,
                                     velocityFESpace.fe(),
                                     M_FESpace.fe(),
                                     velocityFESpace.dof(),
                                     M_FESpace.dof(),
                                     iComponent,
                                     jComponent,
                                     iComponent * velocityTotalDof,
                                     jComponent * meshTotalDof);

                    //assembling the derivative of the convective term
                    if ( convectiveTermDerivative )
                        assembleMatrix ( *M_matrShapeDerVel,
                                         *elementMatrixConvective,
                                         velocityFESpace.fe(),
                                         velocityFESpace.fe(),
                                         velocityFESpace.dof(),
                                         velocityFESpace.dof(),
                                         iComponent,
                                         jComponent,
                                         iComponent * velocityTotalDof,
                                         jComponent * velocityTotalDof );
                }

                assembleMatrix ( *M_matrShapeDerPressure,
                                 *elementMatrixPressure,
                                 pressureFESpace.fe(),
                                 M_FESpace.fe(),
                                 pressureFESpace.dof(),
                                 M_FESpace.dof(),
                                 (UInt) 0,
                                 iComponent,
                                 (UInt) 0,
                                 iComponent * meshTotalDof );
            }
        }

    if (  ! BCh.bcUpdateDone() )
    {
    	// BC boundary information update
    	BCh.bcUpdate ( *velocityFESpace.mesh(), velocityFESpace.feBd(), velocityFESpace.dof() );
    }

    bcManageMatrix ( *M_matrShapeDerVel, *velocityFESpace.mesh(), velocityFESpace.dof(), BCh, velocityFESpace.feBd(), 0.0, 0.0 );

    M_matrShapeDerVel->globalAssemble(M_FESpace.mapPtr(), velocityFESpace.mapPtr());
    M_matrShapeDerPressure->globalAssemble(M_FESpace.mapPtr(), pressureFESpace.mapPtr());


    chrono.stop();
    M_displayer.leaderPrintMax ("done in ", chrono.diff() );
}


} // end namespace LifeV
