
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

    bcManageMatrix ( *M_matrHE, *M_FESpace.mesh(), M_FESpace.dof(), BCh, M_FESpace.feBd(), 1.0, 0. );
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

    chrono.stop();
    M_displayer.leaderPrintMax ("done in " , chrono.diff() );

}

}
