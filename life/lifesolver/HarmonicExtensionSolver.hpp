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
    @file
    @brief Classes to hold algorithms for the mesh motion, for instance, involved in a ALE formulation.
    @author G. Fourestey
    @date 01-10-2007

    @contributor Simone Deparis <simone.deparis@epfl.ch>
    @maintainer Simone Deparis <simone.deparis@epfl.ch>

    This file contains classes which may be used to compute the extension inside the reference domain of a given
    displacement at a specified interface

*/

#ifndef HARMONICEXTENSIONSOLVER_H_
#define HARMONICEXTENSIONSOLVER_H_

#include <life/lifecore/life.hpp>
#include <life/lifefem/DOF.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifefem/AssemblyElemental.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/Assembly.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/FESpace.hpp>


namespace LifeV
{
/*!
  \class HarmonicExtension

  Base class which provides the harmonic extension of a given displacement on a specified part
  of the mesh boundary

  In order to deal with harmonic extensions, we have to provide a mesh (to  be moved), the parameters
  involved in the laplacian discretization: viscosity, quadrature rules and boundary conditions.
  This class contains a   PhysVectUnknown objet wich will hold the extension of the interface displacement.
  The constructor of the class built the global matrix of the discretized laplacian. The extension of the
  displacement is computed by calling the public method update. Finally, this extension can be recovered by
  calling method getDisplacement.

*/

template< typename Mesh,
typename SolverType = LifeV::SolverTrilinos>
class HarmonicExtensionSolver
{
public:

    //! @name Public Types
    //@{

    typedef SolverType                             solver_Type;

    typedef typename solver_Type::matrix_type      matrix_Type;
    typedef typename solver_Type::matrix_ptrtype   matrixPtr_Type;
    typedef typename solver_Type::vector_type      vector_Type;

    // OBSOLETE typedefs
    typedef SolverType                             solver_type;

    typedef typename solver_type::matrix_type      matrix_type;
    typedef typename solver_type::vector_type      vector_type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructors for an harmonics extensions
    /*!
      \param mmFESpace the FEspace that describes the problem
      \param comm  the Epetra_Comm to be used for communication
    */

    HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>&       mmFESpace,
                             boost::shared_ptr<Epetra_Comm>  comm);

    //! Constructors for an harmonics extensions with offset
    /*!
      \param mmFESpace the FEspace that describes the problem
      \param comm  the Epetra_Comm to be used for communication
      \param localMap use localMap instead of M_FESpace.map()
      \param offset use this offset to fill the matrix (both: row and column offset)
    */
    HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>&      mmFESpace,
                             boost::shared_ptr<Epetra_Comm> comm,
                             EpetraMap&                     localMap,
                             UInt                           offset =0
                           );

    //! virtual destructor
    virtual ~HarmonicExtensionSolver() {};

    //@}

    //! @name Methods
    //@{

    //! Set up data from GetPot
    /*!
        @param dataFile GetPot object
     */
    void setUp        ( const GetPot& dataFile );

    //! Initialize with velocityFunction and pressureFunction
    /*!
        @param disp vector of displacements to initialize the simulation
     */
    void initialize   ( const vector_Type& disp )  { setDisplacement(disp); }

    //! Update system
    /*!
      No argument is needd to update the system
     */
    void updateSystem();

    //! Update convective term, boundary condition and solve the linearized ns system
    /*!
        @param bcHandler BC handler
     */
    void iterate(BCHandler& BCh);

    //! returns wheter this processor is the leader.
    bool isLeader() const
    {
        return comm().MyPID() == 0;
    }

    //! prepare to recompute the preconditioner.
    void resetPrec(bool reset = true) { if (reset) M_linearSolver.precReset(); }

    //! manually rescale the system matrix by dt
    void rescaleMatrix(Real& dt) {*M_matrHE *= dt;}

    //! this method is not a setter. It is an operator=+
    void __attribute__ ((__deprecated__)) setMatrix(matrixPtr_Type matr) { addSystemMatrixTo(matr); }

    //! Adds the system matrix to the argument
    void addSystemMatrixTo(matrixPtr_Type matr) const {*matr += *M_matrHE;}

    //! Apply boundary conditions.
    /*!
        @param rightHandSide
        @param bcHandler
     */
    void applyBoundaryConditions(vector_Type& rhs, BCHandler& BCh);

    void computeMatrix();
    void updateDispDiff();

    //@}

    //! @name Set Methods
    //@{
    void setDisplacement(const vector_Type &disp) { M_disp = disp;}
    void setDispOld(const vector_Type &disp)  {    M_dispOld=disp;  }

    //! @name Get Methods
    //@{

    vector_Type const& dispOld() const  {return M_dispOld;}

    vector_Type const& dispDiff() const {return M_dispDiff;}
    vector_Type const& disp()     const {return M_disp;}

    EpetraMap const& getMap() const { return M_localMap; }

    const boost::shared_ptr<Epetra_Comm>& comm() const {return M_Displayer.comm();}
    //@}

private:

    //! Finite Element Space

    FESpace<Mesh, EpetraMap>&      M_FESpace;

    //! local map
    EpetraMap                      M_localMap;

    //! The matrix holding the values
    matrixPtr_Type                 M_matrHE;

    Displayer                      M_Displayer;
    int                            M_me;
    bool                           M_verbose;

    //! Elementary matrix : 3 blocks
    ElemMat                        M_elmat;

    //! The actual extension of the displacement
    vector_Type                    M_disp;
    vector_Type                    M_dispOld;
    vector_Type                    M_dispDiff;

    //! Auxiliary vector holding the second right hand of the system
    vector_Type                    M_secondRHS;

    //! The linear solver
    solver_Type                    M_linearSolver;

    //! Diffusion coefficient for the laplacian operator
    Real                           M_diffusion;

    UInt                           M_offset;
};

// ===================================================
// Constructors & Destructor
// ===================================================


template <typename Mesh, typename SolverType>
HarmonicExtensionSolver<Mesh, SolverType>::
HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>& mmFESpace,
                         boost::shared_ptr<Epetra_Comm>    comm ):
        M_FESpace               ( mmFESpace ),
        M_localMap              ( M_FESpace.map() ),
        M_matrHE                ( new matrix_Type (M_localMap ) ),
        M_Displayer              ( comm ),
        M_me                    ( comm->MyPID() ),
        M_verbose               ( M_me == 0 ),
        M_elmat                 ( M_FESpace.fe().nbFEDof(), nDimensions, nDimensions ),
        M_disp                  ( M_localMap ),
        M_dispOld               ( M_localMap ),
        M_dispDiff              ( M_localMap ),
        M_secondRHS             ( M_localMap ),
        M_linearSolver          ( comm ),
        M_diffusion             ( 1. ),
        M_offset                (0)
{
}

template <typename Mesh, typename SolverType>
HarmonicExtensionSolver<Mesh, SolverType>::
HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>& mmFESpace,
                         boost::shared_ptr<Epetra_Comm>              comm ,
                         EpetraMap& localMap,
                         UInt offset):
        M_FESpace               ( mmFESpace ),
        M_localMap              ( localMap),
        M_matrHE                ( new matrix_Type (M_localMap ) ),
        M_Displayer              ( comm ),
        M_me                    ( comm->MyPID() ),
        M_verbose               ( M_me == 0 ),
        M_elmat                 ( M_FESpace.fe().nbFEDof(), nDimensions, nDimensions ),
        M_disp                  ( mmFESpace.map() ),
        M_dispOld               ( M_disp.map() ),
        M_dispDiff              ( M_disp.map() ),
        M_secondRHS             ( M_disp.map() ),
        M_linearSolver          ( comm ),
        M_diffusion             ( 1. ),
        M_offset                (offset)
{
}

// ===================================================
// Methods
// ===================================================


template <typename Mesh, typename SolverType>
void HarmonicExtensionSolver<Mesh, SolverType>::setUp( const GetPot& dataFile )
{
    M_linearSolver.setDataFromGetPot( dataFile, "mesh_motion/solver" );
    M_linearSolver.setUpPrec(dataFile, "mesh_motion/prec");

    M_diffusion = dataFile("mesh_motion/diffusion",1.0);

    computeMatrix( );
    M_linearSolver.setMatrix( *M_matrHE );
} // end setUp


template <typename Mesh, typename SolverType>
void
HarmonicExtensionSolver<Mesh, SolverType>::updateSystem()
{
    M_Displayer.leaderPrint(" HE-  Updating the system ...                  ");

    M_dispOld = M_disp;

    M_Displayer.leaderPrint("done \n");
}

template <typename Mesh, typename SolverType>
void
HarmonicExtensionSolver<Mesh, SolverType>::iterate( BCHandler& BCh )
{
    Chrono chrono;

    // matrix and vector assembling communication
    M_Displayer.leaderPrint(" HE-  Applying boundary conditions ...         ");
    chrono.start();

    M_secondRHS *= 0.;
    applyBoundaryConditions(M_secondRHS, BCh);

    chrono.stop();
    M_Displayer.leaderPrintMax("done in " , chrono.diff() );

    // solving the system. Note: setMatrix(M_matrHE) done in setUp()
    M_linearSolver.solveSystem( M_secondRHS, M_disp, M_matrHE );
}

template <typename Mesh, typename SolverType>
void
HarmonicExtensionSolver<Mesh, SolverType>::applyBoundaryConditions(vector_Type& rhs, BCHandler& BCh)
{

    // CHANGED BY S. QUINODOZ !
    // "if" exchanged

    if (  ! BCh.bcUpdateDone() )
    {
        // BC boundary information update
        BCh.bcUpdate( *M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );
    }

    if (M_offset)//mans that this is the fullMonolithic case
    {
        BCh.setOffset(M_offset);
    }
    else
    {
        bcManageVector(rhs, *M_FESpace.mesh(), M_FESpace.dof(), BCh, M_FESpace.feBd(), 0., 1.0);
    }

    bcManageMatrix( *M_matrHE, *M_FESpace.mesh(), M_FESpace.dof(), BCh, M_FESpace.feBd(), 1.0, 0. );
}

template <typename Mesh, typename SolverType>
void HarmonicExtensionSolver<Mesh, SolverType>::computeMatrix( )
{
    Chrono chrono;
    chrono.start();
    M_Displayer.leaderPrint(" HE-  Computing constant matrices ...          ");

    M_matrHE.reset( new matrix_Type (M_localMap ) );

    UInt totalDof   = M_FESpace.dof().numTotalDof();
    // Loop on elements
    for ( UInt i = 1; i <= M_FESpace.mesh()->numVolumes(); ++i )
    {
        // Updating derivatives
        M_FESpace.fe().updateFirstDerivQuadPt( M_FESpace.mesh()->volumeList( i ) );
        M_elmat.zero();
        stiff( M_diffusion, M_elmat, M_FESpace.fe(), 0, 0, 3 );
        // Assembling
        for ( UInt j = 0; j < M_FESpace.fieldDim(); ++j )
        {
            assembleMatrix( *M_matrHE, M_elmat, M_FESpace.fe(), M_FESpace.dof(), j, j, j*totalDof+M_offset, j*totalDof+M_offset );
        }
    }

    M_matrHE->globalAssemble();

    chrono.stop();
    M_Displayer.leaderPrintMax("done in " , chrono.diff() );

}


template <typename Mesh, typename SolverType>
void
HarmonicExtensionSolver<Mesh, SolverType>::updateDispDiff()
{
    M_dispDiff =  M_disp ;
    M_dispDiff -= M_dispOld;
}


} // namespace LifeV
#endif //  HARMONICEXTENSIONSOLVER_H_
