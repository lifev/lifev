/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*!
  \file meshMotion.h
  \brief Classes to hold algorithms for the mesh motion, for instance, involved in a ALE formulation.
  \version 1.0
  \author M.A. Fernandez
  \date 11/2002

  This file contains classes which may be used to compute the extension inside the reference domain of a given
  displacement at a specified interface

*/

#ifndef _HARMONICEXTENSIONSOLVER_H_
#define _HARMONICEXTENSIONSOLVER_H_

#include <life/lifecore/life.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifearray/pattern.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifefem/assemb.hpp>
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
          typename SolverType = LifeV::Epetra::SolverTrilinos>
class HarmonicExtensionSolver
{
public:


    typedef SolverType                             solver_type;

    typedef typename solver_type::matrix_type      matrix_type;
    typedef boost::shared_ptr<matrix_type>         matrix_ptrtype;
    typedef typename solver_type::vector_type      vector_type;

    typedef typename solver_type::prec_raw_type    prec_raw_type;
    typedef typename solver_type::prec_type        prec_type;




    //! Constructors for an harmonics extensions
    /*!
      \param mesh the mesh of the reference domain to be moved
      \param mu the "viscosity" in the laplacian operator
      \param Qr the quadrature rule for volumic elementary computations
      \param bdQr the quadrature rule for surface elementary computations
      \param mesh_BCh the list of boundary conditions involved in the harmonic extension

      \note The BCHandler objet (bch) holds the displacement imposed on moving boundary
      in the mesh trhough a BCVetor_Interface objet.

    */

    HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>& mmFESpace,
                             BCHandler&                mesh_BCh,
                             Epetra_Comm&              comm
                             );

    HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>& mmFESpace,
                             Epetra_Comm&              comm
                             );
    //! operators overload



    //! This method updates the extension of the displacement, i.e. it solves the laplacian proglem

    void setUp        ( const GetPot& dataFile );

    void buildSystem();
    void updateSystem();


//     template <typename Mesh>
    void iterate();

//     template <typename Mesh>
//     void updateExtensionTransp( Mesh& mesh, const Real& time = 0.0 );
    //! This method gives a reference to the computed harmonic extension.

    const vector_type& displacement() const {return M_disp;}
    const vector_type& displacementOld() const {return M_dispOld;}

    vector_type& dispDiff() {return M_dispDiff;}


    void setDisplacement(const Vector &disp) { M_disp = disp;}

    //! This method interpolates the mesh velocity when necessary (refFE_u.nbNodes > _mesh.getRefFE().nbNodes)
//     template <typename Mesh>
//     void interpMeshVelocity( Mesh& mesh, const RefFE& refFE_u, const Dof& dof_u, Vector& wInterp );

    //! Returns a reference to the corresponding Dof object
//    const Dof& dofMesh() const;

    //! checking if BC are set
    const bool BCset() const {return M_setBC;}
    //! set the mesh BCs
    void setBC(BCHandler &BCh_harmonicExtension);
    //! returns the BCHandler
    //const BCHandler& BCh_harmonicExtension() const {return *M_BCh_harmonicExtension;}
    const BCHandler& bcHandler() const {return *M_BCh;}

    EpetraMap const& getMap() const { return M_localMap; }
    Epetra_Map const& getRepeatedEpetraMap() const { return *M_localMap.getRepeatedEpetra_Map(); }


private:


    //! Finite Element Space

    FESpace<Mesh, EpetraMap>&      M_FESpace;

    //! local map
    EpetraMap                      M_localMap;

    //! The matrix holding the values
    matrix_ptrtype                 M_matrHE;

    Epetra_Comm*                   M_comm;
    int                            M_me;
    bool                           M_verbose;


//! Current element
//     CurrentFE                      M_fe;

//     //! Current element  on the boundary
//     CurrentBdFE                    M_feBd;

    //! Elementary matrix : 3 blocks
    ElemMat                        M_elmat;

    //! The actual extension of the displacement
    vector_type                    M_disp;
    vector_type                    M_dispOld;
    vector_type                    M_dispDiff;

    //! Auxiliary vector holding the second right hand of the system
    vector_type                    M_f;

    //! The linear solver

    solver_type                    M_linearSolver;

    prec_type                      M_prec;

    //! BC holding the imposed boundary displacement
    BCHandler*                     M_BCh;
    bool                           M_setBC;

    //! Diffusion coefficient for the laplacian operator
    Real                           M_diffusion;


};


//! Constructor for an harmonics extensions
/*!
  \param mesh the mesh of the reference domain to be moved
  \param diffusion the "viscosity" in the laplacian operator
  \param Qr the quadrature rule for volumic elementary computations
  \param bdQr the quadrature rule for surface elementary computations
  \param bch the list of boundary conditions involved in the harmonic extension

  \note The BCHandler objet (bch) holds the displacement imposed on moving boundary
  in the mesh trhough a BCVetor_Interface objet.

*/
template <typename Mesh, typename SolverType>
HarmonicExtensionSolver<Mesh, SolverType>::
HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>& mmFESpace,
                         BCHandler&                bcHandler,
                         Epetra_Comm&              comm ):
    M_FESpace               ( mmFESpace ),
    M_localMap              ( M_FESpace.map() ),
    M_matrHE                ( new matrix_type (M_FESpace.map() ) ),
    M_comm                  ( &comm ),
    M_me                    ( M_comm->MyPID() ),
    M_verbose               ( M_me == 0 ),
    M_disp                  ( M_FESpace.map() ),
    M_dispOld               ( M_FESpace.map() ),
    M_dispDiff              ( M_FESpace.map() ),
    M_f                     ( M_FESpace.map() ),
    M_elmat                 ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_prec                  ( new prec_raw_type() ),
    M_BCh                   ( &bcHandler ),
    M_diffusion             ( 1. )

{
}


template <typename Mesh, typename SolverType>
HarmonicExtensionSolver<Mesh, SolverType>::
HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>& mmFESpace,
                         Epetra_Comm&              comm ):
    M_FESpace               ( mmFESpace ),
    M_localMap              ( M_FESpace.map() ),
    M_matrHE                ( new matrix_type (M_FESpace.map() ) ),
    M_comm                  ( &comm ),
    M_me                    ( M_comm->MyPID() ),
    M_verbose               ( M_me == 0 ),
    M_disp                  ( M_FESpace.map() ),
    M_dispOld               ( M_FESpace.map() ),
    M_dispDiff              ( M_FESpace.map() ),
    M_f                     ( M_FESpace.map() ),
    M_elmat                 ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_prec                   ( new prec_raw_type() ),
    M_diffusion             ( 1. )
{
}



template <typename Mesh, typename SolverType>
void HarmonicExtensionSolver<Mesh, SolverType>::setUp( const GetPot& dataFile )
{

    Chrono chrono;

    if (M_verbose)
        std::cout << " he-  Computing constant matrices ...        ";

    chrono.start();

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
            assembleMatrix( *M_matrHE, M_elmat, M_FESpace.fe(), M_FESpace.dof(), j, j, j*totalDof, j*totalDof );
        }
    }

    M_matrHE->GlobalAssemble();

    chrono.stop();

    if (M_verbose)
        std::cout << " done in " << chrono.diff() << std::endl;


    // Initializations
    //M_disp = ZeroVector( _disp.size() );

    M_linearSolver.setDataFromGetPot( dataFile, "mesh_motion/solver" );
    M_linearSolver.setMatrix( *M_matrHE );

    M_diffusion = dataFile("meshMotion/miscallaneous/diffusion",1.0);

    M_prec->setDataFromGetPot( dataFile, "fluid/prec" );

    // Modify the matrix to handle boundary conditions


}


//! timeadvance method




// This method updates the extension of the displacement, i.e. it solves the laplacian problem

template <typename Mesh, typename SolverType>
void
HarmonicExtensionSolver<Mesh, SolverType>::updateSystem()
{
    M_dispOld = M_disp;
}


    template <typename Mesh, typename SolverType>
void
HarmonicExtensionSolver<Mesh, SolverType>::iterate()
{
    if (M_verbose)
        std::cout << "  HE- Updating boundary conditions ... " << std::flush;

    if (  ! M_BCh->bdUpdateDone() )
    {
        // BC boundary information update
        M_BCh->bdUpdate( *M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );
    }


    // Initializations
    M_f    *= 0.;

    //bcManageMatrix( *M_matrHE, *M_FESpace.mesh(), M_FESpace.dof(), *M_BCh, M_FESpace.feBd(), 1.0 );
    bcManage( *M_matrHE, M_f, *M_FESpace.mesh(), M_FESpace.dof(), *M_BCh, M_FESpace.feBd(), 1.0, 0. );


    M_prec->buildPreconditioner(M_matrHE);

    M_linearSolver.setPreconditioner(M_prec);

    //M_matrHE->spy("H");

    //exit(0);



    // Boundary conditions treatment
    //   bcManageVector( M_f, *M_FESpace.mesh(), M_FESpace.dof(), *M_BCh, M_FESpace.feBd(), 1.0, 1.0 );

    if (M_verbose)
        std::cout << "done.\n" << std::flush;

    if (M_verbose)
        std::cout << "  HE- Solving the system ... " << std::flush;

    M_linearSolver.solve( M_disp, M_f );

    if (M_verbose)
        std::cout << "done.\n" << std::flush;

    M_dispDiff =  M_disp ;
    M_dispDiff -= M_dispOld;


}



// This method updates the extension of the displacement, i.e. it solves the laplacian proglem
// template <typename Mesh, typename SolverType>
// void
// HarmonicExtensionSolver<Mesh, SolverType>::updateExtensionTransp( Mesh& mesh, const Real& 1.0 )
// {
//     // Boundary conditions treatment
//     bcManageVector( M_disp, M_FESpace.mesh(), M_FESpace.dof(), *M_BCh_harmonicExtension, M_FESpace.feBd(), time );
// }


template <typename Mesh, typename SolverType>
void
HarmonicExtensionSolver<Mesh, SolverType>::setBC( BCHandler &BCh_harmonicExtension )
{
    M_BCh    = &BCh_harmonicExtension;
    M_setBC  = true;
}

// template <typename Mesh, typename SolverType>
// const Dof&
// HarmonicExtensionSolver<Mesh, SolverType>::
// dofMesh() const
// {
//     return _dof_mesh;
// }

}
#endif
