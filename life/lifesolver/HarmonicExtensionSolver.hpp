/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
  \author G. Fourestey
  \date 10/2007

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
          typename SolverType = LifeV::SolverTrilinos>
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

    /*HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>& mmFESpace,
                             BCHandler&                mesh_BCh,
                             Epetra_Comm&              comm
                             );*/

    HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>& mmFESpace,
                             Epetra_Comm&              comm);

    HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>& mmFESpace,
                             Epetra_Comm&              comm,
                             EpetraMap&                locMap  ,
                             UInt                      offset =0
                             );
    //! operators overload



    //! This method updates the extension of the displacement, i.e. it solves the laplacian proglem

    void setUp        ( const GetPot& dataFile );

    void initialize   ( const vector_type& disp )  { setDisplacement(disp); }


    void buildSystem();
    void updateSystem();


//     template <typename Mesh>
    void iterate(BCHandler& BCh);

//     template <typename Mesh>
//     void updateExtensionTransp( Mesh& mesh, const Real& time = 0.0 );
    //! This method gives a reference to the computed harmonic extension.

    //const vector_type& displacement() const {return M_disp;}
    vector_type const& dispOld() const  {return M_dispOld;}

    vector_type const& dispDiff() const {return M_dispDiff;}
    vector_type const& disp()     const {return M_disp;}


    void setDisplacement(const vector_type &disp) { M_disp = disp;}

    //! This method interpolates the mesh velocity when necessary (refFE_u.nbNodes > _mesh.getRefFE().nbNodes)
//     template <typename Mesh>
//     void interpMeshVelocity( Mesh& mesh, const RefFE& refFE_u, const Dof& dof_u, Vector& wInterp );

    //! Returns a reference to the corresponding Dof object
//    const Dof& dofMesh() const;

    //! checking if BC are set
    //bool BCset() const {return M_setBC;}
    //! set the mesh BCs
    //void setBC(BCHandler& BCh);
    //! returns the BCHandler
    //const BCHandler& BCh_harmonicExtension() const {return *M_BCh_harmonicExtension;}
    //const BCHandler& bcHandler() const {return *M_BCh;}

    EpetraMap const& getMap() const { return M_localMap; }
    //Epetra_Map const& getRepeatedEpetraMap() const { return *M_localMap.getRepeatedEpetra_Map(); }

    const Epetra_Comm& comm() const {return *M_comm;}

    bool isLeader() const
    {
        return comm().MyPID() == 0;
    }

    void resetPrec() {M_resetPrec = true;}

    void rescaleMatrix(Real& dt){*M_matrHE *= dt;}
    void setMatrix(matrix_ptrtype matr){*matr += *M_matrHE;}
    void applyBoundaryConditions(vector_type& rhs, BCHandler& BCh);
    void computeMatrix();
    void updateDispDiff();
    //vector_type& deltaDisp(){return *M_dispDeltaDiff;}
    //void updateDeltaDisp(vector_type& deltaDisp){M_dispDeltaDiff.reset(new vector_type(deltaDisp));}

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

    //! boolean that indicates if le precond has to be recomputed
    bool                           M_reusePrec;
    int                            M_maxIterForReuse;
    bool                           M_resetPrec;

    //! BC holding the imposed boundary displacement
    //BCHandler*                     M_BCh;
    //bool                           M_setBC;

    //! Diffusion coefficient for the laplacian operator
    Real                           M_diffusion;

    UInt                           M_offset;
    //    boost::shared_ptr<vector_type>                 M_dispDeltaDiff; // useful for shape derivatives and Newton in general
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
/*template <typename Mesh, typename SolverType>
HarmonicExtensionSolver<Mesh, SolverType>::
HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>& mmFESpace,
                         BCHandler&                bcHandler,
                         Epetra_Comm&              comm ):
    M_FESpace               ( mmFESpace ),
    M_localMap              ( M_FESpace.map() ),
    M_matrHE                ( ),
    M_comm                  ( &comm ),
    M_me                    ( M_comm->MyPID() ),
    M_verbose               ( M_me == 0 ),
    M_disp                  ( M_localMap ),
    M_dispOld               ( M_localMap ),
    M_dispDiff              ( M_localMap ),
    M_elmat                 ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_f                     ( M_localMap ),
    M_linearSolver          ( comm ),
    M_prec                  ( ),
    M_reusePrec              ( true ),
    M_maxIterForReuse        ( -1 ),
    M_resetPrec              ( true ),
    //M_BCh                   ( &bcHandler ),
    M_diffusion             ( 1. ),
    M_offset                (0)

{
}*/


template <typename Mesh, typename SolverType>
HarmonicExtensionSolver<Mesh, SolverType>::
HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>& mmFESpace,
                         Epetra_Comm&    comm ):
    M_FESpace               ( mmFESpace ),
    M_localMap              ( M_FESpace.map() ),
    M_matrHE                ( new matrix_type (M_localMap ) ),
    M_comm                  ( &comm ),
    M_me                    ( M_comm->MyPID() ),
    M_verbose               ( M_me == 0 ),
    M_elmat                 ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_disp                  ( M_localMap ),
    M_dispOld               ( M_localMap ),
    M_dispDiff              ( M_localMap ),
    M_f                     ( M_localMap ),
    M_linearSolver          ( comm ),
    M_prec                  (  ),
    M_reusePrec             ( true ),
    M_maxIterForReuse       ( -1 ),
    M_resetPrec             ( true ),
    M_diffusion             ( 1. ),
    M_offset                (0)
{
}

template <typename Mesh, typename SolverType>
HarmonicExtensionSolver<Mesh, SolverType>::
HarmonicExtensionSolver( FESpace<Mesh, EpetraMap>& mmFESpace,
                         Epetra_Comm&              comm ,
                         EpetraMap& localMap,
                         UInt offset):
    M_FESpace               ( mmFESpace ),
    M_localMap              ( localMap),
    M_matrHE                ( new matrix_type (M_localMap ) ),
    M_comm                  ( &comm ),
    M_me                    ( M_comm->MyPID() ),
    M_verbose               ( M_me == 0 ),
    M_elmat                 ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_disp                  ( mmFESpace.map() ),
    M_dispOld               ( M_disp.getMap() ),
    M_dispDiff              ( M_disp.getMap() ),
    M_f                     ( M_disp.getMap() ),
    M_linearSolver          ( ),
    M_prec                   ( ),
    M_reusePrec              ( true ),
    M_maxIterForReuse        ( -1 ),
    M_resetPrec              ( true ),
    M_diffusion             ( 1. ),
    M_offset                (offset)
{
}



template <typename Mesh, typename SolverType>
void HarmonicExtensionSolver<Mesh, SolverType>::setUp( const GetPot& dataFile )
{

    M_linearSolver.setDataFromGetPot( dataFile, "mesh_motion/solver" );
    //    M_prec->setDataFromGetPot( dataFile, "mesh_motion/prec" );
    M_diffusion = dataFile("mesh_motion/diffusion",1.0);

    int maxIterSolver   = dataFile( "mesh_motion/solver/max_iter", -1);
    M_reusePrec       = dataFile( "mesh_motion/prec/reuse", true);
    M_maxIterForReuse = dataFile( "mesh_motion/solver/max_iter_reuse", maxIterSolver*8/10);

    std::string precType = dataFile( "mesh_motion/prec/prectype", "Ifpack");

    M_prec               = prec_type( PRECFactory::instance().createObject( precType ) );

    ASSERT(M_prec.get() != 0, "HE : Preconditioner not set");


    M_prec->setDataFromGetPot( dataFile, "mesh_motion/prec" );

    computeMatrix( );
    M_linearSolver.setMatrix( *M_matrHE );
}


template <typename Mesh, typename SolverType>
void HarmonicExtensionSolver<Mesh, SolverType>::computeMatrix( )
{
    Chrono chrono;
    chrono.start();

    if (M_verbose)
        std::cout << " he-  Computing constant matrices ...        " <<  std::flush;

    M_matrHE.reset( new matrix_type (M_localMap ) );

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

    M_matrHE->GlobalAssemble();

    chrono.stop();

    if (M_verbose)
        std::cout << " done in " << chrono.diff() << std::endl;


    // Initializations
    //M_disp = ZeroVector( _disp.size() );


}

//! timeadvance method
// This method updates the extension of the displacement, i.e. it solves the laplacian problem

template <typename Mesh, typename SolverType>
void
HarmonicExtensionSolver<Mesh, SolverType>::updateSystem()
{
    if (M_verbose)
        std::cout << "  HE- Updating the system      ... " << std::flush;

    M_dispOld = M_disp;

    if (M_verbose)
        std::cout << "  ok." << std::endl;

}


    template <typename Mesh, typename SolverType>
void
HarmonicExtensionSolver<Mesh, SolverType>::iterate( BCHandler& BCh )
{
    if (M_verbose)
        std::cout << "  HE- Updating boundary conditions : " << std::flush;

    // Initializations
    M_f    *= 0.;
    applyBoundaryConditions(M_f, BCh);

    Chrono chrono;


    if ( !M_reusePrec || M_resetPrec || !M_prec->set() )
    {
        chrono.start();

        if (M_verbose)
            std::cout << "  HE-  Computing the precond ...                "  <<  std::flush;

        M_prec->buildPreconditioner(M_matrHE);

        double condest = M_prec->Condest();

        M_linearSolver.setPreconditioner(M_prec);

        chrono.stop();
        if (M_verbose)
        {
            std::cout << "done in " << chrono.diff() << " s.\n";
            std::cout << "  HE-       Estimated condition number = " << condest << "\n" <<  std::flush;
        }


    }
    else
    {
        if (M_verbose)
            std::cout << "  HE-  Reusing  precond ...                \n" <<  std::flush;
    }


    // Real f_norm_inf(M_f.NormInf());

    if (M_verbose)
        std::cout << "  HE- Solving the system ... \n" << std::flush;

    int numIter =  M_linearSolver.solve( M_disp, M_f );

    M_resetPrec = (numIter > M_maxIterForReuse);


    chrono.stop();
    if (M_verbose)
    {
        std::cout << "HE- system solved in " << chrono.diff()
                  << " s. ( " << numIter << "  iterations. ) \n"
                  << std::flush;
    }

    //    M_dispDiff =  M_disp ;
    //    M_dispDiff -= M_dispOld;


}
    template <typename Mesh, typename SolverType>
void
HarmonicExtensionSolver<Mesh, SolverType>::updateDispDiff()
{
    M_dispDiff =  M_disp ;
    M_dispDiff -= M_dispOld;
}


template <typename Mesh, typename SolverType>
void
HarmonicExtensionSolver<Mesh, SolverType>::applyBoundaryConditions(vector_type& rhs, BCHandler& BCh)
{
    if(M_offset)//mans that this is the fullMonolithic case
        {
    	BCh.setOffset(M_offset);
        }
    else
        {
            if (M_verbose) std::cout << "\n  HE- Filling the rhs " << std::flush;
            //    bcManageVector(M_f, *M_FESpace.mesh(), M_FESpace.dof(), *BCh, M_FESpace.feBd(), 0., 1.0 );
            //bcManageVector(rhs, M_FESpace, *BCh, 0., 1.0 );
            bcManageVector(rhs, *M_FESpace.mesh(), M_FESpace.dof(), BCh, M_FESpace.feBd(), 0., 1.0);
            if (M_verbose) std::cout << "\n" << std::flush;
        }
    if (  !BCh.bdUpdateDone() )
    {
        // BC boundary information update
        if (M_verbose) std::cout << "\n      - Updating the BC " << std::flush;
        BCh.bdUpdate( *M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );
        if (M_verbose) std::cout << "\n      - Filling the matrix " ;
        bcManageMatrix( *M_matrHE, *M_FESpace.mesh(), M_FESpace.dof(), BCh, M_FESpace.feBd(), 1.0, 0. );
        if (M_verbose) std::cout << "\n" << std::flush;
    }
}

// This method updates the extension of the displacement, i.e. it solves the laplacian proglem
// template <typename Mesh, typename SolverType>
// void
// HarmonicExtensionSolver<Mesh, SolverType>::updateExtensionTransp( Mesh& mesh, const Real& 1.0 )
// {
//     // Boundary conditions treatment
//     bcManageVector( M_disp, M_FESpace.mesh(), M_FESpace.dof(), *M_BCh_harmonicExtension, M_FESpace.feBd(), time );
// }


//template <typename Mesh, typename SolverType>
//void
//HarmonicExtensionSolver<Mesh, SolverType>::setBC( BCHandler &BCh_harmonicExtension )
//{
//    M_BCh    = &BCh_harmonicExtension;
//    M_setBC  = true;
//}

// template <typename Mesh, typename SolverType>
// const Dof&
// HarmonicExtensionSolver<Mesh, SolverType>::
// dofMesh() const
// {
//     return _dof_mesh;
// }

}
#endif
