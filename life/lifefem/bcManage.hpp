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
  \file bcManage.h
  \brief Functions for imposing boundary condtitions.
  \version 1.0
  \author M.A. Fernandez
  \date 07/2002

  This file contains the functions which may used to impose boundary
  conditions. This file shall be modified to implement new boundary
  coditions (not yet implemented: tangential, normal, etc...)  At the
  moment, we can impose Essential, Natural and Mixte boundary conditions,
  for Scalar or vectors problems.
*/
/*! \note Bug report (V. Martin, 12/02/2003)
  (tested with the mixed hybrid RT0 elements (Darcy))

  The new functions bcManageMatrix and bcManageVector
  are compiling for Essential FUNCTIONS
  but the results are WRONG!

  The new function bcManageMatrix passes for BC VECTORS
  (result: ??)
  The new function bcManageVector fails for BC VECTORS
  (segmentation fault).

*/

#ifndef _BCMANAGE_
#define _BCMANAGE_

#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/dof.hpp>
#include "life/lifefem/FESpace.hpp"


namespace LifeV
{
// ===================================================
// Boundary conditions treatment
// ===================================================
/* bcManage for boundary conditions depending on the solution U too,
   a class PointSolution would be useful */
template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void bcManage( Real (*mu)(Real t,Real x, Real y, Real z, Real u),
               MatrixType& A, VectorType& b, const MeshType& mesh, const Dof& dof,
               const BCHandler& BCh, CurrentBdFE& bdfem, const DataType coef,
               const DataType& t, VectorType& U )
{
    VectorType bRepeated(b.getMap(),Repeated);
    bool globalassemble=false;
    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
        {
            switch ( BCh[ i ].type() )
                {
                case Essential:  // Essential boundary conditions (Dirichlet)
                    globalassemble=true;
                    break;
                case Natural:  // Natural boundary conditions (Neumann)
                    if(BCh[ i ].isUDep())
                        bcNaturalManageUDep(mu, b, mesh, dof, BCh[ i ], bdfem, t,U, BCh.offset());
                    else
                        //in this case mu must be a constant, think about (not still implemented)
                        bcNaturalManage( b, mesh, dof, BCh[ i ], bdfem, t, BCh.offset() );
                    break;
                case Mixte:  // Mixte boundary conditions (Robin)

                    if(BCh[ i ].isUDep())
                        bcMixteManageUDep( A, bRepeated, mesh, dof, BCh[ i ], bdfem, t, U, BCh.offset());    //not implemented yet
                    else
                        bcMixteManage( A, bRepeated, mesh, dof, BCh[ i ], bdfem, t, BCh.offset() );
                    break;
                default:
                    ERROR_MSG( "This BC type is not yet implemented" );
                }
        }

    bRepeated.GlobalAssemble();

    b += bRepeated;
    if(globalassemble)
        A.GlobalAssemble();


    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
        {

            switch ( BCh[ i ].type() )
                {
                case Essential:  // Essential boundary conditions (Dirichlet)
                    if(BCh[ i ].isUDep())
                        bcEssentialManageUDep(A, b, mesh, dof, BCh[ i ], bdfem, coef, t,U, BCh.offset());
                    else
                        bcEssentialManage( A, b, mesh, dof, BCh[ i ], bdfem, coef, t, BCh.offset() );
                    break;
                case Natural:// Natural boundary conditions (Neumann)
                    break;
                case Mixte:  // Mixte boundary conditions (Robin)
                    break;
                default:
                    ERROR_MSG( "This BC type is not yet implemented" );
                }
        }
}
/* needed for ParabolicSolver due to the fact that
   essential bc are handled doing a trick on the matrix */
template <typename MatrixType, typename DataType>
void bcManageMtimeUDep( MatrixType& M, const Dof& dof,
                        const BCHandler& BCh, const DataType coef)
{
    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
        {

            if( BCh[ i ].type()==Essential )
                {
                    const BCBase& BCb=BCh[i];
                    ID idDof;

                    // Number of components involved in this boundary condition
                    UInt nComp = BCb.numberOfComponents();

                    // Number of total scalar Dof
                    UInt totalDof = dof.numTotalDof();

                    if ( BCb.dataVector() )
                        { //! If BC is given under a vectorial form

                            //not possible
                            ERROR_MSG( "This type of BCVector does not exists on bc depentent on solution" );
                        }
                    else
                        { //! If BC is given under a functional form

                            std::vector<ID>   idDofVec(0);
                            idDofVec.reserve(BCb.list_size()*nComp);

                            // Loop on BC identifiers
                            for ( ID i = 1; i <= BCb.list_size(); ++i )
                                {
                                    // Loop on components involved in this boundary condition
                                    for ( ID j = 1; j <= nComp; ++j )
                                        {
                                            // Global Dof
                                            idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof + BCh.offset();
                                            idDofVec.push_back(idDof-1);

                                        }
                                }
                            // Modifying ONLY matrix
                            M.diagonalize( idDofVec, coef );

                        }
                }
        }
}



template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void bcManage( MatrixType& A, VectorType& b, const MeshType& mesh, const Dof& dof,
               const BCHandler& BCh,
               CurrentBdFE& bdfem, const DataType& coef, const DataType& t )
{

    VectorType bRepeated(b.getMap(),Repeated);
    bool globalassemble=false;
    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
        {
//             std::cout << i << " " << BCh[i].type() << " " << Flux << std::endl;
            switch ( BCh[ i ].type() )
                {
                case Essential:  // Essential boundary conditions (Dirichlet)
                    globalassemble=true;
                    break;
                case Natural:  // Natural boundary conditions (Neumann)
                    bcNaturalManage( b, mesh, dof, BCh[ i ], bdfem, t, BCh.offset());
                    break;
                case Mixte:  // Mixte boundary conditions (Robin)
                    bcMixteManage( A, bRepeated, mesh, dof, BCh[ i ], bdfem, t, BCh.offset() );
                    break;
                case Flux:
                    bcFluxManage( A, b, mesh, dof, BCh[ i ], bdfem, t, BCh.offset()+i);
                    break;
                default:
                    ERROR_MSG( "This BC type is not yet implemented" );
                }
        }
    bRepeated.GlobalAssemble();
    b += bRepeated;
    if(globalassemble)
        A.GlobalAssemble();

    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
        {
            switch ( BCh[ i ].type() )
                {
                case Essential:  // Essential boundary conditions (Dirichlet)
                    bcEssentialManage( A, b, mesh, dof, BCh[ i ], bdfem, coef, t, BCh.offset() );
                    break;
                case Natural:  // Natural boundary conditions (Neumann)
                case Mixte:  // Mixte boundary conditions (Robin)
                case Flux:
                    break;
                default:
                    ERROR_MSG( "This BC type is not yet implemented" );
                }
        }
}


//Version that treates only the matrix modifications
//Miguel:10/02  - Mixte : V. Martin: 03/03 // I added the time t for the mixte case. V. Martin
template <typename MatrixType, typename MeshType, typename DataType>
void bcManageMatrix( MatrixType&      A,
                     const MeshType&  mesh,
                     const Dof&       dof,
                     const BCHandler& BCh,
                     CurrentBdFE&     bdfem,
                     const DataType&  coef,
                     const DataType&  t = 0 )
{

    bool globalassemble=false;
    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
        {

            switch ( BCh[ i ].type() )
                {
                case Essential:  // Essential boundary conditions (Dirichlet)
                    globalassemble=true;
                case Natural:  // Natural boundary conditions (Neumann)
                    break;
                case Mixte:  // Mixte boundary conditions (Robin)
                    bcMixteManageMatrix( A, mesh, dof, BCh[ i ], bdfem, t, BCh.offset() );
                    break;
                default:
                    ERROR_MSG( "This BC type is not yet implemented" );
                }
        }
    if(globalassemble)
        A.GlobalAssemble();

    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
        {

            switch ( BCh[ i ].type() )
                {
                case Essential:  // Essential boundary conditions (Dirichlet)
                    bcEssentialManageMatrix( A, dof, BCh[ i ], coef, BCh.offset() );  //! Bug here???
                    break;
                case Natural:  // Natural boundary conditions (Neumann)
                    // Do nothing
                    break;
                case Mixte:  // Mixte boundary conditions (Robin)
                    break;
                default:
                    ERROR_MSG( "This BC type is not yet implemented" );
                }
        }
}


//Version that treates only the vector modifications
//Miguel:10/02  - Mixte : V. Martin: 03/03
template <typename VectorType, typename MeshType, typename DataType>
void bcManageVector( VectorType&      b,
                     const MeshType&  mesh,
                     const Dof&       dof,
                     const BCHandler& BCh,
                     CurrentBdFE&     bdfem,
                     const DataType&  t,
                     const DataType&  coef )
{
    VectorType bRepeated(b.getMap(),Repeated);

    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
        {

            switch ( BCh[ i ].type() )
                {
                case Essential:  // Essential boundary conditions (Dirichlet)
                    bcEssentialManageVector( b, dof, BCh[ i ], t, coef, BCh.offset() );
                    break;
                case Natural:  // Natural boundary conditions (Neumann)
                    bcNaturalManage( b, mesh, dof, BCh[ i ], bdfem, t, BCh.offset() );
                    break;
                case Mixte:  // Mixte boundary conditions (Robin)
                    bcMixteManageVector( bRepeated, mesh, dof, BCh[ i ], bdfem, t, BCh.offset() );
                    break;
                default:
                    ERROR_MSG( "This BC type is not yet implemented" );
                }
        }

    bRepeated.GlobalAssemble();

    b += bRepeated;
}

// FESpace version of bcManageVector

template <typename VectorType, typename DataType, typename Mesh, typename EpetraMap>
void bcManageVector( VectorType&                     b,
                     FESpace<Mesh, EpetraMap>&       fespace,
                     const BCHandler&                BCh,
                     const DataType&                 t,
                     const DataType&                 coef )
{
    VectorType bRepeated(b.getMap(),Repeated);

    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
        {

            switch ( BCh[ i ].type() )
                {
                case Essential:  // Essential boundary conditions (Dirichlet)
                    bcEssentialManageVector( b, fespace.dof(), BCh[ i ], t, coef, BCh.offset() );
                    break;
                case Natural:  // Natural boundary conditions (Neumann)
                    bcNaturalManage( b, *fespace.mesh(), fespace.dof(), BCh[ i ], fespace.feBd(), t, BCh.offset() );
                    break;
                case Mixte:  // Mixte boundary conditions (Robin)
                    bcMixteManageVector( bRepeated, *fespace.mesh(), fespace.dof(), BCh[ i ], fespace.feBd(), t, BCh.offset() );
                    break;
                default:
                    ERROR_MSG( "This BC type is not yet implemented" );
                }
        }

    bRepeated.GlobalAssemble();

    b += bRepeated;
}



//! Version for mixed problem with block matrices
//! Alain, 07/08/02
template <typename MatrixType1, typename MatrixType2, typename VectorType,
          typename MeshType, typename DataType>
void bcManage( MatrixType1& C, MatrixType2& trD, VectorType& b,
               const MeshType& mesh, const Dof& dof, const BCHandler& BCh,
               CurrentBdFE& bdfem, const DataType& coef, const DataType& t )
{
    VectorType bRepeated(b.getMap(),Repeated);

    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
        {

            switch ( BCh[ i ].type() )
                {
                case Essential:  // Essential boundary conditions (Dirichlet)
                    bcEssentialManage( C, trD, b, mesh, dof, BCh[ i ], bdfem, coef, t, BCh.offset() );
                    break;
                case Natural:  // Natural boundary conditions (Neumann)
                    bcNaturalManage( b, mesh, dof, BCh[ i ], bdfem, t, BCh.offset() );
                    break;
                case Mixte:  // Mixte boundary conditions (Robin)
                    bcMixteManage( C, trD, bRepeated, mesh, dof, BCh[ i ], bdfem, t, BCh.offset() );
                    break;
                default:
                    ERROR_MSG( "This BC type is not yet implemented" );
                }
        }
    bRepeated.GlobalAssemble();

    b += bRepeated;
}

//! Alain (nov. 2002)
//! Using row access of trD, it is possible to access to the column of
//! the transpose matrix D. So the global matrix remains symmetric.
template <typename MatrixType1, typename MatrixType2, typename MatrixType3,
          typename VectorType, typename MeshType, typename DataType>
void bcManage( MatrixType1& C, MatrixType2& trD, MatrixType3& D,
               VectorType& b, VectorType& bp, const MeshType& mesh,
               const Dof& dof, const BCHandler& BCh, CurrentBdFE& bdfem,
               const DataType& coef, const DataType& t )
{
    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
        {

            switch ( BCh[ i ].type() )
                {
                case Essential:  // Essential boundary conditions (Dirichlet)
                    bcEssentialManage( C, trD, D, b, bp, mesh, dof, BCh[ i ], bdfem, coef, t, BCh.offset() );
                    break;
                case Natural:  // Natural boundary conditions (Neumann)
                    bcNaturalManage( b, mesh, dof, BCh[ i ], bdfem, t, BCh.offset() );
                    break;
                case Mixte:  // Mixte boundary conditions (Robin)
                    bcMixteManage( C, trD, D, b, bp, mesh, dof, BCh[ i ], bdfem, t, BCh.offset() );
                    break;
                default:
                    ERROR_MSG( "This BC type is not yet implemented" );
                }
        }
}

// ===================================================
// Essential BC
// ===================================================
template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void bcEssentialManageUDep( MatrixType& A, VectorType& b, const MeshType& /*mesh*/, const Dof& dof,
                            const BCBase& BCb, const CurrentBdFE& /*bdfem*/, const DataType& coef,
                            const DataType& t, const VectorType& U , UInt offset=0)
{

    ID idDof;

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    if ( BCb.dataVector() )
        { //! If BC is given under a vectorial form

            //not possible
            ERROR_MSG( "This type of BCVector does not exists on bc depentent on solution" );
        }
    else
        { //! If BC is given under a functional form

            std::vector<ID>   idDofVec(0);
            std::vector<Real> datumVec(0);

            idDofVec.reserve(BCb.list_size()*nComp);
            datumVec.reserve(BCb.list_size()*nComp);

            DataType x, y, z;
            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {
                    // Coordinates of the node where we impose the value
                    x = static_cast< const IdentifierEssential* >( BCb( i ) ) ->x();
                    y = static_cast< const IdentifierEssential* >( BCb( i ) ) ->y();
                    z = static_cast< const IdentifierEssential* >( BCb( i ) ) ->z();

                    // Loop on components involved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                        {
                            // Global Dof
                            idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof + offset;

                            Real datum = BCb( t, x, y, z, BCb.component( j ) ,U[idDof-1]);

                            datumVec.push_back(datum);
                            idDofVec.push_back(idDof-1);

                        }
                }

            // Modifying matrix and right hand side
            A.diagonalize( idDofVec, coef, b, datumVec);
        }
}
template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void bcEssentialManage( MatrixType& A, VectorType& b, const MeshType& /*mesh*/, const Dof& dof, const BCBase& BCb,
                        const CurrentBdFE& /*bdfem*/, const DataType& coef, const DataType& t, UInt offset )
{

    ID idDof;

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();


    std::vector<ID>   idDofVec(0);
    std::vector<Real> datumVec(0);

    idDofVec.reserve(BCb.list_size()*nComp);
    datumVec.reserve(BCb.list_size()*nComp);

    if ( BCb.dataVector() )
        { //! If BC is given under a vectorial form

            const BCVectorInterface* pId = static_cast< const BCVectorInterface* > (BCb.pointerToBCVector());
            assert( pId != 0);

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {

                    if ( !pId->dofInterface().isMyInterfaceDof(BCb( i ) ->id()))
                        {
                            continue;
                        }

                    // Loop on components involved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                        {
                            // Global Dof
                            idDof = BCb( i )->id() + ( BCb.component( j ) - 1 ) * totalDof + offset;

                            datumVec.push_back(BCb( BCb( i ) ->id(), BCb.component( j ) ));
                            idDofVec.push_back(idDof - 1);
                        }
                }
        }
    else
        { //! If BC is given under a functional form

            DataType x, y, z;
            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {
                    // Coordinates of the node where we impose the value
                    x = static_cast< const IdentifierEssential* >( BCb( i ) ) ->x();
                    y = static_cast< const IdentifierEssential* >( BCb( i ) ) ->y();
                    z = static_cast< const IdentifierEssential* >( BCb( i ) ) ->z();

                    // Loop on components involved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                        {
                            // Global Dof
                            idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof + offset;

                            datumVec.push_back(BCb( t, x, y, z, BCb.component( j ) ));
                            idDofVec.push_back(idDof - 1);

                        }
                }
        }

	// Modifying matrix and right hand side
	A.diagonalize( idDofVec, coef, b, datumVec);

}


//Version that treates only the matrix modifications
//Miguel:10/02
template <typename MatrixType, typename DataType>
void bcEssentialManageMatrix( MatrixType& A, const Dof& dof, const BCBase& BCb, const DataType& coef, UInt offset )
{

    ID idDof;
    UInt totalDof;

    // Number of total scalar Dof
    totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    std::vector<ID>   idDofVec(0);
    idDofVec.reserve(BCb.list_size()*nComp);

    // Loop on BC identifiers
    for ( ID i = 1; i <= BCb.list_size(); ++i )
        {
            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
                {
                    // Global Dof
                    idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                    idDofVec.push_back(idDof-1);
                    // Modifying ONLY matrix
                    //A.diagonalize( idDof - 1, coef );
                }
        }
    // Modifying ONLY matrix
    A.diagonalize( idDofVec, coef, offset);

}

//Version that treates only the vector modifications
//Miguel:10/02
template <typename VectorType, typename DataType>
void bcEssentialManageVector( VectorType&     b,
                              const Dof&      dof,
                              const BCBase&   BCb,
                              const DataType& t,
                              const DataType& coef,
                              UInt            offset )
{


    ID idDof;
    UInt totalDof;

    // Number of total scalar Dof
    totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

//     std::cout << " totaldof = " << totalDof << std::endl;
//     std::cout << " nComp    = " << nComp << std::endl;

    std::vector<int>   idDofVec(0);
    idDofVec.reserve(BCb.list_size()*nComp);
    std::vector<Real> datumVec(0);
    datumVec.reserve(BCb.list_size()*nComp);

    if ( BCb.dataVector() )
        {  //! If BC is given under a vectorial form

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {
                    // Loop on components involved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                        {

                            // Global Dof
                            idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof + offset;
                            //                            std::cout << "idof = " << idDof << std::endl;

//                             std::cout <<  "comp = " << BCb.component( j ) << " idof = " << BCb(i)->id() << " -> ";
//                             std::cout <<  BCb( BCb( i ) ->id(), BCb.component( j ) );
//                             std::cout << " ." << std::endl;

                            idDofVec.push_back( idDof );
                            datumVec.push_back( coef * BCb( BCb( i ) ->id(), BCb.component( j ) ) );

                            //                             std::cout << "--" << std::endl;
                        }
                }
        }
    else
        {  //! If BC is given under a functionnal form

            DataType x, y, z;
            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {
                    // Coordinates of the node where we impose the value
                    x = static_cast< const IdentifierEssential* >( BCb( i ) ) ->x();
                    y = static_cast< const IdentifierEssential* >( BCb( i ) ) ->y();
                    z = static_cast< const IdentifierEssential* >( BCb( i ) ) ->z();

                    // Loop on components involved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                        {
                            // Global Dof

                            idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof + offset;
                            // Modifying right hand side
                            idDofVec.push_back(idDof);
                            datumVec.push_back( coef * BCb( t, x, y, z, BCb.component( j ) ) );
                        }
                }
        }

    b.replaceGlobalValues( idDofVec, datumVec);
    //     b.GlobalAssemble(Insert);

}


//! version for mixed problem with block matrices
//! Alain, 07/08/02
template <typename MatrixType1, typename MatrixType2, typename VectorType,
          typename MeshType, typename DataType>
void bcEssentialManage( MatrixType1& A,
                        MatrixType2& trD,
                        VectorType& b,
                        const MeshType& /*mesh*/,
                        const Dof& dof,
                        const BCBase& BCb,
                        const CurrentBdFE& /*bdfem*/,
                        const DataType& coef,
                        const DataType& t ,
                        UInt offset)
{
    ID idDof;
    DataType x, y, z;
    UInt totalDof;

    // Number of total scalar Dof
    totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    std::vector<ID>   idDofVec(0);
    std::vector<Real> datumVec(0);

    idDofVec.reserve(BCb.list_size()*nComp);
    datumVec.reserve(BCb.list_size()*nComp);


    if ( BCb.dataVector() )
        {  //! If BC is given under a vectorial form

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {
                    // Loop on components involved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                        {
                            // Global Dof
                            idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                            // Modifying matrix and right hand side

                            datumVec.push_back(BCb( BCb( i ) ->id(), BCb.component( j ) ));
                            idDofVec.push_back(idDof-1);

                            //A.diagonalize( idDof - 1, coef, b, BCb( BCb( i ) ->id(), BCb.component( j ) ) );
                            trD.zero_row( idDof - 1 );
                        }
                }
        }
    else
        {  //! If BC is given under a functionnal form

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {

                    // Coordinates of the node where we impose the value
                    x = static_cast< const IdentifierEssential* >( BCb( i ) ) ->x();
                    y = static_cast< const IdentifierEssential* >( BCb( i ) ) ->y();
                    z = static_cast< const IdentifierEssential* >( BCb( i ) ) ->z();

                    // Loop on components involved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                        {

                            // Global Dof
                            idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                            // Modifying matrix and right hand side
                            datumVec.push_back(BCb( t, x, y, z, BCb.component( j ) ));
                            idDofVec.push_back(idDof-1);
                            //A.diagonalize( idDof - 1, coef, b, BCb( t, x, y, z, BCb.component( j ) ) );
                            trD.zero_row( idDof - 1 );
                        }
                }
        }

    // Modifying matrix and right hand side
    A.diagonalize( idDofVec, coef, b, datumVec, offset );
}

//! version with column diagonalization of D.
//! Alain. nov 2002.
template <typename MatrixType1, typename MatrixType2, typename MatrixType3,
          typename VectorType, typename MeshType, typename DataType>
void bcEssentialManage( MatrixType1& A, MatrixType2& trD, MatrixType3& D,
                        VectorType& b, VectorType& bp, const MeshType& mesh,
                        const Dof& dof, const BCBase& BCb,
                        const CurrentBdFE& bdfem,
                        const DataType& coef, const DataType& t,
                        UInt offset)
{
    ID idDof;
    DataType x, y, z;
    UInt totalDof;

    // Number of total scalar Dof
    totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    std::vector<ID>   idDofVec(0);
    std::vector<Real> datumVec(0);

    idDofVec.reserve(BCb.list_size()*nComp);
    datumVec.reserve(BCb.list_size()*nComp);

    if ( BCb.dataVector() )
        {  //! If BC is given under a vectorial form

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {
                    // Loop on components involved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                        {
                            // Global Dof
                            idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                            // Modifying matrix and right hand side
                            datumVec.push_back(BCb( BCb( i ) ->id(), BCb.component( j ) ));
                            idDofVec.push_back(idDof-1);

                            //A.diagonalize( idDof - 1, coef, b, BCb( BCb( i ) ->id(), BCb.component( j ) ) );
                            zero_row_col( idDof - 1, trD, D, bp, BCb( t, x, y, z, BCb.component( j ) ) );
                        }
                }
        }
    else
        {  //! If BC is given under a functionnal form

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {

                    // Coordinates of the node where we impose the value
                    x = static_cast< const IdentifierEssential* >( BCb( i ) ) ->x();
                    y = static_cast< const IdentifierEssential* >( BCb( i ) ) ->y();
                    z = static_cast< const IdentifierEssential* >( BCb( i ) ) ->z();

                    // Loop on components involved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                        {

                            // Global Dof
                            idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                            // Modifying matrix and right hand side
                            datumVec.push_back(BCb( t, x, y, z, BCb.component( j ) ));
                            idDofVec.push_back(idDof-1);
                            //A.diagonalize( idDof - 1, coef, b, BCb( t, x, y, z, BCb.component( j ) ) );
                            zero_row_col( idDof - 1, trD, D, bp, BCb( t, x, y, z, BCb.component( j ) ) );
                        }
                }
        }

    // Modifying matrix and right hand side
    A.diagonalize( idDofVec, coef, b, datumVec, offset );

}

// ===================================================
// Natural BC
// ===================================================
template <typename VectorType, typename MeshType, typename DataType>
void bcNaturalManageUDep( Real (*mu)(Real t,Real x, Real y, Real z, Real u),
                          VectorType& b, const MeshType& mesh, const Dof& dof,
                          const BCBase& BCb, CurrentBdFE& bdfem,
                          const DataType& t, const VectorType& U, UInt offset )
{

    // Number of local Dof (i.e. nodes) in this face
    UInt nDofF = bdfem.nbNode;

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    const IdentifierNatural* pId;
    ID ibF, idDof ;

    if ( BCb.dataVector() )
        { //! If BC is given under a vectorial form
            ERROR_MSG( "This type of BCVector does not exists on bc depentent on solution\n" );
        }
    else
        {  //! If BC is given under a functionnal form

            DataType x, y, z;

            if(nComp!=1)
                {
                    ERROR_MSG("For now bcNaturalManageUDep cannot handle non scalar solutions\n");
                }

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {

                    // Pointer to the i-th itdentifier in the list
                    pId = static_cast< const IdentifierNatural* >( BCb( i ) );

                    // Number of the current boundary face
                    ibF = pId->id();

                    // Updating face stuff
                    bdfem.updateMeas( mesh.boundaryFace( ibF ) );

                    std::vector<Real> locU(nDofF+1);    //assumes U is a vec of reals, TODO: deal with more comp
                    Real uPt;            //value in the point
                    for(ID idofLocU=0;idofLocU<nDofF;idofLocU++)
                        {
                            ID idGDofU=pId->bdLocalToGlobal(idofLocU+1)+( BCb.component( 1 ) - 1 ) * totalDof + offset;
                            locU[idofLocU]=U[idGDofU-1];
                        }


                    // Loop on total Dof per Face
                    for ( ID idofF = 1; idofF <= nDofF; ++idofF )
                        {  //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                            //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                            // Loop on components involved in this boundary condition
                            for ( ID j = 1; j <= nComp; ++j )
                                {

                                    //global Dof
                                    idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                    // Loop on quadrature points
                                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                        {
                                            bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                                            uPt=0.0;
                                            for(ID idofLocU=0;idofLocU<nDofF;idofLocU++)
                                                {
                                                    //    Debug(800)<<"debug* naturalManageUDep entering ulocU\n";
                                                    uPt+=locU[idofLocU]*bdfem.phi( int( idofLocU  ),l );
                                                    //    Debug(800)<<"debug* naturalManageUDep exiting ulocU\n";
                                                }

                                            // Adding right hand side contribution
                                            b[ idDof ] += bdfem.phi( int( idofF - 1 ), l ) * BCb( t, x, y, z, BCb.component( j ),uPt ) *
                                                mu(t,x,y,z,uPt)*bdfem.weightMeas( l ); // BASEINDEX + 1
                                            //    Debug(800)<<"debug* naturalManageUDep done one ulocU\n";
                                        }
                                }
                        }
                }
        }
    //    Debug(800)<<"debug* end of naturalManageUDep\n";
}


template <typename VectorType, typename MeshType, typename DataType>
void bcNaturalManage( VectorType& b,
                      const MeshType& mesh,
                      const Dof& dof, const
                      BCBase& BCb,
                      CurrentBdFE& bdfem,
                      const DataType& t,
                      UInt offset )
{

    // Number of local Dof (i.e. nodes) in this face
    UInt nDofF = bdfem.nbNode;

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    const IdentifierNatural* pId;
    ID ibF, idDof, icDof, gDof;
    Real sum;

    if ( BCb.dataVector() )
        {
            //! If BC is given under a vectorial form
            switch ( BCb.pointerToBCVector()->type() )
                {
                case 0:  // if the BC is a vector which values don't need to be integrated
                    {
                        VectorType bUnique(b.getMap(),Unique);

                        std::vector<int>  idDofVec(0);
                        idDofVec.reserve(BCb.list_size()*nComp);
                        std::vector<Real> datumVec(0);
                        datumVec.reserve(BCb.list_size()*nComp);
                        // double datum;

                        // Loop on BC identifiers
                        for ( ID i = 1; i <= BCb.list_size(); ++i )
                            {
                                // Loop on components involved in this boundary condition
                                for ( ID j = 1; j <= nComp; ++j )
                                    {
                                        ID id = BCb(i)->id();

                                        // Global Dof
                                        idDof = id + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                        // we have to avoid adding twice an already integrated quantity:
                                        //  if ( b.getMap().getMap(Unique)->LID(idDof) < 0 ) continue;

                                        idDofVec.push_back( idDof);

                                        // Modifying right hand side (assuming BCvector is a flux)
                                        if(BCb.isgammaVec())  datumVec.push_back( BCb.GammaVec(id, BCb.component(j) )* BCb( id , BCb.component( j ) )) ;
                                        else  datumVec.push_back( BCb.gammaCoef()* BCb( id , BCb.component( j ) ));

                                        // b[ idDof ] += BCb.gammaCoef()* BCb( id , BCb.component( j ) );
                                    }
                            }

                        bUnique.replaceGlobalValues(idDofVec, datumVec);
                        bUnique.GlobalAssemble(Insert);

                        ASSERT( b.getMaptype() == Unique , "b must have unique map, otherwise data will be multiply added on cpu interfaces." );
                        b += bUnique;

                    }
                    break;

                case 1:  // if the BC is a vector of values to be integrated
                    {
                        VectorType bRepeated(b.getMap(),Repeated);

                        // Loop on BC identifiers
                        for ( ID i = 1; i <= BCb.list_size(); ++i )
                            {
                                // Pointer to the i-th itdentifier in the list
                                pId = static_cast< const IdentifierNatural* >( BCb( i ) );

                                // Number of the current boundary face
                                ibF = pId->id();

                                // Updating face stuff
                                bdfem.updateMeasNormalQuadPt( mesh.bElement( ibF ) );

                                // Loop on total Dof per Face
                                for ( ID l = 1; l <= nDofF; ++l )
                                    {

                                        gDof = pId->bdLocalToGlobal( l );

                                        // Loop on components involved in this boundary condition
                                        for ( UInt ic = 0; ic < nComp; ++ic )
                                            {
                                                icDof = gDof + ic * totalDof + offset;

                                                // Loop on quadrature points
                                                for ( int iq = 0; iq < bdfem.nbQuadPt; ++iq )
                                                    {
                                                        sum=0.0;
                                                        // data on quadrature point
                                                        for ( ID m = 1; m <= nDofF; ++m )
                                                            sum +=  BCb( pId->bdLocalToGlobal( m ) , 1 ) * bdfem.phi( int( m - 1 ), iq );
                                                                // Adding right hand side contribution
                                                                bRepeated[ icDof ] += sum * bdfem.phi( int( l - 1 ), iq ) * bdfem.normal( int( ic ), iq )
                                                                * bdfem.weightMeas( iq );

                                                    }
                                            }
                                    }
                            }

                        bRepeated.GlobalAssemble();
                        ASSERT( b.getMaptype() == Unique , "here b should passed as repeated, otherwise not sure of what happens at the cpu interfaces ." );
                        b += bRepeated;
                    }
                    break;
                case 2:  // if the BC is a vector of values with components to be integrated
                    // Loop on BC identifiers
                    {
                        VectorType bRepeated(b.getMap(),Repeated);

                        // Loop on BC identifiers
                        for ( ID i = 1; i <= BCb.list_size(); ++i )
                            {
                                // Pointer to the i-th itdentifier in the list
                                pId = static_cast< const IdentifierNatural* >( BCb( i ) );

                                // Number of the current boundary face
                                ibF = pId->id();

                                // Updating face stuff
                                bdfem.updateMeasNormalQuadPt( mesh.bElement( ibF ) );

                                // Loop on total Dof per Face
                                for ( ID idofF = 1; idofF <= nDofF; ++idofF )
                                    {

                                        gDof = pId->bdLocalToGlobal( idofF );

                                        // Loop on space dimensions
                                        for ( ID ic = 1; ic <= nComp; ++ic ) //modifica Matteo 28/07/08 to make "Component"
                                            //  for ( UInt ic = 0; ic < nDimensions; ++ic )
                                            {
                                                icDof = gDof + ( BCb.component( ic ) - 1 ) * totalDof+ offset;   //Components passed separatedly

                                                // Loop on quadrature points
                                                for ( int iq = 0; iq < bdfem.nbQuadPt; ++iq )
                                                    {
                                                        sum = 0;
                                                        // data on quadrature point
                                                        for ( ID m = 1; m <= nDofF; ++m )
                                                            sum +=  BCb( pId->bdLocalToGlobal( m ) , BCb.component( ic ) ) * bdfem.phi( int( m - 1 ), iq );  //Components passed separatedly

                                                        // Adding right hand side contribution
                                                        bRepeated[ icDof ] += sum *  bdfem.phi( int( idofF - 1 ), iq ) *
                                                            bdfem.weightMeas( iq );
                                                    }
                                            }
                                    }
                            }
                        bRepeated.GlobalAssemble();
                        ASSERT( b.getMaptype() == Unique , "here b should passed as unique, otherwise not sure of what happens at the cpu interfaces ." );
                        b += bRepeated;
                    }
                    break;
                    default:
                    ERROR_MSG( "This type of BCVector does not exist" );
                }
        }

    /*
      MODIFIED by V. Martin : what follows...

      Goal : to remain consistant with the functionnal part,
      and assume that the BCvector is a velocity (= a normal
      derivative and not a flux an integrated normal derivative).
      In the functionnal case, it is a velocity that is considered.

      One has to correctly define what his BCvector represents!!
      -> a flux?
      -> a velocity?

      It can be discussed, because without the measure, the BC vector
      still has a meaning : it's a flux (instead of a velocity), and
      can be used in a simple way in some cases.

      if ( BCb.dataVector() ) {//! If BC is given under a vectorial form

      // Loop on BC identifiers
      for (ID i=1; i<=BCb.list_size(); ++i) {

      // Pointer to the i-th itdentifier in the list
      pId = static_cast< const IdentifierNatural* >( BCb(i) );

      // Number of the current boundary face
      ibF = pId->id();

      // Updating face stuff
      bdfem.updateMeas( mesh.bElement(ibF) );

      // Loop on total Dof per Face
      for (ID idofF=1; idofF<= nDofF; ++idofF) {
      // Loop on components involved in this boundary condition
      for (ID j=1; j<=nComp; ++j) {

      //global Dof
      idDof  =  BCb(i)->id() + (BCb.component(j)-1)*totalDof;

      // Loop on quadrature points
      for(int l=0; l<bdfem.nbQuadPt; ++l) {

      // Adding right hand side contribution (assuming BCvector is a velocity (normal derivative))
      b[idDof-1] += bdfem.phi(int(idofF-1),l) * BCb( BCb(i)->id(),BCb.component(j) ) *
      bdfem.weightMeas(l);
      }
      }
      }
      }
      }
    */
    else
        {  //! If BC is given under a functionnal form

            //std::cout << "BC Natural manage w/ function" << std::endl;
            DataType x, y, z;
            VectorType bRepeated(b.getMap(),Repeated);

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {
                    // Pointer to the i-th itdentifier in the list
                    pId = static_cast< const IdentifierNatural* >( BCb( i ) );
                    // Number of the current boundary face
                    ibF = pId->id();
                    // Updating face stuff
                    bdfem.updateMeasNormalQuadPt( mesh.bElement( ibF ) );
                    // Loop on total Dof per Face
                    for ( ID idofF = 1; idofF <= nDofF; ++idofF )
                        {  //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                            //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)
                            // Loop on components involved in this boundary condition
                            for ( ID j = 1; j <= nComp; ++j )
                                {
                                    //global Dof
                                    idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof + offset;
                                    // Loop on quadrature points
                                    for ( int iq = 0; iq < bdfem.nbQuadPt; ++iq )
                                        {
                                            bdfem.coorQuadPt( x, y, z, iq ); // quadrature point coordinates
                                            // Adding right hand side contribution
                                                                                                    // Adding right hand side contribution
                                            switch (BCb.mode())
                                                {
                                                case Full:
                                                    bRepeated[ idDof ] += bdfem.phi( int( idofF - 1 ), iq ) * BCb( t, x, y, z, BCb.component( j ) ) *
                                                        bdfem.weightMeas( iq ); // BASEINDEX + 1
                                                    break;
                                                case Normal:
                                                    //                                                    std::cout << bdfem.normal(int(j), iq) << " ";
                                                    bRepeated[ idDof ] += BCb( t, x, y, z, BCb.component( j ) )*
                                                        bdfem.phi( int( idofF - 1 ), iq )*
                                                        bdfem.weightMeas( iq )*bdfem.normal( int(j - 1), iq );
                                                    break;
//                                                 case Tangential:
//                                                     bRepeated[ idDof ] += sum*
//                                                         bdfem.phi( int( idofF - 1 ), iq )*
//                                                         bdfem.weightMeas( iq )*bdfem.tangent( 0, int( j ), iq );
//                                                     break;
                                                }
//                                             bRepeated[ idDof ] += bdfem.phi( int( idofF - 1 ), iq ) * BCb( t, x, y, z, BCb.component( j ) ) *
//                                                 bdfem.weightMeas( iq ); // BASEINDEX + 1
                                        }
                                }
                            //                            std::cout << "ok" << std::endl;
                        }
                }
            bRepeated.GlobalAssemble();
            ASSERT( b.getMaptype() == Unique , "here b should passed as unique, otherwise not sure of what happens at the cpu interfaces ." );
            //  b=bRepeated;
            b+= bRepeated;
        }
} // bcNaturalManage

// ===================================================
// Mixte BC
// ===================================================
template <typename MatrixType, typename VectorType, typename DataType, typename MeshType>
void bcMixteManageUDep( MatrixType& A, VectorType& b, const MeshType& mesh, const Dof& dof, const BCBase& BCb,
                        CurrentBdFE& bdfem, const DataType& t ,const VectorType& U, UInt offset)
{
    ERROR_MSG("error bcMixteManageUDep not still implemented\n");
}

template <typename MatrixType, typename VectorType, typename DataType, typename MeshType>
void bcMixteManage( MatrixType& A, VectorType& b, const MeshType& mesh, const Dof& dof, const BCBase& BCb,
                    CurrentBdFE& bdfem, const DataType& t, UInt offset )
{

    // Number of local Dof in this face
    UInt nDofF = bdfem.nbNode;

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    DataType sum;

    const IdentifierNatural* pId;
    ID ibF, idDof, jdDof, kdDof;

    if ( BCb.dataVector() )
        {   //! If BC is given under a vectorial form

            //! for the moment, only one coefficient per BCvector.
            DataType mcoef, mbcb;

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {

                    // Pointer to the i-th identifier in the list
                    pId = static_cast< const IdentifierNatural* >( BCb( i ) );

                    // Number of the current boundary face
                    ibF = pId->id();

                    // Updating face stuff
                    bdfem.updateMeas( mesh.bElement( ibF ) );

                    // Loop on total Dof per Face
                    for ( ID idofF = 1; idofF <= nDofF; ++idofF )
                        { //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                            //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                            // Loop on components involved in this boundary condition
                            for ( ID j = 1; j <= nComp; ++j )
                                {

                                    sum = 0;

                                    // Global Dof
                                    //vincent please check again for your Mixte-FE it doesn't work for Q1:
                                    //   idDof  =  BCb(i)->id() + (BCb.component(j)-1)*totalDof;
                                    idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                    // Loop on quadrature points
                                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                        {

                                            // Compute the mixte coefficients on the quadrature point - vector of mixte
                                            // coefficients is given on the element nodes !
                                            mcoef = 0.0;
                                            mbcb = 0.0;
                                            for( ID n = 1; n <= nDofF; ++n)
                                                {
                                                    kdDof=pId->bdLocalToGlobal( n ); // + ( BCb.component( j ) - 1 ) * totalDof;
                                                    if(BCb.ismixteVec())
                                                        mcoef += BCb.MixteVec( kdDof, BCb.component( j ) ) * bdfem.phi( int( n - 1 ), l );

                                                    else  mcoef += BCb.mixteCoef() * bdfem.phi( int( n - 1 ), l );

                                                    if(BCb.isbetaVec())  mbcb += BCb.BetaVec( kdDof, BCb.component( j ) )
                                                                             * BCb( kdDof, BCb.component( j )) * bdfem.phi( int( n - 1 ), l );

                                                    else  mbcb += BCb.betaCoef() * BCb( kdDof, BCb.component( j )) * bdfem.phi( int( n - 1 ), l );

                                                    // kdDof=pId->bdLocalToGlobal( n ); // + ( BCb.component( j ) - 1 ) * totalDof;
                                                    // mcoef += BCb.MixteVec( kdDof, BCb.component( j ) ) * bdfem.phi( int( n - 1 ), l );
                                                    // mbcb += BCb( kdDof, BCb.component( j ) ) * bdfem.phi( int( n - 1 ), l );
                                                }

                                            // Contribution to the diagonal entry of the elementary boundary mass matrix
                                            //                        sum += mcoef * bdfem.phi( int( idofF - 1 ), l ) * bdfem.phi( int( idofF - 1 ), l ) *
                                            //                               bdfem.weightMeas( l );
                                            sum += mcoef* bdfem.phi( int( idofF - 1 ), l ) *
                                                bdfem.phi( int( idofF - 1 ), l ) *bdfem.weightMeas( l );

                                            // Adding right hand side contribution
                                            //vincent please check again for your Mixte-FE it doesn't work for Q1:
                                            //     b[idDof-1] += bdfem.phi(int(idofF-1),l) * BCb(BCb(i)->id(),BCb.component(j)) *
                                            b[ idDof ] += bdfem.phi( int( idofF - 1 ), l ) * mbcb * // BASEINDEX + 1
                                                bdfem.weightMeas( l );

                                        }

                                    // Assembling diagonal entry
                                    A.set_mat_inc( idDof - 1, idDof - 1, sum );
                                }

                            // Upper diagonal columns of the elementary boundary mass matrix
                            for ( ID k = idofF + 1 ; k <= nDofF ; ++k )
                                {

                                    // Loop on components invoved in this boundary condition
                                    for ( ID j = 1; j <= nComp; ++j )
                                        {

                                            sum = 0;

                                            // Globals Dof: row and columns
                                            //vincent please check again for your Mixte-FE it doesn't work for Q1:
                                            //     idDof  =  BCb(i)->id() + (BCb.component(j)-1)*totalDof;
                                            //     jdDof  =  BCb(k)->id() + (BCb.component(j)-1)*totalDof;
                                            idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof + offset;
                                            jdDof = pId->bdLocalToGlobal( k ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                            // Loop on quadrature points
                                            for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                                {

                                                    // Compute the mixte coefficient on the quadrature point - vector of mixte
                                                    // coefficients is given on the element nodes !
                                                    mcoef = 0.0;
                                                    for( ID n = 1; n <= nDofF; ++n)
                                                        {
                                                            kdDof=pId->bdLocalToGlobal( n ); // + ( BCb.component( j ) - 1 ) * totalDof;
                                                            if(BCb.ismixteVec())
                                                                mcoef += BCb.MixteVec( kdDof, BCb.component( j ) ) * bdfem.phi( int( n - 1 ), l );

                                                            else mcoef += BCb.mixteCoef() * bdfem.phi( int( n - 1 ), l );
                                                            //  mcoef += BCb.MixteVec( kdDof, BCb.component( j ) ) * bdfem.phi( int( n - 1 ), l );
                                                        }


                                                    // Upper diagonal entry of the elementary boundary mass matrix
                                                    //                            sum += mcoef * bdfem.phi( int( idofF - 1 ), l ) * bdfem.phi( int( k - 1 ), l ) *
                                                    //                                   bdfem.weightMeas( l );
                                                    sum += mcoef*bdfem.phi( int( idofF - 1 ), l ) *
                                                        bdfem.phi( int( k - 1 ), l ) * bdfem.weightMeas( l );

                                                }

                                            // Assembling upper entry.  The boundary mass matrix is symetric
                                            A.set_mat_inc( idDof - 1, jdDof - 1, sum );
                                            A.set_mat_inc( jdDof - 1, idDof - 1, sum );
                                        }
                                }
                        }
                }
        }

    else
        {  //! If BC is given under a functionnal form

            DataType x, y, z;

            const BCFunctionMixte* pBcF = static_cast<const BCFunctionMixte*>( BCb.pointerToFunctor() );

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {

                    // Pointer to the i-th identifier in the list
                    pId = static_cast< const IdentifierNatural* >( BCb( i ) );

                    // Number of the current boundary face
                    ibF = pId->id();

                    // Updating face stuff
                    bdfem.updateMeas( mesh.bElement( ibF ) );

                    // Loop on total Dof per Face
                    for ( ID idofF = 1; idofF <= nDofF; ++idofF )
                        { //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                            //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                            // Loop on components invoved in this boundary condition
                            for ( ID j = 1; j <= nComp; ++j )
                                {

                                    sum = 0;

                                    // Global Dof (outside the quad point loop. V. Martin)
                                    idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                    // Loop on quadrature points
                                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                        {

                                            bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                                            // Contribution to the diagonal entry of the elementary boundary mass matrix
                                            sum += pBcF->coef( t, x, y, z, BCb.component( j ) ) * bdfem.phi( int( idofF - 1 ), l ) * bdfem.phi( int( idofF - 1 ), l ) *
                                                bdfem.weightMeas( l );

                                            // Global Dof (Why inside this loop?? V. Martin)
                                            // idDof = pId->bdLocalToGlobal(idofF) + (BCb.component(j)-1)*totalDof;

                                            // Adding right hand side contribution
                                            b[ idDof ] += bdfem.phi( int( idofF - 1 ), l ) * BCb( t, x, y, z, BCb.component( j ) ) * // BASEINDEX + 1
                                                bdfem.weightMeas( l );
                                        }

                                    // Global Dof (Why repeat?? V. Martin)
                                    // idDof = pId->bdLocalToGlobal(idofF) + (BCb.component(j)-1)*totalDof;

                                    // Assembling diagonal entry
                                    A.set_mat_inc( idDof - 1, idDof - 1, sum );
                                }

                            // Upper diagonal columns of the elementary boundary mass matrix
                            for ( ID k = idofF + 1 ; k <= nDofF ; ++k )
                                {

                                    // Loop on components invoved in this boundary condition
                                    for ( ID j = 1; j <= nComp; ++j )
                                        {

                                            sum = 0;

                                            // Loop on quadrature points
                                            for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                                {

                                                    bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                                                    // Upper diagonal entry of the elementary boundary mass matrix
                                                    sum += pBcF->coef( t,  x, y, z, BCb.component( j )   ) * bdfem.phi( int( idofF - 1 ), l ) * bdfem.phi( int( k - 1 ), l ) *
                                                        bdfem.weightMeas( l );
                                                }

                                            // Globals Dof: row and columns
                                            idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof + offset;
                                            jdDof = pId->bdLocalToGlobal( k ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                            // Assembling upper entry.  The boundary mass matrix is symetric
                                            A.set_mat_inc( idDof - 1, jdDof - 1, sum );
                                            A.set_mat_inc( jdDof - 1, idDof - 1, sum );
                                        }
                                }
                        }
                }
        }
}

// ===================================================
// Mixte BC
// ===================================================
// Version that treates only the matrix modifications
// V. Martin 02/03/2003
template <typename MatrixType, typename DataType, typename MeshType>
void bcMixteManageMatrix( MatrixType& A, const MeshType& mesh, const Dof& dof,
                          const BCBase& BCb, CurrentBdFE& bdfem, const DataType& t, UInt offset )
{

    // Number of local Dof in this face
    UInt nDofF = bdfem.nbNode;

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    DataType sum;

    const IdentifierNatural* pId;
    ID ibF, idDof, jdDof;

    if ( BCb.dataVector() )
        {   //! If BC is given under a vectorial form

            //! for the moment, only one coefficient per BCvector.
            DataType mcoef = BCb.mixteCoef();   //!< the mixte coefficient

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {

                    // Pointer to the i-th identifier in the list
                    pId = static_cast< const IdentifierNatural* >( BCb( i ) );

                    // Number of the current boundary face
                    ibF = pId->id();

                    // Updating face stuff
                    bdfem.updateMeas( mesh.bElement( ibF ) );

                    // Loop on total Dof per Face
                    for ( ID idofF = 1; idofF <= nDofF; ++idofF )
                        { //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                            //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                            // Loop on components invoved in this boundary condition
                            for ( ID j = 1; j <= nComp; ++j )
                                {

                                    sum = 0;

                                    // Global Dof
                                    idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                    // Loop on quadrature points
                                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                        {

                                            // Contribution to the diagonal entry of the elementary boundary mass matrix
                                            sum += mcoef * bdfem.phi( int( idofF - 1 ), l ) * bdfem.phi( int( idofF - 1 ), l ) *
                                                bdfem.weightMeas( l );
                                        }

                                    // Assembling diagonal entry
                                    A.set_mat_inc( idDof - 1, idDof - 1, sum );
                                }

                            // Upper diagonal columns of the elementary boundary mass matrix
                            for ( ID k = idofF + 1 ; k <= nDofF ; ++k )
                                {

                                    // Loop on components invoved in this boundary condition
                                    for ( ID j = 1; j <= nComp; ++j )
                                        {

                                            sum = 0;

                                            // Loop on quadrature points
                                            for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                                {

                                                    // Upper diagonal entry of the elementary boundary mass matrix
                                                    sum += mcoef * bdfem.phi( int( idofF - 1 ), l ) * bdfem.phi( int( k - 1 ), l ) *
                                                        bdfem.weightMeas( l );
                                                }

                                            // Globals Dof: row and columns
                                            idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof + offset;
                                            jdDof = BCb( k ) ->id() + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                            // Assembling upper entry.  The boundary mass matrix is symetric
                                            A.set_mat_inc( idDof - 1, jdDof - 1, sum );
                                            A.set_mat_inc( jdDof - 1, idDof - 1, sum );
                                        }
                                }
                        }
                }
        }

    else
        {  //! If BC is given under a functionnal form

            DataType x, y, z;

            const BCFunctionMixte* pBcF = static_cast<const BCFunctionMixte*>( BCb.pointerToFunctor() );

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {

                    // Pointer to the i-th identifier in the list
                    pId = static_cast< const IdentifierNatural* >( BCb( i ) );

                    // Number of the current boundary face
                    ibF = pId->id();

                    // Updating face stuff
                    bdfem.updateMeas( mesh.bElement( ibF ) );

                    // Loop on total Dof per Face
                    for ( ID idofF = 1; idofF <= nDofF; ++idofF )
                        { //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                            //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                            // Loop on components invoved in this boundary condition
                            for ( ID j = 1; j <= nComp; ++j )
                                {

                                    sum = 0;

                                    // Global Dof (outside the quad point loop. V. Martin)
                                    idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                    // Loop on quadrature points
                                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                        {

                                            bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                                            // Contribution to the diagonal entry of the elementary boundary mass matrix
                                            sum += pBcF->coef( t, x, y, z, BCb.component(j) ) * bdfem.phi( int( idofF - 1 ), l ) * bdfem.phi( int( idofF - 1 ), l ) *
                                                bdfem.weightMeas( l );

                                            // Global Dof (Why inside this loop?? V. Martin)
                                            // idDof = pId->bdLocalToGlobal(idofF) + (BCb.component(j)-1)*totalDof;
                                        }

                                    // Global Dof (Why repeat?? V. Martin)
                                    // idDof = pId->bdLocalToGlobal(idofF) + (BCb.component(j)-1)*totalDof;

                                    // Assembling diagonal entry
                                    A.set_mat_inc( idDof - 1, idDof - 1, sum );
                                }

                            // Upper diagonal columns of the elementary boundary mass matrix
                            for ( ID k = idofF + 1 ; k <= nDofF ; ++k )
                                {

                                    // Loop on components invoved in this boundary condition
                                    for ( ID j = 1; j <= nComp; ++j )
                                        {

                                            sum = 0;

                                            // Loop on quadrature points
                                            for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                                {

                                                    bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                                                    // Upper diagonal entry of the elementary boundary mass matrix
                                                    sum += pBcF->coef( t, x, y, z, BCb.component( j )  ) * bdfem.phi( int( idofF - 1 ), l ) * bdfem.phi( int( k - 1 ), l ) *
                                                        bdfem.weightMeas( l );
                                                }

                                            // Globals Dof: row and columns
                                            idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof + offset;
                                            jdDof = pId->bdLocalToGlobal( k ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                            // Assembling upper entry.  The boundary mas matrix is symetric
                                            A.set_mat_inc( idDof - 1, jdDof - 1, sum );
                                            A.set_mat_inc( jdDof - 1, idDof - 1, sum );
                                        }
                                }
                        }
                }
        }
}

// ===================================================
// Mixte BC
// ===================================================
// Version that treates only the vector modifications
// V. Martin 02/03/2003
template <typename VectorType, typename DataType, typename MeshType>
void bcMixteManageVector( VectorType& b,
                          const MeshType& mesh,
                          const Dof& dof,
                          const BCBase& BCb,
                          CurrentBdFE& bdfem,
                          const DataType& t,
                          UInt offset )
{

    // Number of local Dof in this face
    UInt nDofF = bdfem.nbNode;

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    const IdentifierNatural* pId;
    ID ibF, idDof;

    if ( BCb.dataVector() )
        {   //! If BC is given under a vectorial form

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {

                    // Pointer to the i-th identifier in the list
                    pId = static_cast< const IdentifierNatural* >( BCb( i ) );

                    // Number of the current boundary face
                    ibF = pId->id();

                    // Updating face stuff
                    bdfem.updateMeas( mesh.bElement( ibF ) );

                    // Loop on total Dof per Face
                    for ( ID idofF = 1; idofF <= nDofF; ++idofF )
                        { //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                            //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                            // Loop on components invoved in this boundary condition
                            for ( ID j = 1; j <= nComp; ++j )
                                {
                                    // Global Dof
                                    idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                    // Loop on quadrature points
                                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                        {
                                            // Adding right hand side contribution
                                            b[ idDof ] += bdfem.phi( int( idofF - 1 ), l ) * BCb( BCb( i ) ->id(), BCb.component( j ) ) * // BASEINDEX + 1
                                                bdfem.weightMeas( l );
                                        }
                                }
                        }
                }
        }

    else
        {  //! If BC is given under a functionnal form

            DataType x, y, z;

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {

                    // Pointer to the i-th identifier in the list
                    pId = static_cast< const IdentifierNatural* >( BCb( i ) );

                    // Number of the current boundary face
                    ibF = pId->id();

                    // Updating face stuff
                    bdfem.updateMeas( mesh.bElement( ibF ) );

                    // Loop on total Dof per Face
                    for ( ID idofF = 1; idofF <= nDofF; ++idofF )
                        { //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                            //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                            // Loop on components invoved in this boundary condition
                            for ( ID j = 1; j <= nComp; ++j )
                                {

                                    // Global Dof (outside the quad point loop. V. Martin)
                                    idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                    // Loop on quadrature points
                                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                        {

                                            bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                                            // Global Dof (Why inside this loop?? V. Martin)
                                            // idDof = pId->bdLocalToGlobal(idofF) + (BCb.component(j)-1)*totalDof;

                                            // Adding right hand side contribution
                                            b[ idDof ] += bdfem.phi( int( idofF - 1 ), l ) * BCb( t, x, y, z, BCb.component( j ) ) *// BASEINDEX + 1
                                                bdfem.weightMeas( l );
                                        }
                                }
                        }
                }
        }
}

//! version for mixed problem with block matrices
//! Alain, 07/08/02
template <typename MatrixType1, typename MatrixType2, typename VectorType,
          typename DataType, typename MeshType>
void bcMixteManage( MatrixType1& A, MatrixType2 & trD, VectorType& b,
                    const MeshType& mesh, const Dof& dof, const BCBase& BCb,
                    CurrentBdFE& bdfem, const DataType& t, UInt offset )
{

    ASSERT( !BCb.dataVector() , "BC Vector not yet implemented for this particular bcMixteManage." );

    // Number of local Dof in this face
    UInt nDofF = bdfem.nbNode;

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    DataType x, y, z, sum;
    const IdentifierNatural* pId;
    ID ibF, idDof, jdDof;

    const BCFunctionMixte* pBcF = static_cast<const BCFunctionMixte*>( BCb.pointerToFunctor() );

    for ( ID i = 1; i <= BCb.list_size(); ++i )
        {

            // Pointer to the i-th identifier in the list
            pId = static_cast< const IdentifierNatural* >( BCb( i ) );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            bdfem.updateMeas( mesh.bElement( ibF ) );

            // Loop on total Dof per Face
            for ( ID i = 1; i <= nDofF; ++i )
                {

                    // Loop on components involved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                        {

                            sum = 0;

                            // Loop on quadrature points
                            for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                {

                                    bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                                    // Contribution to the diagonal entry of the elementary boundary mass matrix
                                    sum += pBcF->coef( t, x, y, z, BCb.component( j )  ) * bdfem.phi( int( i - 1 ), l ) * bdfem.phi( int( i - 1 ), l ) *
                                        bdfem.weightMeas( l );

                                    // Global Dof
                                    idDof = pId->bdLocalToGlobal( i ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                    // Adding right hand side contribution
                                    b[ idDof ] += bdfem.phi( int( i - 1 ), l ) * BCb( t, x, y, z, BCb.component( j ) ) * bdfem.weightMeas( l ); // BASEINDEX + 1
                                }

                            // Global Dof
                            idDof = pId->bdLocalToGlobal( i ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                            // Assembling diagonal entry
                            A.set_mat_inc( idDof - 1, idDof - 1, sum );
                            trD.zero_row( idDof - 1 );
                        }

                    // Upper diagonal columns of the elementary boundary mass matrix
                    for ( ID k = i + 1 ; k <= nDofF ; ++k )
                        {

                            // Loop on components involved in this boundary condition
                            for ( ID j = 1; j <= nComp; ++j )
                                {

                                    sum = 0;

                                    // Loop on quadrature points
                                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                        {

                                            bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                                            // Upper diagonal entry of the elementary boundary mass matrix
                                            sum += pBcF->coef( t, x, y, z, BCb.component( j )  ) * bdfem.phi( int( i - 1 ), l ) * bdfem.phi( int( k - 1 ), l ) *
                                                bdfem.weightMeas( l );
                                        }

                                    // Globals Dof: row and columns
                                    idDof = pId->bdLocalToGlobal( i ) + ( BCb.component( j ) - 1 ) * totalDof + offset;
                                    jdDof = pId->bdLocalToGlobal( k ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                    // Assembling upper entry.  The boundary mas matrix is symetric
                                    A.set_mat_inc( idDof - 1, jdDof - 1, sum );
                                    A.set_mat_inc( jdDof - 1, idDof - 1, sum );
                                }
                        }
                }
        }
}
//! version with diagonalization of D
//! Alain, 07/08/02
template <typename MatrixType1, typename MatrixType2, typename MatrixType3,
          typename VectorType, typename DataType, typename MeshType>
void bcMixteManage( MatrixType1& A, MatrixType2 & trD, MatrixType3 & D,
                    VectorType& b, VectorType& bp, const MeshType& mesh,
                    const Dof& dof, const BCBase& BCb, CurrentBdFE& bdfem,
                    const DataType& t, UInt offset )
{
    ASSERT( !BCb.dataVector() , "BC Vector not yet implemented for this particular bcMixteManage." );

    // Number of local Dof in this face
    UInt nDofF = bdfem.nbNode;

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    DataType x, y, z, sum;
    const IdentifierNatural* pId;
    ID ibF, idDof, jdDof;

    const BCFunctionMixte* pBcF = static_cast<const BCFunctionMixte*>( BCb.pointerToFunctor() );

    for ( ID i = 1; i <= BCb.list_size(); ++i )
        {

            // Pointer to the i-th identifier in the list
            pId = static_cast< const IdentifierNatural* >( BCb( i ) );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            bdfem.updateMeas( mesh.bElement( ibF ) );

            // Loop on total Dof per Face
            for ( ID i = 1; i <= nDofF; ++i )
                {

                    // Loop on components involved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                        {

                            sum = 0;

                            // Loop on quadrature points
                            for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                {

                                    bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                                    // Contribution to the diagonal entry of the elementary boundary mass matrix
                                    sum += pBcF->coef( t, x, y, z, BCb.component( j )  ) * bdfem.phi( int( i - 1 ), l ) * bdfem.phi( int( i - 1 ), l ) *
                                        bdfem.weightMeas( l );

                                    // Global Dof
                                    idDof = pId->bdLocalToGlobal( i ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                    // Adding right hand side contribution
                                    b[ idDof ] += bdfem.phi( int( i - 1 ), l ) * BCb( t, x, y, z, BCb.component( j ) ) * bdfem.weightMeas( l ); // BASEINDEX + 1
                                }

                            // Global Dof
                            idDof = pId->bdLocalToGlobal( i ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                            // Assembling diagonal entry
                            A.set_mat_inc( idDof - 1, idDof - 1, sum );
                            zero_row_col( idDof - 1, trD, D, bp, sum );
                        }

                    // Upper diagonal columns of the elementary boundary mass matrix
                    for ( ID k = i + 1 ; k <= nDofF ; ++k )
                        {

                            // Loop on components involved in this boundary condition
                            for ( ID j = 1; j <= nComp; ++j )
                                {

                                    sum = 0;

                                    // Loop on quadrature points
                                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                                        {

                                            bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                                            // Upper diagonal entry of the elementary boundary mass matrix
                                            sum += pBcF->coef( t, x, y, z, BCb.component( j ) ) * bdfem.phi( int( i - 1 ), l ) * bdfem.phi( int( k - 1 ), l ) *
                                                bdfem.weightMeas( l );
                                        }

                                    // Global Dof: row and columns
                                    idDof = pId->bdLocalToGlobal( i ) + ( BCb.component( j ) - 1 ) * totalDof + offset;
                                    jdDof = pId->bdLocalToGlobal( k ) + ( BCb.component( j ) - 1 ) * totalDof + offset;

                                    // Assembling upper entry.  The boundary mas matrix is symetric
                                    A.set_mat_inc( idDof - 1, jdDof - 1, sum );
                                    A.set_mat_inc( jdDof - 1, idDof - 1, sum );
                                }
                        }
                }
        }
}



// ===================================================
// Flux BC
// ===================================================

template <typename MatrixType,
          typename VectorType,
          typename MeshType,
          typename DataType>
void bcFluxManage( MatrixType&     A,
                    VectorType&    b,
                   const MeshType& mesh,
                   const Dof&      dof,
                   const BCBase&   BCb,
                   CurrentBdFE&    bdfem,
                   const DataType& t,
                   UInt            offset )

{

    // Number of local Dof in this face
    UInt nDofF = bdfem.nbNode;

    //offset += BCb.fluxFlag();
    if(!offset)//tricky way to understand if I am in the monolithic case ore not.
        offset += BCb.offset();
    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    DataType sum;

    const IdentifierNatural* pId;
    ID ibF, idDof, jdDof, kdDof;

    const BCFunctionMixte* pBcF = static_cast<const BCFunctionMixte*>( BCb.pointerToFunctor() );

    //std::cout<<"offset "<<offset<<std::endl;
    b.checkAndSet(offset + 1,BCb(t, 0., 0., 0., 1));

    //std::cout<<"step1"<<std::endl;
    if ( !BCb.dataVector() )
        {
            DataType x = 0., y = 0., z = 0.;

            for ( ID i = 1; i <= BCb.list_size(); ++i )
                {
                    pId = static_cast< const IdentifierNatural* >( BCb( i ) );

                    // Number of the current boundary face
                    ibF = pId->id();
                    // Updating face stuff
                    bdfem.updateMeasNormalQuadPt( mesh.bElement( ibF ) );

                    for ( ID idofF = 1; idofF <= nDofF; ++idofF )
                        {
                            ID gDof = pId->bdLocalToGlobal( idofF );

                            for ( int ic = 1; ic <= nComp; ++ic)
                                {
                                    idDof = pId->bdLocalToGlobal( idofF ) + (ic - 1)*totalDof;

                                    sum = 0.;
                                    for ( int iq = 0; iq < bdfem.nbQuadPt; ++iq )
                                        {
                                            sum += bdfem.phi( int( idofF - 1 ), iq )*
                                                bdfem.normal(int(ic - 1), iq)*
                                                bdfem.weightMeas(iq);
                                        }

                                    jdDof = offset;

//                                     std::cout << jdDof + 1 << " " << idDof << " " << sum
//                                               << " " << b.checkLID(jdDof+1) << std::endl;
                                    A.set_mat_inc( idDof - 1, jdDof    , sum );
                                    A.set_mat_inc( jdDof    , idDof - 1, sum );
                                }
                        }
                }
        }
}
}
#endif
