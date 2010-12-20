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
    @brief Functions for prescribing boundary conditions

    @author Miguel Fernandez <miguel.fernandez@inria.fr>
    @contributor Mauro Perego <perego.mauro@gmail.com>
    @maintainer Mauro Perego <perego.mauro@gmail.com>

    @date 7-2002
 *///@HEADER

#ifndef BCMANAGE_H
#define BCMANAGE_H 1

#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/BCNormalManager.hpp>


namespace LifeV
{

// ===================================================
//!@name Boundary conditions treatment
//@{
// ===================================================

//! Prescribe boundary conditions.
/*!
  The matrix and the right hand side are modified to take into account the boundary conditions
  @param matrix   The system matrix
  @param rightHandSide   The system right hand side
  @param mesh  The mesh
  @param dof  Container of the local to global map of DOFs
  @param bcHandler The boundary conditions handler
  @param currentBdFE Current finite element on boundary
  @parma diagonalizeCoef The coefficient used during the system diagonalization
  @param time The time
 */
template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void
bcManage( MatrixType& matrix,
          VectorType& rightHandSide,
          MeshType const& mesh,
          Dof const& dof,
          BCHandler const& bcHandler,
          CurrentBdFE& currentBdFE,
          DataType const& diagonalizeCoef,
          DataType const& time = 0 );



//! Prescribe boundary conditions. Case in which the user defined function depends on the FE vector feVec
/*!
   The matrix and the right hand side are modified to take into account the boundary conditions
   @param mu   		User defined function
   @param matrix   	The system matrix
   @param rightHandSide   The system right hand side
   @param mesh  	The mesh
   @param dof  		Container of the local to global map of DOFs
   @param bcHandler The boundary conditions handler
   @param currentBdFE Current finite element on boundary
   @parma diagonalizeCoef The coefficient used during the system diagonalization
   @param time The time
   @param feVec The finite element vector
 */
template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void
bcManage( Real (*mu)(Real time,Real x, Real y, Real z, Real u),
          MatrixType& matrix,
          VectorType& rightHandSide,
          const MeshType& mesh,
          const Dof& dof,
          const BCHandler& bcHandler,
          CurrentBdFE& currentBdFE,
          const DataType diagonalizeCoef,
          const DataType& time,
          VectorType& feVec );



//! Prescribe boundary conditions. Case in which only the matrix is modified
/*!
   The matrix and the right hand side are modified to take into account the boundary conditions
   @param matrix   The system matrix
   @param rightHandSide   The system right hand side
   @param mesh  The mesh
   @param dof  Container of the local to global map of DOFs
   @param bcHandler The boundary conditions handler
   @param currentBdFE Current finite element on boundary
   @parma diagonalizeCoef The coefficient used during the system diagonalization
  @param time The time
 */
template <typename MatrixType, typename MeshType, typename DataType>
void
bcManageMatrix( MatrixType&      matrix,
                const MeshType&  mesh,
                const Dof&       dof,
                const BCHandler& bcHandler,
                CurrentBdFE&     currentBdFE,
                const DataType&  diagonalizeCoef,
                const DataType&  time = 0);



//! Prescribe boundary conditions. Case in which only the right hand side is modified
/*!
 * The right hand side is modified to take into account the boundary conditions
 * @param rightHandSide   The system right hand side
 * @param mesh  The mesh
 * @param dof  Container of the local to global map of DOFs
 * @param bcHandler The boundary conditions handler
 * @param currentBdFE Current finite element on boundary
 * @param time The time
 * @parma diagonalizeCoef The coefficient used during the system diagonalization
 */
template <typename VectorType, typename MeshType, typename DataType>
void
bcManageVector( VectorType&      rightHandSide,
                const MeshType&  mesh,
                const Dof&       dof,
                const BCHandler& bcHandler,
                CurrentBdFE&     currentBdFE,
                const DataType&  time,
                const DataType&  diagonalizeCoef );



//! Prescribe boundary conditions. Case in which only the right hand side is modified
/*!
 * The Right hand side is modified to take into account the boundary conditions
 * @param rightHandSide   The system right hand side
 * @param feSpace  The finite element space
 * @param bcHandler The boundary conditions handler
 * @param time The time
 * @parma diagonalizeCoef The coefficient used during the system diagonalization
 */
template <typename VectorType, typename DataType, typename Mesh, typename EpetraMap>
void
bcManageVector( VectorType&                     rightHandSide,
                FESpace<Mesh, EpetraMap>&       feSpace,
                const BCHandler&                bcHandler,
                const DataType&                 time,
                const DataType&                 diagonalizeCoef );

//@}


// ===================================================
//! @name Essential BC
// @{
// ===================================================

//! Prescribe Essential boundary conditions. Case in which the user defined function depends on the FE vector feVec
/*!
 * The matrix and the right hand side are modified to take into account the Essential boundary conditions
 * @param matrix   The system matrix
 * @param rightHandSide   The system right hand side
 * @param mesh  The mesh
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @param currentBdFE Current finite element on boundary
 * @parma diagonalizeCoef The coefficient used during the system diagonalization
 * @param time The time
 * @param offset The boundary condition offset
 */
template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void
bcEssentialManage( MatrixType& matrix,
                   VectorType& rightHandSide,
                   const MeshType& /*mesh*/,
                   const Dof& dof,
                   const BCBase& boundaryCond,
                   const CurrentBdFE& /*currentBdFE*/,
                   const DataType& diagonalizeCoef,
                   const DataType& time,
                   UInt offset );



//! Prescribe Essential boundary conditions. Case in which the user defined function depends on the FE vector feVec
/*!
 * The matrix and the right hand side are modified to take into account the Essential boundary conditions
 * @param matrix   The system matrix
 * @param rightHandSide   The system right hand side
 * @param mesh  The mesh
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @param currentBdFE Current finite element on boundary
 * @parma diagonalizeCoef The coefficient used during the system diagonalization
 * @param time The time
 * @param feVec The finite element vector
 * @param offset The bcCond offset
 */
template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void
bcEssentialManageUDep( MatrixType& matrix,
                       VectorType& rightHandSide,
                       const MeshType& /*mesh*/,
                       const Dof& dof,
                       const BCBase& boundaryCond,
                       const CurrentBdFE& /*currentBdFE*/,
                       const DataType& diagonalizeCoef,
                       const DataType& time,
                       const VectorType& feVec ,
                       UInt offset=0);



//! Prescribe Essential boundary conditions diagonalizing the matrix
/*!
 * The matrix is modified to take into account the Essential boundary conditions
 * @param matrix   The system matrix
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @parma diagonalizeCoef The coefficient used during the system diagonalization
 * @param offset The boundary condition offset
 */
template <typename MatrixType, typename DataType>
void
bcEssentialManageMatrix( MatrixType& matrix,
                         const Dof& dof,
                         const BCBase& boundaryCond,
                         const DataType& diagonalizeCoef,
                         UInt offset );



//! Prescribe Essential boundary conditions on the right hand side
/*!
 * The right hand side is modified to take into account the Essential boundary conditions
 * @param rightHandSide   The system rightHandSide
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @parma diagonalizeCoef The coefficient used during the system diagonalization
 * @param offset The boundary condition offset
 */
template <typename VectorType, typename DataType>
void
bcEssentialManageVector( VectorType&     rightHandSide,
                         const Dof&      dof,
                         const BCBase&   boundaryCond,
                         const DataType& time,
                         const DataType& diagonalizeCoef,
                         UInt            offset );



///! Prescribe Essential boundary conditions.
/*!
 * The matrix and the right hand side are modified to take into account the Essential boundary conditions
 * @param matrix   The system matrix
 * @param dof  Container of the local to global map of DOFs
 * @param bcHandler The boundary condition handler
 * @parma diagonalizeCoef The coefficient used during the system diagonalization
 */
template <typename MatrixType, typename DataType>
void
bcManageMtimeUDep( MatrixType& matrix,
                   const Dof& dof,
                   const BCHandler& bcHandler,
                   const DataType diagonalizeCoef);

// @}

// ===================================================
//! @name Natural BC
// @{
// ===================================================


//! Prescribe Natural boundary condition
/*!
 * The right hand side is modified to take into account the Natural boundary condition
 * @param rightHandSide   The system right hand side
 * @param mesh  The mesh
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @param currentBdFE Current finite element on boundary
 * @param time The time
  * @param offset The boundary condition offset
 */
template <typename VectorType, typename MeshType, typename DataType>
void
bcNaturalManage( VectorType& rightHandSide,
                 const MeshType& mesh,
                 const Dof& dof, const
                 BCBase& boundaryCond,
                 CurrentBdFE& currentBdFE,
                 const DataType& time,
                 UInt offset );




//! Prescribe Natural boundary condition. Case in which the user defined function depends on the FE vector feVec
/*!
 * The right hand side is modified to take into account the Natural boundary condition
 * @param mu   User defined function
 * @param rightHandSide   The system right hand side
 * @param mesh  The mesh
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @param currentBdFE Current finite element on boundary
 * @param time The time
 * @param offset The boundary condition offset
 */
template <typename VectorType, typename MeshType, typename DataType>
void
bcNaturalManageUDep( Real (*mu)(Real time,Real x, Real y, Real z, Real u),
                     VectorType& rightHandSide,
                     const MeshType& mesh,
                     const Dof& dof,
                     const BCBase& boundaryCond,
                     CurrentBdFE& currentBdFE,
                     const DataType& time,
                     const VectorType& feVec,
                     UInt offset );

// @}


// ===================================================
//! @name Mixte BC
// @{
// ===================================================


//! Prescribe Mixte boundary condition
/*!
 * The matrix and the right hand side are modified to take into account the Mixte boundary condition
 * @param matrix   The system matrix
 * @param rightHandSide   The system right hand side
 * @param mesh  The mesh
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @param currentBdFE Current finite element on boundary
 * @param time The time
 * @param offset The boundary condition offset
 */
template <typename MatrixType, typename VectorType, typename DataType, typename MeshType>
void
bcMixteManage( MatrixType& matrix,
               VectorType& rightHandSide,
               const MeshType& mesh,
               const Dof& dof,
               const BCBase& boundaryCond,
               CurrentBdFE& currentBdFE,
               const DataType& time,
               UInt offset );


//! Prescribe Mixte boundary condition only on the matrix
/*!
 * The matrix is modified to take into account the Mixte boundary condition
 * @param matrix   The system matrix
 * @param mesh  The mesh
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @param currentBdFE Current finite element on boundary
 * @param time The time
 * @param offset The boundary condition offset
 */
template <typename MatrixType, typename DataType, typename MeshType>
void
bcMixteManageMatrix( MatrixType& matrix,
                     const MeshType& mesh,
                     const Dof& dof,
                     const BCBase& boundaryCond,
                     CurrentBdFE& currentBdFE,
                     const DataType& time,
                     UInt offset );


//! Prescribe Mixte boundary condition only on the rightHandSide
/*!
 * The matrix is modified to take into account the Mixte boundary condition
 * @param rightHandSide   The system right hand side
 * @param mesh  The mesh
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @param currentBdFE Current finite element on boundary
 * @param time The time
 * @param offset The boundary condition offset
 */
template <typename VectorType, typename DataType, typename MeshType>
void
bcMixteManageVector( VectorType& rightHandSide,
                     const MeshType& mesh,
                     const Dof& dof,
                     const BCBase& boundaryCond,
                     CurrentBdFE& currentBdFE,
                     const DataType& time,
                     UInt offset );

// @}


// ===================================================
//!@name Flux BC
//@{
// ===================================================


//! Prescribe Flux boundary condition only on the matrix
/*!
 * The matrix is modified to take into account the Flux boundary condition
 * @param matrix   The system matrix
 * @param rightHandSide The system right hand side
 * @param mesh  The mesh
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @param currentBdFE Current finite element on boundary
 * @param time The time
 * @param offset The boundary condition offset
 */
template <typename MatrixType,
typename VectorType,
typename MeshType,
typename DataType>
void
bcFluxManage( MatrixType&     matrix,
              VectorType&    rightHandSide,
              const MeshType& mesh,
              const Dof&      dof,
              const BCBase&   boundaryCond,
              CurrentBdFE&    currentBdFE,
              const DataType& time,
              UInt            offset);


//! Prescribe Flux boundary condition only on the right hand side
/*!
 * The matrix is modified to take into account the Flux boundary condition
 * @param rightHandSide The system right hand side
 * @param mesh  The mesh
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @param currentBdFE Current finite element on boundary
 * @param time The time
 * @param offset The boundary condition offset
 */
template <typename VectorType,
typename DataType>
void
bcFluxManageVector(
    VectorType&    rightHandSide,
    const BCBase&   boundaryCond,
    const DataType& time,
    UInt            offset);


//! Prescribe Flux boundary condition only on the matrix
/*!
 * The matrix is modified to take into account the Flux boundary condition
 * @param matrix   The system matrix
 * @param mesh  The mesh
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @param currentBdFE Current finite element on boundary
 * @param time The time
 * @param offset The boundary condition offset
 */
template <typename MatrixType,
typename MeshType,
typename DataType>
void
bcFluxManageMatrix( MatrixType&     matrix,
                    const MeshType& mesh,
                    const Dof&      dof,
                    const BCBase&   boundaryCond,
                    CurrentBdFE&    currentBdFE,
                    const DataType& /*time*/,
                    UInt            offset );
// @}



// ===================================================
//!@name Resistance BC
// @{
// ===================================================


//! Prescribe Resistance boundary condition
/*!
 * The matrix and the right hand side are modified to take into account the Resistance boundary condition
 * @param rightHandSide The system right hand side
 * @param mesh  The mesh
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @param currentBdFE Current finite element on boundary
 * @param time The time
 * @param offset The boundary condition offset
 */
template <typename MatrixType, typename VectorType, typename DataType, typename MeshType>
void
bcResistanceManage( MatrixType& matrix,
                    VectorType& rightHandSide,
                    const MeshType& mesh,
                    const Dof& dof,
                    const BCBase& boundaryCond,
                    CurrentBdFE& currentBdFE,
                    const DataType& /*time*/,
                    UInt offset );

// @}








// ===================================================
// Implementation
// ===================================================


// ===================================================
// Boundary conditions treatment
// ===================================================

template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void
bcManage( MatrixType& matrix,
          VectorType& rightHandSide,
          MeshType const& mesh,
          Dof const& dof,
          BCHandler const& bcHandler,
          CurrentBdFE& currentBdFE,
          DataType const& diagonalizeCoef,
          DataType const& time = 0 )
{

    VectorType rhsRepeated(rightHandSide.map(), Repeated);
    bool globalassemble=false;


    BCManageNormal<MatrixType> bcManageNormal;


    // Loop on boundary conditions
    for ( ID i = 0; i < bcHandler.size(); ++i )
    {
        switch ( bcHandler[ i ].type() )
        {
        case Essential:  // Essential boundary conditions (Dirichlet)
            //Normal, Tangential or Directional boundary conditions
            if ( (bcHandler[ i ].mode() == Tangential) || (bcHandler[ i ].mode() == Normal) || (bcHandler[ i ].mode() == Directional) )
            {
                bcManageNormal.init(bcHandler[ i ],time); //initialize bcManageNormal
            }
        case EssentialEdges:
        case EssentialVertices:
            globalassemble=true;
            break;
        case Natural:    // Natural boundary conditions (Neumann)
            bcNaturalManage( rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset());
            break;
        case Mixte:      // Mixte boundary conditions (Robin)
            bcMixteManage( matrix, rhsRepeated, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            break;
        case Flux:       // Flux boundary condition
            bcFluxManage( matrix, rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset()+bcHandler[i].offset());
            break;
        case Resistance: // Resistance boundary condition
            bcResistanceManage( matrix, rhsRepeated, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }

    rhsRepeated.globalAssemble();
    rightHandSide += rhsRepeated;
    if (globalassemble)
        matrix.globalAssemble();

    //Build the internal structure, if needed
    bcManageNormal.build(mesh, dof,currentBdFE,matrix,bcHandler.offset(),rightHandSide.mapPtr()->commPtr());
    bcManageNormal.exportToParaview("normalAndTangents");

    //Applying the basis change, if needed
    bcManageNormal.bcShiftToNormalTangentialCoordSystem(matrix, rightHandSide);

    // Loop on boundary conditions
    for ( ID i = 0; i < bcHandler.size(); ++i )
    {
        switch ( bcHandler[ i ].type() )
        {
        case Essential:  // Essential boundary conditions (Dirichlet)
        case EssentialEdges:
        case EssentialVertices:
            bcEssentialManage( matrix, rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, diagonalizeCoef, time, bcHandler.offset() );
            break;
        case Natural:   	// Natural boundary conditions (Neumann)
        case Mixte:  	    // Mixte boundary conditions (Robin)
        case Flux:		    // Flux boundary condition
        case Resistance:    // Resistance boundary conditions
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }

    //Return back to the initial basis
    bcManageNormal.bcShiftToCartesianCoordSystem(matrix, rightHandSide);
}



template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void
bcManage( Real (*mu)(Real time,Real x, Real y, Real z, Real u),
          MatrixType& matrix,
          VectorType& rightHandSide,
          const MeshType& mesh,
          const Dof& dof,
          const BCHandler& bcHandler,
          CurrentBdFE& currentBdFE,
          const DataType diagonalizeCoef,
          const DataType& time,
          VectorType& feVec )
{
    VectorType rhsRepeated(rightHandSide.getMap(),Repeated);

    bool globalassemble=false;
    // Loop on boundary conditions
    for ( ID i = 0; i < bcHandler.size(); ++i )
    {
        switch ( bcHandler[ i ].type() )
        {
        case Essential:  // Essential boundary conditions (Dirichlet)
            if ( (bcHandler[ i ].mode() == Tangential) || (bcHandler[ i ].mode() == Normal) || (bcHandler[ i ].mode() == Directional) )
            {
                ERROR_MSG( "This BC mode is not yet implemented for this setting" );
            }
        case EssentialEdges:
        case EssentialVertices:
            globalassemble=true;
            break;
        case Natural:  // Natural boundary conditions (Neumann)
            if (bcHandler[ i ].isUDep())
                bcNaturalManageUDep(mu, rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, time,feVec, bcHandler.offset());
            else
                //in this case mu must be a constant, think about (not still implemented)
                bcNaturalManage( rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            break;
        case Mixte:  // Mixte boundary conditions (Robin)

            if (bcHandler[ i ].isUDep())
            {
                ERROR_MSG( "This BC mode is not yet implemented for this setting" );    //not implemented yet
            }
            else
            {
                bcMixteManage( matrix, rhsRepeated, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            }
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }

    rhsRepeated.GlobalAssemble();

    rightHandSide += rhsRepeated;
    if (globalassemble)
        matrix.GlobalAssemble();


    // Loop on boundary conditions
    for ( ID i = 0; i < bcHandler.size(); ++i )
    {

        switch ( bcHandler[ i ].type() )
        {
        case Essential:  // Essential boundary conditions (Dirichlet)
        case EssentialEdges:
        case EssentialVertices:
            if (bcHandler[ i ].isUDep())
                bcEssentialManageUDep(matrix, rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, diagonalizeCoef, time,feVec, bcHandler.offset());
            else
                bcEssentialManage( matrix, rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, diagonalizeCoef, time, bcHandler.offset() );
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


template <typename MatrixType, typename MeshType, typename DataType>
void
bcManageMatrix( MatrixType&      matrix,
                const MeshType&  mesh,
                const Dof&       dof,
                const BCHandler& bcHandler,
                CurrentBdFE&     currentBdFE,
                const DataType&  diagonalizeCoef,
                const DataType&  time = 0)
{

    bool globalassemble=false;
    // Loop on boundary conditions
    for ( ID i = 0; i < bcHandler.size(); ++i )
    {

        switch ( bcHandler[ i ].type() )
        {
        case Essential:  // Essential boundary conditions (Dirichlet)
        case EssentialEdges:
        case EssentialVertices:
            globalassemble=true;
            if ( (bcHandler[ i ].mode() == Tangential) || (bcHandler[ i ].mode() == Normal) || (bcHandler[ i ].mode() == Directional) )
            {
                ERROR_MSG( "This BC mode is not yet implemented for this setting" );
            }
            break;
        case Natural:  // Natural boundary conditions (Neumann)
            break;
        case Mixte:  // Mixte boundary conditions (Robin)
            bcMixteManageMatrix( matrix, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            break;
        case Flux:  // Flux boundary conditions
            bcFluxManageMatrix( matrix, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset()+bcHandler[i].offset());
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }
    if (globalassemble)
        matrix.globalAssemble();

    // Loop on boundary conditions
    for ( ID i = 0; i < bcHandler.size(); ++i )
    {

        switch ( bcHandler[ i ].type() )
        {
        case Essential:  // Essential boundary conditions (Dirichlet)
        case EssentialEdges:
        case EssentialVertices:
            bcEssentialManageMatrix( matrix, dof, bcHandler[ i ], diagonalizeCoef, bcHandler.offset() );  //! Bug here???
            break;
        case Natural:  // Natural boundary conditions (Neumann)
            // Do nothing
            break;
        case Mixte:  // Mixte boundary conditions (Robin)
            break;
        case Flux:  // Mixte boundary conditions (Robin)
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }
}


template <typename VectorType, typename MeshType, typename DataType>
void
bcManageVector( VectorType&      rightHandSide,
                const MeshType&  mesh,
                const Dof&       dof,
                const BCHandler& bcHandler,
                CurrentBdFE&     currentBdFE,
                const DataType&  time,
                const DataType&  diagonalizeCoef )
{
    VectorType rhsRepeated(rightHandSide.map(),Repeated);

    // Loop on boundary conditions
    for ( ID i = 0; i < bcHandler.size(); ++i )
    {

        switch ( bcHandler[ i ].type() )
        {
        case Essential:  // Essential boundary conditions (Dirichlet)
        case EssentialEdges:
        case EssentialVertices:
            if ( (bcHandler[ i ].mode() == Tangential) || (bcHandler[ i ].mode() == Normal) || (bcHandler[ i ].mode() == Directional) )
            {
                ERROR_MSG( "This BC mode is not yet implemented for this setting" );
            }
            bcEssentialManageVector( rightHandSide, dof, bcHandler[ i ], time, diagonalizeCoef, bcHandler.offset() );
            break;
        case Natural:  // Natural boundary conditions (Neumann)
            bcNaturalManage( rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            break;
        case Mixte:  // Mixte boundary conditions (Robin)
            bcMixteManageVector( rhsRepeated, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            break;
        case Flux:  // Flux boundary conditions
            bcFluxManageVector( rightHandSide, bcHandler[ i ], time, bcHandler.offset()+bcHandler[i].offset() );
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }

    rhsRepeated.globalAssemble();

    rightHandSide += rhsRepeated;
}


template <typename VectorType, typename DataType, typename Mesh, typename EpetraMap>
void
bcManageVector( VectorType&                     rightHandSide,
                FESpace<Mesh, EpetraMap>&       feSpace,
                const BCHandler&                bcHandler,
                const DataType&                 time,
                const DataType&                 diagonalizeCoef )
{
    VectorType rhsRepeated(rightHandSide.getMap(),Repeated);

    // Loop on boundary conditions
    for ( ID i = 0; i < bcHandler.size(); ++i )
    {

        switch ( bcHandler[ i ].type() )
        {
        case Essential:  // Essential boundary conditions (Dirichlet)
        case EssentialEdges:
        case EssentialVertices:
            if ( (bcHandler[ i ].mode() == Tangential) || (bcHandler[ i ].mode() == Normal) || (bcHandler[ i ].mode() == Directional) )
            {
                ERROR_MSG( "This BC mode is not yet implemented for this setting" );
            }
            bcEssentialManageVector( rightHandSide, feSpace.dof(), bcHandler[ i ], time, diagonalizeCoef, bcHandler.offset() );
            break;
        case Natural:  // Natural boundary conditions (Neumann)
            bcNaturalManage( rightHandSide, *feSpace.mesh(), feSpace.dof(), bcHandler[ i ], feSpace.feBd(), time, bcHandler.offset() );
            break;
        case Mixte:  // Mixte boundary conditions (Robin)
            bcMixteManageVector( rhsRepeated, *feSpace.mesh(), feSpace.dof(), bcHandler[ i ], feSpace.feBd(), time, bcHandler.offset() );
            break;
        case Flux:  // Flux boundary conditions
            bcFluxManageVector( rightHandSide, bcHandler[ i ], time, bcHandler.offset()+bcHandler[i].offset() );
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }

    rhsRepeated.GlobalAssemble();

    rightHandSide += rhsRepeated;
}



// ===================================================
// Essential BC
// ===================================================

template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void
bcEssentialManage( MatrixType& matrix,
                   VectorType& rightHandSide,
                   const MeshType& /*mesh*/,
                   const Dof& dof,
                   const BCBase& boundaryCond,
                   const CurrentBdFE& /*currentBdFE*/,
                   const DataType& diagonalizeCoef,
                   const DataType& time,
                   UInt offset )
{

    ID idDof;

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();


    std::vector<ID>   idDofVec(0);
    std::vector<Real> datumVec(0);

    idDofVec.reserve(boundaryCond.list_size()*nComp);
    datumVec.reserve(boundaryCond.list_size()*nComp);

    if ( boundaryCond.dataVector() )
    { //! If BC is given under a vectorial form

        const BCVectorInterface* pId = static_cast< const BCVectorInterface* > (boundaryCond.pointerToBCVector());
        assert( pId != 0);

        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {

            if ( !pId->dofInterface().isMyInterfaceDof(boundaryCond( i ) ->id()))
            {
                continue;
            }

            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {
                // Global Dof
                idDof = boundaryCond( i )->id() + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                datumVec.push_back(boundaryCond( boundaryCond( i ) ->id(), boundaryCond.component( j ) ));
                idDofVec.push_back(idDof - 1);
            }
        }
    }
    else
    { //! If BC is given under a functional form

        DataType x, y, z;
        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {
            // Coordinates of the node where we impose the value
            x = static_cast< const IdentifierEssential* >( boundaryCond( i ) ) ->x();
            y = static_cast< const IdentifierEssential* >( boundaryCond( i ) ) ->y();
            z = static_cast< const IdentifierEssential* >( boundaryCond( i ) ) ->z();

            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {
                // Global Dof
                idDof = boundaryCond( i ) ->id() + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                datumVec.push_back(boundaryCond( time, x, y, z, boundaryCond.component( j ) ));
                idDofVec.push_back(idDof - 1);

            }
        }
    }

    // Modifying matrix and right hand side
    matrix.diagonalize( idDofVec, diagonalizeCoef, rightHandSide, datumVec);

}


template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void
bcEssentialManageUDep( MatrixType& matrix,
                       VectorType& rightHandSide,
                       const MeshType& /*mesh*/,
                       const Dof& dof,
                       const BCBase& boundaryCond,
                       const CurrentBdFE& /*currentBdFE*/,
                       const DataType& diagonalizeCoef,
                       const DataType& time,
                       const VectorType& feVec ,
                       UInt offset=0)
{

    ID idDof;

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    if ( boundaryCond.dataVector() )
    { //! If BC is given under a vectorial form

        //not possible
        ERROR_MSG( "This type of BCVector does not exists on bc depentent on solution" );
    }
    else
    { //! If BC is given under a functional form

        std::vector<ID>   idDofVec(0);
        std::vector<Real> datumVec(0);

        idDofVec.reserve(boundaryCond.list_size()*nComp);
        datumVec.reserve(boundaryCond.list_size()*nComp);

        DataType x, y, z;
        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {
            // Coordinates of the node where we impose the value
            x = static_cast< const IdentifierEssential* >( boundaryCond( i ) ) ->x();
            y = static_cast< const IdentifierEssential* >( boundaryCond( i ) ) ->y();
            z = static_cast< const IdentifierEssential* >( boundaryCond( i ) ) ->z();

            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {
                // Global Dof
                idDof = boundaryCond( i ) ->id() + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                Real datum = boundaryCond( time, x, y, z, boundaryCond.component( j ) ,feVec[idDof-1]);

                datumVec.push_back(datum);
                idDofVec.push_back(idDof-1);

            }
        }

        // Modifying matrix and right hand side
        matrix.diagonalize( idDofVec, diagonalizeCoef, rightHandSide, datumVec);
    }
}


template <typename MatrixType, typename DataType>
void
bcEssentialManageMatrix( MatrixType& matrix,
                         const Dof& dof,
                         const BCBase& boundaryCond,
                         const DataType& diagonalizeCoef,
                         UInt offset )
{

    ID idDof;
    UInt totalDof;

    // Number of total scalar Dof
    totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    std::vector<ID>   idDofVec(0);
    idDofVec.reserve(boundaryCond.list_size()*nComp);

    // Loop on BC identifiers
    for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
    {
        // Loop on components involved in this boundary condition
        for ( ID j = 1; j <= nComp; ++j )
        {
            // Global Dof
            idDof = boundaryCond( i ) ->id() + ( boundaryCond.component( j ) - 1 ) * totalDof;
            idDofVec.push_back(idDof-1);
            // Modifying ONLY matrix
            //matrix.diagonalize( idDof - 1, diagonalizeCoef );
        }
    }
    // Modifying ONLY matrix
    matrix.diagonalize( idDofVec, diagonalizeCoef, offset);

}


template <typename VectorType, typename DataType>
void
bcEssentialManageVector( VectorType&     rightHandSide,
                         const Dof&      dof,
                         const BCBase&   boundaryCond,
                         const DataType& time,
                         const DataType& diagonalizeCoef,
                         UInt            offset )
{
    ID idDof;
    UInt totalDof;

    // Number of total scalar Dof
    totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    std::vector<int>   idDofVec(0);
    idDofVec.reserve(boundaryCond.list_size()*nComp);
    std::vector<Real> datumVec(0);
    datumVec.reserve(boundaryCond.list_size()*nComp);

    if ( boundaryCond.dataVector() )
    {  //! If BC is given under a vectorial form

        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {
            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {

                // Global Dof
                idDof = boundaryCond( i ) ->id() + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                idDofVec.push_back( idDof );
                datumVec.push_back( diagonalizeCoef * boundaryCond( boundaryCond( i ) ->id(), boundaryCond.component( j ) ) );
            }
        }
    }
    else
    {  //! If BC is given under a functional form

        DataType x, y, z;
        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {
            // Coordinates of the node where we impose the value
            x = static_cast< const IdentifierEssential* >( boundaryCond( i ) ) ->x();
            y = static_cast< const IdentifierEssential* >( boundaryCond( i ) ) ->y();
            z = static_cast< const IdentifierEssential* >( boundaryCond( i ) ) ->z();

            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {
                // Global Dof

                idDof = boundaryCond( i ) ->id() + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;
                // Modifying right hand side
                idDofVec.push_back(idDof);
                datumVec.push_back( diagonalizeCoef * boundaryCond( time, x, y, z, boundaryCond.component( j ) ) );
            }
        }
    }

    rightHandSide.setCoefficients( idDofVec, datumVec);
}

// ===================================================
// Natural BC
// ===================================================


template <typename VectorType, typename MeshType, typename DataType>
void
bcNaturalManage( VectorType& rightHandSide,
                 const MeshType& mesh,
                 const Dof& dof, const
                 BCBase& boundaryCond,
                 CurrentBdFE& currentBdFE,
                 const DataType& time,
                 UInt offset )
{

    // Number of local Dof (i.e. nodes) in this face
    UInt nDofF = currentBdFE.nbNode();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    const IdentifierNatural* pId;
    ID ibF, idDof, icDof, gDof;
    Real sum;

    if ( boundaryCond.dataVector() )
    {
        //! If BC is given under a vectorial form
        switch ( boundaryCond.pointerToBCVector()->type() )
        {
        case 0:  // if the BC is a vector which values don'time need to be integrated
        {
            VectorType bUnique(rightHandSide.map(),Unique);

            std::vector<int>  idDofVec(0);
            idDofVec.reserve(boundaryCond.list_size()*nComp);
            std::vector<Real> datumVec(0);
            datumVec.reserve(boundaryCond.list_size()*nComp);
            // double datum;

            // Loop on BC identifiers
            for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
            {
                // Loop on components involved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {
                    ID id = boundaryCond(i)->id();

                    // Global Dof
                    idDof = id + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                    idDofVec.push_back( idDof);

                    // Modifying right hand side (assuming BCvector is a flux)
                    datumVec.push_back( boundaryCond( id , boundaryCond.component( j ) ));
                }
            }

            bUnique.setCoefficients(idDofVec, datumVec);
            bUnique.globalAssemble(Insert);

            ASSERT( rightHandSide.mapType() == Unique , "rightHandSide must have unique map, otherwise data will be multiply added on cpu interfaces." );
            rightHandSide += bUnique;

        }
        break;

        case 1:  // if the BC is a vector of values to be integrated
        {
            VectorType rhsRepeated(rightHandSide.map(),Repeated);

            // Loop on BC identifiers
            for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
            {
                // Pointer to the i-th itdentifier in the list
                pId = static_cast< const IdentifierNatural* >( boundaryCond( i ) );

                // Number of the current boundary face
                ibF = pId->id();

                // Updating face stuff
                currentBdFE.updateMeasNormalQuadPt( mesh.bElement( ibF ) );

                // Loop on total Dof per Face
                for ( ID l = 1; l <= nDofF; ++l )
                {

                    gDof = pId->localToGlobalMap( l );

                    // Loop on components involved in this boundary condition
                    for ( UInt ic = 0; ic < nComp; ++ic )
                    {
                        icDof = gDof + ic * totalDof + offset;

                        // Loop on quadrature points
                        for ( int iq = 0; iq < (int)currentBdFE.nbQuadPt(); ++iq )
                        {
                            sum=0.0;
                            // data on quadrature point
                            for ( ID m = 1; m <= nDofF; ++m )
                                sum +=  boundaryCond( pId->localToGlobalMap( m ) , 1 ) * currentBdFE.phi( int( m - 1 ), iq );
                            // Adding right hand side contribution
                            rhsRepeated[ icDof ] += sum * currentBdFE.phi( int( l - 1 ), iq ) * currentBdFE.normal( int( ic ), iq )
                                                    * currentBdFE.weightMeas( iq );
                        }
                    }
                }
            }

            rhsRepeated.globalAssemble();
            ASSERT( rightHandSide.mapType() == Unique , "here rightHandSide should passed as repeated, otherwise not sure of what happens at the cpu interfaces ." );
            rightHandSide += rhsRepeated;
        }
        break;
        case 2:  // if the BC is a vector of values with components to be integrated
        {
            VectorType rhsRepeated(rightHandSide.map(),Repeated);

            // Loop on BC identifiers
            for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
            {
                // Pointer to the i-th itdentifier in the list
                pId = static_cast< const IdentifierNatural* >( boundaryCond( i ) );

                // Number of the current boundary face
                ibF = pId->id();

                // Updating face stuff
                currentBdFE.updateMeasNormalQuadPt( mesh.bElement( ibF ) );

                // Loop on total Dof per Face
                for ( ID idofF = 1; idofF <= nDofF; ++idofF )
                {

                    gDof = pId->localToGlobalMap( idofF );

                    // Loop on space dimensions
                    for ( ID ic = 1; ic <= nComp; ++ic )
                        //  for ( UInt ic = 0; ic < nDimensions; ++ic )
                    {
                        icDof = gDof + ( boundaryCond.component( ic ) - 1 ) * totalDof+ offset;   //Components passed separately

                        // Loop on quadrature points
                        for ( int iq = 0; iq < (int)currentBdFE.nbQuadPt(); ++iq )
                        {
                            sum = 0;
                            // data on quadrature point
                            for ( ID m = 1; m <= nDofF; ++m )
                                sum +=  boundaryCond( pId->localToGlobalMap( m ) , boundaryCond.component( ic ) ) * currentBdFE.phi( int( m - 1 ), iq );  //Components passed separatedly

                            // Adding right hand side contribution
                            rhsRepeated[ icDof ] += sum *  currentBdFE.phi( int( idofF - 1 ), iq ) *
                                                    currentBdFE.weightMeas( iq );
                        }
                    }
                }
            }
            rhsRepeated.globalAssemble();
            ASSERT( rightHandSide.mapType() == Unique , "here rightHandSide should passed as unique, otherwise not sure of what happens at the cpu interfaces ." );
            rightHandSide += rhsRepeated;
        }
        break;
        default:
            ERROR_MSG( "This type of BCVector does not exist" );
        }
    }

    else
    {  //! If BC is given under a functional form

        DataType x, y, z;
        VectorType rhsRepeated(rightHandSide.map(),Repeated);

        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {
            // Pointer to the i-th identifier in the list
            pId = static_cast< const IdentifierNatural* >( boundaryCond( i ) );
            // Number of the current boundary face
            ibF = pId->id();
            // Updating face stuff
            currentBdFE.updateMeasNormalQuadPt( mesh.bElement( ibF ) );
            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            {
                for ( ID j = 1; j <= nComp; ++j )
                {
                    //global Dof
                    idDof = pId->localToGlobalMap( idofF ) + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;
                    // Loop on quadrature points
                    for ( int iq = 0; iq < (int)currentBdFE.nbQuadPt(); ++iq )
                    {
                        currentBdFE.coorQuadPt( x, y, z, iq ); // quadrature point coordinates
                        switch (boundaryCond.mode())
                        {
                        case Full:
                            rhsRepeated[ idDof ] += currentBdFE.phi( int( idofF - 1 ), iq ) * boundaryCond( time, x, y, z, boundaryCond.component( j ) ) *
                                                    currentBdFE.weightMeas( iq ); // BASEINDEX + 1
                            break;
                        case Component:
                            rhsRepeated[ idDof ] += currentBdFE.phi( int( idofF - 1 ), iq ) * boundaryCond( time, x, y, z, boundaryCond.component( j ) ) *
                                                    currentBdFE.weightMeas( iq ); // BASEINDEX + 1
                            break;
                        case Normal:
                            rhsRepeated[ idDof ] += boundaryCond( time, x, y, z, boundaryCond.component( j ) )*
                                                    currentBdFE.phi( int( idofF - 1 ), iq )*
                                                    currentBdFE.weightMeas( iq )*currentBdFE.normal( int(j - 1), iq );
                            break;
                        default:
                            ERROR_MSG( "This BC mode is not (yet) implemented" );

                        }
                    }
                }
            }
        }
        rhsRepeated.globalAssemble();
        ASSERT( rightHandSide.mapType() == Unique , "here rightHandSide should passed as unique, otherwise not sure of what happens at the cpu interfaces ." );
        //  rightHandSide=rhsRepeated;
        rightHandSide+= rhsRepeated;
    }
} // bcNaturalManage



template <typename VectorType, typename MeshType, typename DataType>
void
bcNaturalManageUDep( Real (*mu)(Real time,Real x, Real y, Real z, Real u),
                     VectorType& rightHandSide,
                     const MeshType& mesh,
                     const Dof& dof,
                     const BCBase& boundaryCond,
                     CurrentBdFE& currentBdFE,
                     const DataType& time,
                     const VectorType& feVec,
                     UInt offset )
{

    // Number of local Dof (i.e. nodes) in this face
    UInt nDofF = currentBdFE.nbNode;

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    const IdentifierNatural* pId;
    ID ibF, idDof ;

    if ( boundaryCond.dataVector() )
    { //! If BC is given under a vectorial form
        ERROR_MSG( "This type of BCVector does not exists on bc depentent on solution\n" );
    }
    else
    {  //! If BC is given under a functional form

        DataType x, y, z;

        if (nComp!=1)
        {
            ERROR_MSG("For now bcNaturalManageUDep cannot handle non scalar solutions\n");
        }

        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {

            // Pointer to the i-th itdentifier in the list
            pId = static_cast< const IdentifierNatural* >( boundaryCond( i ) );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            currentBdFE.updateMeas( mesh.boundaryFace( ibF ) );

            std::vector<Real> locU(nDofF+1);    //assumes feVec is a vec of reals, TODO: deal with more comp
            Real uPt;            //value in the point
            for (ID idofLocU=0; idofLocU<nDofF; idofLocU++)
            {
                ID idGDofU=pId->bdLocalToGlobal(idofLocU+1)+( boundaryCond.component( 1 ) - 1 ) * totalDof + offset;
                locU[idofLocU]=feVec[idGDofU-1];
            }


            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            {
                // Loop on components involved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {

                    //global Dof
                    idDof = pId->bdLocalToGlobal( idofF ) + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( int l = 0; l < (int)currentBdFE.nbQuadPt; ++l )
                    {
                        currentBdFE.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                        uPt=0.0;
                        for (ID idofLocU=0; idofLocU<nDofF; idofLocU++)
                        {
                            uPt+=locU[idofLocU]*currentBdFE.phi( int( idofLocU  ),l );
                        }

                        // Adding right hand side contribution
                        rightHandSide[ idDof ] += currentBdFE.phi( int( idofF - 1 ), l ) * boundaryCond( time, x, y, z, boundaryCond.component( j ),uPt ) *
                                                  mu(time,x,y,z,uPt)*currentBdFE.weightMeas( l ); // BASEINDEX + 1
                    }
                }
            }
        }
    }
}




// ===================================================
// Mixte BC
// ===================================================


template <typename MatrixType, typename VectorType, typename DataType, typename MeshType>
void
bcMixteManage( MatrixType& matrix,
               VectorType& rightHandSide,
               const MeshType& mesh,
               const Dof& dof,
               const BCBase& boundaryCond,
               CurrentBdFE& currentBdFE,
               const DataType& time,
               UInt offset )
{

    // Number of local Dof in this face
    UInt nDofF = currentBdFE.nbNode();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    DataType sum;

    const IdentifierNatural* pId;
    ID ibF, idDof, jdDof, kdDof;

    if ( boundaryCond.dataVector() )
    {   //! If BC is given under a vectorial form

        //! for the moment, only one coefficient per BCvector.
        DataType mcoef, mbcb;

        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {

            // Pointer to the i-th identifier in the list
            pId = static_cast< const IdentifierNatural* >( boundaryCond( i ) );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            currentBdFE.updateMeas( mesh.bElement( ibF ) );

            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            {
                // Loop on components involved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {

                    sum = 0;

                    idDof = pId->localToGlobalMap( idofF ) + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                    {
                        mcoef = 0.0;
                        mbcb = 0.0;
                        for ( ID n = 1; n <= nDofF; ++n)
                        {
                            kdDof=pId->localToGlobalMap( n ); // + ( boundaryCond.component( j ) - 1 ) * totalDof;
                            if (boundaryCond.ismixteVec())
                                mcoef += boundaryCond.MixteVec( kdDof, boundaryCond.component( j ) ) * currentBdFE.phi( int( n - 1 ), l );
                            else  mcoef += boundaryCond.mixteCoef() * currentBdFE.phi( int( n - 1 ), l );

                            if (boundaryCond.isbetaVec())
                                mbcb += boundaryCond.BetaVec( kdDof, boundaryCond.component( j ) )
                                        * boundaryCond( kdDof, boundaryCond.component( j )) * currentBdFE.phi( int( n - 1 ), l );
                            else  mbcb += boundaryCond.betaCoef() * boundaryCond( kdDof, boundaryCond.component( j )) * currentBdFE.phi( int( n - 1 ), l );
                        }

                        sum += mcoef* currentBdFE.phi( int( idofF - 1 ), l ) *
                               currentBdFE.phi( int( idofF - 1 ), l ) *currentBdFE.weightMeas( l );

                        rightHandSide[ idDof ] += currentBdFE.phi( int( idofF - 1 ), l ) * mbcb * // BASEINDEX + 1
                                                  currentBdFE.weightMeas( l );

                    }

                    // Assembling diagonal entry
                    matrix.addToCoefficient( idDof - 1, idDof - 1, sum );
                }

                // Upper diagonal columns of the elementary boundary mass matrix
                for ( ID k = idofF + 1 ; k <= nDofF ; ++k )
                {

                    // Loop on components invoved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                    {

                        sum = 0;

                        idDof = pId->localToGlobalMap( idofF ) + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;
                        jdDof = pId->localToGlobalMap( k ) + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                        // Loop on quadrature points
                        for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                        {
                            mcoef = 0.0;
                            for ( ID n = 1; n <= nDofF; ++n)
                            {
                                kdDof=pId->localToGlobalMap( n ); // + ( boundaryCond.component( j ) - 1 ) * totalDof;
                                if (boundaryCond.ismixteVec())
                                    mcoef += boundaryCond.MixteVec( kdDof, boundaryCond.component( j ) ) * currentBdFE.phi( int( n - 1 ), l );

                                else mcoef += boundaryCond.mixteCoef() * currentBdFE.phi( int( n - 1 ), l );
                            }

                            currentBdFE.weightMeas( l );
                            sum += mcoef*currentBdFE.phi( int( idofF - 1 ), l ) *
                                   currentBdFE.phi( int( k - 1 ), l ) * currentBdFE.weightMeas( l );

                        }

                        // Assembling upper entry.  The boundary mass matrix is symetric
                        matrix.addToCoefficient( idDof - 1, jdDof - 1, sum );
                        matrix.addToCoefficient( jdDof - 1, idDof - 1, sum );
                    }
                }
            }
        }
    }

    else
    {  //! If BC is given under a functional form

        DataType x, y, z;

        const BCFunctionMixte* pBcF = static_cast<const BCFunctionMixte*>( boundaryCond.pointerToFunctor() );

        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {

            // Pointer to the i-th identifier in the list
            pId = static_cast< const IdentifierNatural* >( boundaryCond( i ) );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            currentBdFE.updateMeas( mesh.bElement( ibF ) );

            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            {
                // Loop on components invoved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {

                    sum = 0;

                    // Global Dof (outside the quad point loop. V. Martin)
                    idDof = pId->localToGlobalMap( idofF ) + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                    {

                        currentBdFE.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                        // Contribution to the diagonal entry of the elementary boundary mass matrix
                        sum += pBcF->coef( time, x, y, z, boundaryCond.component( j ) ) * currentBdFE.phi( int( idofF - 1 ), l ) * currentBdFE.phi( int( idofF - 1 ), l ) *
                               currentBdFE.weightMeas( l );

                        // Adding right hand side contribution
                        rightHandSide[ idDof ] += currentBdFE.phi( int( idofF - 1 ), l ) * boundaryCond( time, x, y, z, boundaryCond.component( j ) ) * // BASEINDEX + 1
                                                  currentBdFE.weightMeas( l );
                    }

                    // Assembling diagonal entry
                    matrix.addToCoefficient( idDof - 1, idDof - 1, sum );
                }

                // Upper diagonal columns of the elementary boundary mass matrix
                for ( ID k = idofF + 1 ; k <= nDofF ; ++k )
                {

                    // Loop on components invoved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                    {

                        sum = 0;

                        // Loop on quadrature points
                        for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                        {

                            currentBdFE.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                            // Upper diagonal entry of the elementary boundary mass matrix
                            sum += pBcF->coef( time,  x, y, z, boundaryCond.component( j )   ) * currentBdFE.phi( int( idofF - 1 ), l ) * currentBdFE.phi( int( k - 1 ), l ) *
                                   currentBdFE.weightMeas( l );
                        }

                        // Globals Dof: row and columns
                        idDof = pId->localToGlobalMap( idofF ) + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;
                        jdDof = pId->localToGlobalMap( k ) + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                        // Assembling upper entry.  The boundary mass matrix is symetric
                        matrix.addToCoefficient( idDof - 1, jdDof - 1, sum );
                        matrix.addToCoefficient( jdDof - 1, idDof - 1, sum );
                    }
                }
            }
        }
    }
}  //bcMixteManage


template <typename MatrixType, typename DataType, typename MeshType>
void
bcMixteManageMatrix( MatrixType& matrix,
                     const MeshType& mesh,
                     const Dof& dof,
                     const BCBase& boundaryCond,
                     CurrentBdFE& currentBdFE,
                     const DataType& time,
                     UInt offset )
{
    if ( matrix.matrixPtr()->Filled() )
        matrix.openCrsMatrix();
    // Number of local Dof in this face
    UInt nDofF = currentBdFE.nbNode();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    DataType sum;

    const IdentifierNatural* pId;
    ID ibF, idDof, jdDof;

    if ( boundaryCond.dataVector() )
    {   //! If BC is given under a vectorial form

        //! for the moment, only one coefficient per BCvector.
        DataType mcoef = boundaryCond.mixteCoef();   //!< the mixte coefficient

        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {

            // Pointer to the i-th identifier in the list
            pId = static_cast< const IdentifierNatural* >( boundaryCond( i ) );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            currentBdFE.updateMeas( mesh.bElement( ibF ) );

            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            {
                // Loop on components invoved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {

                    sum = 0;

                    // Global Dof
                    idDof = boundaryCond( i ) ->id() + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                    {

                        // Contribution to the diagonal entry of the elementary boundary mass matrix
                        sum += mcoef * currentBdFE.phi( int( idofF - 1 ), l ) * currentBdFE.phi( int( idofF - 1 ), l ) *
                               currentBdFE.weightMeas( l );
                    }

                    // Assembling diagonal entry
                    matrix.addToCoefficient( idDof - 1, idDof - 1, sum );
                }

                // Upper diagonal columns of the elementary boundary mass matrix
                for ( ID k = idofF + 1 ; k <= nDofF ; ++k )
                {

                    // Loop on components invoved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                    {

                        sum = 0;

                        // Loop on quadrature points
                        for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                        {

                            // Upper diagonal entry of the elementary boundary mass matrix
                            sum += mcoef * currentBdFE.phi( int( idofF - 1 ), l ) * currentBdFE.phi( int( k - 1 ), l ) *
                                   currentBdFE.weightMeas( l );
                        }

                        // Globals Dof: row and columns
                        idDof = boundaryCond( i ) ->id() + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;
                        jdDof = boundaryCond( k ) ->id() + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                        // Assembling upper entry.  The boundary mass matrix is symetric
                        matrix.addToCoefficient( idDof - 1, jdDof - 1, sum );
                        matrix.addToCoefficient( jdDof - 1, idDof - 1, sum );
                    }
                }
            }
        }
    }

    else
    {  //! If BC is given under a functional form

        DataType x, y, z;

        const BCFunctionMixte* pBcF = static_cast<const BCFunctionMixte*>( boundaryCond.pointerToFunctor() );

        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {

            // Pointer to the i-th identifier in the list
            pId = static_cast< const IdentifierNatural* >( boundaryCond( i ) );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            currentBdFE.updateMeas( mesh.bElement( ibF ) );

            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            {
                // Loop on components invoved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {

                    sum = 0;

                    // Global Dof (outside the quad point loop. V. Martin)
                    idDof = pId->localToGlobalMap( idofF ) + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                    {

                        currentBdFE.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                        // Contribution to the diagonal entry of the elementary boundary mass matrix
                        sum += pBcF->coef( time, x, y, z, boundaryCond.component(j) ) * currentBdFE.phi( int( idofF - 1 ), l ) * currentBdFE.phi( int( idofF - 1 ), l ) *
                               currentBdFE.weightMeas( l );
                    }


                    // Assembling diagonal entry
                    matrix.addToCoefficient( idDof - 1, idDof - 1, sum );
                }

                // Upper diagonal columns of the elementary boundary mass matrix
                for ( ID k = idofF + 1 ; k <= nDofF ; ++k )
                {

                    // Loop on components invoved in this boundary condition
                    for ( ID j = 1; j <= nComp; ++j )
                    {

                        sum = 0;

                        // Loop on quadrature points
                        for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                        {

                            currentBdFE.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                            // Upper diagonal entry of the elementary boundary mass matrix
                            sum += pBcF->coef( time, x, y, z, boundaryCond.component( j )  ) * currentBdFE.phi( int( idofF - 1 ), l ) * currentBdFE.phi( int( k - 1 ), l ) *
                                   currentBdFE.weightMeas( l );
                        }

                        // Globals Dof: row and columns
                        idDof = pId->localToGlobalMap( idofF ) + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;
                        jdDof = pId->localToGlobalMap( k ) + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                        // Assembling upper entry.  The boundary mas matrix is symetric
                        matrix.addToCoefficient( idDof - 1, jdDof - 1, sum );
                        matrix.addToCoefficient( jdDof - 1, idDof - 1, sum );
                    }
                }
            }
        }
    }
}   //bcMixteManageMatrix



template <typename VectorType, typename DataType, typename MeshType>
void
bcMixteManageVector( VectorType& rightHandSide,
                     const MeshType& mesh,
                     const Dof& dof,
                     const BCBase& boundaryCond,
                     CurrentBdFE& currentBdFE,
                     const DataType& time,
                     UInt offset )
{

    // Number of local Dof in this face
    UInt nDofF = currentBdFE.nbNode();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    const IdentifierNatural* pId;
    ID ibF, idDof;

    if ( boundaryCond.dataVector() )
    {   //! If BC is given under a vectorial form

        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {

            // Pointer to the i-th identifier in the list
            pId = static_cast< const IdentifierNatural* >( boundaryCond( i ) );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            currentBdFE.updateMeas( mesh.bElement( ibF ) );

            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            {
                // Loop on components invoved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {
                    // Global Dof
                    idDof = boundaryCond( i ) ->id() + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                    {
                        // Adding right hand side contribution
                        rightHandSide[ idDof ] += currentBdFE.phi( int( idofF - 1 ), l ) * boundaryCond( boundaryCond( i ) ->id(), boundaryCond.component( j ) ) * // BASEINDEX + 1
                                                  currentBdFE.weightMeas( l );
                    }
                }
            }
        }
    }

    else
    {  //! If BC is given under a functional form

        DataType x, y, z;

        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {

            // Pointer to the i-th identifier in the list
            pId = static_cast< const IdentifierNatural* >( boundaryCond( i ) );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            currentBdFE.updateMeas( mesh.bElement( ibF ) );

            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            {
                // Loop on components invoved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {

                    // Global Dof (outside the quad point loop. V. Martin)
                    idDof = pId->localToGlobalMap( idofF ) + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                    {

                        currentBdFE.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                        // Adding right hand side contribution
                        rightHandSide[ idDof ] += currentBdFE.phi( int( idofF - 1 ), l ) * boundaryCond( time, x, y, z, boundaryCond.component( j ) ) *// BASEINDEX + 1
                                                  currentBdFE.weightMeas( l );
                    }
                }
            }
        }
    }
} //bcMixteManageVector



// ===================================================
// Flux BC
// ===================================================


template <typename MatrixType,
typename VectorType,
typename MeshType,
typename DataType>
void
bcFluxManage( MatrixType&     matrix,
              VectorType&    rightHandSide,
              const MeshType& mesh,
              const Dof&      dof,
              const BCBase&   boundaryCond,
              CurrentBdFE&    currentBdFE,
              const DataType& time,
              UInt            offset)

{
    bcFluxManageVector(rightHandSide,boundaryCond,time,offset);
    bcFluxManageMatrix(matrix, mesh, dof, boundaryCond, currentBdFE, time, offset);
} //bcFluxManage


template <typename VectorType,
typename DataType>
void
bcFluxManageVector(
    VectorType&    rightHandSide,
    const BCBase&   boundaryCond,
    const DataType& time,
    UInt            offset)

{
    rightHandSide.setCoefficient(offset + 1,boundaryCond(time, 0., 0., 0., 1));
}  //bcFluxManageVector


template <typename MatrixType,
typename MeshType,
typename DataType>
void
bcFluxManageMatrix( MatrixType&     matrix,
                    const MeshType& mesh,
                    const Dof&      dof,
                    const BCBase&   boundaryCond,
                    CurrentBdFE&    currentBdFE,
                    const DataType& /*time*/,
                    UInt            offset )
{
    if ( matrix.matrixPtr()->Filled() )
        matrix.openCrsMatrix();

    // Number of local Dof in this face
    UInt nDofF = currentBdFE.nbNode();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    DataType sum;

    const IdentifierNatural* pId;
    ID ibF, idDof, jdDof/*, kdDof*/;

    if ( !boundaryCond.dataVector() )
    {
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {
            pId = static_cast< const IdentifierNatural* >( boundaryCond( i ) );

            // Number of the current boundary face
            ibF = pId->id();
            // Updating face stuff
            currentBdFE.updateMeasNormalQuadPt( mesh.bElement( ibF ) );

            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            {
                for ( int ic = 1; ic <= (int)nComp; ++ic)
                {
                    idDof = pId->localToGlobalMap( idofF ) + (ic - 1)*totalDof;

                    sum = 0.;
                    for ( int iq = 0; iq < (int)currentBdFE.nbQuadPt(); ++iq )
                    {
                        sum += currentBdFE.phi( int( idofF - 1 ), iq )*
                               currentBdFE.normal(int(ic - 1), iq)*
                               currentBdFE.weightMeas(iq);
                    }

                    jdDof = offset;

                    matrix.addToCoefficient( idDof - 1, jdDof    , sum );
                    matrix.addToCoefficient( jdDof    , idDof - 1, sum );
                }
            }
        }
    }
} // bcFluxManageMatrix



template <typename MatrixType, typename VectorType, typename DataType, typename MeshType>
void
bcResistanceManage( MatrixType& matrix,
                    VectorType& rightHandSide,
                    const MeshType& mesh,
                    const Dof& dof,
                    const BCBase& boundaryCond,
                    CurrentBdFE& currentBdFE,
                    const DataType& /*time*/,
                    UInt offset )
{
    // Number of local Dof in this face
    UInt nDofF = currentBdFE.nbNode();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    std::set<ID> resistanceDofs;

    const IdentifierNatural* pId;
    ID ibF, idDof, jdDof, kdDof;

    if ( boundaryCond.dataVector() )
    {
        //auxiliary vector
        VectorType vv(rightHandSide.map(), Repeated);

        DataType  mbcb;

        // Loop on BC identifiers
        for ( ID i = 1; i <= boundaryCond.list_size(); ++i )
        {
            // Pointer to the i-th identifier in the list
            pId = static_cast< const IdentifierNatural* >( boundaryCond( i ) );

            // Number of the current boundary face
            ibF = pId->id();

            currentBdFE.updateMeasNormalQuadPt( mesh.boundaryFace( ibF ) );

            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            {
                resistanceDofs.insert( pId->localToGlobalMap( idofF ) );

                // Loop on components involved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {
                    idDof = pId->localToGlobalMap( idofF ) + ( boundaryCond.component( j ) - 1 ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                    {
                        vv[idDof] += currentBdFE.phi( int( idofF-1 ), l ) *  currentBdFE.normal( int( j-1 ), l ) * currentBdFE.weightMeas( l );

                        mbcb=0;

                        // data on quadrature point
                        for ( ID n = 1; n <= nDofF; ++n)
                        {
                            kdDof=pId->localToGlobalMap( n );
                            mbcb += boundaryCond( kdDof, boundaryCond.component( j ) )* currentBdFE.phi( int( n - 1 ), l ) ;
                        }

                        rightHandSide[ idDof ] +=  mbcb* currentBdFE.phi( int( idofF - 1 ), l ) *  currentBdFE.normal( int( j-1 ), l ) * // BASEINDEX + 1
                                                   currentBdFE.weightMeas( l );

                    }
                }
            }
        }


        for ( std::set<ID>::iterator iDofIt = resistanceDofs.begin();
                iDofIt != resistanceDofs.end(); ++iDofIt )
        {
            for ( UInt iComp = 1; iComp <= nComp; ++iComp )
            {
                idDof = *iDofIt + ( boundaryCond.component( iComp ) - 1 ) * totalDof + offset;
                for ( std::set<ID>::iterator jDofIt = resistanceDofs.begin();
                        jDofIt != resistanceDofs.end(); ++jDofIt )
                {
                    for ( UInt jComp = 1; jComp <= nComp; ++jComp )
                    {
                        jdDof = *jDofIt + ( boundaryCond.component( jComp ) - 1 ) * totalDof + offset;

                        matrix.addToCoefficient( idDof-1,  jdDof-1, boundaryCond.resistanceCoef() * vv[idDof] * vv[jdDof] );
                    }
                }
            }
        }

    }
    else
        ERROR_MSG( "This BC type is not yet implemented" );
} //bcResistanceManage

} // end of namespace LifeV
#endif
