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

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManageNormal.hpp>


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
          DOF const& dof,
          BCHandler const& bcHandler,
          CurrentBoundaryFE& currentBdFE,
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
          const DOF& dof,
          const BCHandler& bcHandler,
          CurrentBoundaryFE& currentBdFE,
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
                const DOF&       dof,
                const BCHandler& bcHandler,
                CurrentBoundaryFE&     currentBdFE,
                const DataType&  diagonalizeCoef,
                const DataType&  time = 0 );



//! Prescribe boundary conditions. Case in which only the right hand side is modified
/*! This method is deprecated since the order of diagonalizeCoef and time are switched wrt to bcManage.
 *  Use instead bcManageRhs ad be careful to use the correct order.
 *
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
LIFEV_DEPRECATED ( void )
bcManageVector( VectorType&      rightHandSide,
                const MeshType&  mesh,
                const DOF&       dof,
                const BCHandler& bcHandler,
                CurrentBoundaryFE&     currentBdFE,
                const DataType&  time,
                const DataType&  diagonalizeCoef );


//! Prescribe boundary conditions. Case in which only the right hand side is modified
/*!
 * The right hand side is modified to take into account the boundary conditions
 * @param rightHandSide   The system right hand side
 * @param mesh  The mesh
 * @param dof  Container of the local to global map of DOFs
 * @param bcHandler The boundary conditions handler
 * @param currentBdFE Current finite element on boundary
 * @parma diagonalizeCoef The coefficient used during the system diagonalization
 * @param time The time
 */
template <typename VectorType, typename MeshType, typename DataType>
void
bcManageRhs( VectorType&      rightHandSide,
             const MeshType&  mesh,
             const DOF&       dof,
             const BCHandler& bcHandler,
             CurrentBoundaryFE&     currentBdFE,
             const DataType&  diagonalizeCoef,
             const DataType&  time );



//! Prescribe boundary conditions. Case in which only the right hand side is modified
/*! This method is deprecated since the order of diagonalizeCoef and time are switched wrt to bcManage.
 *  Use instead bcManageRhs ad be careful to use the correct order.
 *
 * The Right hand side is modified to take into account the boundary conditions
 * @param rightHandSide   The system right hand side
 * @param feSpace  The finite element space
 * @param bcHandler The boundary conditions handler
 * @param time The time
 * @parma diagonalizeCoef The coefficient used during the system diagonalization
 */
template <typename VectorType, typename DataType, typename Mesh, typename MapEpetra>
LIFEV_DEPRECATED ( void )
bcManageVector( VectorType&                     rightHandSide,
                FESpace<Mesh, MapEpetra>&       feSpace,
                const BCHandler&                bcHandler,
                const DataType&                 time,
                const DataType&                 diagonalizeCoef );

//! Prescribe boundary conditions. Case in which only the residual is available
/*
The residual and the right hand side are modified to take into account the boundary conditions
@param res residual vector
@param rhs right hand side
@param sol the solution vector used to compute the residual. It must be a Unique VectorEpetra
@param feSpace the finite element space
@param bcHandler the boundary condition handler
@param time the current time level
@param diagonalizeCoef the coefficient put in the diagonal entry (of a matrix) when applying Dirichlet boundary conditions
 */
template <typename VectorType, typename DataType, typename Mesh, typename MapEpetra>
void
bcManageResidual( VectorType&                     res,
                  VectorType&                     rhs,
                  const VectorType&               sol,
                  FESpace<Mesh, MapEpetra>&       feSpace,
                  const BCHandler&                bcHandler,
                  const DataType&                 time,
                  const DataType&                 diagonalizeCoef );


//! Prescribe boundary conditions. Case in which only the right hand side is modified
/*!
 * The Right hand side is modified to take into account the boundary conditions
 * @param rightHandSide   The system right hand side
 * @param feSpace  The finite element space
 * @param bcHandler The boundary conditions handler
 * @param time The time
 * @parma diagonalizeCoef The coefficient used during the system diagonalization
 */
template <typename VectorType, typename DataType, typename Mesh, typename MapEpetra>
void
bcManageRhs( VectorType&                     rightHandSide,
             FESpace<Mesh, MapEpetra>&       feSpace,
             const BCHandler&                bcHandler,
             const DataType&                 diagonalizeCoef,
             const DataType&                 time );

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
                   const DOF& dof,
                   const BCBase& boundaryCond,
                   const CurrentBoundaryFE& /*currentBdFE*/,
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
                       const DOF& dof,
                       const BCBase& boundaryCond,
                       const CurrentBoundaryFE& /*currentBdFE*/,
                       const DataType& diagonalizeCoef,
                       const DataType& time,
                       const VectorType& feVec ,
                       UInt offset = 0 );



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
                         const DOF& dof,
                         const BCBase& boundaryCond,
                         const DataType& diagonalizeCoef,
                         UInt offset );



//! Prescribe Essential boundary conditions on the right hand side
/*! This method is deprecated since the order of diagonalizeCoef and time are switched wrt to bcManage.
 *  Use instead bcManageRhs ad be careful to use the correct order.
 *
 * The right hand side is modified to take into account the Essential boundary conditions
 * @param rightHandSide   The system rightHandSide
 * @param dof  Container of the local to global map of DOFs
 * @param boundaryCond The boundary condition (@c BCBase)
 * @parma diagonalizeCoef The coefficient used during the system diagonalization
 * @param offset The boundary condition offset
 */
template <typename VectorType, typename DataType>
LIFEV_DEPRECATED ( void )
bcEssentialManageVector( VectorType&     rightHandSide,
                         const DOF&      dof,
                         const BCBase&   boundaryCond,
                         const DataType& time,
                         const DataType& diagonalizeCoef,
                         UInt            offset );

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
bcEssentialManageRhs( VectorType&     rightHandSide,
                      const DOF&      dof,
                      const BCBase&   boundaryCond,
                      const DataType& diagonalizeCoef,
                      const DataType& time,
                      UInt            offset );

//! Prescribe all the Essential boundary conditions on the right hand side and forgetting about the other BCs.
/*!
 * The right hand side is modified to take into account the Essential boundary conditions
 * This is useful when imposing homogeneous BC, in conjuction with coeff = 0.
 * @param rightHandSide   The system rightHandSide
 * @param dof  Container of the local to global map of DOFs
 * @param bcHandler The boundary conditions handler
 * @parma diagonalizeCoef The coefficient used during the system diagonalization
 * @param offset The boundary condition offset
 * Remark: another possible name would be bcManageHomogeneousRhs and set diagonalizeCoef = 0.
 */
template <typename VectorType, typename DataType>
void
bcEssentialManageRhs( VectorType&     rightHandSide,
                      const DOF&      dof,
                      const BCHandler& bcHandler,
                      const DataType& diagonalizeCoef,
                      const DataType& time);



//! Prescribe essential boundary conditions. Case in which only the residual is available
/*
The residual and the right hand side are modified to take into account the boundary conditions
@param res residual vector
@param rhs right hand side
@param sol the solution vector used to compute the residual. It must be a Unique VectorEpetra
@param dof the DOF instance
@param boundaryCond the specific boundary condition (BCBase)
@param time the current time level
@param diagonalizeCoef the coefficient put in the diagonal entry (of a matrix) when applying Dirichlet boundary conditions
@param offset the UInt offset for the boundary condition
 */
template <typename VectorType, typename DataType>
void
bcEssentialManageResidual(VectorType&     res,
                          VectorType&     rhs,
                          const VectorType&     sol,
                          const DOF&      dof,
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
                   const DOF& dof,
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
                 const DOF& dof, const
                 BCBase& boundaryCond,
                 CurrentBoundaryFE& currentBdFE,
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
                     const DOF& dof,
                     const BCBase& boundaryCond,
                     CurrentBoundaryFE& currentBdFE,
                     const DataType& time,
                     const VectorType& feVec,
                     UInt offset );

// @}


// ===================================================
//! @name Robin BC
// @{
// ===================================================


//! Prescribe Robin boundary condition
/*!
 * The matrix and the right hand side are modified to take into account the Robin boundary condition
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
bcRobinManage( MatrixType& matrix,
               VectorType& rightHandSide,
               const MeshType& mesh,
               const DOF& dof,
               const BCBase& boundaryCond,
               CurrentBoundaryFE& currentBdFE,
               const DataType& time,
               UInt offset );


//! Prescribe Robin boundary condition only on the matrix
/*!
 * The matrix is modified to take into account the Robin boundary condition
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
bcRobinManageMatrix( MatrixType& matrix,
                     const MeshType& mesh,
                     const DOF& dof,
                     const BCBase& boundaryCond,
                     CurrentBoundaryFE& currentBdFE,
                     const DataType& time,
                     UInt offset );




//! Prescribe Robin boundary conditions. Case in which only the residual is available
/*

The residual and the right hand side are modified to take into account the boundary conditions
@param res residual vector
@param rhs right hand side
@param sol the solution vector used to compute the residual. It must be a Unique VectorEpetra
@param dof the DOF instance
@param currentBdFE the current boundary finite element
@param boundaryCond the specific boundary condition (BCBase)
@param time the current time level
@param diagonalizeCoef the coefficient put in the diagonal entry (of a matrix) when applying Dirichlet boundary conditions
@param offset the UInt offset for the boundary condition
 */
template <typename VectorType, typename DataType, typename MeshType>
void
bcRobinManageResidual( VectorType& residual,
                       VectorType& rightHandSide,
                       const VectorType& solution,
                       const MeshType& mesh,
                       const DOF& dof,
                       const BCBase& boundaryCond,
                       CurrentBoundaryFE& currentBdFE,
                       const DataType& time,
                       UInt offset );


//! Prescribe Robin boundary condition only on the rightHandSide
/*!
 * The matrix is modified to take into account the Robin boundary condition
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
bcRobinManageVector( VectorType& rightHandSide,
                     const MeshType& mesh,
                     const DOF& dof,
                     const BCBase& boundaryCond,
                     CurrentBoundaryFE& currentBdFE,
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
              const DOF&      dof,
              const BCBase&   boundaryCond,
              CurrentBoundaryFE&    currentBdFE,
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
                    const DOF&      dof,
                    const BCBase&   boundaryCond,
                    CurrentBoundaryFE&    currentBdFE,
                    const DataType& /*time*/,
                    UInt            offset );


//! Prescribe Flux boundary conditions. Case in which only the residual is available
/*
The residual and the right hand side are modified to take into account the boundary conditions
@param res residual vector
@param rhs right hand side
@param sol the solution vector used to compute the residual. It must be a Unique VectorEpetra
@param dof the DOF instance
@param currentBdFE the current boundary finite element
@param boundaryCond the specific boundary condition (BCBase)
@param time the current time level
@param diagonalizeCoef the coefficient put in the diagonal entry (of a matrix) when applying Dirichlet boundary conditions
@param offset the UInt offset for the boundary condition
 */
template <typename VectorType,
typename MeshType,
typename DataType>
void
bcFluxManageResidual( VectorType&     residual,
                      VectorType&     rightHandSide,
                      const VectorType&     solution,
                      const MeshType& mesh,
                      const DOF&      dof,
                      const BCBase&   boundaryCond,
                      CurrentBoundaryFE&    currentBdFE,
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
                    const DOF& dof,
                    const BCBase& boundaryCond,
                    CurrentBoundaryFE& currentBdFE,
                    const DataType& /*time*/,
                    UInt offset );
template <typename VectorType, typename DataType, typename MeshType>
void
bcResistanceManageVector( VectorType& rightHandSide,
                    const MeshType& mesh,
                    const DOF& dof,
                    const BCBase& boundaryCond,
                    CurrentBoundaryFE& currentBdFE,
                    const DataType& /*time*/,
                    UInt offset );
template <typename MatrixType, typename DataType, typename MeshType>
void
bcResistanceManageMatrix( MatrixType& matrix,
                    const MeshType& mesh,
                    const DOF& dof,
                    const BCBase& boundaryCond,
                    CurrentBoundaryFE& currentBdFE,
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
          DOF const& dof,
          BCHandler const& bcHandler,
          CurrentBoundaryFE& currentBdFE,
          DataType const& diagonalizeCoef,
          DataType const& time )
{

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
        case Robin:      // Robin boundary conditions (Robin)
            bcRobinManage( matrix, rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            globalassemble=true;
            break;
        case Flux:       // Flux boundary condition
            bcFluxManage( matrix, rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset()+bcHandler[i].offset());
            globalassemble=true;
            break;
        case Resistance: // Resistance boundary condition
            bcResistanceManage( matrix, rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            globalassemble=true;
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }

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
        case Robin:  	    // Robin boundary conditions (Robin)
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
          const DOF& dof,
          const BCHandler& bcHandler,
          CurrentBoundaryFE& currentBdFE,
          const DataType diagonalizeCoef,
          const DataType& time,
          VectorType& feVec )
{

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
        case Robin:  // Robin boundary conditions (Robin)

            if (bcHandler[ i ].isUDep())
            {
                ERROR_MSG( "This BC mode is not yet implemented for this setting" );    //not implemented yet
            }
            else
            {
                bcRobinManage( matrix, rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            }
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }


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
        case Robin:  // Robin boundary conditions (Robin)
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
                const DOF&       dof,
                const BCHandler& bcHandler,
                CurrentBoundaryFE&     currentBdFE,
                const DataType&  diagonalizeCoef,
                const DataType&  time )
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
        case Robin:  // Robin boundary conditions (Robin)
            bcRobinManageMatrix( matrix, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            globalassemble=true;
            break;
        case Flux:  // Flux boundary conditions
            bcFluxManageMatrix( matrix, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset()+bcHandler[i].offset());
            globalassemble=true;
            break;
        case Resistance:
            bcResistanceManageMatrix( matrix, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            globalassemble=true;
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
        case Robin:  // Robin boundary conditions (Robin)
            break;
        case Flux:  // Robin boundary conditions (Robin)
            break;
        case Resistance:
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
                const DOF&       dof,
                const BCHandler& bcHandler,
                CurrentBoundaryFE&     currentBdFE,
                const DataType&  time,
                const DataType&  diagonalizeCoef )
{
    bcManageRhs( rightHandSide,
                    mesh,
                    dof,
                    bcHandler,
                    currentBdFE,
                    diagonalizeCoef,
                    time );
}

template <typename VectorType, typename MeshType, typename DataType>
void
bcManageRhs( VectorType&      rightHandSide,
                const MeshType&  mesh,
                const DOF&       dof,
                const BCHandler& bcHandler,
                CurrentBoundaryFE&     currentBdFE,
                const DataType&  diagonalizeCoef,
                const DataType&  time )
{

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
            bcEssentialManageRhs( rightHandSide, dof, bcHandler[ i ], diagonalizeCoef, time, bcHandler.offset() );
            break;
        case Natural:  // Natural boundary conditions (Neumann)
            bcNaturalManage( rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            break;
        case Robin:  // Robin boundary conditions (Robin)
            bcRobinManageVector( rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            break;
        case Flux:  // Flux boundary conditions
            bcFluxManageVector( rightHandSide, bcHandler[ i ], time, bcHandler.offset()+bcHandler[i].offset() );
            break;
        case Resistance:
            bcResistanceManageVector( rightHandSide, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }

}


template <typename VectorType, typename MeshType, typename DataType>
void
bcManageResidual( VectorType&                     res,
                  VectorType&                     rhs,
                  const VectorType&                     sol,
                  const MeshType&  mesh,
                  const DOF&       dof,
                  const BCHandler& bcHandler,
                  CurrentBoundaryFE&     currentBdFE,
                  const DataType&  time,
                  const DataType&  diagonalizeCoef )
{
    VectorType rhsRepeated(rhs.map(),Repeated);

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
            bcEssentialManageResidual( res, rhs, sol, dof, bcHandler[ i ], time, diagonalizeCoef, bcHandler.offset() );
            break;
        case Natural:  // Natural boundary conditions (Neumann)
            bcNaturalManage( rhs, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset() );
            break;
        case Robin:  // Robin boundary conditions (Robin) to be implemented
            bcRobinManageResidual( res,   rhs, sol, mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset());
            break;
        case Flux:  // Flux boundary conditions to be implemented
            bcFluxManageResidual( res, rhs, sol,  mesh, dof, bcHandler[ i ], currentBdFE, time, bcHandler.offset()+bcHandler[i].offset() );
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }

    //    rhsRepeated.globalAssemble();

    //rhs += rhsRepeated;
}


template <typename VectorType, typename DataType, typename Mesh, typename MapEpetra>
void
bcManageVector( VectorType&                     rightHandSide,
                FESpace<Mesh, MapEpetra>&       feSpace,
                const BCHandler&                bcHandler,
                const DataType&                 time,
                const DataType&                 diagonalizeCoef )
{
    bcManageRhs( rightHandSide,
                 feSpace,
                 bcHandler,
                 diagonalizeCoef,
                 time );
}

template <typename VectorType, typename DataType, typename Mesh, typename MapEpetra>
void
bcManageRhs( VectorType&                     rightHandSide,
                FESpace<Mesh, MapEpetra>&       feSpace,
                const BCHandler&                bcHandler,
                const DataType&                 diagonalizeCoef,
                const DataType&                 time )
{
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
            bcEssentialManageRhs( rightHandSide, feSpace.dof(), bcHandler[ i ], diagonalizeCoef, time, bcHandler.offset() );
            break;
        case Natural:  // Natural boundary conditions (Neumann)
            bcNaturalManage( rightHandSide, *feSpace.mesh(), feSpace.dof(), bcHandler[ i ], feSpace.feBd(), time, bcHandler.offset() );
            break;
        case Robin:  // Robin boundary conditions (Robin)
            bcRobinManageVector( rightHandSide, *feSpace.mesh(), feSpace.dof(), bcHandler[ i ], feSpace.feBd(), time, bcHandler.offset() );
            break;
        case Flux:  // Flux boundary conditions
            bcFluxManageVector( rightHandSide, bcHandler[ i ], time, bcHandler.offset()+bcHandler[i].offset() );
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
void
bcEssentialManage( MatrixType& matrix,
                   VectorType& rightHandSide,
                   const MeshType& /*mesh*/,
                   const DOF& dof,
                   const BCBase& boundaryCond,
                   const CurrentBoundaryFE& /*currentBdFE*/,
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

    if ( boundaryCond.isDataAVector() )
    { //! If BC is given under a vectorial form

        assert( static_cast< const BCVectorInterface* > (boundaryCond.pointerToBCVector()) != 0 );

        // Loop on BC identifiers
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {

            // Loop on components involved in this boundary condition
            for ( ID j = 0; j < nComp; ++j )
            {
                // Global Dof
                idDof = boundaryCond[ i ]->id() + boundaryCond.component( j ) * totalDof + offset;
                datumVec.push_back(boundaryCond( boundaryCond[ i ] ->id(), boundaryCond.component( j ) ));
                idDofVec.push_back(idDof);
            }
        }
    }
    else
    { //! If BC is given under a functional form

        DataType x, y, z;
        // Loop on BC identifiers
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {
            // Coordinates of the node where we impose the value
            x = static_cast< const BCIdentifierEssential* >( boundaryCond[ i ] ) ->x();
            y = static_cast< const BCIdentifierEssential* >( boundaryCond[ i ] ) ->y();
            z = static_cast< const BCIdentifierEssential* >( boundaryCond[ i ] ) ->z();

            // Loop on components involved in this boundary condition
            for ( ID j = 0; j < nComp; ++j )
            {
                // Global Dof
                idDof = boundaryCond[ i ] ->id() + boundaryCond.component( j ) * totalDof + offset;

                datumVec.push_back(boundaryCond( time, x, y, z, boundaryCond.component( j ) ));
                idDofVec.push_back(idDof);

            }
        }
    }

    // If there is an offset than there is a Lagrange multiplier (flux BC)
    if (boundaryCond.offset() > 0)
    {
        // bcType has been changed Flux -> Essential, need to diagonalize also the Lagrange multiplier
       idDofVec.push_back(offset + boundaryCond.offset());
       datumVec.push_back( 0. );
    }

    // Modifying matrix and right hand side
    matrix.diagonalize( idDofVec, diagonalizeCoef, rightHandSide, datumVec);

}


template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void
bcEssentialManageUDep( MatrixType& matrix,
                       VectorType& rightHandSide,
                       const MeshType& /*mesh*/,
                       const DOF& dof,
                       const BCBase& boundaryCond,
                       const CurrentBoundaryFE& /*currentBdFE*/,
                       const DataType& diagonalizeCoef,
                       const DataType& time,
                       const VectorType& feVec ,
                       UInt offset )
{

    ID idDof;

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    if ( boundaryCond.isDataAVector() )
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
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {
            // Coordinates of the node where we impose the value
            x = static_cast< const BCIdentifierEssential* >( boundaryCond[ i ] ) ->x();
            y = static_cast< const BCIdentifierEssential* >( boundaryCond[ i ] ) ->y();
            z = static_cast< const BCIdentifierEssential* >( boundaryCond[ i ] ) ->z();

            // Loop on components involved in this boundary condition
            for ( ID j = 0; j < nComp; ++j )
            {
                // Global Dof
                idDof = boundaryCond[ i ] ->id() + boundaryCond.component( j ) * totalDof + offset;

                Real datum = boundaryCond( time, x, y, z, boundaryCond.component( j ) ,feVec[idDof]);

                datumVec.push_back(datum);
                idDofVec.push_back(idDof);

            }
        }

    // If there is an offset than there is a Lagrange multiplier (flux BC)
        if (boundaryCond.offset() > 0)
        {
            // bcType has been changed Flux -> Essential, need to diagonalize also the Lagrange multiplier
            idDofVec.push_back(offset + boundaryCond.offset());
        }

        // Modifying matrix and right hand side
        matrix.diagonalize( idDofVec, diagonalizeCoef, rightHandSide, datumVec);
    }
}


template <typename MatrixType, typename DataType>
void
bcEssentialManageMatrix( MatrixType& matrix,
                         const DOF& dof,
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
    for ( ID i = 0; i < boundaryCond.list_size(); ++i )
    {
        // Loop on components involved in this boundary condition
        for ( ID j = 0; j < nComp; ++j )
        {
            // Global Dof
            idDof = boundaryCond[ i ] ->id() + boundaryCond.component( j ) * totalDof;
            idDofVec.push_back(idDof);
        }
    }

    // If there is an offset than there is a Lagrange multiplier (flux BC)
    if (boundaryCond.offset() > 0)
      {
        // bcType has been changed Flux -> Essential, need to diagonalize also the Lagrange multiplier
        idDofVec.push_back(offset + boundaryCond.offset());
      }

    // Modifying ONLY matrix
    matrix.diagonalize( idDofVec, diagonalizeCoef, offset);
}


template <typename VectorType, typename DataType>
void
bcEssentialManageVector( VectorType&     rightHandSide,
                         const DOF&      dof,
                         const BCBase&   boundaryCond,
                         const DataType& time,
                         const DataType& diagonalizeCoef,
                         UInt            offset )
{
    bcEssentialManageRhs( rightHandSide,
                          dof,
                          boundaryCond,
                          diagonalizeCoef,
                          time,
                          offset );

}

template <typename VectorType, typename DataType>
void
bcEssentialManageRhs( VectorType&     rightHandSide,
                      const DOF&      dof,
                      const BCHandler& bcHandler,
                      const DataType& diagonalizeCoef,
                      const DataType& time)
{
    // Loop on boundary conditions
    for ( ID i = 0; i < bcHandler.size(); ++i )
    {

        switch ( bcHandler[ i ].type() )
        {
        case Essential:  // Essential boundary conditions (Dirichlet)
        case EssentialEdges:
        case EssentialVertices:
            if ( (bcHandler[ i ].mode() == Tangential) ||
                 (bcHandler[ i ].mode() == Normal) || (bcHandler[ i ].mode() == Directional) )
            {
                ERROR_MSG( "This BC mode is not yet implemented for this setting" );
            }
            bcEssentialManageRhs( rightHandSide, dof, bcHandler[ i ], diagonalizeCoef, time, bcHandler.offset() );
            break;
            // Not considering the other cases.
        case Natural:  // Natural boundary conditions (Neumann)
        case Robin:  // Robin boundary conditions (Robin)
        case Flux:  // Flux boundary conditions
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }

}


template <typename VectorType, typename DataType>
void
bcEssentialManageRhs( VectorType&     rightHandSide,
                      const DOF&      dof,
                      const BCBase&   boundaryCond,
                      const DataType& diagonalizeCoef,
                      const DataType& time,
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

    if ( boundaryCond.isDataAVector() )
    {  //! If BC is given under a vectorial form

        // Loop on BC identifiers
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {
            // Loop on components involved in this boundary condition
            for ( ID j = 0; j < nComp; ++j )
            {

                // Global Dof
                idDof = boundaryCond[ i ] ->id() + boundaryCond.component( j ) * totalDof + offset;

                idDofVec.push_back( idDof );
                datumVec.push_back( diagonalizeCoef * boundaryCond( boundaryCond[ i ] ->id(), boundaryCond.component( j ) ) );
            }
        }
    }
    else
    {  //! If BC is given under a functional form

        DataType x, y, z;
        // Loop on BC identifiers
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {
            // Coordinates of the node where we impose the value
            x = static_cast< const BCIdentifierEssential* >( boundaryCond[ i ] ) ->x();
            y = static_cast< const BCIdentifierEssential* >( boundaryCond[ i ] ) ->y();
            z = static_cast< const BCIdentifierEssential* >( boundaryCond[ i ] ) ->z();

            // Loop on components involved in this boundary condition
            for ( ID j = 0; j < nComp; ++j )
            {
                // Global Dof

                idDof = boundaryCond[ i ] ->id() + boundaryCond.component( j ) * totalDof + offset;
                // Modifying right hand side
                idDofVec.push_back(idDof);
                datumVec.push_back( diagonalizeCoef * boundaryCond( time, x, y, z, boundaryCond.component( j ) ) );
            }
        }
    }

    rightHandSide.setCoefficients( idDofVec, datumVec);
}


template <typename VectorType, typename DataType>
void
bcEssentialManageResidual(VectorType&     res,
                          VectorType&     rhs,
                          const VectorType&     sol,
                          const DOF&      dof,
                          const BCBase&   boundaryCond,
                          const DataType& time,
                          const DataType& diagonalizeCoef,
                          UInt            offset )
{

    if(sol.mapType()==Unique)
    {
        std::cout<<"pass me a repeated solution"<<std::endl;
        VectorType repeatedSolution(sol, Repeated);
        bcEssentialManageResidual(  res,
				    rhs,
				    repeatedSolution,
				    dof,
				    boundaryCond,
				    time,
				    diagonalizeCoef,
				    offset );
        return;
    }

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
    std::vector<Real> rhsVec(0);
    rhsVec.reserve(boundaryCond.list_size()*nComp);

    if ( boundaryCond.isDataAVector() )
    {  //! If BC is given under a vectorial form
        // Loop on BC identifiers
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {
            // Loop on components involved in this boundary condition
            for ( ID j = 0; j < nComp; ++j )
            {
                // Global Dof
                idDof = boundaryCond[ i ] ->id() + boundaryCond.component( j ) * totalDof + offset;
                idDofVec.push_back( idDof );
                datumVec.push_back( diagonalizeCoef*sol(idDof) );
                rhsVec.push_back(diagonalizeCoef*boundaryCond( boundaryCond[ i ] ->id(), boundaryCond.component( j ) ));
            }
        }
    }
    else
    {  //! If BC is given under a functional form
        DataType x, y, z;
        // Loop on BC identifiers
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {
            // Coordinates of the node where we impose the value
            x = static_cast< const BCIdentifierEssential* >( boundaryCond[ i ] ) ->x();
            y = static_cast< const BCIdentifierEssential* >( boundaryCond[ i ] ) ->y();
            z = static_cast< const BCIdentifierEssential* >( boundaryCond[ i ] ) ->z();

            // Loop on components involved in this boundary condition
            for ( ID j = 0; j < nComp; ++j )
            {
                // Global Dof

                idDof = boundaryCond[ i ] ->id() + boundaryCond.component( j ) * totalDof + offset;
                // Modifying right hand side
                idDofVec.push_back(idDof);
                datumVec.push_back( diagonalizeCoef*sol(idDof) );
                rhsVec.push_back(diagonalizeCoef*boundaryCond( time, x, y, z, boundaryCond.component( j ) ));
            }
        }
    }

    res.setCoefficients( idDofVec, datumVec);
    rhs.setCoefficients( idDofVec, rhsVec);
}

// ===================================================
// Natural BC
// ===================================================


template <typename VectorType, typename MeshType, typename DataType>
void
bcNaturalManage( VectorType& rightHandSide,
                 const MeshType& mesh,
                 const DOF& dof, const
                 BCBase& boundaryCond,
                 CurrentBoundaryFE& currentBdFE,
                 const DataType& time,
                 UInt offset )
{

    // Number of local DOF (i.e. nodes) in this face
    UInt nDofF = currentBdFE.nbNode();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    const BCIdentifierNatural* pId;
    ID ibF, idDof, icDof, gDof;
    Real sum;

    if ( boundaryCond.isDataAVector() )
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
            for ( ID i = 0; i < boundaryCond.list_size(); ++i )
            {
                // Loop on components involved in this boundary condition
                for ( ID j = 0; j < nComp; ++j )
                {
                    ID id = boundaryCond[i]->id();

                    // Global Dof
                    idDof = id + boundaryCond.component( j ) * totalDof + offset;

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
            for ( ID i = 0; i < boundaryCond.list_size(); ++i )
            {
                // Pointer to the i-th itdentifier in the list
                pId = static_cast< const BCIdentifierNatural* >( boundaryCond[ i ] );

                // Number of the current boundary face
                ibF = pId->id();

                // Updating face stuff
                currentBdFE.updateMeasNormal( mesh.boundaryFacet( ibF ) );

                // Loop on total DOF per Face
                for ( ID l = 0; l < nDofF; ++l )
                {

                    gDof = pId->boundaryLocalToGlobalMap( l );

                    // Loop on components involved in this boundary condition
                    for ( UInt ic = 0; ic < nComp; ++ic )
                    {
                        icDof = gDof + ic * totalDof + offset;

                        // Loop on quadrature points
                        for ( int iq = 0; iq < (int)currentBdFE.nbQuadPt(); ++iq )
                        {
                            sum=0.0;
                            // data on quadrature point
                            for ( ID m = 0; m < nDofF; ++m )
                                sum +=  boundaryCond( pId->boundaryLocalToGlobalMap( m ) , 0 ) * currentBdFE.phi( int( m ), iq );
                            // Adding right hand side contribution
                            rhsRepeated[ icDof ] += sum * currentBdFE.phi( int( l ), iq ) * currentBdFE.normal( int( ic ), iq )
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
            for ( ID i = 0; i < boundaryCond.list_size(); ++i )
            {
                // Pointer to the i-th itdentifier in the list
                pId = static_cast< const BCIdentifierNatural* >( boundaryCond[ i ] );

                // Number of the current boundary face
                ibF = pId->id();

                // Updating face stuff
                currentBdFE.updateMeas( mesh.boundaryFacet( ibF ) );

                // Loop on total DOF per Face
                for ( ID idofF = 0; idofF < nDofF; ++idofF )
                {

                    gDof = pId->boundaryLocalToGlobalMap( idofF );

                    // Loop on space dimensions
                    for ( ID ic = 0; ic < nComp; ++ic )
                        //  for ( UInt ic = 0; ic < nDimensions; ++ic )
                    {
                        icDof = gDof +  boundaryCond.component( ic ) * totalDof+ offset;   //Components passed separately

                        // Loop on quadrature points
                        for ( int iq = 0; iq < (int)currentBdFE.nbQuadPt(); ++iq )
                        {
                            sum = 0;
                            // data on quadrature point
                            for ( ID m = 0; m < nDofF; ++m )
                                sum +=  boundaryCond( pId->boundaryLocalToGlobalMap( m ) , boundaryCond.component( ic ) ) * currentBdFE.phi( int( m ), iq );  //Components passed separatedly

                            // Adding right hand side contribution
                            rhsRepeated[ icDof ] += sum *  currentBdFE.phi( int( idofF ), iq ) *
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
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {
            // Pointer to the i-th identifier in the list
            pId = static_cast< const BCIdentifierNatural* >( boundaryCond[ i ] );
            // Number of the current boundary face
            ibF = pId->id();
            // Updating face stuff
            currentBdFE.updateMeasNormalQuadPt( mesh.boundaryFacet( ibF ) );
            // Loop on total DOF per Face
            for ( ID idofF = 0; idofF < nDofF; ++idofF )
            {
                for ( ID j = 0; j < nComp; ++j )
                {
                    //global Dof
                    idDof = pId->boundaryLocalToGlobalMap( idofF ) + boundaryCond.component( j ) * totalDof + offset;
                    // Loop on quadrature points
                    for ( int iq = 0; iq < (int)currentBdFE.nbQuadPt(); ++iq )
                    {
                        // quadrature point coordinates
                        x = currentBdFE.quadPt(iq, 0);
                        y = currentBdFE.quadPt(iq, 1);
                        z = currentBdFE.quadPt(iq, 2);

                        switch (boundaryCond.mode())
                        {
                        case Full:
                            rhsRepeated[ idDof ] += currentBdFE.phi( int( idofF ), iq ) * boundaryCond( time, x, y, z, boundaryCond.component( j ) ) *
                                                    currentBdFE.weightMeas( iq );
                            break;
                        case Component:
                            rhsRepeated[ idDof ] += currentBdFE.phi( int( idofF ), iq ) * boundaryCond( time, x, y, z, boundaryCond.component( j ) ) *
                                                    currentBdFE.weightMeas( iq );
                            break;
                        case Normal:
                            rhsRepeated[ idDof ] += boundaryCond( time, x, y, z, boundaryCond.component( j ) )*
                                                    currentBdFE.phi( int( idofF ), iq )*
                                                    currentBdFE.weightMeas( iq )*currentBdFE.normal( int(j), iq );
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
        rightHandSide += rhsRepeated;
    }
} // bcNaturalManage



template <typename VectorType, typename MeshType, typename DataType>
void
bcNaturalManageUDep( Real (*mu)(Real time,Real x, Real y, Real z, Real u),
                     VectorType& rightHandSide,
                     const MeshType& mesh,
                     const DOF& dof,
                     const BCBase& boundaryCond,
                     CurrentBoundaryFE& currentBdFE,
                     const DataType& time,
                     const VectorType& feVec,
                     UInt offset )
{

    // Number of local DOF (i.e. nodes) in this face
    UInt nDofF = currentBdFE.nbNode();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    const BCIdentifierNatural* pId;
    ID ibF, idDof ;

    if ( boundaryCond.isDataAVector() )
    { //! If BC is given under a vectorial form
        ERROR_MSG( "This type of BCVector does not exists on bc depentent on solution\n" );
    }
    else
    {  //! If BC is given under a functional form

        DataType x, y, z;
        VectorType rhsRepeated(rightHandSide.map(),Repeated);

        if (nComp!=0)
        {
            ERROR_MSG("For now bcNaturalManageUDep cannot handle non scalar solutions\n");
        }

        // Loop on BC identifiers
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {

            // Pointer to the i-th itdentifier in the list
            pId = static_cast< const BCIdentifierNatural* >( boundaryCond[ i ] );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            currentBdFE.updateMeasQuadPt( mesh.boundaryFacet( ibF ) );

            std::vector<Real> locU(nDofF);    //assumes feVec is a vec of reals, TODO: deal with more comp
            Real uPt;            //value in the point
            for (ID idofLocU=0; idofLocU<nDofF; idofLocU++)
            {
                ID idGDofU=pId->boundaryLocalToGlobalMap(idofLocU)+ boundaryCond.component( 0 ) * totalDof + offset;
                locU[idofLocU]=feVec[idGDofU];
            }


            // Loop on total Dof per Face
            for ( ID idofF = 0; idofF < nDofF; ++idofF )
            {
                // Loop on components involved in this boundary condition
                for ( ID j = 0; j < nComp; ++j )
                {

                    //global Dof
                    idDof = pId->boundaryLocalToGlobalMap( idofF ) + boundaryCond.component( j ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                    {
                    	// quadrature point coordinates
						x = currentBdFE.quadPt(l, 0);
						y = currentBdFE.quadPt(l, 1);
						z = currentBdFE.quadPt(l, 2);

                        uPt=0.0;
                        for (ID idofLocU=0; idofLocU<nDofF; idofLocU++)
                        {
                            uPt+=locU[idofLocU]*currentBdFE.phi( int( idofLocU  ),l );
                        }

                        // Adding right hand side contribution
                        rhsRepeated[ idDof ] += currentBdFE.phi( int( idofF ), l ) * boundaryCond( time, x, y, z, boundaryCond.component( j ),uPt ) *
                                                  mu(time,x,y,z,uPt)*currentBdFE.weightMeas( l );
                    }
                }
            }
        }

        rhsRepeated.globalAssemble();
        ASSERT( rightHandSide.mapType() == Unique , "here rightHandSide should passed as unique, otherwise not sure of what happens at the cpu interfaces ." );
        rightHandSide += rhsRepeated;
    }
}




// ===================================================
// Robin BC
// ===================================================


template <typename MatrixType, typename VectorType, typename DataType, typename MeshType>
void
bcRobinManage( MatrixType& matrix,
               VectorType& rightHandSide,
               const MeshType& mesh,
               const DOF& dof,
               const BCBase& boundaryCond,
               CurrentBoundaryFE& currentBdFE,
               const DataType& time,
               UInt offset )
{
	bcRobinManageMatrix( matrix, mesh, dof, boundaryCond, currentBdFE, time, offset );
    bcRobinManageVector( rightHandSide, mesh, dof, boundaryCond, currentBdFE, time, offset );
}  //bcRobinManage


template <typename MatrixType, typename DataType, typename MeshType>
void
bcRobinManageMatrix( MatrixType& matrix,
                     const MeshType& mesh,
                     const DOF& dof,
                     const BCBase& boundaryCond,
                     CurrentBoundaryFE& currentBdFE,
                     const DataType& time,
                     UInt offset )
{
    // Open the matrix if it is closed
    if ( matrix.matrixPtr()->Filled() )
        matrix.openCrsMatrix();

    // Number of local DOF in this face
    UInt nDofF = currentBdFE.nbNode();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    DataType sum;

    const BCIdentifierNatural* pId;
    ID ibF, idDof, jdDof, kdDof;

    if ( boundaryCond.isDataAVector() )
    {   //! If BC is given under a vectorial form

        //! for the moment, only one coefficient per BCvector.
        DataType mcoef;

        // Loop on BC identifiers
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {

            // Pointer to the i-th identifier in the list
            pId = static_cast< const BCIdentifierNatural* >( boundaryCond[ i ] );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            currentBdFE.updateMeas( mesh.boundaryFacet( ibF ) );

            // Loop on total DOF per Face
            for ( ID idofF = 0; idofF < nDofF; ++idofF )
            {
                // Loop on components invoved in this boundary condition
                for ( ID j = 0; j < nComp; ++j )
                {

                    sum = 0;

                    // Global Dof
                    idDof = pId->boundaryLocalToGlobalMap( idofF ) + boundaryCond.component( j ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( UInt l = 0; l < currentBdFE.nbQuadPt(); ++l )
                    {
                        mcoef = 0.0;
                        for ( UInt n = 0; n < nDofF; ++n )
                        {
                            kdDof=pId->boundaryLocalToGlobalMap( n );
                            if (boundaryCond.isRobinCoeffAVector())
                                mcoef += boundaryCond.robinCoeffVector( kdDof, boundaryCond.component( j ) ) * currentBdFE.phi( n, l );
                            else
                                mcoef += boundaryCond.robinCoeff() * currentBdFE.phi( n, l );
                        }

                        // Contribution to the diagonal entry of the elementary boundary mass matrix
                        sum += mcoef * currentBdFE.phi( idofF, l ) * currentBdFE.phi( idofF, l ) * currentBdFE.weightMeas( l );
                    }

                    // Assembling diagonal entry
                    matrix.addToCoefficient( idDof, idDof, sum );
                }

                // Upper diagonal columns of the elementary boundary mass matrix
                for ( ID k = idofF + 1 ; k <nDofF ; ++k )
                {

                    // Loop on components invoved in this boundary condition
                    for ( ID j = 0; j < nComp; ++j )
                    {

                        sum = 0;

                        // Loop on quadrature points
                        for ( UInt l = 0; l < currentBdFE.nbQuadPt(); ++l )
                        {
                            mcoef = 0.0;
                            for ( UInt n = 0; n < nDofF; ++n)
                            {
                                kdDof = pId->boundaryLocalToGlobalMap( n );
                                if (boundaryCond.isRobinCoeffAVector())
                                    mcoef += boundaryCond.robinCoeffVector( kdDof, boundaryCond.component( j ) ) * currentBdFE.phi( n, l );

                                else
                                    mcoef += boundaryCond.robinCoeff() * currentBdFE.phi( n, l );
                            }

                            // Upper diagonal entry of the elementary boundary mass matrix
                            sum += mcoef * currentBdFE.phi( idofF, l ) * currentBdFE.phi( k, l ) *
                                   currentBdFE.weightMeas( l );
                        }

                        // Globals DOF: row and columns
                        idDof = pId->boundaryLocalToGlobalMap( idofF ) + boundaryCond.component( j ) * totalDof + offset;
                        jdDof = pId->boundaryLocalToGlobalMap( k ) + boundaryCond.component( j ) * totalDof + offset;

                        // Assembling upper entry.  The boundary mass matrix is symetric
                        matrix.addToCoefficient( idDof, jdDof, sum );
                        matrix.addToCoefficient( jdDof, idDof, sum );
                    }
                }
            }
        }
    }

    else
    {  //! If BC is given under a functional form

        DataType x, y, z;

        const BCFunctionRobin* pBcF = static_cast<const BCFunctionRobin*>( boundaryCond.pointerToFunctor() );

        // Loop on BC identifiers
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {

            // Pointer to the i-th identifier in the list
            pId = static_cast< const BCIdentifierNatural* >( boundaryCond[ i ] );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            currentBdFE.updateMeasQuadPt( mesh.boundaryFacet( ibF ) );

            // Loop on total DOF per Face
            for ( ID idofF = 0; idofF < nDofF; ++idofF )
            {
                // Loop on components involved in this boundary condition
                for ( ID j = 0; j < nComp; ++j )
                {
                    sum = 0;

                    // Global DOF (outside the quad point loop. V. Martin)
                    idDof = pId->boundaryLocalToGlobalMap( idofF ) + boundaryCond.component( j ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( UInt l = 0; l < currentBdFE.nbQuadPt(); ++l )
                    {
                    	// quadrature point coordinates
						x = currentBdFE.quadPt(l, 0);
						y = currentBdFE.quadPt(l, 1);
						z = currentBdFE.quadPt(l, 2);

                        // Contribution to the diagonal entry of the elementary boundary mass matrix
                        sum += pBcF->coef( time, x, y, z, boundaryCond.component(j) ) * currentBdFE.phi( idofF, l ) * currentBdFE.phi( idofF, l ) *
                               currentBdFE.weightMeas( l );
                    }


                    // Assembling diagonal entry
                    matrix.addToCoefficient( idDof, idDof, sum );
                }

                // Upper diagonal columns of the elementary boundary mass matrix
                for ( ID k = idofF + 1 ; k < nDofF ; ++k )
                {

                    // Loop on components involved in this boundary condition
                    for ( ID j = 0; j < nComp; ++j )
                    {
                        sum = 0;

                        // Loop on quadrature points
                        for ( UInt l = 0; l < currentBdFE.nbQuadPt(); ++l )
                        {

                        	// quadrature point coordinates
							x = currentBdFE.quadPt(l, 0);
							y = currentBdFE.quadPt(l, 1);
							z = currentBdFE.quadPt(l, 2);

                            // Upper diagonal entry of the elementary boundary mass matrix
                            sum += pBcF->coef( time, x, y, z, boundaryCond.component( j )  ) * currentBdFE.phi( idofF, l ) * currentBdFE.phi( k, l ) *
                                   currentBdFE.weightMeas( l );
                        }

                        // Globals DOF: row and columns
                        idDof = pId->boundaryLocalToGlobalMap( idofF ) + boundaryCond.component( j ) * totalDof + offset;
                        jdDof = pId->boundaryLocalToGlobalMap( k ) + boundaryCond.component( j ) * totalDof + offset;

                        // Assembling upper entry.  The boundary mas matrix is symetric
                        matrix.addToCoefficient( idDof, jdDof, sum );
                        matrix.addToCoefficient( jdDof, idDof, sum );
                    }
                }
            }
        }
    }
}   //bcRobinManageMatrix



template <typename VectorType, typename DataType, typename MeshType>
void
bcRobinManageVector( VectorType& rightHandSide,
                     const MeshType& mesh,
                     const DOF& dof,
                     const BCBase& boundaryCond,
                     CurrentBoundaryFE& currentBdFE,
                     const DataType& time,
                     UInt offset )
{

    // Number of local DOF in this face
    UInt nDofF = currentBdFE.nbNode();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    const BCIdentifierNatural* pId;
    ID ibF, idDof, kdDof;

    if ( boundaryCond.isDataAVector() )
    {   //! If BC is given under a vectorial form

        VectorType rhsRepeated(rightHandSide.map(),Repeated);

        //! Defining the coefficients
        DataType mbcb;

        // Loop on BC identifiers
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {

            // Pointer to the i-th identifier in the list
            pId = static_cast< const BCIdentifierNatural* >( boundaryCond[ i ] );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            currentBdFE.updateMeas( mesh.boundaryFacet( ibF ) );

            // Loop on total DOF per Face
            for ( ID idofF = 0; idofF < nDofF; ++idofF )
            {
                // Loop on components invoved in this boundary condition
                for ( ID j = 0; j < nComp; ++j )
                {
                    // Global Dof
                    idDof = pId->boundaryLocalToGlobalMap( idofF ) + boundaryCond.component( j ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( UInt l = 0; l < currentBdFE.nbQuadPt(); ++l )
                    {
                        mbcb = 0.0;
                        for ( UInt n = 0; n < nDofF; ++n )
                        {
                            kdDof=pId->boundaryLocalToGlobalMap( n );

                            if ( boundaryCond.isBetaCoeffAVector() )
                                mbcb += boundaryCond.betaCoeffVector( kdDof, boundaryCond.component( j ) )
                                        * boundaryCond( kdDof, boundaryCond.component( j )) * currentBdFE.phi( n, l );
                            else
                                mbcb += boundaryCond.betaCoeff() * boundaryCond( kdDof, boundaryCond.component( j )) * currentBdFE.phi( n, l );
                        }

                        // Adding right hand side contribution
                        rhsRepeated[ idDof ] += currentBdFE.phi( idofF, l ) * mbcb * currentBdFE.weightMeas( l );
                    }
                }
            }
        }
        rhsRepeated.globalAssemble();
        ASSERT( rightHandSide.mapType() == Unique, "here rightHandSide should passed as unique, otherwise not sure of what happens at the cpu interfaces ." );
        rightHandSide += rhsRepeated;
    }

    else
    {  //! If BC is given under a functional form

        VectorType rhsRepeated(rightHandSide.map(),Repeated);

        DataType x, y, z;

        // Loop on BC identifiers
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {

            // Pointer to the i-th identifier in the list
            pId = static_cast< const BCIdentifierNatural* >( boundaryCond[ i ] );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            currentBdFE.updateMeasQuadPt( mesh.boundaryFacet( ibF ) );

            // Loop on total DOF per Face
            for ( ID idofF = 0; idofF < nDofF; ++idofF )
            {
                // Loop on components involved in this boundary condition
                for ( ID j = 0; j < nComp; ++j )
                {

                    // Global DOF (outside the quad point loop. V. Martin)
                    idDof = pId->boundaryLocalToGlobalMap( idofF ) + boundaryCond.component( j ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( UInt l = 0; l < currentBdFE.nbQuadPt(); ++l )
                    {
                    	// quadrature point coordinates
						x = currentBdFE.quadPt(l, 0);
						y = currentBdFE.quadPt(l, 1);
						z = currentBdFE.quadPt(l, 2);

                        // Adding right hand side contribution
                        rhsRepeated[ idDof ] += currentBdFE.phi( idofF, l ) * boundaryCond( time, x, y, z, boundaryCond.component( j ) ) *
                                                  currentBdFE.weightMeas( l );
                    }
                }
            }
        }
        rhsRepeated.globalAssemble();
        ASSERT( rightHandSide.mapType() == Unique, "here rightHandSide should passed as unique, otherwise not sure of what happens at the cpu interfaces ." );
        rightHandSide += rhsRepeated;
    }
} //bcRobinManageVector




template <typename VectorType, typename DataType, typename MeshType>
void bcRobinManageResidual( VectorType& residual,
                       VectorType& rightHandSide,
                       const VectorType& solution, //solution must not be repeated
                       const MeshType& mesh,
                       const DOF& dof,
                       const BCBase& boundaryCond,
                       CurrentBoundaryFE& currentBdFE,
                       const DataType& time,
                       UInt offset )
{

    if(solution.mapType()==Repeated)
    {
        std::cout<<"pass me a non-repeated solution"<<std::endl;
        VectorType uniqueSolution(solution, Unique, Zero);
        bcRobinManageResidual(  residual,
                                rightHandSide,
                                uniqueSolution,
                                mesh,
                                dof,
                                boundaryCond,
                                currentBdFE,
                                time,
                                offset );
        return;
    }

    MatrixEpetra<Real> matrix(solution.map());
    bcRobinManage( matrix, rightHandSide, mesh, dof, boundaryCond, currentBdFE, time, offset );
    matrix.globalAssemble();
    residual += matrix*solution;
}   //bcRobinManageResidual






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
              const DOF&      dof,
              const BCBase&   boundaryCond,
              CurrentBoundaryFE&    currentBdFE,
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
    rightHandSide.setCoefficient(offset, boundaryCond(time, 0., 0., 0., 0));
}  //bcFluxManageVector


template <typename MatrixType,
typename MeshType,
typename DataType>
void
bcFluxManageMatrix( MatrixType&     matrix,
                    const MeshType& mesh,
                    const DOF&      dof,
                    const BCBase&   boundaryCond,
                    CurrentBoundaryFE&    currentBdFE,
                    const DataType& /*time*/,
                    UInt            offset )
{
    if ( matrix.matrixPtr()->Filled() )
        matrix.openCrsMatrix();

    // Number of local DOF in this face
    UInt nDofF = currentBdFE.nbNode();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    DataType sum;

    const BCIdentifierNatural* pId;
    ID ibF, idDof, jdDof/*, kdDof*/;

    if ( !boundaryCond.isDataAVector() )
    {
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {
            pId = static_cast< const BCIdentifierNatural* >( boundaryCond[ i ] );

            // Number of the current boundary face
            ibF = pId->id();
            // Updating face stuff
            currentBdFE.updateMeasNormal( mesh.boundaryFacet( ibF ) );

            for ( ID idofF = 0; idofF < nDofF; ++idofF )
            {
                for ( int ic = 0; ic < (int)nComp; ++ic)
                {
                    idDof = pId->boundaryLocalToGlobalMap( idofF ) + ic * totalDof;

                    sum = 0.;
                    for ( int iq = 0; iq < (int)currentBdFE.nbQuadPt(); ++iq )
                    {
                        sum += currentBdFE.phi( int( idofF ), iq )*
                               currentBdFE.normal(ic , iq)*
                               currentBdFE.weightMeas(iq);
                    }

                    jdDof = offset;

                    matrix.addToCoefficient( idDof    , jdDof    , sum );
                    matrix.addToCoefficient( jdDof    , idDof    , sum );
                }
            }
        }
    }
} // bcFluxManageMatrix

template <typename VectorType,
typename MeshType,
typename DataType>
void
bcFluxManageResidual( VectorType&      residual,
                      VectorType&      rightHandSide,
                      const VectorType&     solution,
                      const MeshType&  mesh,
                      const DOF&       dof,
                      const BCBase&    boundaryCond,
                      CurrentBoundaryFE&    currentBdFE,
                      const DataType&  time,
                      UInt             offset )
{
    if(solution.mapType()==Repeated)
    {
        std::cout<<"pass me a non-repeated solution"<<std::endl;
        VectorType uniqueSolution(solution, Unique, Zero);
        bcFluxManageResidual(  residual,
                               rightHandSide,
                                uniqueSolution,
                                mesh,
                                dof,
                                boundaryCond,
                                currentBdFE,
                                time,
                                offset );
        return;
    }
    MatrixEpetra<Real> matrix(solution.map());
    bcFluxManage(matrix, rightHandSide, mesh, dof, boundaryCond, currentBdFE, time, offset );
    matrix.globalAssemble();
    residual += matrix*solution;
}



template <typename MatrixType, typename VectorType, typename DataType, typename MeshType>
void
bcResistanceManage( MatrixType& matrix,
                    VectorType& rightHandSide,
                    const MeshType& mesh,
                    const DOF& dof,
                    const BCBase& boundaryCond,
                    CurrentBoundaryFE& currentBdFE,
                    const DataType& time,
                    UInt offset )
{
    bcResistanceManageMatrix( matrix, mesh, dof, boundaryCond, currentBdFE, time, offset );
    bcResistanceManageVector( rightHandSide, mesh, dof, boundaryCond, currentBdFE, time, offset );
} //bcResistanceManage


template <typename VectorType, typename DataType, typename MeshType>
void
bcResistanceManageVector( VectorType& rightHandSide,
                    const MeshType& mesh,
                    const DOF& dof,
                    const BCBase& boundaryCond,
                    CurrentBoundaryFE& currentBdFE,
                    const DataType& /*time*/,
                    UInt offset )
{
    // Number of local DOF in this face
    UInt nDofF = currentBdFE.nbNode();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    std::set<ID> resistanceDofs;

    const BCIdentifierNatural* pId;
    ID ibF, idDof, kdDof;

    if ( boundaryCond.isDataAVector() )
    {
        VectorType rhsRepeated(rightHandSide.map(),Repeated);

        DataType  mbcb;

        // Loop on BC identifiers
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {
            // Pointer to the i-th identifier in the list
            pId = static_cast< const BCIdentifierNatural* >( boundaryCond[ i ] );

            // Number of the current boundary face
            ibF = pId->id();

            currentBdFE.updateMeasNormal( mesh.boundaryFacet( ibF ) );

            // Loop on total DOF per Face
            for ( ID idofF = 0; idofF < nDofF; ++idofF )
            {
                resistanceDofs.insert( pId->boundaryLocalToGlobalMap( idofF ) );

                // Loop on components involved in this boundary condition
                for ( ID j = 0; j < nComp; ++j )
                {
                    idDof = pId->boundaryLocalToGlobalMap( idofF ) + boundaryCond.component( j ) * totalDof + offset;

                    // Loop on quadrature points
                    for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                    {
                        mbcb=0;

                        // data on quadrature point
                        for ( ID n = 0; n < nDofF; ++n)
                        {
                            kdDof=pId->boundaryLocalToGlobalMap( n );
                            mbcb += boundaryCond( kdDof, boundaryCond.component( j ) )* currentBdFE.phi( int( n ), l ) ;
                        }

                        rhsRepeated[ idDof ] +=  mbcb* currentBdFE.phi( int( idofF ), l ) *  currentBdFE.normal( int( j ), l ) *
                                                   currentBdFE.weightMeas( l );

                    }
                }
            }
        }
        rhsRepeated.globalAssemble();
        ASSERT( rightHandSide.mapType() == Unique, "here rightHandSide should passed as unique, otherwise not sure of what happens at the cpu interfaces ." );
        rightHandSide += rhsRepeated;
    }
    else
        ERROR_MSG( "This BC type is not yet implemented" );
} //bcResistanceManageVector


template <typename MatrixType, typename DataType, typename MeshType>
void
bcResistanceManageMatrix( MatrixType& matrix,
                    const MeshType& mesh,
                    const DOF& dof,
                    const BCBase& boundaryCond,
                    CurrentBoundaryFE& currentBdFE,
                    const DataType& /*time*/,
                    UInt offset )
{
    // Open the matrix if it is closed:
    if ( matrix.matrixPtr()->Filled() )
        matrix.openCrsMatrix();

    // Number of local DOF in this face
    UInt nDofF = currentBdFE.nbNode();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = boundaryCond.numberOfComponents();

    std::set<ID> resistanceDofs;

    const BCIdentifierNatural* pId;
    ID ibF, idDof, jdDof;

    if ( boundaryCond.isDataAVector() )
    {
        //auxiliary vector
        VectorEpetra vv(boundaryCond.pointerToBCVector()->rhsVector().map(), Repeated);
        vv *= 0.;

        // Loop on BC identifiers
        for ( ID i = 0; i < boundaryCond.list_size(); ++i )
        {
            // Pointer to the i-th identifier in the list
            pId = static_cast< const BCIdentifierNatural* >( boundaryCond[ i ] );

            // Number of the current boundary face
            ibF = pId->id();

            currentBdFE.updateMeasNormal( mesh.boundaryFacet( ibF ) );

            // Loop on total DOF per Face
            for ( ID idofF = 0; idofF < nDofF; ++idofF )
            {
                resistanceDofs.insert( pId->boundaryLocalToGlobalMap( idofF ) );

                // Loop on components involved in this boundary condition
                for ( ID j = 0; j < nComp; ++j )
                {
                    idDof = pId->boundaryLocalToGlobalMap( idofF ) + boundaryCond.component( j ) * totalDof + offset;

                    // std::cout << "\nDOF " << idDof << " is involved in Resistance BC" << std::endl;

                    // Loop on quadrature points
                    for ( int l = 0; l < (int)currentBdFE.nbQuadPt(); ++l )
                    {
                        vv[idDof] += currentBdFE.phi( int( idofF), l ) *  currentBdFE.normal( int( j ), l ) * currentBdFE.weightMeas( l );
                    }
                }
            }
        }

        // this vector now contains the values needed to modify the matrix
        vv.globalAssemble();

        // I want it on a unique processor (the root processor) that will take care of modifying
        // the matrix on the LHS
        VectorEpetra vvReduced( vv, 0 );

        // I need to tell the root processor what are the IDs of the DOFs on the "resistance" boundary.
        // Each processor finds numMyResistanceDofs in its own portion of the mesh
        Int numMyResistanceDofs( resistanceDofs.size() );
        // Summing together the numMyResistanceDofs values, I obtain the global number of
        // resistance DOFs, **including repeated DOFs**
        Int numGlobalResistanceDofs(0);
        vv.map().comm().SumAll(&numMyResistanceDofs, &numGlobalResistanceDofs, 1);

        // I need to share the list of resistance DOFs, via MPI calls. I will need vectors (arrays), not sets.
        // Each process will store its resistance DOFs in a vector
        std::vector<Int> myResistanceDofs( numGlobalResistanceDofs, -1 );
        // And each processor will receive the other processors' resistance DOFs in a gathered vector
        // i.e. an "array of arrays", the i-th sub-array being a copy of myResistanceDofs from processor i
        std::vector<Int> globalResistanceDofs( numGlobalResistanceDofs*vv.map().comm().NumProc(), 0 );

        // Fill myResistanceDofs with the actual DOF IDs (only the first numMyResistanceDofs will be nonzero)
        UInt iCount(0);
        for ( std::set<ID>::iterator iDofIt = resistanceDofs.begin();
                        iDofIt != resistanceDofs.end(); ++iDofIt, ++iCount )
        {
            myResistanceDofs[iCount] = *iDofIt;
        }

        // Gather the lists of resistance DOF IDs from all processors
        vv.map().comm().GatherAll(&myResistanceDofs[0], &globalResistanceDofs[0], numGlobalResistanceDofs);

        // Create a unique list of IDs: here I make use of a set
        std::set<ID> globalResistanceDofSet;
        for( Int iDof = 0; iDof < numGlobalResistanceDofs*vv.map().comm().NumProc(); ++iDof )
        {
            if( globalResistanceDofs[iDof] > -1 )
            {
                globalResistanceDofSet.insert( globalResistanceDofs[iDof] );
                // std::cout << "\n(after gathering) DOF " << globalResistanceDofs[iDof] << std::endl;
            }
        }
        // std::cout << "\n(after gathering) number of DOFs = " << globalResistanceDofs.size() << std::endl;

        // Only the root processor has the needed values to modify the matrix
        if( ! vv.map().comm().MyPID() )
        {
            for ( std::set<ID>::iterator iDofIt = globalResistanceDofSet.begin();
                            iDofIt != globalResistanceDofSet.end(); ++iDofIt )
            {
                for ( UInt iComp = 0; iComp < nComp; ++iComp )
                {
                    idDof = *iDofIt + boundaryCond.component( iComp ) * totalDof + offset;
                    for ( std::set<ID>::iterator jDofIt = globalResistanceDofSet.begin();
                                    jDofIt != globalResistanceDofSet.end(); ++jDofIt )
                    {
                        for ( UInt jComp = 0; jComp < nComp; ++jComp )
                        {
                            jdDof = *jDofIt + boundaryCond.component( jComp ) * totalDof + offset;

                            matrix.addToCoefficient( idDof,  jdDof, boundaryCond.resistanceCoeff() *
                                                     vvReduced[idDof] * vvReduced[jdDof] );
                        }
                    }
                }
            }
        }
    }
    else
        ERROR_MSG( "This BC type is not yet implemented" );

} //bcResistanceManageMatrix

} // end of namespace LifeV
#endif
