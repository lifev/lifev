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


    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
    {

        switch ( BCh[ i ].type() )
        {
        case Essential:  // Essential boundary conditions (Dirichlet)
        if(BCh[ i ].isUDep())
          bcEssentialManageUDep(A, b, mesh, dof, BCh[ i ], bdfem, coef, t,U);
        else
              bcEssentialManage( A, b, mesh, dof, BCh[ i ], bdfem, coef, t );
            break;
        case Natural:  // Natural boundary conditions (Neumann)
        if(BCh[ i ].isUDep())
          bcNaturalManageUDep(mu, b, mesh, dof, BCh[ i ], bdfem, t,U);
            else
          //in this case mu must be a constant, think about (not still implemented)
              bcNaturalManage( b, mesh, dof, BCh[ i ], bdfem, t );
            break;
        case Mixte:  // Mixte boundary conditions (Robin)
        if(BCh[ i ].isUDep())
          bcMixteManageUDep( A, b, mesh, dof, BCh[ i ], bdfem, t, U);    //not still implemented
        else
              bcMixteManage( A, b, mesh, dof, BCh[ i ], bdfem, t );
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

            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
            {
              // Loop on components involved in this boundary condition
              for ( ID j = 1; j <= nComp; ++j )
              {
                // Glogal Dof
                idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;

#if USE_BOOST_MATRIX
                using namespace boost::numeric::ublas;
                matrix_row<MatrixType> mr (M, idDof-1);
                mr *= 0;
                M( idDof-1, idDof-1 ) = coef;
#else
                // Modifying matrix
                M.diagonalize_row( idDof - 1, coef);
#endif
              }
            }
          }
        }
    }
}



template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void bcManage( MatrixType& A, VectorType& b, const MeshType& mesh, const Dof& dof,
                const BCHandler& BCh,
                CurrentBdFE& bdfem, const DataType& coef, const DataType& t )
{


    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
    {

        switch ( BCh[ i ].type() )
        {
        case Essential:  // Essential boundary conditions (Dirichlet)
            bcEssentialManage( A, b, mesh, dof, BCh[ i ], bdfem, coef, t );
            break;
        case Natural:  // Natural boundary conditions (Neumann)
            bcNaturalManage( b, mesh, dof, BCh[ i ], bdfem, t );
            break;
        case Mixte:  // Mixte boundary conditions (Robin)
            bcMixteManage( A, b, mesh, dof, BCh[ i ], bdfem, t );
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }
}


//Version that treates only the matrix modifications
//Miguel:10/02  - Mixte : V. Martin: 03/03 // I added the time t for the mixte case. V. Martin
template <typename MatrixType, typename MeshType, typename DataType>
void bcManageMatrix( MatrixType& A, const MeshType& mesh, const Dof& dof,
                       const BCHandler& BCh,
                       CurrentBdFE& bdfem, const DataType& coef, const DataType& t = 0 )
{

    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
    {

        switch ( BCh[ i ].type() )
        {
        case Essential:  // Essential boundary conditions (Dirichlet)
            bcEssentialManageMatrix( A, dof, BCh[ i ], coef );  //! Bug here???
            break;
        case Natural:  // Natural boundary conditions (Neumann)
            // Do nothing
            break;
        case Mixte:  // Mixte boundary conditions (Robin)
            bcMixteManageMatrix( A, mesh, dof, BCh[ i ], bdfem, t );
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }
}


//Version that treates only the vector modifications
//Miguel:10/02  - Mixte : V. Martin: 03/03
template <typename VectorType, typename MeshType, typename DataType>
void bcManageVector( VectorType& b, const MeshType& mesh, const Dof& dof,
                       const BCHandler& BCh, CurrentBdFE& bdfem, const DataType& t, const DataType& coef )
{

    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
    {

        switch ( BCh[ i ].type() )
        {
        case Essential:  // Essential boundary conditions (Dirichlet)
            bcEssentialManageVector( b, dof, BCh[ i ], t, coef );
            break;
        case Natural:  // Natural boundary conditions (Neumann)
            bcNaturalManage( b, mesh, dof, BCh[ i ], bdfem, t );
            break;
        case Mixte:  // Mixte boundary conditions (Robin)
            bcMixteManageVector( b, mesh, dof, BCh[ i ], bdfem, t );
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }
}



//! Version for mixed problem with block matrices
//! Alain, 07/08/02
template <typename MatrixType1, typename MatrixType2, typename VectorType,
typename MeshType, typename DataType>
void bcManage( MatrixType1& C, MatrixType2& trD, VectorType& b,
                const MeshType& mesh, const Dof& dof, const BCHandler& BCh,
                CurrentBdFE& bdfem, const DataType& coef, const DataType& t )
{
    // Loop on boundary conditions
    for ( Index_t i = 0; i < BCh.size(); ++i )
    {

        switch ( BCh[ i ].type() )
        {
        case Essential:  // Essentila boundary conditions (Dirichlet)
            bcEssentialManage( C, trD, b, mesh, dof, BCh[ i ], bdfem, coef, t );
            break;
        case Natural:  // Natural boundary conditions (Neumann)
            bcNaturalManage( b, mesh, dof, BCh[ i ], bdfem, t );
            break;
        case Mixte:  // Mixte boundary conditions (Robin)
            bcMixteManage( C, trD, b, mesh, dof, BCh[ i ], bdfem, t );
            break;
        default:
            ERROR_MSG( "This BC type is not yet implemented" );
        }
    }
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
            bcEssentialManage( C, trD, D, b, bp, mesh, dof, BCh[ i ], bdfem, coef, t );
            break;
        case Natural:  // Natural boundary conditions (Neumann)
            bcNaturalManage( b, mesh, dof, BCh[ i ], bdfem, t );
            break;
        case Mixte:  // Mixte boundary conditions (Robin)
            bcMixteManage( C, trD, D, b, bp, mesh, dof, BCh[ i ], bdfem, t );
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
void bcEssentialManageUDep( MatrixType& A, VectorType& b, const MeshType& mesh, const Dof& dof,
    const BCBase& BCb, const CurrentBdFE& bdfem, const DataType& coef,
    const DataType& t, const VectorType& U )
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

        DataType x, y, z;
        // Loop on BC identifiers
        for ( ID i = 1; i <= BCb.list_size(); ++i )
        {
            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {
                // Glogal Dof
                idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                // Coordinates of the node where we impose the value
                x = static_cast< const IdentifierEssential* >( BCb( i ) ) ->x();
                y = static_cast< const IdentifierEssential* >( BCb( i ) ) ->y();
                z = static_cast< const IdentifierEssential* >( BCb( i ) ) ->z();

                Real datum = BCb( t, x, y, z, BCb.component( j ) ,U[idDof-1]);


#if USE_BOOST_MATRIX
                using namespace boost::numeric::ublas;
                matrix_row<MatrixType> mr (A, idDof-1);
                mr *= 0;
                matrix_column<MatrixType> mc (A, idDof-1);
                b -= mc*datum; // correct rhs
                mc *= 0;
                A( idDof-1, idDof-1 ) = coef;
                b( idDof-1 ) = coef*datum;
#else
                // Modifying matrix and right hand side
                A.diagonalize( idDof - 1, coef, b, datum );
#endif
            }
        }
    }
}
template <typename MatrixType, typename VectorType, typename MeshType, typename DataType>
void bcEssentialManage( MatrixType& A, VectorType& b, const MeshType& /*mesh*/, const Dof& dof, const BCBase& BCb,
                        const CurrentBdFE& /*bdfem*/, const DataType& coef, const DataType& t )
{

    ID idDof;

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    if ( BCb.dataVector() )
    { //! If BC is given under a vectorial form

        // Loop on BC identifiers
        for ( ID i = 1; i <= BCb.list_size(); ++i )
        {
            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {
                // Global Dof
                idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
#if USE_BOOST_MATRIX
                Real datum = BCb( BCb( i ) ->id(), BCb.component( j ) );

                using namespace boost::numeric::ublas;
                matrix_row<MatrixType> mr (A, idDof-1);
                mr *= 0;
                matrix_column<MatrixType> mc (A, idDof-1);
                b -= mc*datum; // correct rhs
                mc *= 0;
                A( idDof-1, idDof-1 ) = coef;
                b( idDof-1 ) = coef*datum;
#else
                // Modifying matrix and right hand side
                A.diagonalize( idDof - 1, coef, b, BCb( BCb( i ) ->id(), BCb.component( j ) ) );
#endif
            }
        }
    }
    else
    { //! If BC is given under a functional form

        DataType x, y, z;
        // Loop on BC identifiers
        for ( ID i = 1; i <= BCb.list_size(); ++i )
        {
            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {
                // Glogal Dof
                idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                // Coordinates of the node where we impose the value
                x = static_cast< const IdentifierEssential* >( BCb( i ) ) ->x();
                y = static_cast< const IdentifierEssential* >( BCb( i ) ) ->y();
                z = static_cast< const IdentifierEssential* >( BCb( i ) ) ->z();

#if USE_BOOST_MATRIX
                Real datum = BCb( t, x, y, z, BCb.component( j ) );


                using namespace boost::numeric::ublas;
                matrix_row<MatrixType> mr (A, idDof-1);
                mr *= 0;
                matrix_column<MatrixType> mc (A, idDof-1);
                b -= mc*datum; // correct rhs
                mc *= 0;
                A( idDof-1, idDof-1 ) = coef;
                b( idDof-1 ) = coef*datum;
#else
                // Modifying matrix and right hand side
                A.diagonalize( idDof - 1, coef, b, BCb( t, x, y, z, BCb.component( j ) ) );
#endif
            }
        }
    }
}


//Version that treates only the matrix modifications
//Miguel:10/02
template <typename MatrixType, typename DataType>
void bcEssentialManageMatrix( MatrixType& A, const Dof& dof, const BCBase& BCb, const DataType& coef )
{

    ID idDof;
    UInt totalDof;

    // Number of total scalar Dof
    totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    // Loop on BC identifiers
    for ( ID i = 1; i <= BCb.list_size(); ++i )
    {
        // Loop on components involved in this boundary condition
        for ( ID j = 1; j <= nComp; ++j )
        {
            // Glogal Dof
            idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
            // Modifying ONLY matrix
            A.diagonalize_row( idDof - 1, coef );
        }
    }
}

//Version that treates only the vector modifications
//Miguel:10/02
template <typename VectorType, typename DataType>
void bcEssentialManageVector( VectorType& b, const Dof& dof, const BCBase& BCb, const DataType& t, const DataType& coef )
{


    ID idDof;
    UInt totalDof;

    // Number of total scalar Dof
    totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    if ( BCb.dataVector() )
    {  //! If BC is given under a vectorial form

        // Loop on BC identifiers
        for ( ID i = 1; i <= BCb.list_size(); ++i )
        {
            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {

                // Glogal Dof
                idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                // Modifying right hand side
                b( idDof - 1 ) = coef * BCb( BCb( i ) ->id(), BCb.component( j ) );
            }
        }
    }
    else
    {  //! If BC is given under a functionnal form

        DataType x, y, z;
        // Loop on BC identifiers
        for ( ID i = 1; i <= BCb.list_size(); ++i )
        {
            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {
                // Glogal Dof
                idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                // Coordinates of the node where we impose the value
                x = static_cast< const IdentifierEssential* >( BCb( i ) ) ->x();
                y = static_cast< const IdentifierEssential* >( BCb( i ) ) ->y();
                z = static_cast< const IdentifierEssential* >( BCb( i ) ) ->z();
                // Modifying right hand side
                b( idDof - 1 ) = coef * BCb( t, x, y, z, BCb.component( j ) );
            }
        }
    }
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
                        const DataType& t )
{
    ID idDof;
    DataType x, y, z;
    UInt totalDof;

    // Number of total scalar Dof
    totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    if ( BCb.dataVector() )
    {  //! If BC is given under a vectorial form

        // Loop on BC identifiers
        for ( ID i = 1; i <= BCb.list_size(); ++i )
        {
            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {
                // Glogal Dof
                idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                // Modifying matrix and right hand side
                A.diagonalize( idDof - 1, coef, b, BCb( BCb( i ) ->id(), BCb.component( j ) ) );
                trD.zero_row( idDof - 1 );
            }
        }
    }
    else
    {  //! If BC is given under a functionnal form

        // Loop on BC identifiers
        for ( ID i = 1; i <= BCb.list_size(); ++i )
        {

            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {

                // Glogal Dof
                idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                // Coordinates of the node where we impose the value
                x = static_cast< const IdentifierEssential* >( BCb( i ) ) ->x();
                y = static_cast< const IdentifierEssential* >( BCb( i ) ) ->y();
                z = static_cast< const IdentifierEssential* >( BCb( i ) ) ->z();
                // Modifying matrix and right hand side
                A.diagonalize( idDof - 1, coef, b, BCb( t, x, y, z, BCb.component( j ) ) );
                trD.zero_row( idDof - 1 );
            }
        }
    }
}

//! version with column diagonalization of D.
//! Alain. nov 2002.
template <typename MatrixType1, typename MatrixType2, typename MatrixType3,
typename VectorType, typename MeshType, typename DataType>
void bcEssentialManage( MatrixType1& A, MatrixType2& trD, MatrixType3& D,
                        VectorType& b, VectorType& bp, const MeshType& mesh,
                        const Dof& dof, const BCBase& BCb,
                        const CurrentBdFE& bdfem,
                        const DataType& coef, const DataType& t )
{
    ID idDof;
    DataType x, y, z;
    UInt totalDof;

    // Number of total scalar Dof
    totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    if ( BCb.dataVector() )
    {  //! If BC is given under a vectorial form

        // Loop on BC identifiers
        for ( ID i = 1; i <= BCb.list_size(); ++i )
        {
            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {
                // Glogal Dof
                idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                // Modifying matrix and right hand side
                A.diagonalize( idDof - 1, coef, b, BCb( BCb( i ) ->id(), BCb.component( j ) ) );
                zero_row_col( idDof - 1, trD, D, bp, BCb( t, x, y, z, BCb.component( j ) ) );
            }
        }
    }
    else
    {  //! If BC is given under a functionnal form

        // Loop on BC identifiers
        for ( ID i = 1; i <= BCb.list_size(); ++i )
        {

            // Loop on components involved in this boundary condition
            for ( ID j = 1; j <= nComp; ++j )
            {

                // Glogal Dof
                idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                // Coordinates of the node where we impose the value
                x = static_cast< const IdentifierEssential* >( BCb( i ) ) ->x();
                y = static_cast< const IdentifierEssential* >( BCb( i ) ) ->y();
                z = static_cast< const IdentifierEssential* >( BCb( i ) ) ->z();
                // Modifying matrix and right hand side
                A.diagonalize( idDof - 1, coef, b, BCb( t, x, y, z, BCb.component( j ) ) );
                zero_row_col( idDof - 1, trD, D, bp, BCb( t, x, y, z, BCb.component( j ) ) );
            }
        }
    }
}

// ===================================================
// Natural BC
// ===================================================
template <typename VectorType, typename MeshType, typename DataType>
void bcNaturalManageUDep( Real (*mu)(Real t,Real x, Real y, Real z, Real u),
            VectorType& b, const MeshType& mesh, const Dof& dof,
            const BCBase& BCb, CurrentBdFE& bdfem,
            const DataType& t, const VectorType& U )
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
            ID idGDofU=pId->bdLocalToGlobal(idofLocU+1)+( BCb.component( 1 ) - 1 ) * totalDof;
        locU[idofLocU]=U[idGDofU-1];
            }


            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            {  //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                // Loop on components involved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {

                    //glogal Dof
                    idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof;

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
                        b[ idDof - 1 ] += bdfem.phi( int( idofF - 1 ), l ) * BCb( t, x, y, z, BCb.component( j ),uPt ) *
                    mu(t,x,y,z,uPt)*bdfem.weightMeas( l );
//    Debug(800)<<"debug* naturalManageUDep done one ulocU\n";
                    }
                }
            }
        }
    }
//    Debug(800)<<"debug* end of naturalManageUDep\n";
}


template <typename VectorType, typename MeshType, typename DataType>
void bcNaturalManage( VectorType& b, const MeshType& mesh, const Dof& dof, const BCBase& BCb,
                      CurrentBdFE& bdfem, const DataType& t )
{

    // Number of local Dof (i.e. nodes) in this face
    UInt nDofF = bdfem.nbNode;

    // Number of total scalar Dof
    UInt totalDof = dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCb.numberOfComponents();

    const IdentifierNatural* pId;
    ID ibF, idDof, icDof, gDof;

    if ( BCb.dataVector() )
    { //! If BC is given under a vectorial form
        switch ( BCb.pointerToBCVector() ->type() )
        {
        case 0:  // if the BC is a vector which values don't need to be integrated
            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
            {
                // Loop on components involved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {
                    ID id = BCb(i)->id();
                    // Glogal Dof
                    idDof = id + ( BCb.component( j ) - 1 ) * totalDof;

                    // Modifying right hand side (assuming BCvector is a flux)
                    b[ idDof - 1 ] += BCb( id , BCb.component( j ) );
                }
            }
            break;
        case 1:  // if the BC is a vector of values to be integrated
            // Loop on BC identifiers
            for ( ID i = 1; i <= BCb.list_size(); ++i )
            {

                // Pointer to the i-th itdentifier in the list
                pId = static_cast< const IdentifierNatural* >( BCb( i ) );

                // Number of the current boundary face
                ibF = pId->id();

                // Updating face stuff
                bdfem.updateMeasNormalQuadPt( mesh.boundaryFace( ibF ) );

                // Loop on total Dof per Face
                for ( ID l = 1; l <= nDofF; ++l )
                {

                    gDof = pId->bdLocalToGlobal( l );

                    // Loop on components involved in this boundary condition
                    for ( UInt ic = 0; ic < nComp; ++ic )
                    {
                        icDof = gDof + ic * totalDof;

                        // Loop on quadrature points
                        for ( int iq = 0; iq < bdfem.nbQuadPt; ++iq )
                        {
                            // Adding right hand side contribution
                            b[ icDof - 1 ] += BCb( gDof , 1 ) * bdfem.phi( int( l - 1 ), iq ) * bdfem.normal( int( ic ), iq ) * bdfem.weightMeas( iq );
                        }
                    }
                }
            }
            break;
    case 2:  // if the BC is a vector of values with components to be integrated
      // Loop on BC identifiers
      for ( ID i = 1; i <= BCb.list_size(); ++i )
            {

          // Pointer to the i-th itdentifier in the list
          pId = static_cast< const IdentifierNatural* >( BCb( i ) );

          // Number of the current boundary face
          ibF = pId->id();

          // Updating face stuff
          bdfem.updateMeasNormalQuadPt( mesh.boundaryFace( ibF ) );

          // Loop on total Dof per Face
          for ( ID l = 1; l <= nDofF; ++l )
                {

          gDof = pId->bdLocalToGlobal( l );

          // Loop on space dimensions condition
          for ( UInt ic = 0; ic < nDimensions; ++ic )
                    {
              // Loop on quadrature points
              for ( int iq = 0; iq < bdfem.nbQuadPt; ++iq )
                        {
              // Adding right hand side contribution
              b[ gDof - 1 ] += BCb( gDof , ic+1 ) * bdfem.phi( int( l - 1 ), iq ) * bdfem.normal( int( ic ), iq ) * bdfem.weightMeas( iq );
                        }
                    }
                }
            }
      break;
        default:
      ERROR_MSG( "This type of BCVector does not exists" );
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
        bdfem.updateMeas( mesh.boundaryFace(ibF) );

        // Loop on total Dof per Face
        for (ID idofF=1; idofF<= nDofF; ++idofF) {
    // Loop on components involved in this boundary condition
    for (ID j=1; j<=nComp; ++j) {

     //glogal Dof
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

        DataType x, y, z;

        // Loop on BC identifiers
        for ( ID i = 1; i <= BCb.list_size(); ++i )
        {

            // Pointer to the i-th itdentifier in the list
            pId = static_cast< const IdentifierNatural* >( BCb( i ) );

            // Number of the current boundary face
            ibF = pId->id();

            // Updating face stuff
            bdfem.updateMeas( mesh.boundaryFace( ibF ) );

            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            {  //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                // Loop on components involved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {

                    //glogal Dof
                    idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof;

                    // Loop on quadrature points
                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                    {

                        bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                        // Adding right hand side contribution
                        b[ idDof - 1 ] += bdfem.phi( int( idofF - 1 ), l ) * BCb( t, x, y, z, BCb.component( j ) ) *
                                          bdfem.weightMeas( l );
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
void bcMixteManageUDep( MatrixType& A, VectorType& b, const MeshType& mesh, const Dof& dof, const BCBase& BCb,
                    CurrentBdFE& bdfem, const DataType& t ,const VectorType& U)
{
  ERROR_MSG("error bcMixteManageUDep not still implemented\n");
}
template <typename MatrixType, typename VectorType, typename DataType, typename MeshType>
void bcMixteManage( MatrixType& A, VectorType& b, const MeshType& mesh, const Dof& dof, const BCBase& BCb,
                    CurrentBdFE& bdfem, const DataType& t )
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
            bdfem.updateMeas( mesh.boundaryFace( ibF ) );

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
                    idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof;

                    // Loop on quadrature points
                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                    {

                       // Compute the mixte coefficients on the quadrature point - vector of mixte
                       // coefficients is given on the element nodes !
		       mcoef = 0.0;
		       mbcb = 0.0;
		       for( ID n = 1; n <= nDofF; ++n)
		       {
			  kdDof=pId->bdLocalToGlobal( n ) + ( BCb.component( j ) - 1 ) * totalDof;
			  mcoef += BCb.MixteVec( kdDof, BCb.component( j ) ) * bdfem.phi( int( n - 1 ), l );
			  mbcb += BCb( kdDof, BCb.component( j ) ) * bdfem.phi( int( n - 1 ), l );
		       }

                        // Contribution to the diagonal entry of the elementary boundary mass matrix
//                        sum += mcoef * bdfem.phi( int( idofF - 1 ), l ) * bdfem.phi( int( idofF - 1 ), l ) *
//                               bdfem.weightMeas( l );
	                sum += mcoef* bdfem.phi( int( idofF - 1 ), l ) *
                               bdfem.phi( int( idofF - 1 ), l ) *bdfem.weightMeas( l );

                        // Adding right hand side contribution
                        //vincent please check again for your Mixte-FE it doesn't work for Q1:
                        //     b[idDof-1] += bdfem.phi(int(idofF-1),l) * BCb(BCb(i)->id(),BCb.component(j)) *
                        b[ idDof - 1 ] += bdfem.phi( int( idofF - 1 ), l ) * mbcb *
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

                        // Glogals Dof: row and columns
                        //vincent please check again for your Mixte-FE it doesn't work for Q1:
                        //     idDof  =  BCb(i)->id() + (BCb.component(j)-1)*totalDof;
                        //     jdDof  =  BCb(k)->id() + (BCb.component(j)-1)*totalDof;
                        idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof;
                        jdDof = pId->bdLocalToGlobal( k ) + ( BCb.component( j ) - 1 ) * totalDof;

                        // Loop on quadrature points
                        for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                        {

                       // Compute the mixte coefficient on the quadrature point - vector of mixte
                       // coefficients is given on the element nodes !
		       mcoef = 0.0;
		       for( ID n = 1; n <= nDofF; ++n)
		       {
			  kdDof=pId->bdLocalToGlobal( n ) + ( BCb.component( j ) - 1 ) * totalDof;
			  mcoef += BCb.MixteVec( kdDof, BCb.component( j ) ) * bdfem.phi( int( n - 1 ), l );
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
            bdfem.updateMeas( mesh.boundaryFace( ibF ) );

            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            { //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                // Loop on components invoved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {

                    sum = 0;

                    // Glogal Dof (outside the quad point loop. V. Martin)
                    idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof;

                    // Loop on quadrature points
                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                    {

                        bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                        // Contribution to the diagonal entry of the elementary boundary mass matrix
                        sum += pBcF->coef( t, x, y, z, j ) * bdfem.phi( int( idofF - 1 ), l ) * bdfem.phi( int( idofF - 1 ), l ) *
                               bdfem.weightMeas( l );

                        // Glogal Dof (Why inside this loop?? V. Martin)
                        // idDof = pId->bdLocalToGlobal(idofF) + (BCb.component(j)-1)*totalDof;

                        // Adding right hand side contribution
                        b[ idDof - 1 ] += bdfem.phi( int( idofF - 1 ), l ) * BCb( t, x, y, z, BCb.component( j ) ) *
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
                            sum += pBcF->coef( x, y, z, t, j ) * bdfem.phi( int( idofF - 1 ), l ) * bdfem.phi( int( k - 1 ), l ) *
                                   bdfem.weightMeas( l );
                        }

                        // Glogals Dof: row and columns
                        idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof;
                        jdDof = pId->bdLocalToGlobal( k ) + ( BCb.component( j ) - 1 ) * totalDof;

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
// Version that treates only the matrix modifications
// V. Martin 02/03/2003
template <typename MatrixType, typename DataType, typename MeshType>
void bcMixteManageMatrix( MatrixType& A, const MeshType& mesh, const Dof& dof,
                          const BCBase& BCb, CurrentBdFE& bdfem, const DataType& t )
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
            bdfem.updateMeas( mesh.boundaryFace( ibF ) );

            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            { //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                // Loop on components invoved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {

                    sum = 0;

                    // Glogal Dof
                    idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;

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

                        // Glogals Dof: row and columns
                        idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
                        jdDof = BCb( k ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;

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
            bdfem.updateMeas( mesh.boundaryFace( ibF ) );

            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            { //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                // Loop on components invoved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {

                    sum = 0;

                    // Glogal Dof (outside the quad point loop. V. Martin)
                    idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof;

                    // Loop on quadrature points
                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                    {

                        bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                        // Contribution to the diagonal entry of the elementary boundary mass matrix
                        sum += pBcF->coef( t, x, y, z, j ) * bdfem.phi( int( idofF - 1 ), l ) * bdfem.phi( int( idofF - 1 ), l ) *
                               bdfem.weightMeas( l );

                        // Glogal Dof (Why inside this loop?? V. Martin)
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
                            sum += pBcF->coef( x, y, z, t, j ) * bdfem.phi( int( idofF - 1 ), l ) * bdfem.phi( int( k - 1 ), l ) *
                                   bdfem.weightMeas( l );
                        }

                        // Glogals Dof: row and columns
                        idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof;
                        jdDof = pId->bdLocalToGlobal( k ) + ( BCb.component( j ) - 1 ) * totalDof;

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
void bcMixteManageVector( VectorType& b, const MeshType& mesh, const Dof& dof,
                          const BCBase& BCb, CurrentBdFE& bdfem, const DataType& t )
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
            bdfem.updateMeas( mesh.boundaryFace( ibF ) );

            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            { //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                // Loop on components invoved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {
                    // Glogal Dof
                    idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;

                    // Loop on quadrature points
                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                    {
                        // Adding right hand side contribution
                        b[ idDof - 1 ] += bdfem.phi( int( idofF - 1 ), l ) * BCb( BCb( i ) ->id(), BCb.component( j ) ) *
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
            bdfem.updateMeas( mesh.boundaryFace( ibF ) );

            // Loop on total Dof per Face
            for ( ID idofF = 1; idofF <= nDofF; ++idofF )
            { //! fixed a possible BUG(??): it was the same variable : i for list and nDofF! (V. Martin)
                //! Checked only for RT0-Q0 fe. (to be tested with Q1 or Q2...)

                // Loop on components invoved in this boundary condition
                for ( ID j = 1; j <= nComp; ++j )
                {

                    // Glogal Dof (outside the quad point loop. V. Martin)
                    idDof = pId->bdLocalToGlobal( idofF ) + ( BCb.component( j ) - 1 ) * totalDof;

                    // Loop on quadrature points
                    for ( int l = 0; l < bdfem.nbQuadPt; ++l )
                    {

                        bdfem.coorQuadPt( x, y, z, l ); // quadrature point coordinates

                        // Glogal Dof (Why inside this loop?? V. Martin)
                        // idDof = pId->bdLocalToGlobal(idofF) + (BCb.component(j)-1)*totalDof;

                        // Adding right hand side contribution
                        b[ idDof - 1 ] += bdfem.phi( int( idofF - 1 ), l ) * BCb( t, x, y, z, BCb.component( j ) ) *
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
                    CurrentBdFE& bdfem, const DataType& t )
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
        bdfem.updateMeas( mesh.boundaryFace( ibF ) );

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
                    sum += pBcF->coef( t, x, y, z, j ) * bdfem.phi( int( i - 1 ), l ) * bdfem.phi( int( i - 1 ), l ) *
                           bdfem.weightMeas( l );

                    // Glogal Dof
                    idDof = pId->bdLocalToGlobal( i ) + ( BCb.component( j ) - 1 ) * totalDof;

                    // Adding right hand side contribution
                    b[ idDof - 1 ] += bdfem.phi( int( i - 1 ), l ) * BCb( t, x, y, z, BCb.component( j ) ) * bdfem.weightMeas( l );
                }

                // Global Dof
                idDof = pId->bdLocalToGlobal( i ) + ( BCb.component( j ) - 1 ) * totalDof;

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
                        sum += pBcF->coef( x, y, z, t, j ) * bdfem.phi( int( i - 1 ), l ) * bdfem.phi( int( k - 1 ), l ) *
                               bdfem.weightMeas( l );
                    }

                    // Glogals Dof: row and columns
                    idDof = pId->bdLocalToGlobal( i ) + ( BCb.component( j ) - 1 ) * totalDof;
                    jdDof = pId->bdLocalToGlobal( k ) + ( BCb.component( j ) - 1 ) * totalDof;

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
                    const DataType& t )
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
        bdfem.updateMeas( mesh.boundaryFace( ibF ) );

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
                    sum += pBcF->coef( t, x, y, z, j ) * bdfem.phi( int( i - 1 ), l ) * bdfem.phi( int( i - 1 ), l ) *
                           bdfem.weightMeas( l );

                    // Glogal Dof
                    idDof = pId->bdLocalToGlobal( i ) + ( BCb.component( j ) - 1 ) * totalDof;

                    // Adding right hand side contribution
                    b[ idDof - 1 ] += bdfem.phi( int( i - 1 ), l ) * BCb( t, x, y, z, BCb.component( j ) ) * bdfem.weightMeas( l );
                }

                // Global Dof
                idDof = pId->bdLocalToGlobal( i ) + ( BCb.component( j ) - 1 ) * totalDof;

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
                        sum += pBcF->coef( x, y, z, t, j ) * bdfem.phi( int( i - 1 ), l ) * bdfem.phi( int( k - 1 ), l ) *
                               bdfem.weightMeas( l );
                    }

                    // Glogals Dof: row and columns
                    idDof = pId->bdLocalToGlobal( i ) + ( BCb.component( j ) - 1 ) * totalDof;
                    jdDof = pId->bdLocalToGlobal( k ) + ( BCb.component( j ) - 1 ) * totalDof;

                    // Assembling upper entry.  The boundary mas matrix is symetric
                    A.set_mat_inc( idDof - 1, jdDof - 1, sum );
                    A.set_mat_inc( jdDof - 1, idDof - 1, sum );
                }
            }
        }
    }
}
}

#endif
