//===========================================================================//
//===========================================================================//
// Obsolete not compiling Templates functions                                //
//===========================================================================//
//===========================================================================//

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
                    if ( (BCh[ i ].mode() == Tangential) || (BCh[ i ].mode() == Normal) || (BCh[ i ].mode() == Directional) )
                    {
                        ERROR_MSG( "This BC mode is not yet implemented for this setting" );
                    }
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
                    if ( (BCh[ i ].mode() == Tangential) || (BCh[ i ].mode() == Normal) || (BCh[ i ].mode() == Directional) )
                    {
                        ERROR_MSG( "This BC mode is not yet implemented for this setting" );
                    }
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
