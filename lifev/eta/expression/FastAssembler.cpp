#include <lifev/eta/expression/FastAssembler.hpp>
#include <chrono>
#include <omp.h>
using namespace std::chrono;

using namespace std;

using namespace LifeV;

//=========================================================================
// Constructor //
FastAssembler::FastAssembler( const meshPtr_Type& mesh, const commPtr_Type& comm, const ReferenceFE* refFE, const qr_Type* qr ) :
		M_mesh ( mesh ),
		M_comm ( comm ),
		M_referenceFE ( refFE ),
		M_qr ( qr )
{
}
//=========================================================================
// Destructor //
FastAssembler::~FastAssembler()
{
	//---------------------

	delete[] M_detJacobian;

	//---------------------

	for( int i = 0 ; i < M_numElements ; i++ )
	{
		for ( int j = 0; j < 3; j++ )
		{
			delete [] M_invJacobian[i][j];
		}
		delete [] M_invJacobian[i];
	}
	delete [] M_invJacobian;

	//---------------------

	for ( int i = 0; i < M_referenceFE->nbDof(); i++ )
	{
		delete [] M_phi[i];
	}

	delete [] M_phi;

	//---------------------

	for( int i = 0 ; i < M_referenceFE->nbDof() ; i++ )
	{
		for ( int j = 0; j < M_qr->nbQuadPt(); j++ )
		{
			delete [] M_dphi[i][j];
		}
		delete [] M_dphi[i];
	}
	delete [] M_dphi;

	//---------------------

	for( int i = 0 ; i < M_numElements ; i++ )
	{
		for ( int j = 0; j < M_referenceFE->nbDof(); j++ )
		{
			delete [] M_vals[i][j];
		}
		delete [] M_vals[i];
	}
	delete [] M_vals;


	for( int i = 0 ; i < M_numElements ; i++ )
	{
		delete [] M_rows[i];
		delete [] M_cols[i];
	}
	delete [] M_rows;
	delete [] M_cols;
}
//=========================================================================
void
FastAssembler::allocateSpace ( const int& numElements, CurrentFE* fe, const fespacePtr_Type& fespace )
{
	M_numElements = numElements;
	M_numScalarDofs = fespace->dof().numTotalDof();

	M_detJacobian = new double[ M_numElements ];

	M_invJacobian = new double** [ M_numElements ];

	//-------------------------------------------------------------------------------------------------

	M_invJacobian = new double** [ M_numElements];

	for ( int i = 0; i <  M_numElements; i++ )
	{
		M_invJacobian[i] = new double* [ 3 ];
		for ( int j = 0; j < 3 ; j++ )
		{
			M_invJacobian[i][j] = new double [ 3 ];
		}
	}

	for ( int i = 0; i < M_numElements; i++ )
	{
		fe->update( M_mesh->element (i), UPDATE_DPHI );
		M_detJacobian[i] = fe->detJacobian(0);
		for ( int j = 0; j < 3; j++ )
		{
			for ( int k = 0; k < 3; k++ )
			{
				M_invJacobian[i][j][k] = fe->tInverseJacobian(j,k,2);
			}
		}
	}

	//-------------------------------------------------------------------------------------------------

	M_phi = new double* [ M_referenceFE->nbDof() ];
	for ( int i = 0; i < M_referenceFE->nbDof(); i++ )
	{
		M_phi[i] = new double [ M_qr->nbQuadPt() ];
	}

    // PHI REF
    for (UInt q (0); q < M_qr->nbQuadPt(); ++q)
    {
        for (UInt j (0); j < M_referenceFE->nbDof(); ++j)
        {
        	M_phi[j][q] =  M_referenceFE->phi (j, M_qr->quadPointCoor (q) );
        }
    }

    //-------------------------------------------------------------------------------------------------

    M_dphi = new double** [ M_referenceFE->nbDof() ];

    for ( int i = 0; i <  M_referenceFE->nbDof(); i++ )
    {
    	M_dphi[i] = new double* [ M_qr->nbQuadPt() ];
    	for ( int j = 0; j < M_qr->nbQuadPt() ; j++ )
    	{
    		M_dphi[i][j] = new double [ 3 ];
    	}
    }

    //DPHI REF
    for (UInt q (0); q < M_qr->nbQuadPt(); ++q)
    {
    	for (UInt i (0); i < M_referenceFE->nbDof(); ++i)
    	{
    		for (UInt j (0); j < 3; ++j)
    		{
    			M_dphi[i][q][j] = M_referenceFE->dPhi (i, j, M_qr->quadPointCoor (q) );
    		}
    	}
    }

    //-------------------------------------------------------------------------------------------------

    M_elements = new double* [ M_numElements ];

    for ( int i = 0; i <  M_numElements; i++ )
    {
        M_elements[i] = new double [ M_referenceFE->nbDof() ];
    }

    for ( int i = 0; i <  M_numElements; i++ )
    {
        for ( int j = 0; j <  M_referenceFE->nbDof(); j++ )
        {
            M_elements[i][j] =     fespace->dof().localToGlobalMap (i, j);
        }
    }

    //-------------------------------------------------------------------------------------------------

    // Allocate space for M_rows, M_cols and M_vals

    M_vals = new double** [M_numElements];

    for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
    {
    	M_vals[i_elem] = new double* [ M_referenceFE->nbDof() ];
    	for ( int i = 0; i <  M_referenceFE->nbDof(); i++ )
    	{
    		M_vals[i_elem][i] = new double [ M_referenceFE->nbDof() ];
    	}
    }

    M_rows = new int* [M_numElements];

    for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
    {
    	M_rows[i_elem] = new int [ M_referenceFE->nbDof() ];
    }

    M_cols = new int* [M_numElements];

    for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
    {
    	M_cols[i_elem] = new int [ M_referenceFE->nbDof() ];
    }

    int M_nbRow = M_referenceFE->nbDof();
    int M_nbColumn = M_referenceFE->nbDof();
}
//=========================================================================
void
FastAssembler::allocateSpace( const int& numElements, CurrentFE* fe, const fespacePtr_Type& fespace, const UInt* meshSub_elements )
{
	M_numElements = numElements;
	M_numScalarDofs = fespace->dof().numTotalDof();

	M_detJacobian = new double[ M_numElements ];

	M_invJacobian = new double** [ M_numElements ];

	//-------------------------------------------------------------------------------------------------

	M_invJacobian = new double** [ M_numElements];

	for ( int i = 0; i <  M_numElements; i++ )
	{
		M_invJacobian[i] = new double* [ 3 ];
		for ( int j = 0; j < 3 ; j++ )
		{
			M_invJacobian[i][j] = new double [ 3 ];
		}
	}

	for ( int i = 0; i < M_numElements; i++ )
	{
		fe->update( M_mesh->element (meshSub_elements[i]), UPDATE_DPHI );
		M_detJacobian[i] = fe->detJacobian(0);
		for ( int j = 0; j < 3; j++ )
		{
			for ( int k = 0; k < 3; k++ )
			{
				M_invJacobian[i][j][k] = fe->tInverseJacobian(j,k,0);
			}
		}
	}

	//-------------------------------------------------------------------------------------------------

	M_phi = new double* [ M_referenceFE->nbDof() ];
	for ( int i = 0; i < M_referenceFE->nbDof(); i++ )
	{
		M_phi[i] = new double [ M_qr->nbQuadPt() ];
	}

    // PHI REF
    for (UInt q (0); q < M_qr->nbQuadPt(); ++q)
    {
        for (UInt j (0); j < M_referenceFE->nbDof(); ++j)
        {
        	M_phi[j][q] =  M_referenceFE->phi (j, M_qr->quadPointCoor (q) );
        }
    }

    //-------------------------------------------------------------------------------------------------

    M_dphi = new double** [ M_referenceFE->nbDof() ];

    for ( int i = 0; i <  M_referenceFE->nbDof(); i++ )
    {
    	M_dphi[i] = new double* [ M_qr->nbQuadPt() ];
    	for ( int j = 0; j < M_qr->nbQuadPt() ; j++ )
    	{
    		M_dphi[i][j] = new double [ 3 ];
    	}
    }

    //DPHI REF
    for (UInt q (0); q < M_qr->nbQuadPt(); ++q)
    {
    	for (UInt i (0); i < M_referenceFE->nbDof(); ++i)
    	{
    		for (UInt j (0); j < 3; ++j)
    		{
    			M_dphi[i][q][j] = M_referenceFE->dPhi (i, j, M_qr->quadPointCoor (q) );
    		}
    	}
    }

    //-------------------------------------------------------------------------------------------------

    M_elements = new double* [ M_numElements ];

    for ( int i = 0; i <  M_numElements; i++ )
    {
        M_elements[i] = new double [ M_referenceFE->nbDof() ];
    }

    for ( int i = 0; i <  M_numElements; i++ )
    {
        for ( int j = 0; j <  M_referenceFE->nbDof(); j++ )
        {
            M_elements[i][j] =     fespace->dof().localToGlobalMap (meshSub_elements[i], j);
        }
    }

}

//=========================================================================
void
FastAssembler::assembleGradGrad_scalar( matrixPtr_Type& matrix )
{
    int ndof = M_referenceFE->nbDof();
    int NumQuadPoints = M_qr->nbQuadPt();

    double w_quad[NumQuadPoints];
    for ( int q = 0; q < NumQuadPoints ; q++ )
    {
    	w_quad[q] = M_qr->weight(q);
    }

    double dtime = omp_get_wtime();
    double dtime2 = dtime;

    //#pragma omp parallel shared(M_rows,M_cols,M_vals) firstprivate( w_quad, ndof, NumQuadPoints)
    #pragma omp parallel firstprivate( w_quad, ndof, NumQuadPoints)
    {
        int i_el, i_elem, i_dof, q, d1, d2, i_test, i_trial;
        double integral;

        double dphi_phys[ndof][NumQuadPoints][3];

        // ELEMENTI
        #pragma omp for
        for ( i_el = 0; i_el <  M_numElements ; i_el++ )
        {
        	i_elem = i_el;

            // DOF
            for ( i_dof = 0; i_dof <  ndof; i_dof++ )
            {
                // QUAD
                for (  q = 0; q < NumQuadPoints ; q++ )
                {
                    // DIM 1
                    for ( d1 = 0; d1 < 3 ; d1++ )
                    {
                        dphi_phys[i_dof][q][d1] = 0.0;

                        // DIM 2
                        for ( d2 = 0; d2 < 3 ; d2++ )
                        {
                            dphi_phys[i_dof][q][d1] += M_invJacobian[i_elem][d1][d2] * M_dphi[i_dof][q][d2];
                        }
                    }
                }
            }

            // DOF - test
            for ( i_test = 0; i_test <  ndof; i_test++ )
            {
                M_rows[i_elem][i_test] = M_elements[i_elem][i_test];

                // DOF - trial
                for ( i_trial = 0; i_trial <  ndof; i_trial++ )
                {
                    M_cols[i_elem][i_trial] = M_elements[i_elem][i_trial];

                    integral = 0.0;
                    // QUAD
                    for ( q = 0; q < NumQuadPoints ; q++ )
                    {
                        // DIM 1
                        for ( d1 = 0; d1 < 3 ; d1++ )
                        {
                            integral += dphi_phys[i_test][q][d1] * dphi_phys[i_trial][q][d1]*w_quad[q];
                        }
                    }
                    M_vals[i_elem][i_test][i_trial] = integral * M_detJacobian[i_elem];
                }
            }
//			#pragma omp critical
//            matrix->matrixPtr()->InsertGlobalValues ( ndof, M_rows[i_elem], ndof, M_cols[i_elem], M_vals[i_elem], Epetra_FECrsMatrix::ROW_MAJOR);
        }
    }

    dtime = omp_get_wtime() - dtime;
    printf("\n\nelapsed time partial: %f\n\n", dtime);

    for ( int k = 0; k < M_numElements; ++k )
    {
    	matrix->matrixPtr()->InsertGlobalValues ( ndof, M_rows[k], ndof, M_cols[k], M_vals[k], Epetra_FECrsMatrix::ROW_MAJOR);
    }

    dtime2 = omp_get_wtime() - dtime2;
    printf("\n\nelapsed time total: %f\n\n", dtime2);

}
//=========================================================================
void
FastAssembler::assembleGradGrad_vectorial( matrixPtr_Type& matrix )
{
    int ndof = M_referenceFE->nbDof();
    int NumQuadPoints = M_qr->nbQuadPt();

    double w_quad[NumQuadPoints];
    for ( int q = 0; q < NumQuadPoints ; q++ )
    {
    	w_quad[q] = M_qr->weight(q);
    }

    #pragma omp parallel firstprivate( w_quad, ndof, NumQuadPoints)
    {
        int i_el, i_elem, i_dof, q, d1, d2, i_test, i_trial;
        double integral;

        double dphi_phys[ndof][NumQuadPoints][3];

        // ELEMENTI
        #pragma omp for
        for ( i_el = 0; i_el <  M_numElements ; i_el++ )
        {
        	i_elem = i_el;

            // DOF
            for ( i_dof = 0; i_dof <  ndof; i_dof++ )
            {
                // QUAD
                for (  q = 0; q < NumQuadPoints ; q++ )
                {
                    // DIM 1
                    for ( d1 = 0; d1 < 3 ; d1++ )
                    {
                        dphi_phys[i_dof][q][d1] = 0.0;

                        // DIM 2
                        for ( d2 = 0; d2 < 3 ; d2++ )
                        {
                            dphi_phys[i_dof][q][d1] += M_invJacobian[i_elem][d1][d2] * M_dphi[i_dof][q][d2];
                        }
                    }
                }
            }

            // DOF - test
            for ( i_test = 0; i_test <  ndof; i_test++ )
            {
                M_rows[i_elem][i_test] = M_elements[i_elem][i_test];

                // DOF - trial
                for ( i_trial = 0; i_trial <  ndof; i_trial++ )
                {
                    M_cols[i_elem][i_trial] = M_elements[i_elem][i_trial];

                    integral = 0.0;
                    // QUAD
                    for ( q = 0; q < NumQuadPoints ; q++ )
                    {
                        // DIM 1
                        for ( d1 = 0; d1 < 3 ; d1++ )
                        {
                            integral += dphi_phys[i_test][q][d1] * dphi_phys[i_trial][q][d1]*w_quad[q];
                        }
                    }
                    M_vals[i_elem][i_test][i_trial] = integral * M_detJacobian[i_elem];
                }
            }
        }
    }

    for ( int k = 0; k < M_numElements; ++k )
    {
    	matrix->matrixPtr()->InsertGlobalValues ( ndof, M_rows[k], ndof, M_cols[k], M_vals[k], Epetra_FECrsMatrix::ROW_MAJOR);
    	for ( UInt d1 = 1; d1 < 3 ; d1++ )
    	{
    		for ( UInt i = 0; i <  ndof; i++ )
    		{
    			M_rows[k][i] += M_numScalarDofs;
    			M_cols[k][i] += M_numScalarDofs;
    		}
    		matrix->matrixPtr()->InsertGlobalValues ( ndof, M_rows[k], ndof, M_cols[k], M_vals[k], Epetra_FECrsMatrix::ROW_MAJOR);
    	}
    }

}
//=========================================================================
void
FastAssembler::assembleMass_vectorial( matrixPtr_Type& matrix )
{
    int ndof = M_referenceFE->nbDof();
    int NumQuadPoints = M_qr->nbQuadPt();

    double w_quad[NumQuadPoints];
    for ( int q = 0; q < NumQuadPoints ; q++ )
    {
    	w_quad[q] = M_qr->weight(q);
    }

    #pragma omp parallel firstprivate( w_quad, ndof, NumQuadPoints)
    {
        int i_el, i_elem, i_dof, q, d1, d2, i_test, i_trial;
        double integral;

        double dphi_phys[ndof][NumQuadPoints][3];

        // ELEMENTI
        #pragma omp for
        for ( i_el = 0; i_el <  M_numElements ; i_el++ )
        {
        	i_elem = i_el;

            // DOF - test
            for ( i_test = 0; i_test <  ndof; i_test++ )
            {
                M_rows[i_elem][i_test] = M_elements[i_elem][i_test];

                // DOF - trial
                for ( i_trial = 0; i_trial <  ndof; i_trial++ )
                {
                    M_cols[i_elem][i_trial] = M_elements[i_elem][i_trial];

                    integral = 0.0;
                    // QUAD
                    for ( q = 0; q < NumQuadPoints ; q++ )
                    {
                    	integral += M_phi[i_test][q] * M_phi[i_trial][q]*w_quad[q];
                    }
                    M_vals[i_elem][i_test][i_trial] = integral * M_detJacobian[i_elem];
                }
            }
        }
    }

    for ( int k = 0; k < M_numElements; ++k )
    {
    	matrix->matrixPtr()->InsertGlobalValues ( ndof, M_rows[k], ndof, M_cols[k], M_vals[k], Epetra_FECrsMatrix::ROW_MAJOR);
    	for ( UInt d1 = 1; d1 < 3 ; d1++ )
    	{
    		for ( UInt i = 0; i <  ndof; i++ )
    		{
    			M_rows[k][i] += M_numScalarDofs;
    			M_cols[k][i] += M_numScalarDofs;
    		}
    		matrix->matrixPtr()->InsertGlobalValues ( ndof, M_rows[k], ndof, M_cols[k], M_vals[k], Epetra_FECrsMatrix::ROW_MAJOR);
    	}
    }

}
//=========================================================================
void
FastAssembler::assembleMass_scalar( matrixPtr_Type& matrix )
{
    int ndof = M_referenceFE->nbDof();
    int NumQuadPoints = M_qr->nbQuadPt();

    double w_quad[NumQuadPoints];
    for ( int q = 0; q < NumQuadPoints ; q++ )
    {
    	w_quad[q] = M_qr->weight(q);
    }

    #pragma omp parallel firstprivate( w_quad, ndof, NumQuadPoints)
    {
        int i_el, i_elem, i_dof, q, d1, d2, i_test, i_trial;
        double integral;

        double dphi_phys[ndof][NumQuadPoints][3];

        // ELEMENTI
        #pragma omp for
        for ( i_el = 0; i_el <  M_numElements ; i_el++ )
        {
        	i_elem = i_el;

            // DOF - test
            for ( i_test = 0; i_test <  ndof; i_test++ )
            {
                M_rows[i_elem][i_test] = M_elements[i_elem][i_test];

                // DOF - trial
                for ( i_trial = 0; i_trial <  ndof; i_trial++ )
                {
                    M_cols[i_elem][i_trial] = M_elements[i_elem][i_trial];

                    integral = 0.0;
                    // QUAD
                    for ( q = 0; q < NumQuadPoints ; q++ )
                    {
                    	integral += M_phi[i_test][q] * M_phi[i_trial][q]*w_quad[q];
                    }
                    M_vals[i_elem][i_test][i_trial] = integral * M_detJacobian[i_elem];
                }
            }
        }
    }

    for ( int k = 0; k < M_numElements; ++k )
    {
    	matrix->matrixPtr()->InsertGlobalValues ( ndof, M_rows[k], ndof, M_cols[k], M_vals[k], Epetra_FECrsMatrix::ROW_MAJOR);
    }

}
//=========================================================================
void
FastAssembler::assembleConvective( matrixPtr_Type& matrix, const vector_Type& u_h )
{
	assembleConvective( *matrix,  u_h );
}
//=========================================================================
void
FastAssembler::assembleConvective( matrix_Type& matrix, const vector_Type& u_h )
{
    double ** local_matrix;

    int ndof = M_referenceFE->nbDof();
    int NumQuadPoints = M_qr->nbQuadPt();
    int ndof_vec = M_referenceFE->nbDof()*3;

    double*** M_vals = new double** [M_numElements];

    for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
    {
        M_vals[i_elem] = new double* [ ndof ];
        for ( int i = 0; i <  ndof; i++ )
        {
            M_vals[i_elem][i] = new double [ ndof ];
        }
    }

    int** M_rows = new int* [M_numElements];

    for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
    {
    	M_rows[i_elem] = new int [ ndof ];
    }

    int** M_cols = new int* [M_numElements];

    for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
    {
    	M_cols[i_elem] = new int [ ndof ];
    }

    double w_quad[NumQuadPoints];
    for ( int q = 0; q < NumQuadPoints ; q++ )
	{
        w_quad[q] = M_qr->weight(q);
    }

    //double dtime = omp_get_wtime();
    //double dtime2 = dtime;

    //#pragma omp parallel shared(M_rows,M_cols,M_vals) firstprivate( w_quad, ndof, NumQuadPoints)
    #pragma omp parallel firstprivate( w_quad, ndof, NumQuadPoints)
    {
        int i_elem, i_dof, q, d1, d2, i_test, i_trial, e_idof;;
        double integral;

        double dphi_phys[ndof][NumQuadPoints][3];

        double uhq[3][NumQuadPoints];

        // ELEMENTI,
        #pragma omp for
        for ( i_elem = 0; i_elem <  M_numElements; i_elem++ )
        {
            // DOF
            for ( i_dof = 0; i_dof <  ndof; i_dof++ )
            {
                // QUAD
                for (  q = 0; q < NumQuadPoints ; q++ )
                {
                    // DIM 1
                    for ( d1 = 0; d1 < 3 ; d1++ )
                    {
                        dphi_phys[i_dof][q][d1] = 0.0;

                        // DIM 2
                        for ( d2 = 0; d2 < 3 ; d2++ )
                        {
                            dphi_phys[i_dof][q][d1] += M_invJacobian[i_elem][d1][d2] * M_dphi[i_dof][q][d2];
                        }
                    }
                }
            }

            // QUAD
            for (  q = 0; q < NumQuadPoints ; q++ )
            {
            	for ( d1 = 0; d1 < 3 ; d1++ )
            	{
            		uhq[d1][q] = 0.0;
            		for ( i_dof = 0; i_dof <  ndof; i_dof++ )
            		{
            			e_idof =  M_elements[i_elem][i_dof] + d1*M_numScalarDofs  ;
            			uhq[d1][q] += u_h[e_idof] * M_phi[i_dof][q];
            			//printf("\n u_h[%d] = %f",  e_idof, u_h[e_idof]);
            		}
            	}
            }

            // DOF - test
            for ( i_test = 0; i_test <  ndof; i_test++ )
            {
                M_rows[i_elem][i_test] = M_elements[i_elem][i_test];

                // DOF - trial
                for ( i_trial = 0; i_trial <  ndof; i_trial++ )
                {
                    M_cols[i_elem][i_trial] = M_elements[i_elem][i_trial];

                    integral = 0.0;
                    // QUAD
                    for ( q = 0; q < NumQuadPoints ; q++ )
                    {
                        // DIM 1
                        for ( d1 = 0; d1 < 3 ; d1++ )
                        {
                            integral += uhq[d1][q] * dphi_phys[i_trial][q][d1] * M_phi[i_test][q] * w_quad[q];
                        }
                    }
                    M_vals[i_elem][i_test][i_trial] = integral *  M_detJacobian[i_elem];
                }
            }
        }
    }

    //dtime = omp_get_wtime() - dtime;
    //printf("\n\nelapsed time partial: %f\n\n", dtime);

    for ( int k = 0; k < M_numElements; ++k )
    {
    	matrix.matrixPtr()->InsertGlobalValues ( ndof, M_rows[k], ndof, M_cols[k], M_vals[k], Epetra_FECrsMatrix::ROW_MAJOR);
    	for ( int d1 = 1; d1 < 3 ; d1++ )
    	{
    		for ( int i = 0; i <  ndof; i++ )
    		{
                M_rows[k][i] += M_numScalarDofs;
                M_cols[k][i] += M_numScalarDofs;
    		}
        	matrix.matrixPtr()->InsertGlobalValues ( ndof, M_rows[k], ndof, M_cols[k], M_vals[k], Epetra_FECrsMatrix::ROW_MAJOR);
    	}
    }

    //dtime2 = omp_get_wtime() - dtime2;
    //printf("\n\nelapsed time total: %f\n\n", dtime2);

    // cleaning memory
    for( int i = 0 ; i < M_numElements ; i++ )
    {
        for ( int j = 0; j < ndof; j++ )
        {
            delete [] M_vals[i][j];
        }
        delete [] M_vals[i];
    }
    delete [] M_vals;


    for( int i = 0 ; i < M_numElements ; i++ )
    {
    	delete [] M_rows[i];
    	delete [] M_cols[i];
    }
    delete [] M_rows;
    delete [] M_cols;

}
//=========================================================================
