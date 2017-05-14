#include <lifev/core/fem/FastAssemblerMixed.hpp>
#include <chrono>
#ifdef EPETRA_HAVE_OMP
#include <omp.h>
#endif
using namespace std::chrono;

using namespace std;

using namespace LifeV;

//=========================================================================
// Constructor //
FastAssemblerMixed::FastAssemblerMixed( const meshPtr_Type& mesh, const commPtr_Type& comm,
							            const ReferenceFE* refFE_test, const ReferenceFE* refFE_trial,
							            const qr_Type* qr_integration ) :
		M_mesh ( mesh ),
		M_comm ( comm ),
		M_referenceFE_test ( refFE_test ),
		M_referenceFE_trial ( refFE_trial ),
		M_qr_integration ( qr_integration )
{
}
//=========================================================================
// Destructor //
FastAssemblerMixed::~FastAssemblerMixed()
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

	for ( int i = 0; i < M_referenceFE_test->nbDof(); i++ )
	{
		delete [] M_phi_test[i];
	}

	delete [] M_phi_test;

	//---------------------

	for ( int i = 0; i < M_referenceFE_trial->nbDof(); i++ )
	{
		delete [] M_phi_trial[i];
	}

	delete [] M_phi_trial;

	//---------------------

	for( int i = 0 ; i < M_referenceFE_test->nbDof() ; i++ )
	{
		for ( int j = 0; j < M_qr_integration->nbQuadPt(); j++ )
		{
			delete [] M_dphi_test[i][j];
		}
		delete [] M_dphi_test[i];
	}
	delete [] M_dphi_test;

	//---------------------

	for( int i = 0 ; i < M_referenceFE_trial->nbDof() ; i++ )
	{
		for ( int j = 0; j < M_qr_integration->nbQuadPt(); j++ )
		{
			delete [] M_dphi_trial[i][j];
		}
		delete [] M_dphi_trial[i];
	}
	delete [] M_dphi_trial;

	//---------------------

	for( int i = 0 ; i < M_numElements ; i++ )
	{
		delete [] M_elements_test[i];
	}
	delete [] M_elements_test;

	//---------------------

	for( int i = 0 ; i < M_numElements ; i++ )
	{
		delete [] M_elements_trial[i];
	}
	delete [] M_elements_trial;

	//---------------------

	for( int i = 0 ; i < M_numElements ; i++ )
	{
		for ( int k = 0; k < 3; k++ )
		{
			for ( int j = 0; j <  M_referenceFE_test->nbDof(); j++ )
			{
				delete [] M_vals[i][k][j];
			}
			delete [] M_vals[i][k];
		}
		delete [] M_rows[i];
		delete [] M_cols[i];
		delete [] M_vals[i];
	}
	delete [] M_rows;
	delete [] M_cols;
	delete [] M_vals;
}
//=========================================================================
void
FastAssemblerMixed::allocateSpace ( const int& numElements,
		                            CurrentFE* fe_test, const fespacePtr_Type& fespace_test,
		                            CurrentFE* fe_trial, const fespacePtr_Type& fespace_trial )
{
	M_numElements = numElements;

	M_numScalarDofs_test = fespace_test->dof().numTotalDof();

	M_numScalarDofs_trial = fespace_trial->dof().numTotalDof();

	//-------------------------------------------------------------------------------------------------

	M_detJacobian = new double[ M_numElements ];

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
		fe_test->update( M_mesh->element (i), UPDATE_DPHI );
		M_detJacobian[i] = fe_test->detJacobian(0);
		for ( int j = 0; j < 3; j++ )
		{
			for ( int k = 0; k < 3; k++ )
			{
				M_invJacobian[i][j][k] = fe_test->tInverseJacobian(j,k,2);
			}
		}
	}

	//-------------------------------------------------------------------------------------------------
	// TEST FUNCTIONS - phi_i and d_phi_i
	//-------------------------------------------------------------------------------------------------

	M_phi_test = new double* [ M_referenceFE_test->nbDof() ];
	for ( int i = 0; i < M_referenceFE_test->nbDof(); i++ )
	{
		M_phi_test[i] = new double [ M_qr_integration->nbQuadPt() ];
	}

	// PHI REF
	for (UInt q (0); q < M_qr_integration->nbQuadPt(); ++q)
	{
		for (UInt j (0); j < M_referenceFE_test->nbDof(); ++j)
		{
			M_phi_test[j][q] =  M_referenceFE_test->phi (j, M_qr_integration->quadPointCoor (q) );
		}
	}

	//-------------------------------------------------------------------------------------------------

	M_dphi_test = new double** [ M_referenceFE_test->nbDof() ];

	for ( int i = 0; i <  M_referenceFE_test->nbDof(); i++ )
	{
		M_dphi_test[i] = new double* [ M_qr_integration->nbQuadPt() ];
		for ( int j = 0; j < M_qr_integration->nbQuadPt() ; j++ )
		{
			M_dphi_test[i][j] = new double [ 3 ];
		}
	}

	//DPHI REF
	for (UInt q (0); q < M_qr_integration->nbQuadPt(); ++q)
	{
		for (UInt i (0); i < M_referenceFE_test->nbDof(); ++i)
		{
			for (UInt j (0); j < 3; ++j)
			{
				M_dphi_test[i][q][j] = M_referenceFE_test->dPhi (i, j, M_qr_integration->quadPointCoor (q) );
			}
		}
	}

	//-------------------------------------------------------------------------------------------------
	// TRIAL FUNCTIONS - phi_j and d_phi_j
	//-------------------------------------------------------------------------------------------------

	M_phi_trial = new double* [ M_referenceFE_trial->nbDof() ];
	for ( int i = 0; i < M_referenceFE_trial->nbDof(); i++ )
	{
		M_phi_trial[i] = new double [ M_qr_integration->nbQuadPt() ];
	}

	// PHI REF
	for (UInt q (0); q < M_qr_integration->nbQuadPt(); ++q)
	{
		for (UInt j (0); j < M_referenceFE_trial->nbDof(); ++j)
		{
			M_phi_trial[j][q] =  M_referenceFE_trial->phi (j, M_qr_integration->quadPointCoor (q) );
		}
	}

	//-------------------------------------------------------------------------------------------------

	M_dphi_trial = new double** [ M_referenceFE_trial->nbDof() ];

	for ( int i = 0; i <  M_referenceFE_trial->nbDof(); i++ )
	{
		M_dphi_trial[i] = new double* [ M_qr_integration->nbQuadPt() ];
		for ( int j = 0; j < M_qr_integration->nbQuadPt() ; j++ )
		{
			M_dphi_trial[i][j] = new double [ 3 ];
		}
	}

	//DPHI REF
	for (UInt q (0); q < M_qr_integration->nbQuadPt(); ++q)
	{
		for (UInt i (0); i < M_referenceFE_trial->nbDof(); ++i)
		{
			for (UInt j (0); j < 3; ++j)
			{
				M_dphi_trial[i][q][j] = M_referenceFE_trial->dPhi (i, j, M_qr_integration->quadPointCoor (q) );
			}
		}
	}

	//-------------------------------------------------------------------------------------------------
	// CONNECTIVITY TEST
	//-------------------------------------------------------------------------------------------------

	M_elements_test = new double* [ M_numElements ];

	for ( int i = 0; i <  M_numElements; i++ )
	{
		M_elements_test[i] = new double [ M_referenceFE_test->nbDof() ];
	}

	for ( int i = 0; i <  M_numElements; i++ )
	{
		for ( int j = 0; j <  M_referenceFE_test->nbDof(); j++ )
		{
			M_elements_test[i][j] = fespace_test->dof().localToGlobalMap (i, j);
		}
	}

	//-------------------------------------------------------------------------------------------------
	// CONNECTIVITY TRIAL
	//-------------------------------------------------------------------------------------------------

	M_elements_trial = new double* [ M_numElements ];

	for ( int i = 0; i <  M_numElements; i++ )
	{
		M_elements_trial[i] = new double [ M_referenceFE_trial->nbDof() ];
	}

	for ( int i = 0; i <  M_numElements; i++ )
	{
		for ( int j = 0; j <  M_referenceFE_trial->nbDof(); j++ )
		{
			M_elements_trial[i][j] = fespace_trial->dof().localToGlobalMap (i, j);
		}
	}

	//-------------------------------------------------------------------------------------------------
	// Allocate space for M_rows, M_cols and M_vals
	//-------------------------------------------------------------------------------------------------

	M_vals = new double*** [M_numElements];

	for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
	{
		M_vals[i_elem] = new double** [ 3 ];
		for ( int k = 0; k < 3; k++ )
		{
			M_vals[i_elem][k] = new double* [ M_referenceFE_test->nbDof() ];
			for ( int i = 0; i <  M_referenceFE_test->nbDof(); i++ )
			{
				M_vals[i_elem][k][i] = new double [ M_referenceFE_trial->nbDof() ];
			}
		}
	}

	M_rows = new int* [M_numElements];

	for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
	{
		M_rows[i_elem] = new int [ M_referenceFE_test->nbDof() ];
	}

	M_cols = new int* [M_numElements];

	for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
	{
		M_cols[i_elem] = new int [ M_referenceFE_trial->nbDof() ];
	}
}
//=====================================================================================================
void
FastAssemblerMixed::assemble_NS_block01 ( matrixPtr_Type& matrix )
{

	// In this case we need to assemble: -1.0 * div(phi_i)*phi_j

	int ndof_test = M_referenceFE_test->nbDof();
	int ndof_trial = M_referenceFE_trial->nbDof();
	int NumQuadPoints = M_qr_integration->nbQuadPt(); // WARNING!!

	double w_quad[NumQuadPoints]; // Note: suppose that test has more precise qr than trial
	for ( int q = 0; q < NumQuadPoints ; q++ )
	{
		w_quad[q] = M_qr_integration->weight(q);
	}

	#pragma omp parallel firstprivate ( w_quad, ndof_test, ndof_trial, NumQuadPoints )
	{
		int i_el, i_elem, i_dof, q, d1, d2, i_test, i_trial, dim_mat;
		double integral;

		double dphi_phys[ndof_test][NumQuadPoints][3]; // Gradient in the physical domain

		// ELEMENTI
		#pragma omp for
		for ( i_el = 0; i_el <  M_numElements ; i_el++ )
		{
			i_elem = i_el;

			// DOF
			for ( i_dof = 0; i_dof <  ndof_test; i_dof++ )
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
							dphi_phys[i_dof][q][d1] += M_invJacobian[i_elem][d1][d2] * M_dphi_test[i_dof][q][d2];
						}
					}
				}
			}

			for ( dim_mat = 0; dim_mat < 3 ; dim_mat++ )
			{
				// DOF - test
				for ( i_test = 0; i_test <  ndof_test; i_test++ )
				{
					M_rows[i_elem][i_test] = M_elements_test[i_elem][i_test];

					// DOF - trial
					for ( i_trial = 0; i_trial < ndof_trial; i_trial++ )
					{
						M_cols[i_elem][i_trial] = M_elements_trial[i_elem][i_trial];

						integral = 0.0;
						// QUAD
						for ( q = 0; q < NumQuadPoints ; q++ )
						{
							integral += dphi_phys[i_test][q][dim_mat]*M_phi_trial[i_trial][q]*w_quad[q];
						}

						M_vals[i_elem][dim_mat][i_test][i_trial] = -1.0 * integral * M_detJacobian[i_elem];
					}
				}
			}
		}
	}

	for ( int k = 0; k < M_numElements; ++k )
	{
		matrix->matrixPtr()->InsertGlobalValues ( ndof_test, M_rows[k], ndof_trial, M_cols[k], M_vals[k][0], Epetra_FECrsMatrix::ROW_MAJOR);
	}

	for ( UInt d1 = 1; d1 < 3 ; d1++ )
	{
		for ( int k = 0; k < M_numElements; ++k )
		{
			for ( UInt i = 0; i <  ndof_test; i++ )
			{
				M_rows[k][i] += M_numScalarDofs_test;
			}
			matrix->matrixPtr()->InsertGlobalValues ( ndof_test, M_rows[k], ndof_trial, M_cols[k], M_vals[k][d1], Epetra_FECrsMatrix::ROW_MAJOR);
		}
	}
}
//=====================================================================================================
void
FastAssemblerMixed::assemble_NS_block10 ( matrixPtr_Type& matrix )
{
	// In this case we need to assemble: phi_i * div(phi_j)

	int ndof_test = M_referenceFE_test->nbDof();
	int ndof_trial = M_referenceFE_trial->nbDof();
	int NumQuadPoints = M_qr_integration->nbQuadPt();

	double w_quad[NumQuadPoints];
	for ( int q = 0; q < NumQuadPoints ; q++ )
	{
		w_quad[q] = M_qr_integration->weight(q);
	}

	#pragma omp parallel firstprivate ( w_quad, ndof_test, ndof_trial, NumQuadPoints )
	{
		int i_el, i_elem, i_dof, q, d1, d2, i_test, i_trial, dim_mat;
		double integral;

		double dphi_phys[ndof_trial][NumQuadPoints][3]; // Gradient in the physical domain

		// ELEMENTI
		#pragma omp for
		for ( i_el = 0; i_el <  M_numElements ; i_el++ )
		{
			i_elem = i_el;

			// DOF
			for ( i_dof = 0; i_dof <  ndof_trial; i_dof++ )
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
							dphi_phys[i_dof][q][d1] += M_invJacobian[i_elem][d1][d2] * M_dphi_trial[i_dof][q][d2];
						}
					}
				}
			}

			for ( dim_mat = 0; dim_mat < 3 ; dim_mat++ )
			{
				// DOF - test
				for ( i_test = 0; i_test <  ndof_test; i_test++ )
				{
					M_rows[i_elem][i_test] = M_elements_test[i_elem][i_test];

					// DOF - trial
					for ( i_trial = 0; i_trial < ndof_trial; i_trial++ )
					{
						M_cols[i_elem][i_trial] = M_elements_trial[i_elem][i_trial];

						integral = 0.0;
						// QUAD
						for ( q = 0; q < NumQuadPoints ; q++ )
						{
							integral += dphi_phys[i_trial][q][dim_mat]*M_phi_test[i_test][q]*w_quad[q];
						}

						M_vals[i_elem][dim_mat][i_test][i_trial] = integral * M_detJacobian[i_elem];
					}
				}
			}
		}
	}

	for ( int k = 0; k < M_numElements; ++k )
	{
		matrix->matrixPtr()->InsertGlobalValues ( ndof_test, M_rows[k], ndof_trial, M_cols[k], M_vals[k][0], Epetra_FECrsMatrix::ROW_MAJOR);
	}

	for ( UInt d1 = 1; d1 < 3 ; d1++ )
	{
		for ( int k = 0; k < M_numElements; ++k )
		{
			for ( UInt i = 0; i <  ndof_trial; i++ )
			{
				M_cols[k][i] += M_numScalarDofs_trial;
			}
			matrix->matrixPtr()->InsertGlobalValues ( ndof_test, M_rows[k], ndof_trial, M_cols[k], M_vals[k][d1], Epetra_FECrsMatrix::ROW_MAJOR);
		}
	}
}
//=====================================================================================================
void
FastAssemblerMixed::assemble_SUPG_block10 ( matrixPtr_Type& matrix, const vector_Type& u_h )
{
    int ndof_test = M_referenceFE_test->nbDof();
    int ndof_trial = M_referenceFE_trial->nbDof();
    int NumQuadPoints = M_qr_integration->nbQuadPt();
    
    double w_quad[NumQuadPoints];
    for ( int q = 0; q < NumQuadPoints ; q++ )
    {
        w_quad[q] = M_qr_integration->weight(q);
    }
    
    #pragma omp parallel firstprivate ( w_quad, ndof_test, ndof_trial, NumQuadPoints )
    {
        int i_elem, i_dof, q, d1, d2, i_test, i_trial, dim_mat, e_idof;
        double integral, integral_partial;
        
        double dphi_phys_test[ndof_test][NumQuadPoints][3]; // Gradient in the physical domain
        double dphi_phys_trial[ndof_trial][NumQuadPoints][3];
        
        double uhq[3][NumQuadPoints];
        
        // ELEMENTI
        #pragma omp for
        for ( i_elem = 0; i_elem <  M_numElements ; i_elem++ )
        {
            // DOF
            for ( i_dof = 0; i_dof <  ndof_test; i_dof++ )
            {
                // QUAD
                for (  q = 0; q < NumQuadPoints ; q++ )
                {
                    // DIM 1
                    for ( d1 = 0; d1 < 3 ; d1++ )
                    {
                        dphi_phys_test[i_dof][q][d1] = 0.0;
                        
                        // DIM 2
                        for ( d2 = 0; d2 < 3 ; d2++ )
                        {
                            dphi_phys_test[i_dof][q][d1] += M_invJacobian[i_elem][d1][d2] * M_dphi_test[i_dof][q][d2];
                        }
                    }
                }
            }
            
            // DOF -- nota che può essere ottimizzato rishapando il gradiente fisico. Prima q poi d1 poi d2 e poi i_dof
            for ( i_dof = 0; i_dof <  ndof_trial; i_dof++ )
            {
                // QUAD
                for (  q = 0; q < NumQuadPoints ; q++ )
                {
                    // DIM 1
                    for ( d1 = 0; d1 < 3 ; d1++ )
                    {
                        dphi_phys_trial[i_dof][q][d1] = 0.0;

                        // DIM 2
                        for ( d2 = 0; d2 < 3 ; d2++ )
                        {
                            dphi_phys_trial[i_dof][q][d1] += M_invJacobian[i_elem][d1][d2] * M_dphi_trial[i_dof][q][d2];
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
                    for ( i_dof = 0; i_dof < ndof_trial; i_dof++ )
                    {
                        e_idof =  M_elements_trial[i_elem][i_dof] + d1*M_numScalarDofs_trial  ;
                        uhq[d1][q] += u_h[e_idof] * M_phi_trial[i_dof][q];
                    }
                }
            }
            
            for ( dim_mat = 0; dim_mat < 3 ; dim_mat++ )
            {
                // DOF - test
                for ( i_test = 0; i_test <  ndof_test; i_test++ )
                {
                    M_rows[i_elem][i_test] = M_elements_test[i_elem][i_test];
                    
                    // DOF - trial
                    for ( i_trial = 0; i_trial < ndof_trial; i_trial++ )
                    {
                        M_cols[i_elem][i_trial] = M_elements_trial[i_elem][i_trial];
                        
                        integral = 0.0;
                        
                        // QUAD
                        for ( q = 0; q < NumQuadPoints ; q++ )
                        {
                            integral_partial = 0.0;
                            for ( d1 = 0; d1 < 3 ; d1++ )
                            {
                                integral_partial += dphi_phys_trial[i_trial][q][d1] * uhq[d1][q];
                            }
                            integral += ( dphi_phys_test[i_test][q][dim_mat] * M_phi_trial[i_trial][q] +
                                          dphi_phys_test[i_test][q][dim_mat] * integral_partial ) * w_quad[q];
                        }
                        
                        M_vals[i_elem][dim_mat][i_test][i_trial] = integral * M_detJacobian[i_elem];
                    }
                }
            }
        }
    }
    
    for ( int k = 0; k < M_numElements; ++k )
    {
        matrix->matrixPtr()->InsertGlobalValues ( ndof_test, M_rows[k], ndof_trial, M_cols[k], M_vals[k][0], Epetra_FECrsMatrix::ROW_MAJOR);
    }
    
    for ( UInt d1 = 1; d1 < 3 ; d1++ )
    {
        for ( int k = 0; k < M_numElements; ++k )
        {
            for ( UInt i = 0; i <  ndof_trial; i++ )
            {
                M_cols[k][i] += M_numScalarDofs_trial;
            }
            matrix->matrixPtr()->InsertGlobalValues ( ndof_test, M_rows[k], ndof_trial, M_cols[k], M_vals[k][d1], Epetra_FECrsMatrix::ROW_MAJOR);
        }
    }
}
//=====================================================================================================
void
FastAssemblerMixed::assemble_SUPG_block01 ( matrixPtr_Type& matrix, const vector_Type& u_h )
{
    int ndof_test = M_referenceFE_test->nbDof();
    int ndof_trial = M_referenceFE_trial->nbDof();
    int NumQuadPoints = M_qr_integration->nbQuadPt();
    
    double w_quad[NumQuadPoints];
    for ( int q = 0; q < NumQuadPoints ; q++ )
    {
        w_quad[q] = M_qr_integration->weight(q);
    }
    
    #pragma omp parallel firstprivate ( w_quad, ndof_test, ndof_trial, NumQuadPoints )
    {
        int i_elem, i_dof, q, d1, d2, i_test, i_trial, dim_mat, e_idof;
        double integral, integral_partial;
        
        double dphi_phys_test[ndof_test][NumQuadPoints][3]; // Gradient in the physical domain
        double dphi_phys_trial[ndof_trial][NumQuadPoints][3];
        
        double uhq[3][NumQuadPoints];
        
        // ELEMENTI
        #pragma omp for
        for ( i_elem = 0; i_elem <  M_numElements ; i_elem++ )
        {
            // DOF
            for ( i_dof = 0; i_dof <  ndof_test; i_dof++ )
            {
                // QUAD
                for (  q = 0; q < NumQuadPoints ; q++ )
                {
                    // DIM 1
                    for ( d1 = 0; d1 < 3 ; d1++ )
                    {
                        dphi_phys_test[i_dof][q][d1] = 0.0;
                        
                        // DIM 2
                        for ( d2 = 0; d2 < 3 ; d2++ )
                        {
                            dphi_phys_test[i_dof][q][d1] += M_invJacobian[i_elem][d1][d2] * M_dphi_test[i_dof][q][d2];
                        }
                    }
                }
            }
            
            // DOF -- nota che può essere ottimizzato rishapando il gradiente fisico. Prima q poi d1 poi d2 e poi i_dof
            for ( i_dof = 0; i_dof <  ndof_trial; i_dof++ )
            {
                // QUAD
                for (  q = 0; q < NumQuadPoints ; q++ )
                {
                    // DIM 1
                    for ( d1 = 0; d1 < 3 ; d1++ )
                    {
                        dphi_phys_trial[i_dof][q][d1] = 0.0;
                        
                        // DIM 2
                        for ( d2 = 0; d2 < 3 ; d2++ )
                        {
                            dphi_phys_trial[i_dof][q][d1] += M_invJacobian[i_elem][d1][d2] * M_dphi_trial[i_dof][q][d2];
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
                    for ( i_dof = 0; i_dof < ndof_test; i_dof++ )
                    {
                        e_idof =  M_elements_test[i_elem][i_dof] + d1*M_numScalarDofs_test  ;
                        uhq[d1][q] += u_h[e_idof] * M_phi_test[i_dof][q];
                    }
                }
            }
            
            for ( dim_mat = 0; dim_mat < 3 ; dim_mat++ )
            {
                // DOF - test
                for ( i_test = 0; i_test <  ndof_test; i_test++ )
                {
                    M_rows[i_elem][i_test] = M_elements_test[i_elem][i_test];
                    
                    // DOF - trial
                    for ( i_trial = 0; i_trial < ndof_trial; i_trial++ )
                    {
                        M_cols[i_elem][i_trial] = M_elements_trial[i_elem][i_trial];
                        
                        integral = 0.0;
                        
                        // QUAD
                        for ( q = 0; q < NumQuadPoints ; q++ )
                        {
                            integral_partial = 0.0;
                            for ( d1 = 0; d1 < 3 ; d1++ )
                            {
                                integral_partial += dphi_phys_test[i_test][q][d1] * uhq[d1][q];
                            }
                            integral += dphi_phys_trial[i_trial][q][dim_mat] * integral_partial * w_quad[q];
                        }
                        
                        M_vals[i_elem][dim_mat][i_test][i_trial] = integral * M_detJacobian[i_elem];
                    }
                }
            }
        }
    }
    
    for ( int k = 0; k < M_numElements; ++k )
    {
        matrix->matrixPtr()->InsertGlobalValues ( ndof_test, M_rows[k], ndof_trial, M_cols[k], M_vals[k][0], Epetra_FECrsMatrix::ROW_MAJOR);
    }
    
    for ( UInt d1 = 1; d1 < 3 ; d1++ )
    {
        for ( int k = 0; k < M_numElements; ++k )
        {
            for ( UInt i = 0; i <  ndof_test; i++ )
            {
                M_rows[k][i] += M_numScalarDofs_test;
            }
            matrix->matrixPtr()->InsertGlobalValues ( ndof_test, M_rows[k], ndof_trial, M_cols[k], M_vals[k][d1], Epetra_FECrsMatrix::ROW_MAJOR);
        }
    }
}




