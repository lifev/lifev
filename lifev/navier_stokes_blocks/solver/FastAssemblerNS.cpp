#include <lifev/navier_stokes_blocks/solver/FastAssemblerNS.hpp>
#include <chrono>
#include <omp.h>
using namespace std::chrono;

using namespace std;

using namespace LifeV;

//=========================================================================
// Constructor //
FastAssemblerNS::FastAssemblerNS( const meshPtr_Type& mesh, const commPtr_Type& comm,
        						  const ReferenceFE* refFE_velocity, const ReferenceFE* refFE_pressure,
        						  const fespacePtr_Type& fespace_velocity, const fespacePtr_Type& fespace_pressure,
        						  const qr_Type* qr ) :
		M_mesh ( mesh ),
		M_comm ( comm ),
		M_referenceFE_velocity ( refFE_velocity ),
		M_referenceFE_pressure ( refFE_pressure ),
		M_fespace_velocity ( fespace_velocity ),
		M_fespace_pressure ( fespace_pressure ),
		M_qr ( qr ),
		M_useSUPG ( false )
{
}
//=========================================================================
// Destructor //
FastAssemblerNS::~FastAssemblerNS()
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

	for ( int i = 0; i < M_referenceFE_velocity->nbDof(); i++ )
	{
		delete [] M_phi_velocity[i];
	}

	delete [] M_phi_velocity;

	//---------------------

	for( int i = 0 ; i < M_referenceFE_velocity->nbDof() ; i++ )
	{
		for ( int j = 0; j < M_qr->nbQuadPt(); j++ )
		{
			delete [] M_dphi_velocity[i][j];
		}
		delete [] M_dphi_velocity[i];
	}
	delete [] M_dphi_velocity;

	//---------------------

	for( int i = 0 ; i < M_referenceFE_velocity->nbDof() ; i++ )
	{
		for ( int j = 0; j < M_qr->nbQuadPt(); j++ )
		{
			for ( int k = 0; k < 3; k++ )
			{
				delete [] M_d2phi_velocity[i][j][k];
			}
			delete [] M_d2phi_velocity[i][j];
		}
		delete [] M_d2phi_velocity[i];
	}
	delete [] M_d2phi_velocity;

	//---------------------

	for ( int i = 0; i < M_referenceFE_pressure->nbDof(); i++ )
	{
		delete [] M_phi_pressure[i];
	}

	delete [] M_phi_pressure;

	//---------------------

	for( int i = 0 ; i < M_referenceFE_pressure->nbDof() ; i++ )
	{
		for ( int j = 0; j < M_qr->nbQuadPt(); j++ )
		{
			delete [] M_dphi_pressure[i][j];
		}
		delete [] M_dphi_pressure[i];
	}
	delete [] M_dphi_pressure;

	//---------------------

	for( int i = 0 ; i < M_numElements ; i++ )
	{
		delete [] M_elements_velocity[i];
	}
	delete [] M_elements_velocity;

	//---------------------

	for( int i = 0 ; i < M_numElements ; i++ )
	{
		delete [] M_elements_pressure[i];
	}
	delete [] M_elements_pressure;

	//---------------------

	for( int i = 0 ; i < M_numElements ; i++ )
	{
		delete [] M_rows_velocity[i];
		delete [] M_cols_velocity[i];
		delete [] M_rows_pressure[i];
		delete [] M_cols_pressure[i];
	}
	delete [] M_rows_velocity;
	delete [] M_cols_velocity;
	delete [] M_rows_pressure;
	delete [] M_cols_pressure;

	//---------------------

	for( int i = 0 ; i < M_numElements ; i++ )
	{
		for ( int j = 0; j < M_referenceFE_velocity->nbDof(); j++ )
		{
			delete [] M_vals_00[i][j];
		}
		delete [] M_vals_00[i];
	}
	delete [] M_vals_00;

	//---------------------

	for( int i = 0 ; i < M_numElements ; i++ )
	{
		for ( int j = 0; j < M_referenceFE_velocity->nbDof(); j++ )
		{
			delete [] M_vals_01[i][j];
		}
		delete [] M_vals_01[i];
	}
	delete [] M_vals_01;

	//---------------------

	for( int i = 0 ; i < M_numElements ; i++ )
	{
		for ( int j = 0; j < M_referenceFE_pressure->nbDof(); j++ )
		{
			delete [] M_vals_10[i][j];
		}
		delete [] M_vals_10[i];
	}
	delete [] M_vals_10;

	//---------------------

	for( int i = 0 ; i < M_numElements ; i++ )
	{
		for ( int j = 0; j < M_referenceFE_pressure->nbDof(); j++ )
		{
			delete [] M_vals_11[i][j];
		}
		delete [] M_vals_11[i];
	}
	delete [] M_vals_11;

	if ( M_useSUPG )
	{
		for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
		{
			for ( int k = 0; k <  3; k++ )
			{
				for ( int z = 0; z <  3; z++ )
				{
					for ( int i = 0; i <  M_referenceFE_velocity->nbDof(); i++ )
					{
						delete [] M_vals_supg[i_elem][k][z][i];
					}
					delete [] M_vals_supg[i_elem][k][z];
				}
				delete [] M_vals_supg[i_elem][k];
			}
			delete [] M_vals_supg[i_elem];
		}
		delete [] M_vals_supg;

		for( int i = 0 ; i < M_numElements ; i++ )
		{
			for ( int k = 0; k < 3; k++ )
			{
				for ( int j = 0; j <  M_referenceFE_velocity->nbDof(); j++ )
				{
					delete [] M_vals_supg_01[i][k][j];
				}
				delete [] M_vals_supg_01[i][k];
			}
			delete [] M_vals_supg_01[i];
		}
		delete [] M_vals_supg_01;

		for( int i = 0 ; i < M_numElements ; i++ )
		{
			for ( int k = 0; k < 3; k++ )
			{
				for ( int j = 0; j <  M_referenceFE_pressure->nbDof(); j++ )
				{
					delete [] M_vals_supg_10[i][k][j];
				}
				delete [] M_vals_supg_10[i][k];
			}
			delete [] M_vals_supg_10[i];
		}
		delete [] M_vals_supg_10;

		for( int i = 0 ; i < M_numElements ; i++ )
		{
			delete [] M_rows_tmp[i];
			delete [] M_cols_tmp[i];
		}
		delete [] M_rows_tmp;
		delete [] M_cols_tmp;


		for( int i = 0 ; i < M_numElements ; i++ )
		{
			for ( int j = 0; j < 3; j++ )
			{
				delete [] M_G[i][j];
			}
			delete [] M_G[i];
			delete [] M_g[i];
			delete [] M_Tau_M[i];
			delete [] M_Tau_C[i];
            delete [] M_Tau_M_hat[i];
		}
		delete [] M_G;
		delete [] M_g;
		delete [] M_Tau_M;
		delete [] M_Tau_C;
        delete [] M_Tau_M_hat;
	}

}
//=========================================================================
void
FastAssemblerNS::setConstants_NavierStokes( const Real& density, const Real& viscosity, const Real& timestep, const Real& orderBDF, const Real& C_I )
{
	M_density = density;
	M_viscosity = viscosity;
	M_timestep = timestep;
	M_orderBDF = orderBDF;
	M_C_I = C_I;
}
//=========================================================================
void
FastAssemblerNS::setConstants_NavierStokes( const Real& density, const Real& viscosity, const Real& timestep, const Real& orderBDF, const Real& C_I, const Real& alpha )
{
    M_density = density;
    M_viscosity = viscosity;
    M_timestep = timestep;
    M_orderBDF = orderBDF;
    M_C_I = C_I;
    M_alpha = alpha;
}
//=========================================================================
void
FastAssemblerNS::allocateSpace ( CurrentFE* current_fe_velocity, const bool& use_supg )
{
	M_useSUPG = use_supg;

	M_numElements = M_mesh->numVolumes();

	M_numScalarDofs = M_fespace_velocity->dof().numTotalDof();

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
		current_fe_velocity->update( M_mesh->element (i), UPDATE_D2PHI );
		M_detJacobian[i] = current_fe_velocity->detJacobian(0);
		for ( int j = 0; j < 3; j++ )
		{
			for ( int k = 0; k < 3; k++ )
			{
				M_invJacobian[i][j][k] = current_fe_velocity->tInverseJacobian(j,k,2);
			}
		}
	}

	//-------------------------------------------------------------------------------------------------
	// PHI VELOCITY
	//-------------------------------------------------------------------------------------------------

	M_phi_velocity = new double* [ M_referenceFE_velocity->nbDof() ];
	for ( int i = 0; i < M_referenceFE_velocity->nbDof(); i++ )
	{
		M_phi_velocity[i] = new double [ M_qr->nbQuadPt() ];
	}

	// PHI REF
	for (UInt q (0); q < M_qr->nbQuadPt(); ++q)
	{
		for (UInt j (0); j < M_referenceFE_velocity->nbDof(); ++j)
		{
			M_phi_velocity[j][q] =  M_referenceFE_velocity->phi (j, M_qr->quadPointCoor (q) );
		}
	}

	//-------------------------------------------------------------------------------------------------
	// DPHI VELOCITY
	//-------------------------------------------------------------------------------------------------

	M_dphi_velocity = new double** [ M_referenceFE_velocity->nbDof() ];

	for ( int i = 0; i <  M_referenceFE_velocity->nbDof(); i++ )
	{
		M_dphi_velocity[i] = new double* [ M_qr->nbQuadPt() ];
		for ( int j = 0; j < M_qr->nbQuadPt() ; j++ )
		{
			M_dphi_velocity[i][j] = new double [ 3 ];
		}
	}

	//DPHI REF
	for (UInt q (0); q < M_qr->nbQuadPt(); ++q)
	{
		for (UInt i (0); i < M_referenceFE_velocity->nbDof(); ++i)
		{
			for (UInt j (0); j < 3; ++j)
			{
				M_dphi_velocity[i][q][j] = M_referenceFE_velocity->dPhi (i, j, M_qr->quadPointCoor (q) );
			}
		}
	}

	//-------------------------------------------------------------------------------------------------
	// D2PHI VELOCITY
	//-------------------------------------------------------------------------------------------------

	M_d2phi_velocity = new double*** [ M_referenceFE_velocity->nbDof() ];

	for ( int i = 0; i <  M_referenceFE_velocity->nbDof(); i++ )
	{
		M_d2phi_velocity[i] = new double** [ M_qr->nbQuadPt() ];
		for ( int j = 0; j < M_qr->nbQuadPt() ; j++ )
		{
			M_d2phi_velocity[i][j] = new double* [ 3 ];
			for ( int k = 0; k < 3 ; k++ )
			{
				M_d2phi_velocity[i][j][k] = new double [ 3 ];
			}
		}
	}

	//D2PHI REF
	for (UInt q (0); q < M_qr->nbQuadPt(); ++q)
	{
		for (UInt i (0); i < M_referenceFE_velocity->nbDof(); ++i)
		{
			for (UInt j (0); j < 3; ++j)
			{
				for (UInt k (0); k < 3; ++k)
				{
					M_d2phi_velocity[i][q][j][k] = M_referenceFE_velocity->d2Phi (i, j, k, M_qr->quadPointCoor (q) );
				}
			}
		}
	}

	//-------------------------------------------------------------------------------------------------
	// PHI PRESSURE
	//-------------------------------------------------------------------------------------------------

	M_phi_pressure = new double* [ M_referenceFE_pressure->nbDof() ];
	for ( int i = 0; i < M_referenceFE_pressure->nbDof(); i++ )
	{
		M_phi_pressure[i] = new double [ M_qr->nbQuadPt() ];
	}

	// PHI REF
	for (UInt q (0); q < M_qr->nbQuadPt(); ++q)
	{
		for (UInt j (0); j < M_referenceFE_pressure->nbDof(); ++j)
		{
			M_phi_pressure[j][q] =  M_referenceFE_pressure->phi (j, M_qr->quadPointCoor (q) );
		}
	}

	//-------------------------------------------------------------------------------------------------
	// DPHI PRESSURE
	//-------------------------------------------------------------------------------------------------

	M_dphi_pressure = new double** [ M_referenceFE_pressure->nbDof() ];

	for ( int i = 0; i <  M_referenceFE_pressure->nbDof(); i++ )
	{
		M_dphi_pressure[i] = new double* [ M_qr->nbQuadPt() ];
		for ( int j = 0; j < M_qr->nbQuadPt() ; j++ )
		{
			M_dphi_pressure[i][j] = new double [ 3 ];
		}
	}

	//DPHI REF
	for (UInt q (0); q < M_qr->nbQuadPt(); ++q)
	{
		for (UInt i (0); i < M_referenceFE_pressure->nbDof(); ++i)
		{
			for (UInt j (0); j < 3; ++j)
			{
				M_dphi_pressure[i][q][j] = M_referenceFE_pressure->dPhi (i, j, M_qr->quadPointCoor (q) );
			}
		}
	}

	//-------------------------------------------------------------------------------------------------
	// CONNECTIVITY TEST
	//-------------------------------------------------------------------------------------------------

	M_elements_velocity = new double* [ M_numElements ];

	for ( int i = 0; i <  M_numElements; i++ )
	{
		M_elements_velocity[i] = new double [ M_referenceFE_velocity->nbDof() ];
	}

	for ( int i = 0; i <  M_numElements; i++ )
	{
		for ( int j = 0; j <  M_referenceFE_velocity->nbDof(); j++ )
		{
			M_elements_velocity[i][j] = M_fespace_velocity->dof().localToGlobalMap (i, j);
		}
	}

	//-------------------------------------------------------------------------------------------------
	// CONNECTIVITY TRIAL
	//-------------------------------------------------------------------------------------------------

	M_elements_pressure = new double* [ M_numElements ];

	for ( int i = 0; i <  M_numElements; i++ )
	{
		M_elements_pressure[i] = new double [ M_referenceFE_pressure->nbDof() ];
	}

	for ( int i = 0; i <  M_numElements; i++ )
	{
		for ( int j = 0; j <  M_referenceFE_pressure->nbDof(); j++ )
		{
			M_elements_pressure[i][j] = M_fespace_pressure->dof().localToGlobalMap (i, j);
		}
	}

	//-------------------------------------------------------------------------------------------------
	// ROWS VELOCITY
	//-------------------------------------------------------------------------------------------------

	M_rows_velocity = new int* [M_numElements];

	for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
	{
		M_rows_velocity[i_elem] = new int [ M_referenceFE_velocity->nbDof() ];
	}

	//-------------------------------------------------------------------------------------------------
	// ROWS PRESSURE
	//-------------------------------------------------------------------------------------------------

	M_rows_pressure = new int* [M_numElements];

	for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
	{
		M_rows_pressure[i_elem] = new int [ M_referenceFE_pressure->nbDof() ];
	}

	//-------------------------------------------------------------------------------------------------
	// COLS VELOCITY
	//-------------------------------------------------------------------------------------------------

	M_cols_velocity = new int* [M_numElements];

	for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
	{
		M_cols_velocity[i_elem] = new int [ M_referenceFE_velocity->nbDof() ];
	}

	//-------------------------------------------------------------------------------------------------
	// COLS PRESSURE
	//-------------------------------------------------------------------------------------------------

	M_cols_pressure = new int* [M_numElements];

	for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
	{
		M_cols_pressure[i_elem] = new int [ M_referenceFE_pressure->nbDof() ];
	}

	//-------------------------------------------------------------------------------------------------
	// VALUES BLOCK (0,0)
	//-------------------------------------------------------------------------------------------------

	M_vals_00 = new double** [M_numElements];

	for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
	{
		M_vals_00[i_elem] = new double* [ M_referenceFE_velocity->nbDof() ];
		for ( int i = 0; i <  M_referenceFE_velocity->nbDof(); i++ )
		{
			M_vals_00[i_elem][i] = new double [ M_referenceFE_velocity->nbDof() ];
		}
	}

	//-------------------------------------------------------------------------------------------------
	// VALUES BLOCK (0,1)
	//-------------------------------------------------------------------------------------------------

	M_vals_01 = new double** [M_numElements];

	for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
	{
		M_vals_01[i_elem] = new double* [ M_referenceFE_velocity->nbDof() ];
		for ( int i = 0; i <  M_referenceFE_velocity->nbDof(); i++ )
		{
			M_vals_01[i_elem][i] = new double [ M_referenceFE_pressure->nbDof() ];
		}
	}

	//-------------------------------------------------------------------------------------------------
	// VALUES BLOCK (1,0)
	//-------------------------------------------------------------------------------------------------

	M_vals_10 = new double** [M_numElements];

	for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
	{
		M_vals_10[i_elem] = new double* [ M_referenceFE_pressure->nbDof() ];
		for ( int i = 0; i <  M_referenceFE_pressure->nbDof(); i++ )
		{
			M_vals_10[i_elem][i] = new double [ M_referenceFE_velocity->nbDof() ];
		}
	}

	//-------------------------------------------------------------------------------------------------
	// VALUES BLOCK (1,1)
	//-------------------------------------------------------------------------------------------------

	M_vals_11 = new double** [M_numElements];

	for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
	{
		M_vals_11[i_elem] = new double* [ M_referenceFE_pressure->nbDof() ];
		for ( int i = 0; i <  M_referenceFE_pressure->nbDof(); i++ )
		{
			M_vals_11[i_elem][i] = new double [ M_referenceFE_pressure->nbDof() ];
		}
	}

	if ( M_useSUPG )
	{
		//-------------------------------------------------------------------------------------------------
		// SUPG MATRIX BLOCK (0,0) - because term with Tau_C leads to blocks on mixed components
		//-------------------------------------------------------------------------------------------------

		M_vals_supg = new double**** [M_numElements];

		for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
		{
			M_vals_supg[i_elem] = new double*** [ 3 ];
			for ( int k = 0; k <  3; k++ )
			{
				M_vals_supg[i_elem][k] = new double** [ 3 ];
				for ( int z = 0; z < 3; z++ )
				{
					M_vals_supg[i_elem][k][z] = new double* [ M_referenceFE_velocity->nbDof() ];
					for ( int i = 0; i <  M_referenceFE_velocity->nbDof(); i++ )
					{
						M_vals_supg[i_elem][k][z][i] = new double [ M_referenceFE_velocity->nbDof() ];
						for ( int j = 0; j <  M_referenceFE_velocity->nbDof(); j++ )
						{
							M_vals_supg[i_elem][k][z][i][j] = 0.0;
						}
					}
				}
			}
		}

		M_vals_supg_01 = new double*** [M_numElements];
		for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
		{
			M_vals_supg_01[i_elem] = new double** [ 3 ];
			for ( int k = 0; k <  3; k++ )
			{
				M_vals_supg_01[i_elem][k] = new double* [ M_referenceFE_velocity->nbDof() ];
				for ( int i = 0; i <  M_referenceFE_velocity->nbDof(); i++ )
				{
					M_vals_supg_01[i_elem][k][i] = new double [ M_referenceFE_pressure->nbDof() ];
				}
			}
		}

		M_vals_supg_10 = new double*** [M_numElements];
		for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
		{
			M_vals_supg_10[i_elem] = new double** [ 3 ];
			for ( int k = 0; k <  3; k++ )
			{
				M_vals_supg_10[i_elem][k] = new double* [ M_referenceFE_pressure->nbDof() ];
				for ( int i = 0; i <  M_referenceFE_pressure->nbDof(); i++ )
				{
					M_vals_supg_10[i_elem][k][i] = new double [ M_referenceFE_velocity->nbDof() ];
				}
			}
		}

		M_rows_tmp = new int* [M_numElements];

		for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
		{
			M_rows_tmp[i_elem] = new int [ M_referenceFE_velocity->nbDof() ];
		}

		M_cols_tmp = new int* [M_numElements];

		for ( int i_elem = 0; i_elem <  M_numElements; i_elem++ )
		{
			M_cols_tmp[i_elem] = new int [ M_referenceFE_velocity->nbDof() ];
		}

		//-------------------------------------------------------------------------------------------------
		// METRIC TENSOR AND VECTOR
		//-------------------------------------------------------------------------------------------------

		M_G = new double** [ M_numElements];
		M_g = new double* [ M_numElements];

		for ( int i = 0; i <  M_numElements; i++ )
		{
			M_G[i] = new double* [ 3 ];
			M_g[i] = new double [ 3 ];
			for ( int j = 0; j < 3 ; j++ )
			{
				M_G[i][j] = new double [ 3 ];
			}
		}

		//-------------------------------------------------------------------------------------------------
		// STABILIZATION PARAMETERS
		//-------------------------------------------------------------------------------------------------

		M_Tau_M = new double* [ M_numElements];
		M_Tau_C = new double* [ M_numElements];
        M_Tau_M_hat = new double* [ M_numElements];
		for ( int i = 0; i <  M_numElements; i++ )
		{
			M_Tau_M[i] = new double [ M_qr->nbQuadPt() ];
			M_Tau_C[i] = new double [ M_qr->nbQuadPt() ];
            M_Tau_M_hat[i] = new double [ M_qr->nbQuadPt() ];
		}

		//-------------------------------------------------------------------------------------------------

		for ( int i = 0; i < M_numElements; i++ )
		{
			M_G[i][0][0] = M_invJacobian[i][0][0]*M_invJacobian[i][0][0] + M_invJacobian[i][0][1]*M_invJacobian[i][0][1] + M_invJacobian[i][0][2]*M_invJacobian[i][0][2];
			M_G[i][0][1] = M_invJacobian[i][0][0]*M_invJacobian[i][1][0] + M_invJacobian[i][0][1]*M_invJacobian[i][1][1] + M_invJacobian[i][0][2]*M_invJacobian[i][1][2];
			M_G[i][0][2] = M_invJacobian[i][0][0]*M_invJacobian[i][2][0] + M_invJacobian[i][0][1]*M_invJacobian[i][2][1] + M_invJacobian[i][0][2]*M_invJacobian[i][2][2];

			M_G[i][1][0] = M_invJacobian[i][1][0]*M_invJacobian[i][0][0] + M_invJacobian[i][1][1]*M_invJacobian[i][0][1] + M_invJacobian[i][1][2]*M_invJacobian[i][0][2];
			M_G[i][1][1] = M_invJacobian[i][1][0]*M_invJacobian[i][1][0] + M_invJacobian[i][1][1]*M_invJacobian[i][1][1] + M_invJacobian[i][1][2]*M_invJacobian[i][1][2];
			M_G[i][1][2] = M_invJacobian[i][1][0]*M_invJacobian[i][2][0] + M_invJacobian[i][1][1]*M_invJacobian[i][2][1] + M_invJacobian[i][1][2]*M_invJacobian[i][2][2];

			M_G[i][2][0] = M_invJacobian[i][2][0]*M_invJacobian[i][0][0] + M_invJacobian[i][2][1]*M_invJacobian[i][0][1] + M_invJacobian[i][2][2]*M_invJacobian[i][0][2];
			M_G[i][2][1] = M_invJacobian[i][2][0]*M_invJacobian[i][1][0] + M_invJacobian[i][2][1]*M_invJacobian[i][1][1] + M_invJacobian[i][2][2]*M_invJacobian[i][1][2];
			M_G[i][2][2] = M_invJacobian[i][2][0]*M_invJacobian[i][2][0] + M_invJacobian[i][2][1]*M_invJacobian[i][2][1] + M_invJacobian[i][2][2]*M_invJacobian[i][2][2];

			M_g[i][0] = M_invJacobian[i][0][0] + M_invJacobian[i][0][1] + M_invJacobian[i][0][2];
			M_g[i][1] = M_invJacobian[i][1][0] + M_invJacobian[i][1][1] + M_invJacobian[i][1][2];
			M_g[i][2] = M_invJacobian[i][2][0] + M_invJacobian[i][2][1] + M_invJacobian[i][2][2];
		}
	}
}
//=========================================================================
void
FastAssemblerNS::assemble_supg_terms( matrixPtr_Type& block00, matrixPtr_Type& block01, matrixPtr_Type& block10, matrixPtr_Type& block11, const vector_Type& u_h )
{
	int ndof_velocity = M_referenceFE_velocity->nbDof();
	int ndof_pressure = M_referenceFE_pressure->nbDof();
	int NumQuadPoints = M_qr->nbQuadPt();
	int ndof_vec = M_referenceFE_velocity->nbDof()*3;
	M_numScalarDofs =  M_fespace_velocity->dof().numTotalDof();

	double w_quad[NumQuadPoints];
	for ( int q = 0; q < NumQuadPoints ; q++ )
	{
		w_quad[q] = M_qr->weight(q);
	}

	#pragma omp parallel firstprivate( w_quad, ndof_velocity, ndof_pressure, NumQuadPoints)
	{
		int i_elem, i_dof, q, d1, d2, i_test, i_trial, e_idof, dim_mat, iCoor, jCoor, k1, k2;
		double integral, integral_test, integral_trial, integral_partial, integral_lapl, partialSum;

		double dphi_phys_velocity[ndof_velocity][NumQuadPoints][3];
		double dphi_phys_pressure[ndof_pressure][NumQuadPoints][3];
		double d2phi_phys_velocity[ndof_velocity][NumQuadPoints][3][3];
		double uhq[3][NumQuadPoints];

		// ELEMENTI,
		#pragma omp for
		for ( i_elem = 0; i_elem <  M_numElements; i_elem++ )
		{
			// DOF
			for ( i_dof = 0; i_dof <  ndof_velocity; i_dof++ )
			{
				// QUAD
				for (  q = 0; q < NumQuadPoints ; q++ )
				{
					// DIM 1
					for ( d1 = 0; d1 < 3 ; d1++ )
					{
						dphi_phys_velocity[i_dof][q][d1] = 0.0;

						// DIM 2
						for ( d2 = 0; d2 < 3 ; d2++ )
						{
							dphi_phys_velocity[i_dof][q][d1] += M_invJacobian[i_elem][d1][d2] * M_dphi_velocity[i_dof][q][d2];
						}
					}
				}
			}

			// DOF -- nota che può essere ottimizzato rishapando il gradiente fisico. Prima q poi d1 poi d2 e poi i_dof
			for ( i_dof = 0; i_dof <  ndof_pressure; i_dof++ )
			{
				// QUAD
				for (  q = 0; q < NumQuadPoints ; q++ )
				{
					// DIM 1
					for ( d1 = 0; d1 < 3 ; d1++ )
					{
						dphi_phys_pressure[i_dof][q][d1] = 0.0;

						// DIM 2
						for ( d2 = 0; d2 < 3 ; d2++ )
						{
							dphi_phys_pressure[i_dof][q][d1] += M_invJacobian[i_elem][d1][d2] * M_dphi_pressure[i_dof][q][d2];
						}
					}
				}
			}

			// QUAD
			for (  q = 0; q < NumQuadPoints ; q++ )
			{
				// DOF
				for ( i_dof = 0; i_dof <  ndof_velocity; i_dof++ )
				{
					// DIM 1
					for ( iCoor = 0; iCoor < 3 ; iCoor++ )
					{
						// DIM 2
						for ( jCoor = 0; jCoor < 3 ; jCoor++ )
						{
							partialSum = 0.0;
							// DIM 1
							for ( k1 = 0; k1 < 3 ; k1++ )
							{
								// DIM 2
								for ( k2 = 0; k2 < 3 ; k2++ )
								{
									partialSum += M_invJacobian[i_elem][iCoor][k1]
									              * M_d2phi_velocity[i_dof][q][k1][k2]
									              * M_invJacobian[i_elem][jCoor][k2];
								}
							}
							d2phi_phys_velocity[i_dof][q][iCoor][jCoor] = partialSum;
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
					for ( i_dof = 0; i_dof <  ndof_velocity; i_dof++ )
					{
						e_idof =  M_elements_velocity[i_elem][i_dof] + d1*M_numScalarDofs  ;
						uhq[d1][q] += u_h[e_idof] * M_phi_velocity[i_dof][q];
					}
				}
			}

			// STABILIZZAZIONE - coefficienti Tau_M e Tau_C
			for (  q = 0; q < NumQuadPoints ; q++ )
			{
				M_Tau_M[i_elem][q] = 1.0/std::sqrt(
						M_density*M_density*M_orderBDF*M_orderBDF/(M_timestep*M_timestep) // TAU_M_DEN_DT
						+M_density*M_density*( uhq[0][q]*( M_G[i_elem][0][0]* uhq[0][q] +  M_G[i_elem][0][1] * uhq[1][q] + M_G[i_elem][0][2] * uhq[2][q] ) +
											   uhq[1][q]*( M_G[i_elem][1][0]* uhq[0][q] +  M_G[i_elem][1][1] * uhq[1][q] + M_G[i_elem][1][2] * uhq[2][q] ) +
											   uhq[2][q]*( M_G[i_elem][2][0]* uhq[0][q] +  M_G[i_elem][2][1] * uhq[1][q] + M_G[i_elem][2][2] * uhq[2][q]   )
											 ) 	// TAU_M_DEN_VEL
						+M_C_I*M_viscosity*M_viscosity*(
								M_G[i_elem][0][0]*M_G[i_elem][0][0] + M_G[i_elem][0][1]*M_G[i_elem][0][1] + M_G[i_elem][0][2]*M_G[i_elem][0][2] +
								M_G[i_elem][1][0]*M_G[i_elem][1][0] + M_G[i_elem][1][1]*M_G[i_elem][1][1] + M_G[i_elem][1][2]*M_G[i_elem][1][2] +
								M_G[i_elem][2][0]*M_G[i_elem][2][0] + M_G[i_elem][2][1]*M_G[i_elem][2][1] + M_G[i_elem][2][2]*M_G[i_elem][2][2]
											 )  // TAU_M_DEN_VISC*/
				);

				M_Tau_C[i_elem][q] = 1.0/( M_g[i_elem][0]*M_Tau_M[i_elem][q]*M_g[i_elem][0] +
										   M_g[i_elem][1]*M_Tau_M[i_elem][q]*M_g[i_elem][1] +
										   M_g[i_elem][2]*M_Tau_M[i_elem][q]*M_g[i_elem][2]
										  );

			}

			// DOF - test
			for ( i_test = 0; i_test <  ndof_velocity; i_test++ )
			{
				M_rows_velocity[i_elem][i_test] = M_elements_velocity[i_elem][i_test];

				// DOF - trial
				for ( i_trial = 0; i_trial <  ndof_velocity; i_trial++ )
				{
					M_cols_velocity[i_elem][i_trial] = M_elements_velocity[i_elem][i_trial];

					integral = 0.0;
					// QUAD
					for ( q = 0; q < NumQuadPoints ; q++ )
					{
						integral_test = 0.0;
						integral_trial = 0.0;
						integral_lapl = 0.0;

						// DIM 1
						for ( d1 = 0; d1 < 3 ; d1++ )
						{
							integral_lapl += d2phi_phys_velocity[i_trial][q][d1][d1];
							integral_test += uhq[d1][q] * dphi_phys_velocity[i_test][q][d1];   // w grad(phi_i)
							integral_trial += uhq[d1][q] * dphi_phys_velocity[i_trial][q][d1]; // w grad(phi_j) terzo indice dphi_phys e' derivata in x,y,z
						}
						integral += M_Tau_M[i_elem][q] * (integral_test * integral_trial + integral_test * M_phi_velocity[i_trial][q]
						                                  - integral_test * integral_lapl ) * w_quad[q];
						// above, the term "integral_test * M_phi[i_trial][q]" is the one which comes from (w grad(phi_i), phi_j)
					}
					M_vals_supg[i_elem][0][0][i_test][i_trial] = integral *  M_detJacobian[i_elem]; // (w grad(phi_i), w grad(phi_j) )
					M_vals_supg[i_elem][1][1][i_test][i_trial] = integral *  M_detJacobian[i_elem];
					M_vals_supg[i_elem][2][2][i_test][i_trial] = integral *  M_detJacobian[i_elem];

					for ( d1 = 0; d1 < 3; d1++ )
					{
						for ( d2 = 0; d2 < 3; d2++ )
						{
							integral = 0.0;
							// QUAD
							for ( q = 0; q < NumQuadPoints ; q++ )
							{
								integral += M_Tau_C[i_elem][q] * dphi_phys_velocity[i_test][q][d1] * dphi_phys_velocity[i_trial][q][d2] * w_quad[q]; // div(phi_i) * div(phi_j)
							}
							M_vals_supg[i_elem][d1][d2][i_test][i_trial] += integral *  M_detJacobian[i_elem];
						}
					}

				}

				for ( i_trial = 0; i_trial <  ndof_pressure; i_trial++ )
				{
					for ( dim_mat = 0; dim_mat < 3 ; dim_mat++ )
					{
						M_cols_pressure[i_elem][i_trial] = M_elements_pressure[i_elem][i_trial];

						integral = 0.0;

						// QUAD
						for ( q = 0; q < NumQuadPoints ; q++ )
						{
							integral_partial = 0.0;
							for ( d1 = 0; d1 < 3 ; d1++ )
							{
								integral_partial += ( dphi_phys_velocity[i_test][q][d1] * uhq[d1][q] );
							}
							integral += M_Tau_M[i_elem][q] * dphi_phys_pressure[i_trial][q][dim_mat] * integral_partial * w_quad[q];
						}
						M_vals_supg_01[i_elem][dim_mat][i_test][i_trial] = integral * M_detJacobian[i_elem];
					}
				}

			}

			// DOF - test - assemblo blocchi (1,0) e (1,1)
			for ( i_test = 0; i_test <  ndof_pressure; i_test++ )
			{
				M_rows_pressure[i_elem][i_test] = M_elements_pressure[i_elem][i_test];

				// DOF - trial
				for ( i_trial = 0; i_trial <  ndof_pressure; i_trial++ )
				{
					integral = 0.0;
					// QUAD
					for ( q = 0; q < NumQuadPoints ; q++ )
					{
						// DIM 1
						for ( d1 = 0; d1 < 3 ; d1++ )
						{
							integral += M_Tau_M[i_elem][q] * dphi_phys_pressure[i_test][q][d1] * dphi_phys_pressure[i_trial][q][d1]*w_quad[q];
						}
					}
					M_vals_11[i_elem][i_test][i_trial] = 0.0;
					M_vals_11[i_elem][i_test][i_trial] = integral * M_detJacobian[i_elem];
				}
			}

			for ( dim_mat = 0; dim_mat < 3 ; dim_mat++ )
			{
				// DOF - test
				for ( i_test = 0; i_test <  ndof_pressure; i_test++ )
				{
					// DOF - trial
					for ( i_trial = 0; i_trial < ndof_velocity; i_trial++ )
					{
						integral = 0.0;

						// QUAD
						for ( q = 0; q < NumQuadPoints ; q++ )
						{
							integral_partial = 0.0;
							for ( d1 = 0; d1 < 3 ; d1++ )
							{
								integral_partial += ( dphi_phys_velocity[i_trial][q][d1] * uhq[d1][q] - d2phi_phys_velocity[i_trial][q][d1][d1] );
							}
							integral += M_Tau_M[i_elem][q] * ( dphi_phys_pressure[i_test][q][dim_mat] * M_phi_velocity[i_trial][q] +
										                       dphi_phys_pressure[i_test][q][dim_mat] * integral_partial ) * w_quad[q];
						}

						M_vals_supg_10[i_elem][dim_mat][i_test][i_trial] = integral * M_detJacobian[i_elem];
					}
				}
			}
		}
	}

	for ( UInt d1 = 0; d1 < 3 ; d1++ ) // row index
	{
		for ( int k = 0; k < M_numElements; ++k )
		{
			for ( UInt i = 0; i <  ndof_velocity; i++ )
			{
				M_rows_tmp[k][i] = M_rows_velocity[k][i] + d1 * M_numScalarDofs;
			}
			block00->matrixPtr()->InsertGlobalValues ( ndof_velocity, M_rows_tmp[k], ndof_velocity, M_cols_velocity[k], M_vals_supg[k][d1][0], Epetra_FECrsMatrix::ROW_MAJOR);
			block01->matrixPtr()->InsertGlobalValues ( ndof_velocity, M_rows_tmp[k], ndof_pressure, M_cols_pressure[k], M_vals_supg_01[k][d1], Epetra_FECrsMatrix::ROW_MAJOR);
		}

		for ( int k = 0; k < M_numElements; ++k )
		{
			for ( UInt i = 0; i <  ndof_velocity; i++ )
			{
				M_cols_tmp[k][i] = M_cols_velocity[k][i] + M_numScalarDofs;
			}
			block00->matrixPtr()->InsertGlobalValues ( ndof_velocity, M_rows_tmp[k], ndof_velocity, M_cols_tmp[k], M_vals_supg[k][d1][1], Epetra_FECrsMatrix::ROW_MAJOR);
		}

		for ( int k = 0; k < M_numElements; ++k )
		{
			for ( UInt i = 0; i <  ndof_velocity; i++ )
			{
				M_cols_tmp[k][i] = M_cols_velocity[k][i] + 2 * M_numScalarDofs;
			}
			block00->matrixPtr()->InsertGlobalValues ( ndof_velocity, M_rows_tmp[k], ndof_velocity, M_cols_tmp[k], M_vals_supg[k][d1][2], Epetra_FECrsMatrix::ROW_MAJOR);
		}
	}

	for ( int k = 0; k < M_numElements; ++k )
	{
		block11->matrixPtr()->InsertGlobalValues ( ndof_pressure, M_rows_pressure[k], ndof_pressure, M_cols_pressure[k], M_vals_11[k], Epetra_FECrsMatrix::ROW_MAJOR);
	}

	for ( UInt d1 = 0; d1 < 3 ; d1++ )
	{
		for ( int k = 0; k < M_numElements; ++k )
		{
			for ( UInt i = 0; i <  ndof_velocity; i++ )
			{
				M_cols_tmp[k][i] = M_cols_velocity[k][i] + d1 * M_numScalarDofs;
			}
			block10->matrixPtr()->InsertGlobalValues ( ndof_pressure, M_rows_pressure[k], ndof_velocity, M_cols_tmp[k], M_vals_supg_10[k][d1], Epetra_FECrsMatrix::ROW_MAJOR);
		}
	}
}
//=========================================================================
void
FastAssemblerNS::assemble_constant_terms( matrixPtr_Type& mass, matrixPtr_Type& stiffness, matrixPtr_Type& b01, matrixPtr_Type& b10 )
{
    int ndof_velocity = M_referenceFE_velocity->nbDof();
    int ndof_pressure = M_referenceFE_pressure->nbDof();
    int NumQuadPoints = M_qr->nbQuadPt();
    int ndof_vec = M_referenceFE_velocity->nbDof()*3;
    M_numScalarDofs =  M_fespace_velocity->dof().numTotalDof();
    
    double w_quad[NumQuadPoints];
    for ( int q = 0; q < NumQuadPoints ; q++ )
    {
        w_quad[q] = M_qr->weight(q);
    }
    
    #pragma omp parallel firstprivate( w_quad, ndof_velocity, ndof_pressure, NumQuadPoints)
    {
        int i_elem, i_dof, q, d1, d2, i_test, i_trial, dim_mat, dim_mat1, dim_mat2;
        double integral, integral_00, integral_01, integral_02, integral_10, integral_11, integral_12;
        double integral_20, integral_21, integral_22;
        
        double dphi_phys_velocity[ndof_velocity][NumQuadPoints][3]; // Gradient in the physical domain
        
        // ELEMENTI
        #pragma omp for
        for ( i_elem = 0; i_elem <  M_numElements ; i_elem++ )
        {
            // DOF
            for ( i_dof = 0; i_dof <  ndof_velocity; i_dof++ )
            {
                // QUAD
                for (  q = 0; q < NumQuadPoints ; q++ )
                {
                    // DIM 1
                    for ( d1 = 0; d1 < 3 ; d1++ )
                    {
                        dphi_phys_velocity[i_dof][q][d1] = 0.0;
                        
                        // DIM 2
                        for ( d2 = 0; d2 < 3 ; d2++ )
                        {
                            dphi_phys_velocity[i_dof][q][d1] += M_invJacobian[i_elem][d1][d2] * M_dphi_velocity[i_dof][q][d2];
                        }
                    }
                }
            }
            
            // DOF - test
            for ( i_test = 0; i_test <  ndof_velocity; i_test++ )
            {
                M_rows_velocity[i_elem][i_test] = M_elements_velocity[i_elem][i_test];
                
                // DOF - trial
                for ( i_trial = 0; i_trial <  ndof_velocity; i_trial++ )
                {
                    M_cols_velocity[i_elem][i_trial] = M_elements_velocity[i_elem][i_trial];
                    M_vals_00[i_elem][i_test][i_trial] = 0.0;
                    integral = 0.0;
                    // QUAD
                    for ( q = 0; q < NumQuadPoints ; q++ )
                    {
                        integral += M_phi_velocity[i_test][q] * M_phi_velocity[i_trial][q]*w_quad[q];
                    }
                    M_vals_00[i_elem][i_test][i_trial] = M_density * integral * M_detJacobian[i_elem];
                }
            }
            
            /*
            for ( dim_mat1 = 0; dim_mat1 < 3 ; dim_mat1++ )
            {
                for ( dim_mat2 = 0; dim_mat2 < 3 ; dim_mat2++ )
                {
                    // DOF - test
                    for ( i_test = 0; i_test <  ndof_velocity; i_test++ )
                    {
                        // DOF - trial
                        for ( i_trial = 0; i_trial <  ndof_velocity; i_trial++ )
                        {
                            M_vals_supg[i_elem][dim_mat1][dim_mat2][i_test][i_trial] = 0.0;
                            integral = 0.0;
                            // QUAD
                            for ( q = 0; q < NumQuadPoints ; q++ )
                            {
                                integral += dphi_phys_velocity[i_test][q][dim_mat2] *
                                            ( dphi_phys_velocity[i_trial][q][dim_mat2] +
                                              dphi_phys_velocity[i_trial][q][dim_mat1]
                                            ) * w_quad[q];
                            }
                            M_vals_supg[i_elem][dim_mat1][dim_mat2][i_test][i_trial] =
                                                                        integral * M_detJacobian[i_elem];
                        }
                    }
                }
            }
            */

            // DOF - test
            for ( i_test = 0; i_test <  ndof_velocity; i_test++ )
            {
                // DOF - trial
                for ( i_trial = 0; i_trial <  ndof_velocity; i_trial++ )
                {
                    M_vals_supg[i_elem][0][0][i_test][i_trial] = 0.0;
                    M_vals_supg[i_elem][0][1][i_test][i_trial] = 0.0;
                    M_vals_supg[i_elem][0][2][i_test][i_trial] = 0.0;
                    M_vals_supg[i_elem][1][0][i_test][i_trial] = 0.0;
                    M_vals_supg[i_elem][1][1][i_test][i_trial] = 0.0;
                    M_vals_supg[i_elem][1][2][i_test][i_trial] = 0.0;
                    M_vals_supg[i_elem][2][0][i_test][i_trial] = 0.0;
                    M_vals_supg[i_elem][2][1][i_test][i_trial] = 0.0;
                    M_vals_supg[i_elem][2][2][i_test][i_trial] = 0.0;
                    integral_00 = 0.0;
                    integral_01 = 0.0;
                    integral_02 = 0.0;
                    integral_10 = 0.0;
                    integral_11 = 0.0;
                    integral_12 = 0.0;
                    integral_20 = 0.0;
                    integral_21 = 0.0;
                    integral_22 = 0.0;
                    // QUAD
                    for ( q = 0; q < NumQuadPoints ; q++ )
                    {
                        integral_00 += ( 2*dphi_phys_velocity[i_test][q][0] * dphi_phys_velocity[i_trial][q][0] +
                                         dphi_phys_velocity[i_test][q][1] * dphi_phys_velocity[i_trial][q][1]   +
                                         dphi_phys_velocity[i_test][q][2] * dphi_phys_velocity[i_trial][q][2] ) * w_quad[q];
                        
                        integral_01 += dphi_phys_velocity[i_test][q][1] * dphi_phys_velocity[i_trial][q][0] * w_quad[q];
                        
                        integral_02 += dphi_phys_velocity[i_test][q][2] * dphi_phys_velocity[i_trial][q][0] * w_quad[q];
                        
                        integral_10 += dphi_phys_velocity[i_test][q][0] * dphi_phys_velocity[i_trial][q][1] * w_quad[q];
                        
                        integral_11 += ( 2*dphi_phys_velocity[i_test][q][1] * dphi_phys_velocity[i_trial][q][1] +
                                        dphi_phys_velocity[i_test][q][0] * dphi_phys_velocity[i_trial][q][0]   +
                                        dphi_phys_velocity[i_test][q][2] * dphi_phys_velocity[i_trial][q][2] ) * w_quad[q];
                        
                        integral_12 += dphi_phys_velocity[i_test][q][2] * dphi_phys_velocity[i_trial][q][1] * w_quad[q];
                        
                        integral_20 += dphi_phys_velocity[i_test][q][0] * dphi_phys_velocity[i_trial][q][2] * w_quad[q];
                        
                        integral_21 += dphi_phys_velocity[i_test][q][1] * dphi_phys_velocity[i_trial][q][2] * w_quad[q];
                        
                        integral_22 += ( 2*dphi_phys_velocity[i_test][q][2] * dphi_phys_velocity[i_trial][q][2] +
                                        dphi_phys_velocity[i_test][q][0] * dphi_phys_velocity[i_trial][q][0]   +
                                        dphi_phys_velocity[i_test][q][1] * dphi_phys_velocity[i_trial][q][1] ) * w_quad[q];
                        
                    }
                    M_vals_supg[i_elem][0][0][i_test][i_trial] = integral_00 * M_detJacobian[i_elem];
                    M_vals_supg[i_elem][0][1][i_test][i_trial] = integral_01 * M_detJacobian[i_elem];
                    M_vals_supg[i_elem][0][2][i_test][i_trial] = integral_02 * M_detJacobian[i_elem];
                    M_vals_supg[i_elem][1][0][i_test][i_trial] = integral_10 * M_detJacobian[i_elem];
                    M_vals_supg[i_elem][1][1][i_test][i_trial] = integral_11 * M_detJacobian[i_elem];
                    M_vals_supg[i_elem][1][2][i_test][i_trial] = integral_12 * M_detJacobian[i_elem];
                    M_vals_supg[i_elem][2][0][i_test][i_trial] = integral_20 * M_detJacobian[i_elem];
                    M_vals_supg[i_elem][2][1][i_test][i_trial] = integral_21 * M_detJacobian[i_elem];
                    M_vals_supg[i_elem][2][2][i_test][i_trial] = integral_22 * M_detJacobian[i_elem];
                }
            }
            
            for ( dim_mat = 0; dim_mat < 3 ; dim_mat++ )
            {
                // DOF - test
                for ( i_test = 0; i_test <  ndof_velocity; i_test++ )
                {
                    // DOF - trial
                    for ( i_trial = 0; i_trial < ndof_pressure; i_trial++ )
                    {
                        M_vals_supg_01[i_elem][dim_mat][i_test][i_trial] = 0.0;
                        M_cols_pressure[i_elem][i_trial] = M_elements_pressure[i_elem][i_trial];
                        integral = 0.0;
                        // QUAD
                        for ( q = 0; q < NumQuadPoints ; q++ )
                        {
                            integral += dphi_phys_velocity[i_test][q][dim_mat]*M_phi_pressure[i_trial][q]*w_quad[q];
                        }
                        M_vals_supg_01[i_elem][dim_mat][i_test][i_trial] = -1.0 * integral * M_detJacobian[i_elem];
                    }
                }
            }
            
            for ( dim_mat = 0; dim_mat < 3 ; dim_mat++ )
            {
                // DOF - test
                for ( i_test = 0; i_test <  ndof_pressure; i_test++ )
                {
                    M_rows_pressure[i_elem][i_test] = M_elements_pressure[i_elem][i_test];
                    
                    // DOF - trial
                    for ( i_trial = 0; i_trial < ndof_velocity; i_trial++ )
                    {
                        M_vals_supg_10[i_elem][dim_mat][i_test][i_trial] = 0.0;
                        integral = 0.0;
                        // QUAD
                        for ( q = 0; q < NumQuadPoints ; q++ )
                        {
                            integral += dphi_phys_velocity[i_trial][q][dim_mat]*M_phi_pressure[i_test][q]*w_quad[q];
                        }
                        
                        M_vals_supg_10[i_elem][dim_mat][i_test][i_trial] = integral * M_detJacobian[i_elem];
                    }
                }
            }
        }
    }
    
    for ( int k = 0; k < M_numElements; ++k )
    {
        mass->matrixPtr()->InsertGlobalValues ( ndof_velocity, M_rows_velocity[k],
                                                ndof_velocity, M_cols_velocity[k],
                                                M_vals_00[k],
                                                Epetra_FECrsMatrix::ROW_MAJOR);
        
        b01->matrixPtr()->InsertGlobalValues ( ndof_velocity,
                                               M_rows_velocity[k],
                                               ndof_pressure,
                                               M_cols_pressure[k],
                                               M_vals_supg_01[k][0],
                                               Epetra_FECrsMatrix::ROW_MAJOR);
        
        b10->matrixPtr()->InsertGlobalValues ( ndof_pressure, M_rows_pressure[k],
                                               ndof_velocity, M_cols_velocity[k],
                                               M_vals_supg_10[k][0],
                                               Epetra_FECrsMatrix::ROW_MAJOR);
        
    }
    
    for ( UInt d1 = 0; d1 < 3 ; d1++ ) // row index
    {
        for ( int k = 0; k < M_numElements; ++k )
        {
            for ( UInt i = 0; i <  ndof_velocity; i++ )
            {
                M_rows_tmp[k][i] = M_rows_velocity[k][i] + d1 * M_numScalarDofs;
            }
            stiffness->matrixPtr()->InsertGlobalValues ( ndof_velocity, M_rows_tmp[k],
                                                         ndof_velocity, M_cols_velocity[k],
                                                         M_vals_supg[k][d1][0],
                                                         Epetra_FECrsMatrix::ROW_MAJOR);

        }
        
        for ( int k = 0; k < M_numElements; ++k )
        {
            for ( UInt i = 0; i <  ndof_velocity; i++ )
            {
                M_cols_tmp[k][i] = M_cols_velocity[k][i] + M_numScalarDofs;
            }
            stiffness->matrixPtr()->InsertGlobalValues ( ndof_velocity, M_rows_tmp[k],
                                                         ndof_velocity, M_cols_tmp[k],
                                                         M_vals_supg[k][d1][1],
                                                         Epetra_FECrsMatrix::ROW_MAJOR);
        }
        
        for ( int k = 0; k < M_numElements; ++k )
        {
            for ( UInt i = 0; i <  ndof_velocity; i++ )
            {
                M_cols_tmp[k][i] = M_cols_velocity[k][i] + 2 * M_numScalarDofs;
            }
            stiffness->matrixPtr()->InsertGlobalValues ( ndof_velocity, M_rows_tmp[k],
                                                         ndof_velocity, M_cols_tmp[k],
                                                         M_vals_supg[k][d1][2],
                                                         Epetra_FECrsMatrix::ROW_MAJOR);
        }
    }
    
    for ( UInt d1 = 1; d1 < 3 ; d1++ )
    {
        for ( int k = 0; k < M_numElements; ++k )
        {
            for ( UInt i = 0; i <  ndof_velocity; i++ )
            {
                M_rows_velocity[k][i] += M_numScalarDofs;
                M_cols_velocity[k][i] += M_numScalarDofs;
            }
            mass->matrixPtr()->InsertGlobalValues ( ndof_velocity, M_rows_velocity[k],
                                                    ndof_velocity, M_cols_velocity[k],
                                                    M_vals_00[k],
                                                    Epetra_FECrsMatrix::ROW_MAJOR);
            
            b01->matrixPtr()->InsertGlobalValues ( ndof_velocity, M_rows_velocity[k],
                                                   ndof_pressure, M_cols_pressure[k],
                                                   M_vals_supg_01[k][d1],
                                                   Epetra_FECrsMatrix::ROW_MAJOR);
            
            b10->matrixPtr()->InsertGlobalValues ( ndof_pressure, M_rows_pressure[k],
                                                   ndof_velocity, M_cols_velocity[k],
                                                   M_vals_supg_10[k][d1],
                                                   Epetra_FECrsMatrix::ROW_MAJOR);
            
        }
    }
}
