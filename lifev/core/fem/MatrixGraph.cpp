#include <lifev/core/fem/MatrixGraph.hpp>
#include <chrono>
using namespace std::chrono;

using namespace std;

using namespace LifeV;

//=========================================================================
// Constructor
MatrixGraph::MatrixGraph( const meshPtr_Type& mesh, const commPtr_Type& comm, const ReferenceFE* refFE ) :
		M_mesh ( mesh ),
		M_comm ( comm ),
		M_referenceFE ( refFE )
{
}
//=========================================================================
// Destructor //
MatrixGraph::~MatrixGraph()
{
	//---------------------

	for( int i = 0 ; i < M_numElements ; i++ )
	{
		delete [] M_elements[i];
	}
	delete [] M_elements;

	//---------------------

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
MatrixGraph::buildGraph ( const int& numElements, CurrentFE* fe, const fespacePtr_Type& fespace, graphPtr_Type& matrix_graph )
{
	M_numElements = numElements;

	M_numScalarDofs = fespace->dof().numTotalDof();

	//-------------------------------------------------------------------------------------------------
	// ALLOCATE MEMORY
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
            M_elements[i][j] = fespace->dof().localToGlobalMap (i, j);
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

    //-------------------------------------------------------------------------------------------------
    // BUILD GRAPH
    //-------------------------------------------------------------------------------------------------

    matrix_graph.reset (new graph_Type (Copy, *(fespace->map().map (Unique) ), 0 ) );

    int ndof = M_referenceFE->nbDof();

    for ( int i_elem = 0; i_elem <  M_numElements ; i_elem++ )
    {
    	// DOF - test
    	for ( int i_test = 0; i_test <  ndof; i_test++ )
    	{
    		M_rows[i_elem][i_test] = M_elements[i_elem][i_test];

    		// DOF - trial
    		for ( int i_trial = 0; i_trial <  ndof; i_trial++ )
    		{
    			M_cols[i_elem][i_trial] = M_elements[i_elem][i_trial];
    		}
    	}

    	matrix_graph->InsertGlobalIndices ( ndof, M_rows[i_elem], ndof, M_cols[i_elem] );

    	for ( UInt d1 = 1; d1 < 3 ; d1++ )
    	{
    		for ( UInt i = 0; i <  ndof; i++ )
    		{
    			M_rows[i_elem][i] += M_numScalarDofs;
    			M_cols[i_elem][i] += M_numScalarDofs;
    		}
    		matrix_graph->InsertGlobalIndices ( ndof, M_rows[i_elem], ndof, M_cols[i_elem] );
    	}
    }

    matrix_graph->GlobalAssemble();

}
