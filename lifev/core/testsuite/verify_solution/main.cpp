/*
 * VerifySolutions_test.cpp
 *
 *  Created on: May 10, 2014
 *      Author: Simone Deparis
 */

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/LifeV.hpp>

#include "lifev/core/util/VerifySolutions.hpp"
#include "lifev/core/array/MapEpetra.hpp"
#include "lifev/core/array/VectorEpetra.hpp"

#include <cmath>

// ===================================================
//! Namespaces
// ===================================================
using namespace LifeV;

/* Matlab/Octave code:
  map = [0:99];
  v1 = cos(map);
  v2 = sin(map);
  v3 = 1./(map+1);
  v = (v1+v2+v3)/3;
  norm(v) % => 3.4180
  u1 = v1-v; u2 = v2-v; u3 = v3-v;
  C = [v1*v1' v2*v1' v3*v1';
  v1*v2' v2*v2' v3*v2';
  v1*v3' v2*v3' v3*v3';];
  % => C =
    49.9880    0.3006    0.9148
    0.3006   50.0120    0.5379
    0.9148    0.5379    1.6350
 */
// ===================================================
//! Main
// ===================================================
int main (int argc, char** argv)
{

#ifdef HAVE_MPI

    MPI_Init ( &argc, &argv );

    std::cout << "MPI Initialization" << std::endl;

    bool testSuccess(true);

#endif
    {
#ifdef EPETRA_MPI
    boost::shared_ptr<Epetra_Comm>   comm( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    boost::shared_ptr<Epetra_Comm>   comm( new Epetra_SerialComm() );
#endif

    // Build Distributed linear map based on [0:99]

    Int numGlobalElements(100);
    MapEpetra linearMap(numGlobalElements, 0, comm);

    // Build Distributed vectors [cos(index)] [sin(index)] [1/(index+1)]
    VectorEpetra v1(linearMap), v2(linearMap), v3(linearMap);

    Int numLocalEntries(linearMap.map(Unique)->NumMyElements () );

    Int* globalElements (linearMap.map(Unique)->MyGlobalElements () );
    Int globalIndex(0);
    for (int i(0); i < numLocalEntries; ++i)
    {
        globalIndex = globalElements[i];
        v1.setCoefficient(globalIndex, cos(globalIndex));
        v2.setCoefficient(globalIndex, sin(globalIndex));
        v3.setCoefficient(globalIndex, 1./(globalIndex+1));
    }

    Epetra_SerialDenseMatrix refM(3,3);
    refM(0,0) = 4.998801864473413e+01; refM(0,1) = 3.006425761130278e-01; refM(0,2) = 9.147988589055555e-01;
    refM(1,0) = 3.006425761130278e-01; refM(1,1) = 5.001198135526587e+01; refM(1,2) = 5.379437373113811e-01;
    refM(2,0) = 9.147988589055555e-01; refM(2,1) = 5.379437373113811e-01; refM(2,2) = 1.634983900184893e+00;

    Real referenceMeanNorm (3.417955110635026e+00);

    // Insert vectors into VerifySolution
    VerifySolutions verify;
    // Tolerance between the error and the errorKnown
    const LifeV::Real tolerance ( 1e-8 );

    verify.PushBack(v1);
    verify.PushBack(v2);
    verify.PushBack(v3);

    // Compute the Mean and Correlation matrix
    verify.ComputeCorrelation();

    // Verify that it is the same as pre-computed.
    bool isMeanOk   ( verify.Check( referenceMeanNorm, tolerance ) );
    bool isMatrixOk ( verify.Check( refM, tolerance ) );

    testSuccess = isMeanOk && isMatrixOk;

    }

#ifdef HAVE_MPI

    MPI_Finalize();

    std::cout << "MPI Finalization" << std::endl;

#endif

    if ( testSuccess )
    {
        return ( EXIT_SUCCESS );
    }
    else
    {
        return ( EXIT_FAILURE );
    }
}




