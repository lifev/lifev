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
    @brief Tutorial for laplacian of trial functions.

    @author Davide Forti <davide.forti@epfl.ch>
    @date 2013-05

    Test to check (and to show) how the laplacian pf the test
    and trial functions is computed

 */

// ---------------------------------------------------------------
// We reuse the same files as in the previous tutorial. We add the
// chrono to measure the timings.
// ---------------------------------------------------------------


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <lifev/core/util/LifeChrono.hpp>

#include <boost/shared_ptr.hpp>

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;

Real uOnes ( const Real& /* t */, const Real&  x , const Real&  y , const Real& z , const ID& i )
{
    return 1.0;
}

Real uFunction ( const Real& /* t */, const Real&  x , const Real&  y , const Real& z , const ID& i )
{
    if ( i == 0 )
    	return x*x*0.5;
    else if ( i == 1)
    	return y*y*0.5;
    else
    	return z*z*0.5;
}

int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    const bool verbose (Comm->MyPID() == 0);

    const UInt Nelements (10);

    boost::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);

    Real length = 3.0;

    regularMesh3D ( *fullMeshPtr, 1, Nelements, Nelements, Nelements, false,
    		        length,   length,   length,
                    0.0,  0.0,  0.0);

    MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, Comm);
    boost::shared_ptr< mesh_Type > meshPtr (meshPart.meshPartition() );

    fullMeshPtr.reset();

    std::string uOrder ("P2");

    ///////////////////////////////////
    // Testing the scalar field case //
    ///////////////////////////////////

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace
    ( new FESpace< mesh_Type, MapEpetra > (meshPtr, uOrder, 1, Comm) );

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > ETuSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPtr, & (uSpace->refFE() ), & (uSpace->fe().geoMap() ), Comm) );

    vector_Type myFEvector (uSpace->map(), Unique);
    uSpace->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (uFunction), myFEvector, 0.0);
    vector_Type myFEvectorRepeated (myFEvector, Repeated);

    vector_Type vectorOnes (uSpace->map(), Unique);
    uSpace->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (uOnes), vectorOnes, 0.0);

    vector_Type integral (uSpace->map() );

    {
    	using namespace ExpressionAssembly;

    	integrate ( elements (ETuSpace->mesh() ),
    				uSpace->qr(),
    				ETuSpace,
    				laplacian(ETuSpace, myFEvectorRepeated) * phi_i
    			  )
    		>> integral;
    }

    integral.globalAssemble();

    Real result = 0.0;

    result = integral.dot(vectorOnes);

    std::cout << "\n\nSCALAR CASE " << std::endl;
    std::cout << "\nThe volume is = " << length*length*length << std::endl;
    std::cout << "\nThe result is = " << result << std::endl;
    std::cout << "\nThe error is = " << result-(length*length*length) << std::endl;

    ///////////////////////////////////
    // Testing the vector field case //
    ///////////////////////////////////

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpaceVec
    ( new FESpace< mesh_Type, MapEpetra > (meshPtr, uOrder, 3, Comm) );

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 3 > > ETuSpaceVec
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 3 > (meshPtr, & (uSpaceVec->refFE() ), & (uSpaceVec->fe().geoMap() ), Comm) );

    vector_Type myFEvectorVec (uSpaceVec->map(), Unique);
    uSpaceVec->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (uFunction), myFEvectorVec, 0.0);
    vector_Type myFEvectorRepeatedVec (myFEvectorVec, Repeated);

    vector_Type vectorOnesVec (uSpaceVec->map(), Unique);
    uSpaceVec->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (uOnes), vectorOnesVec, 0.0);

    vector_Type integralVec (uSpaceVec->map() );

    {
    	using namespace ExpressionAssembly;

    	integrate ( elements (ETuSpaceVec->mesh() ),
    			    uSpaceVec->qr(),
    			    ETuSpaceVec,
    			    dot ( laplacian(ETuSpaceVec, myFEvectorVec), phi_i )
    			  )
    		>> integralVec;
    }

    integralVec.globalAssemble();

    result = 0.0;

    result = integralVec.dot(vectorOnesVec);

    std::cout << "\n\nVECTORIAL CASE " << std::endl;
    std::cout << "\nThe volume is = " << length*length*length << std::endl;
    std::cout << "\nThe result is = " << result << std::endl;
    std::cout << "\nThe error is = " << result-(length*length*length) << "\n\n";

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}


