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
    @brief Tutorial for gradient interpolation.

    @author Luca Pasquale <luca.pasquale@mail.polimi.it>
    @date 2013-05

    In this tutorial we use two expressions available when
    integrating with Expression Templates:
    - interpoaltion of the gradient of a given function
    - multiplication by a MatrixSmall
    We do this by multiplying the gradient by the identity matrix
    (i.e. we're integrating the divergence)

    Tutorials that should be read before: 1,3

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

Real uOne ( const Real& /* t */, const Real&  x , const Real&  y , const Real& z , const ID& i )
{
	return 1.0;
}

Real uTestFunctions ( const Real& /* t */, const Real&  x , const Real&  y , const Real& z , const ID& i )
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

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace ( new FESpace< mesh_Type, MapEpetra > (meshPtr, uOrder, 1, Comm) );

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > ETuSpace ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPtr, & (uSpace->refFE() ), & (uSpace->fe().geoMap() ), Comm) );

    vector_Type vectorTestFunctions (uSpace->map(), Unique);
    uSpace->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (uTestFunctions), vectorTestFunctions, 0.0);
    vector_Type vectorTestFunctionsRepeated (vectorTestFunctions, Repeated);

    vector_Type vectorOnes (uSpace->map(), Repeated);
    vectorOnes.zero();
    vectorOnes += 1;

    vector_Type integral (uSpace->map() );

    {
    	using namespace ExpressionAssembly;

    	integrate ( elements (ETuSpace->mesh() ),
    				uSpace->qr(),
    				ETuSpace,
    				value(ETuSpace,vectorOnes) * laplacian(phi_i)
    			  )
    		>> integral;
    }

    integral.globalAssemble();

    Real result = 0.0;

    result = integral.dot(vectorTestFunctions);

    if ( Comm->MyPID() == 0 )
    {
		std::cout << "\n\nSCALAR CASE " << std::endl;
		std::cout << "\nThe volume is = " << length*length*length << std::endl;
		std::cout << "\nThe result is = " << result << std::endl;
		std::cout << "\nThe error is = " << result-(length*length*length) << std::endl;
    }

    Real error_scalar = result-length*length*length;

    ///////////////////////////////////
    // Testing the vector field case //
    ///////////////////////////////////

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpaceVec
    ( new FESpace< mesh_Type, MapEpetra > (meshPtr, uOrder, 3, Comm) );

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 3 > > ETuSpaceVec
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 3 > (meshPtr, & (uSpaceVec->refFE() ), & (uSpaceVec->fe().geoMap() ), Comm) );

    vector_Type vectorTestFunctionsVec (uSpaceVec->map(), Unique);
    uSpaceVec->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (uTestFunctions), vectorTestFunctionsVec, 0.0);
    vector_Type vectorTestFunctionsRepeatedVec (vectorTestFunctionsVec, Repeated);

    vector_Type vectorOnesVec (uSpaceVec->map(), Repeated );
    vectorOnesVec.zero();
    vectorOnesVec += 1;

    vector_Type integralVec (uSpaceVec->map() );

    {
    	using namespace ExpressionAssembly;

    	integrate ( elements (ETuSpaceVec->mesh() ),
    			    uSpaceVec->qr(),
    			    ETuSpaceVec,
    			    dot ( value(ETuSpaceVec,vectorOnesVec), laplacian(phi_i) )
    			  )
    		>> integralVec;
    }

    integralVec.globalAssemble();

    result = 0.0;
    result = integralVec.dot(vectorTestFunctionsVec);

    if ( Comm->MyPID() == 0 )
    {
    	std::cout << "\n\nVECTORIAL CASE " << std::endl;
    	std::cout << "\nThe volume is = " << (length*length*length)*3 << std::endl;
    	std::cout << "\nThe result is = " << result << std::endl;
    	std::cout << "\nThe error is = " << result-((length*length*length)*3) << "\n\n";
    }

    Real error_vectorial = result-3*length*length*length;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    if ( std::fabs(error_scalar) < 1e-9 && std::fabs(error_vectorial) < 1e-9 )
    {
    	return ( EXIT_SUCCESS );
    }

    return ( EXIT_FAILURE );
}


