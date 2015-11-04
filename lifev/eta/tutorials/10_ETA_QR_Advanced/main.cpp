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
    @brief Advanced tutorial about the quadrature in the ETA.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 08-08-2012

    This tutorial shows an advanced feature of the ETA framework:
    adaptative quadratures.

    Indeed, the ETA framework allows the user to define different
    ways to integrate, that have possibly not been forseen before,
    without having to enter the whole machinery.

    In this tutorial, we show how to create an "adaptative"
    quadrature rule which consists in a standard quadrature
    rule excepted for some precise elements where the integral
    is known to be discontinuous. In the elements crossed by this
    discontinuity, we implement a Monte-Carlo type integral
    (which is probably not the best way to compute the integral,
    but for this example, this is fine).

    Please, read tutorial 1, 3, 6 and 9 before reading this tutorial.
 */

// ---------------------------------------------------------------
// We include here exactly the same files that in tutorial 1,
// excepted for the QRAdapterBase.hpp file and the cstdlib, which
// is needed for the generation of the random numbers.
// ---------------------------------------------------------------

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/eta/fem/QRAdapterBase.hpp>

#include <boost/shared_ptr.hpp>

#include <cstdlib>

using namespace LifeV;

// ---------------------------------------------------------------
// We start directly with the definition of the class that we
// will need for the adaptative quadrature.
//
// Remark that, for the sake of this tutorial, we implement all
// the methods "inline", while it should be located in a separated
// .hpp file.
//
// First of all, the class is template on the mesh type,
// which is quite standard for the LifeV class, to avoid
// depending on a specific type of mesh (even if we suppose that
// it contains only tetrahedra).
//
// Then, we make the class derive from the QRAdapterBase class,
// which takes in argument exactly the name of the class. This is
// required by the ETA framework to recognize that this object is
// used to define a quadrature.
// ---------------------------------------------------------------

template< typename MeshType >
class MCQuadrature : public QRAdapterBase< MCQuadrature <MeshType> >
{

private:

    // ---------------------------------------------------------------
    // A small typedef for convenience.
    // ---------------------------------------------------------------

    typedef QRAdapterBase< MCQuadrature<MeshType> > base_Type;

    // ---------------------------------------------------------------
    // We start the description of the class by the members that it
    // will contains.
    //
    // This class contains:
    // - a quadrature rule : the unadapted quadrature
    // - a pointer to a quadrature rule :  the adapted quadrature
    // - a pointer to the mesh
    // - a boolean indicating whether to use the adapted quadrature
    //   or not.
    // - an integer representing the number of Monte-Carlo points
    // that we want for the quadrature.
    // ---------------------------------------------------------------

    QuadratureRule M_stdQR;

    std::shared_ptr<QuadratureRule> M_adaptedQR;

    std::shared_ptr<MeshType> M_mesh;

    bool M_isAdaptedElement;

    UInt M_nbPoints;

public:

    // ---------------------------------------------------------------
    // There is a fixed interface that must be present for the class
    // to work properly. They are:
    // - a copy constructor
    // - the update method
    // - the isAdaptedElement getter
    // - the standardQR getter
    // - the adaptedQR getter
    //
    // These elements are mandatory, other interfaces can of course
    // be added! We here add another constructor.
    // ---------------------------------------------------------------

    // ---------------------------------------------------------------
    // Simple copy constructor
    // ---------------------------------------------------------------

    MCQuadrature ( const MCQuadrature<MeshType>& qr)
        : base_Type(),
          M_stdQR (qr.M_stdQR),
          M_adaptedQR (new QuadratureRule (* (qr.M_adaptedQR) ) ),
          M_mesh (qr.M_mesh),
          M_isAdaptedElement (qr.M_isAdaptedElement),
          M_nbPoints (qr.M_nbPoints)
    {}

    // ---------------------------------------------------------------
    // This method should compute whether the element needs an
    // adapted quadrature or not, and if yes, it compute the
    // quadrature.
    //
    // We delegate here the work to a private method.
    // ---------------------------------------------------------------

    void update (UInt elementID)
    {
        computeQuadRule (elementID);
    }


    // ---------------------------------------------------------------
    // Says if the element (for which the update has been called)
    // needs an adapted quadrature.
    // ---------------------------------------------------------------

    bool isAdaptedElement() const
    {
        return M_isAdaptedElement;
    }


    // ---------------------------------------------------------------
    // Getter for the non adapted QR.
    // ---------------------------------------------------------------

    const QuadratureRule& standardQR() const
    {
        return M_stdQR;
    }


    // ---------------------------------------------------------------
    // Getter for the adapted QR.
    // ---------------------------------------------------------------

    const QuadratureRule& adaptedQR() const
    {
        return *M_adaptedQR;
    }


    // ---------------------------------------------------------------
    // Our special constructor, which uses the mesh, the non-adapted
    // quadrature rule and the number of Monte-Carlo points.
    // ---------------------------------------------------------------

    MCQuadrature (std::shared_ptr<MeshType> mesh, const QuadratureRule& stdQR, const UInt& nbPoints)
        : base_Type(),
          M_stdQR (stdQR),
          M_adaptedQR (new QuadratureRule (stdQR) ),
          M_mesh (mesh),
          M_isAdaptedElement (false),
          M_nbPoints (nbPoints)
    {}

private:

    // ---------------------------------------------------------------
    // Here is the core of the class: the computation of the
    // quadrature rule if needed.
    //
    // The principle is to check whether the element is crossed by
    // the discontinuity.
    // - if yes, generate a Monte-Carlo quadrature
    // - if no, use a standard quadrature rule.
    // ---------------------------------------------------------------

    void computeQuadRule (UInt elementID)
    {

        // ---------------------------------------------------------------
        // We know that we want to adapt the quadrature if the element
        // crosses the plan given by x-2y=0. To know
        // we an element is crossed, we check the sign of x-2y (which
        // is positive on one side and negative on the other). If we
        // get both signes, we need to adapt the quadrature.
        // ---------------------------------------------------------------

        bool positiveVertices (false);
        bool negativeVertices (false);

        for (UInt iterVertex (0); iterVertex < 4; ++iterVertex)
        {
            Real x (M_mesh->element (elementID).point (iterVertex).coordinate (0) );
            Real y (M_mesh->element (elementID).point (iterVertex).coordinate (1) );

            if (x - 2 * y < 0)
            {
                negativeVertices = true;
            }
            else if (x - 2 * y > 0)
            {
                positiveVertices = true;
            };
            // If a vertex gives zero, we do not care about it!
        }

        M_isAdaptedElement = ( positiveVertices && negativeVertices );

        // ---------------------------------------------------------------
        // Now, if the element is to be adapted, we generate the Monte-
        // Carlo quadrature.
        // ---------------------------------------------------------------

        if (M_isAdaptedElement)
        {
            M_adaptedQR.reset (new QuadratureRule ("Monte-Carlo", TETRA, 3, 0, 0) );

            UInt nbPts (0);
            while ( nbPts < M_nbPoints)
            {
                Real x ( std::rand() / Real (RAND_MAX) );
                Real y ( std::rand() / Real (RAND_MAX) );
                Real z ( std::rand() / Real (RAND_MAX) );

                if ( x + y + z <= 1)
                {
                    M_adaptedQR->addPoint (QuadraturePoint (x, y, z, 1.0 / (6.0 * M_nbPoints) ) );
                    nbPts += 1;
                }
            } // end of while
        } // end of if
    }
};


// ---------------------------------------------------------------
// We define then the function to integrate: it is the indicator
// function of (x-2y>0).
// ---------------------------------------------------------------


class fFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() (const VectorSmall<3>& position)
    {
        if (position[0] - 2 * position[1] > 0)
        {
            return 1;
        }
        return 0;
    }

    fFunctor() {}
    fFunctor (const fFunctor&) {}
    ~fFunctor() {}
};


// ---------------------------------------------------------------
// We proceed then exactly has described in tutorial 1 until the
// assembly procedure is reached.
// ---------------------------------------------------------------


typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;


int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    const bool verbose (Comm->MyPID() == 0);


    if (verbose)
    {
        std::cout << " -- Building and partitioning the mesh ... " << std::flush;
    }

    const UInt Nelements (10);

    std::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);

    regularMesh3D ( *fullMeshPtr, 1, Nelements, Nelements, Nelements, false,
                    2.0,   2.0,   2.0,
                    -1.0,  -1.0,  -1.0);

    MeshPartitioner< mesh_Type >  meshPart (fullMeshPtr, Comm);

    fullMeshPtr.reset();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // First of all, we compute the integral using a standard
    // quadrature rule. We see that we do not get the correct answer
    // (which is 4)
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Computing the volume (std) ... " << std::flush;
    }

    Real localIntegralStd (0.0);

    {
        using namespace ExpressionAssembly;

        std::shared_ptr<fFunctor> fFct (new fFunctor);

        integrate (  elements (meshPart.meshPartition() ),
                     quadRuleTetra1pt,
                     eval (fFct, X)
                  )
                >> localIntegralStd;
    }

    Real globalIntegralStd (0.0);

    Comm->Barrier();
    Comm->SumAll (&localIntegralStd, &globalIntegralStd, 1);


    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " Integral: " << globalIntegralStd << std::endl;
    }


    // ---------------------------------------------------------------
    // Now we use our adapted quadrature. We simply replace the
    // quadrature by an instance of the quadrature adapter that we
    // created.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Computing the volume (adapted) ... " << std::flush;
    }

    Real localIntegralAdapted (0.0);

    {
        using namespace ExpressionAssembly;

        std::shared_ptr<fFunctor> fFct (new fFunctor);

        MCQuadrature<mesh_Type> myQuad (meshPart.meshPartition(), quadRuleTetra1pt, 1e4);

        integrate (  elements (meshPart.meshPartition() ),
                     myQuad,
                     eval (fFct, X)
                  )
                >> localIntegralAdapted;
    }

    Real globalIntegralAdapted (0.0);

    Comm->Barrier();
    Comm->SumAll (&localIntegralAdapted, &globalIntegralAdapted, 1);


    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " Integral: " << globalIntegralAdapted << std::endl;
    }

    // ---------------------------------------------------------------
    // Finally, we check that the errors are not too high.
    // ---------------------------------------------------------------


#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    if ( std::abs (4 - globalIntegralAdapted) < 1e-2 )
    {
        return ( EXIT_SUCCESS );
    }
    return ( EXIT_FAILURE );

}


