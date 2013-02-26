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
    @brief Example usage of MeshEntityContainer features

    @author Antonio Cervone <ant.cervone@gmail.com>
    @contributor
    @maintainer Antonio Cervone <ant.cervone@gmail.com>

    @date 07-08-2012

    Extract elements and facets from a mesh based on a functor

 */

// ===================================================
//! Includes
// ===================================================
// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

#include <lifev/core/filter/GetPot.hpp>

#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/RegionMesh2DStructured.hpp>
#include <lifev/core/mesh/MeshEntityContainer.hpp>

using namespace LifeV;

// put all local stuff in a local namespace
namespace
{

typedef RegionMesh<LinearTriangle>      mesh_Type;
typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

// interrogator that checks if the barycenter of the mesh entity satisfies the condition
// given in the constructor wrt the given circle
// (default: barycenter inside the circle)
template < typename MeshEntityType,
         typename ComparisonPolicyType = boost::function2 < bool,
         Real const&,
         Real const& > >
class CircleInterrogator
{
public:
    typedef MeshEntityType       meshEntity_Type;
    typedef ComparisonPolicyType comparisonPolicy_Type;

    CircleInterrogator ( Vector3D const& center,
                         Real radius,
                         comparisonPolicy_Type const& policy = std::less<Real>() )
        : M_center ( center ),
          M_radius ( radius ),
          M_policy ( policy ) {}

    bool operator() ( const meshEntity_Type& entity ) const
    {
        // compute the barycenter of the entity
        // (this should be a method of the object)
        Vector3D barycenter;
        for ( UInt k = 0; k < meshEntity_Type::S_numPoints; k++ )
        {
            barycenter += entity.point ( k ).coordinates();
        }
        barycenter /= meshEntity_Type::S_numPoints;

        // check if the distance between the barycenter and the center of the circle
        // satisfies the policy (default: distance less than the radius)
        return M_policy ( ( barycenter - M_center ).norm(), M_radius );
    }

private:
    const Vector3D M_center;
    const Real M_radius;
    const comparisonPolicy_Type M_policy;

}; // CircleInterrogator

// interrogator that checks if there is a change in sign at the points
// of the entity of the level set function that is passed to the constructor
template < typename MeshEntityType >
class LevelSetInterrogator
{
public:
    typedef MeshEntityType                          meshEntity_Type;
    typedef boost::function1<Real, Vector3D const&> distanceFunction_Type;

    explicit LevelSetInterrogator ( distanceFunction_Type const& distFun ) : M_distanceFunction ( distFun ) {}

    bool operator() ( meshEntity_Type const& entity ) const
    {
        // get the value of the level set in each point of the entity
        std::vector<Real> sign ( meshEntity_Type::S_numPoints );
        for ( UInt k = 0; k < meshEntity_Type::S_numPoints; k++ )
        {
            sign[ k ] = M_distanceFunction ( entity.point ( k ).coordinates() );
        }

        // if there is a change in sign the level set is crossing
        // try all combinations between points
        for ( UInt i = 0; i < meshEntity_Type::S_numPoints - 1; i++ )
            for ( UInt j = i + 1; j < meshEntity_Type::S_numPoints; j++ )
                if ( sign[ i ] * sign[ j ] < 0. )
                {
                    return true;
                }

        return false;
    }

private:
    distanceFunction_Type M_distanceFunction;

}; // LevelSetInterrogator

// the simplest form of level set is the signed distance to a line in the 2D Cartesian plane
class lineDistance
{
public:
    lineDistance ( Real a, Real b, Real c ) :
        M_a ( a ),
        M_b ( b ),
        M_c ( c ) {}

    // this is the distance of a generic point from the line: ( a*x + b*y - c ) / sqrt( a^2 + b^2 )
    Real operator() ( Vector3D const& point )
    {
        return ( M_a * point[ 0 ] + M_b * point[ 1 ] + M_c ) / std::sqrt ( M_a * M_a + M_b * M_b );
    }

private:
    Real const M_a;
    Real const M_b;
    Real const M_c;

}; // lineDistance

} // anonymous namespace

int main ( int argc, char** argv )
{
    // verbosity
    bool verbose = 1;

    // communicator
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    verbose = comm->MyPID() == 0;
#else
    boost::shared_ptr<Epetra_Comm> comm ( new Epetra_SerialComm );
#endif

    // number of elements for each direction in the mesh
    const UInt numMeshElem = 3;

    // init mesh
    meshPtr_Type mesh ( new mesh_Type ( comm ) );

    // build 2D mesh on the unit square
    regularMesh2D ( *mesh,
                    1,
                    numMeshElem, numMeshElem,
                    false,
                    1.0, 1.0,
                    0.0, 0.0 );

    // ===== TEST 1 =====
    {
        // get all elements that are inside the circle of center (0.5,0.5) and radius 0.25
        CircleInterrogator<mesh_Type::element_Type> myCircleInterrogator ( Vector3D ( 0.5, 0.5, 0.0 ), 0.25 );

        // get the number of entities that satisfy the predicate
        UInt numExtractedElements = mesh->elementList().countAccordingToPredicate ( myCircleInterrogator );
        if ( verbose ) std::cout << "the number of elements that stay inside the circle is "
                                     << numExtractedElements << std::endl;

        std::vector<mesh_Type::element_Type const*> extractedElements ( numExtractedElements ); // i love c++11 auto...

        // get a const pointer to those entities
        extractedElements = mesh->elementList().extractAccordingToPredicate ( myCircleInterrogator );
        if ( verbose )
        {
            std::cout << "the elements ids are: ";
        }
        for ( UInt k = 0; k < numExtractedElements; k++ )
        {
            if ( verbose )
            {
                std::cout << extractedElements[ k ]->id() << " ";
            }
        }
        if ( verbose )
        {
            std::cout << std::endl;
        }
    }
    // ==================

    // ===== TEST 2 =====
    {
        // define a line that crosses the mesh
        // ( the case when the level set passes exactly through a point is left to the reader )
        lineDistance myLine ( 1., 1., -1.1 );

        // create the interrogator that uses the line defined above
        LevelSetInterrogator<mesh_Type::facet_Type> myLevelSetInterrogator ( myLine ); // here a lambda function would reduce code

        // get the number of entities that satisfy the predicate
        UInt numExtractedFacets = mesh->facetList().countAccordingToPredicate ( myLevelSetInterrogator );
        if ( verbose ) std::cout << "the number of facets that are crossed by the level set is "
                                     << numExtractedFacets << std::endl;

        std::vector<mesh_Type::facet_Type const*> extractedFacets ( numExtractedFacets ); // i love c++11 auto...

        // get a const pointer to those entities
        extractedFacets = mesh->facetList().extractAccordingToPredicate ( myLevelSetInterrogator );
        if ( verbose )
        {
            std::cout << "the facets ids are: ";
        }
        for ( UInt k = 0; k < numExtractedFacets; k++ )
        {
            if ( verbose )
            {
                std::cout << extractedFacets[ k ]->id() << " ";
            }
        }
        if ( verbose )
        {
            std::cout << std::endl;
        }
    }
    // ==================

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}
