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
    @brief Utilities

    @contributor Simone Palamara <palamara.simone@gmail.com>
    @maintainer Simone Palamara <palamara.simone@gmail.com>

    This file contains a set of base utilities used to applied current on a specified point.
 */

#ifndef ELECTROPHYSIOLOGYUTILITY_H
#define ELECTROPHYSIOLOGYUTILITY_H 1

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <Teuchos_ScalarTraitsDecl.hpp>
#include <time.h>       /* time */
#include <math.h>       /* floor */
#include <lifev/core/mesh/NeighborMarker.hpp>


namespace LifeV
{

// Predeclaration

namespace ElectrophysiologyUtility
{

typedef boost::unordered_set<ID> neighbors_Type;
typedef std::vector<neighbors_Type> neighborList_Type;

//! HeartUtility - A string parser grammar based on \c boost::spirit::qi
/*!
 *  @author(s) Simone Palamara
 *
 *  \c ElectrophysiologyUtility contains methods for applied current on a specified point.
 *
 */

//! @name Methods
//@{

//! Find closest point within radius and applied a constant current
/*!
 * @param point Vector of real containing the coordinates of the point within the radius
 * @param radius Radius used to find point.
 * @param appliedCurrentVector    Vector epetra containing the applied current.
 * @param valueAppliedCurrent    Value of the current to apply at the specified point.
 * @param fullMesh   Pointer to the mesh.
 * @param Comm   EpetraMpi comunicator.
 */
template<typename Mesh> inline void appliedCurrentClosestPointWithinRadius (std::vector<Real>& point, Real Radius, boost::shared_ptr<VectorEpetra> appliedCurrentVector, Real valueAppliedCurrent,  boost::shared_ptr< Mesh > fullMesh, boost::shared_ptr<Epetra_Comm> Comm )
{
    Int n = appliedCurrentVector -> epetraVector().MyLength();

    Int ids;
    for ( UInt i (0); i < n; i++)
    {
        int iGID = appliedCurrentVector -> blockMap().GID ( static_cast<EpetraInt_Type>(i) );
        Real px = fullMesh -> point ( iGID ).x();
        Real py = fullMesh -> point ( iGID ).y();
        Real pz = fullMesh -> point ( iGID ).z();

        Real distance = std::sqrt ( ( point[0] - px) * (point[0] - px)
                                    + ( point[1] - py) * (point[1] - py)
                                    + ( point[2] - pz) * (point[2] - pz) );
        if (distance <= Radius)
        {
            ids = iGID;
            Radius = distance;
        }


    }
    Real localRadius = Radius;
    Real globalRadius (0);
    Comm->Barrier();
    Comm->MinAll (&localRadius, &globalRadius, 1);

    if (globalRadius == localRadius)
    {
        appliedCurrentVector-> operator [] ( ids ) = valueAppliedCurrent;
    }

}

//! Find all the points within radius and applied a constant current
/*!
 * @param point Vector of real containing the coordinates of the point within the radius
 * @param radius Radius used to find point.
 * @param appliedCurrentVector    Vector epetra containing the applied current.
 * @param valueAppliedCurrent    Value of the current to apply at the specified point.
 * @param fullMesh   Pointer to the mesh.
 */
template<typename Mesh> inline void appliedCurrentPointsWithinRadius (std::vector<Real>& point, Real Radius, boost::shared_ptr<VectorEpetra> appliedCurrentVector, Real valueAppliedCurrent,  boost::shared_ptr< Mesh > fullMesh )
{
    Int n = appliedCurrentVector -> epetraVector().MyLength();

    std::vector<UInt> ids;
    for ( UInt i (0); i < n; i++)
    {
        int iGID = appliedCurrentVector -> blockMap().GID ( static_cast<EpetraInt_Type>(i) );
        Real px = fullMesh -> point ( iGID ).x();
        Real py = fullMesh -> point ( iGID ).y();
        Real pz = fullMesh -> point ( iGID ).z();

        Real distance = std::sqrt ( ( point[0] - px) * (point[0] - px)
                                    + ( point[1] - py) * (point[1] - py)
                                    + ( point[2] - pz) * (point[2] - pz) );
        if (distance <= Radius)
        {
            ids.push_back (iGID);
        }


    }

    for (int i (0); i < ids.size(); i++)
    {
        appliedCurrentVector -> operator [] ( ids.at (i) ) = valueAppliedCurrent;
    }
}

//! Find closest point to a given one in the mesh
/*!
 * @param point Vector of real containing the coordinates of the point within the radius
 * @param vec    Vector epetra used to recover Global id.
 * @param fullMesh   Pointer to the mesh.
 * @param Comm   EpetraMpi comunicator.
 * @return The ids of the closest point
 */
template<typename Mesh> inline UInt findClosestPoint (std::vector<Real>& point, boost::shared_ptr<VectorEpetra> vec,  boost::shared_ptr< Mesh > fullMesh, boost::shared_ptr<Epetra_Comm> Comm )
{
    Int n = vec -> epetraVector().MyLength();
    Real Radius = 100000.0;
    Int ids;
    for ( UInt i (0); i < n; i++)
    {
        Int iGID = vec -> blockMap().GID ( static_cast<EpetraInt_Type>(i) );
        Real px = fullMesh -> point ( iGID ).x();
        Real py = fullMesh -> point ( iGID ).y();
        Real pz = fullMesh -> point ( iGID ).z();

        Real distance = std::sqrt ( ( point[0] - px) * (point[0] - px)
                                    + ( point[1] - py) * (point[1] - py)
                                    + ( point[2] - pz) * (point[2] - pz) );
        if (distance <= Radius)
        {
            ids = iGID;
            Radius = distance;
        }


    }
    Real localRadius = Radius;
    Real globalRadius (0);
    Comm->Barrier();
    Comm->MinAll (&localRadius, &globalRadius, 1);
    if (globalRadius == localRadius)
    {
        return ids;
    }
    else
    {
        return -1;
    }

}

//! Collect all the ids of points with a given flag in a local vector for each processor
/*!
 * @param containerIds Vector containig all the ids of the points with a given flag
 * @param flag    Flag.
 * @param vec    Vector epetra used to recover Global id.
 * @param fullMesh   Pointer to the mesh.
 */
template<typename Mesh> inline void allIdsPointsWithGivenFlag (std::vector<ID>& containerIds, UInt flag, boost::shared_ptr<VectorEpetra> vec,  boost::shared_ptr< Mesh > fullMesh )
{
    for ( Int j (0); j < vec->epetraVector().MyLength() ; ++j )
    {
        if ( fullMesh -> point ( vec->blockMap().GID ( static_cast<EpetraInt_Type>(j) ) ).markerID() == flag )
        {
            containerIds.push_back (vec->blockMap().GID ( static_cast<EpetraInt_Type>(j) ) );
        }
    }
}

//! Select randomly a value in a given set and update the set by deleting the chosen value from the set
/*!
 * @param container vector containig all the value
 * @return Selected value
 */
template<typename Type> inline Type randomPointInSet (std::vector<Type>& container )
{
    /* generate uniform random number: */
    UInt number = std::floor ( (Teuchos::ScalarTraits<Real>::random() + 1.0) / 2.0 * container.size() );
    Type value = container[number];
    container.erase (container.begin() + number);
    return value;
}


//! Select randomly a point in a given set and its neighborhood and update the set by deleting the chosen value from the set
/*!
 * @param container vector containig all the value
 * @param selectedPoints selected points
 * @param neighbors cointainer of the neighbors
 */
inline void randomPointInSetAndNeighborhood (std::vector<ID>& container, std::vector<ID>& selectedPoints, neighborList_Type& neighbors )
{
    /* generate uniform random number: */
    UInt number = std::floor ( (Teuchos::ScalarTraits<Real>::random() + 1.0) / 2.0 * container.size() );
    selectedPoints.push_back (container[number]);
    neighbors_Type::iterator it = neighbors[container[number]].begin();
    for (int i = 0; i < neighbors[container[number]].size(); i++)
    {
        selectedPoints.push_back (*it);
        it++;
    }
    container.erase (container.begin() + number);
}



//! Select randomly N values in a given set and update the set by deleting the chosen N values from the set
/*!
 * @param container Set containig all the value
 * @param selectedValue vector containg the selected values
 * @param N number of values we want to select
 */
template<typename Type> inline void randomNPointsInSet (std::vector<Type>& container, std::vector<Type>& selectedValue, UInt N )
{
    /* initialize random seed: */
    Teuchos::ScalarTraits<Real>::seedrandom (time (NULL) );
    for (int i = 0; i < N; i++)
    {
        selectedValue.push_back (randomPointInSet<Type> (container) );
    }
}

//! Select randomly N values in a given set and update the set by deleting the chosen N values from the set
/*!
 * @param container Set containig all the value
 * @param selectedValue vector containg the selected values
 * @param N number of values we want to select
 * @param fullMesh mesh used to determine the neighborhood
 */
template<typename Mesh> inline void randomNPointsInSetAndNeighborhood (std::vector<ID>& container, std::vector<ID>& selectedValue, std::vector<Real>& activationTime, Real deltaT, UInt N,  Mesh& fullMesh, boost::shared_ptr<Epetra_Comm> Comm)
{
    /* initialize random seed: */
    Teuchos::ScalarTraits<Real>::seedrandom (time (NULL) );
    neighborList_Type listNeighborhood;
    createPointNeighbors<Mesh> (fullMesh, listNeighborhood);
    Real timePmj (Comm->MyPID() *deltaT);
    for (int i = 0; i < N; i++)
    {
        std::vector<ID> localValue;
        randomPointInSetAndNeighborhood (container, localValue, listNeighborhood);
        for (int j = 0; j < localValue.size(); j++)
        {
            selectedValue.push_back (localValue[j]);
            activationTime.push_back (timePmj);
        }
    }
}



//! Apply given current to a set of points
/*!
 * @param container Set containig all the globalIDs of points to which we want to apply the current
 * @param appliedCurrentVector vector containg the current
 * @param valueAppliedCurrent value of the current we want to apply
 */
template<typename Type> inline void applyCurrentGivenSetOfPoints (std::vector<Type>& container, boost::shared_ptr<VectorEpetra> appliedCurrentVector, Real valueAppliedCurrent )
{
    for (int i (0); i < container.size(); i++)
    {
        Int lrow = appliedCurrentVector ->blockMap().LID ( static_cast<EpetraInt_Type>(container[i]) );
        if (lrow >= 0)
        {
            appliedCurrentVector -> operator [] ( container[i] ) = valueAppliedCurrent;
        }
    }
}


//! Apply given current with a set of PMJs.
/*!
 * @param container vector containig all the globalIDs of the PMJs to which we want to apply the current
 * @param timeActivation vector containig the activation time of the PMJs
 * @param appliedCurrentVector vector containg the applied current
 * @param valueAppliedCurrent vector containig the apply current during the all simulation time
 * @param shiftVector vector containig the index that refers to the value of applied current in each PMJ
 * @param currentTime current time
 */
template<typename Type> inline void applyCurrentGivenSetOfPoints (std::vector<Type>& container, std::vector<Real> activationTime, boost::shared_ptr<VectorEpetra> appliedCurrentVector, std::vector<Real>& valueAppliedCurrent, std::vector<UInt>& shiftVector, Real currentTime )
{
    for (int i (0); i < container.size(); i++)
    {
        if (activationTime[i] <= currentTime)
        {
            shiftVector[i] += 1;
            Int lrow = appliedCurrentVector ->blockMap().LID ( static_cast<EpetraInt_Type>(container[i]) );
            if (lrow >= 0)
            {
                appliedCurrentVector -> operator [] ( container[i] ) = valueAppliedCurrent[shiftVector[i]] / 85.7;
            }
        }
        else
        {
            Int lrow = appliedCurrentVector ->blockMap().LID ( static_cast<EpetraInt_Type>(container[i]) );
            if (lrow >= 0)
            {
                appliedCurrentVector -> operator [] ( container[i] ) = 0;
            }
        }


    }
}

//@}

} // namespace ElectrophysiologyUtility

} // namespace LifeV

#endif /* HEARTUTILITY_H */
