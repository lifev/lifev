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
    @brief Class for interfacing dofs between two 3D meshes

    @author M.A. Fernandez
    @date 00-11-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @contributor Mauro Perego <mperego@fsu.edu>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    This file contains the class which may be used to update and hold the connections between the dof
    on two matching meshes.
 */

#ifndef _DOFINTERFACE3DTO3D_HH
#define _DOFINTERFACE3DTO3D_HH

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/MarkerDefinitions.hpp>

#include <lifev/core/fem/DOFInterface.hpp>
#include <lifev/core/fem/ReferenceFE.hpp>
#include <lifev/core/fem/DOF.hpp>
#include <lifev/core/fem/CurrentBoundaryFE.hpp>
#include <lifev/core/util/LifeChrono.hpp>

namespace LifeV
{
/*!
  \class DOFInterface3Dto3D

  Base class which builds and holds the connections of the DOF associated with matching facets with given flags.
  The connections can be between two different matching meshes, or within the same mesh (having matching boundary facets).

  In order to hold the interface connections the user must give the
  ReferenceFE elements and DOF used in both meshes.  The connections may be
  built by calling the update method. An interpolate method has been
  provided in order to interpolate data at the interface. Finally
  method getInterfaceDof gives the connections.

*/
class DOFInterface3Dto3D:
    public DOFInterface
{
public:
    //! @name Public typedefs
    //@{
    typedef boost::numeric::ublas::vector<Real> Vector;
    typedef boost::function<bool ( const std::vector<Real>&, const std::vector<Real>&, const Real& ) > fct;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! default constructor
    DOFInterface3Dto3D()
        :
        M_refFE1 ( 0 ),
        M_dof1 ( 0 ),
        M_refFE2 ( 0 ),
        M_dof2 ( 0 )
    {}

    //! Constructor for interfacing DOF of the same type (ReferenceFE)
    /*!
     \param refFe the reference FE used in both meshes
     \param dof1 the DOF object of the mesh in which we want to make the computations
     \param dof2 the DOF object of the mesh which provides the data at the interface
    */
    DOFInterface3Dto3D ( const ReferenceFE& refFE, const DOF& dof1, const DOF& dof2 );

    //! Constructor for interfacing DOF of different type (ReferenceFE)
    /*!
      \param refFe1 the reference FE used in the mesh in which we want to make the computations
      \param dof1 the DOF object of the mesh in which we want to make the computations
      \param refFe2 the reference FE used in the mesh which provides the data at the interface
      \param dof2 the DOF object of the mesh which provides the data at the interface
     */
    DOFInterface3Dto3D ( const ReferenceFE& refFE1, const DOF& dof1, const ReferenceFE& refFE2, const DOF& dof2 );

    //! Constructor for interfacing DOF of the same type (ReferenceFE)
    /*!
     \param dof the DOF object of the mesh in which we want to make the computations
    */
    DOFInterface3Dto3D ( const ReferenceFE& refFE, const DOF& dof);


    //@}


    //! @name Methods
    //@{

    //! Setup method for common finite element
    void setup (const ReferenceFE& refFE, const DOF& dof1, const DOF& dof2 );

    //! Setup method for possibly different finite element
    void setup (const ReferenceFE& refFE1, const DOF& dof1, const ReferenceFE& refFE2, const DOF& dof2 );

    //! This method builds the DOF connections at the interface
    /*!
      \param mesh1 the mesh in which we want to make the computations
      \param flag1 the marker of the interface in the mesh1
      \param mesh2 the mesh which provides de data at the interface
      \param flag2 the marker of the interface in the mesh2
      \param tol tolerance for connecting points of both meshes at the interface
      \param flag3 the marker of a region of interface in the mesh1
      \brief{The parameter flag3 is used in test_meshReorder to export the part of interface determined by flag3 on mesh2.}
     */

    template <typename MeshType>
    void update ( MeshType& mesh1,
                  const markerID_Type& flag1,
                  MeshType& mesh2,
                  const markerID_Type& flag2,
                  const Real& tol,
                  Int const* const flag3 = 0 );

    //! This method builds the DOF connections between two matching surfaces belonging to the same mesh
    /*!
      \param mesh the mesh in which we want to make the computations
      \param flag1 the marker of the first set of facets in the mesh
      \param flag2 the marker of facets to be connected with facets marked with flag1
      \param tol tolerance for connecting points of both meshes at the interface
      \param coupled a function of points p1 and p2, that returns true if the points are coupled, false otherwise.
             points p1 and p2 are constituted of std::vector<Real>(3) containing the coordinates.
    */
    template <typename MeshType>
    void update ( MeshType& mesh,
                  const markerID_Type& flag1,
                  const markerID_Type& flag2,
                  const Real& tol, const fct& coupled );


    //! This method interpolate data when using different FE
    /*!
      \param mesh2 the mesh which provides de data at the interface
      \param v the data vector on mesh2 holding dofs of type refFE1
      \param vI the interpolated data vector on mesh2 holding dofs of type refFE2

      \note We should use this method ONLY when the accuracy of the data is less that
            the accuracy of the unknowns on mesh1, i.e., when refFE1 is more accurate
     than refFE2
     */
    template <typename MeshType>
    void interpolate ( MeshType& mesh2, const UInt nbComp, const Vector& v, Vector& vI );

    //@}


    //! @name Get Methods
    //@{

    //! This method returns the corresponding dof object when interpolation is used
    const UInt& numTotalDof() const
    {
        return M_dof->numTotalDof();
    }

    //! reference to the list of connected facets.
    const std::list< std::pair<ID, ID> >& connectedFacetMap() const
    {
        return M_facetToFacetConnectionList;
    }

    //@}

private:

    //!  STL iterator type for the lists
    typedef std::list< std::pair<ID, ID> >::iterator Iterator;


    //! @name Private Methods
    //@{

    //! This method builds the connections between facets at the interface (M_facetToFacetConnectionList container)
    /*!
      \param mesh1 the mesh in which we want to make the computations
      \param flag1 the marker of the interface in the mesh1
      \param mesh2 the mesh which provides de data at the interface
      \param flag2 the marker of the interface in the mesh2
      \param tol tolerance for connecting points of both meshes at the interface
      \param coupled a function of points p1 and p2, that returns true if the points are coupled, false otherwise.
             points p1 and p2 are constituted of std::vector<Real>(3) containing the coordinates.
     */
    template <typename MeshType>
    void updateFacetConnections ( const MeshType& mesh1, const markerID_Type& flag1,
                                  const MeshType& mesh2, const markerID_Type& flag2, const Real& tol, const fct& coupled );
    //! This method builds the connections between DOF at the interface (M_dofToDofConnectionList container)
    /*!
      \param mesh1 the mesh in which we want to make the computations
      \param dof1 the DOF object of the mesh in which we want to make the computations
      \param mesh2 the mesh which provides de data at the interface
      \param dof2 the DOF object of the mesh which provides de data at the interface
      \param tol tolerance for connecting points of both meshes at the interface
      \param coupled a function of points p1 and p2, that returns true if the points are coupled, false otherwise.
             points p1 and p2 are constituted of std::vector<Real>(3) containing the coordinates.
      \param flag3 the marker of a region of interface in the mesh1
      \brief{The parameter flag3 is used in test_meshReorder to export the part of interface determined by flag3 on mesh2.}
      */
    template <typename MeshType>
    void updateDofConnections ( const MeshType& mesh1, const DOF& dof1,
                                const MeshType& mesh2, const DOF& dof2, const Real& tol, const fct& coupled,
                                Int const* const flag3 = 0 );

    //@}


    //! ReferenceFE object used in the mesh in which we want to make the computations
    const ReferenceFE* M_refFE1;

    //! DOF object of the mesh in which we want to make the computations
    const DOF* M_dof1;

    //! ReferenceFE object used in the mesh which provides de data at the interface
    const ReferenceFE* M_refFE2;

    //! DOF object of the mesh which provides the data at the interface
    const DOF* M_dof2;

    //! Auxiliary DOF object of the mesh which provides the local to global table
    //! when interpolation is used
    boost::shared_ptr<DOF> M_dof;

    //! STL list which holds the connections between facets at the interface
    std::list< std::pair<ID, ID> > M_facetToFacetConnectionList;

    //!  Auxiliary STL list which holds the connections between DOF at the interface
    //! Empty after calling update
    std::list< std::pair<ID, ID> > M_dofToDofConnectionList;

};

// ===================================================
// Helpers
// ===================================================

//! Returns true if the Points p1 and p2 are equal with respect to the tolerance tol (in norm 1)
bool coincide ( const std::vector<Real>& p1, const std::vector<Real>& p2, const Real& tol );

// ===================================================
// Methods
// ===================================================

template <typename MeshType>
void DOFInterface3Dto3D::update ( MeshType& mesh1, const markerID_Type& flag1,
                                  MeshType& mesh2, const markerID_Type& flag2,
                                  const Real& tol, Int const* const flag3 )
{

    // Updating facet connections at the interface
    updateFacetConnections ( mesh1, flag1, mesh2, flag2, tol, coincide );

    if ( M_refFE1->nbDof() > M_refFE2->nbDof() )
    {
        // Update of the DOF connections when we need interpolation
        M_dof->update ( mesh2 ); // Building auxiliary dof
        updateDofConnections ( mesh1, *M_dof1, mesh2, *M_dof, tol, coincide, flag3 ); // Update of the DOF connections
    }
    else
        // Update of the DOF connections without interpolation
    {
        updateDofConnections ( mesh1, *M_dof1, mesh2, *M_dof2, tol, coincide, flag3 );
    }
}


template <typename MeshType>
void DOFInterface3Dto3D::update ( MeshType& mesh, const markerID_Type& flag1,
                                  const markerID_Type& flag2, const Real& tol, const fct& coupled)
{
    // Updating facet connections at the interface
    updateFacetConnections ( mesh, flag1, mesh, flag2, tol, coupled );

    // Update of the DOF connections without interpolation
    updateDofConnections ( mesh, *M_dof1, mesh, *M_dof2, tol, coupled );
}


template <typename MeshType>
void DOFInterface3Dto3D::interpolate ( MeshType& mesh2, const UInt nbComp, const Vector& v, Vector& vI )
{

    typedef typename MeshType::elementShape_Type GeoShape; // Element shape
    typedef typename GeoShape::GeoBShape GeoBShape;  // Facet Shape

    UInt nbVertexPerFacet = GeoBShape::S_numVertices; // Number of facet's vertices
    UInt nbRidgePerFacet = GeoBShape::S_numRidges;    // Number of facet's ridges

    UInt nbDofPerVertex = M_refFE1->nbDofPerVertex(); // number of DOF per vertices
    UInt nbDofPerRidge = M_refFE1->nbDofPerRidge();   // number of DOF per ridges
    UInt nbDofPerFacet = M_refFE1->nbDofPerFacet();   // number of DOF per facets

    UInt nbVertexPerElement = GeoShape::S_numVertices; // Number of element's vertices
    UInt nbRidgePerElement = GeoShape::S_numRidges;    // Number of element's ridges

    UInt nbDofPerElement = M_refFE2->nbDof(); // Number of DOF per element in the lowDof mesh

    UInt nbVertexDofPerElement = nbVertexPerElement * nbDofPerVertex; // number of vertex's DOF on a Element
    UInt nbRidgeDofPerElement = nbRidgePerElement * nbDofPerRidge; // number of ridge's DOF on a Element

    ID ibF, iElAd, iFaEl, iVeEl, lDof, iEdEl;

    Real x, y, z, sum;

    std::vector<Real> vLoc ( nbDofPerElement * nbComp );


    // Loop on facets at the interface (matching facets)
    for ( Iterator i = M_facetToFacetConnectionList.begin(); i != M_facetToFacetConnectionList.end(); ++i )
    {

        ibF = i->second; // Facet number at the interface

        iElAd = mesh2.boundaryFacet ( ibF ).firstAdjacentElementIdentity(); // id of the element adjacent to the facet
        iFaEl = mesh2.boundaryFacet ( ibF ).firstAdjacentElementPosition(); // local id of the facet in its adjacent element

        // Updating the local dof of the data vector in the adjacent element
        for ( UInt icmp = 0; icmp < nbComp; ++icmp )
            for ( ID idof = 0; idof < nbDofPerElement; ++idof )
            {
                vLoc[ icmp * nbDofPerElement + idof ] = v [ icmp * M_dof2->numTotalDof() + M_dof2->localToGlobalMap ( iElAd, idof ) ];
            }

        // Vertex based Dof
        if ( nbDofPerVertex )
        {

            // loop on facet vertices
            for ( ID iVeFa = 0; iVeFa < nbVertexPerFacet; ++iVeFa )
            {

                iVeEl = GeoShape::facetToPoint ( iFaEl, iVeFa ); // local vertex number (in element)

                // Loop number of DOF per vertex
                for ( ID l = 0; l < nbDofPerVertex; ++l )
                {
                    lDof = iVeEl * nbDofPerVertex + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    x = M_refFE1->xi ( lDof );
                    y = M_refFE1->eta ( lDof );
                    z = M_refFE1->zeta ( lDof );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        sum = 0;
                        for ( ID idof = 0; idof < nbDofPerElement; ++idof )  // Loop on local DOF on the adjacent element
                        {
                            sum += vLoc[ icmp * nbDofPerElement + idof ] * M_refFE2->phi ( idof, x, y, z );
                        }

                        // Updating interpolating vector
                        vI ( icmp * M_dof->numTotalDof() + M_dof->localToGlobalMap ( iElAd, lDof ) ) = sum;
                    }
                }
            }
        }

        // Ridge based Dof
        if ( nbDofPerRidge )
        {

            // loop on facet ridges
            for ( ID iEdFa = 0; iEdFa < nbRidgePerFacet; ++iEdFa )
            {

                iEdEl = GeoShape::faceToRidge ( iFaEl, iEdFa ).first; // local ridge number (in element)

                // Loop number of DOF per ridge
                for ( ID l = 0; l < nbDofPerRidge; ++l )
                {
                    lDof = nbVertexDofPerElement + iEdEl * nbDofPerRidge + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    x = M_refFE1->xi ( lDof );
                    y = M_refFE1->eta ( lDof );
                    z = M_refFE1->zeta ( lDof );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        sum = 0;
                        for ( ID idof = 0; idof < nbDofPerElement; ++idof )   // Loop on local DOF on the adjacent element
                        {
                            sum += vLoc[ icmp * nbDofPerElement + idof ] * M_refFE2->phi ( idof, x, y, z );
                        }

                        // Updating interpolating vector
                        vI ( icmp * M_dof->numTotalDof() + M_dof->localToGlobalMap ( iElAd, lDof ) ) = sum;
                    }
                }
            }
        }

        // Loop on number of DOF per facet
        for ( ID l = 0; l < nbDofPerFacet; ++l )
        {
            lDof = nbRidgeDofPerElement + nbVertexDofPerElement + iFaEl * nbDofPerFacet + l; // Local dof in the adjacent Element

            // Nodal coordinates
            x = M_refFE1->xi ( lDof );
            y = M_refFE1->eta ( lDof );
            z = M_refFE1->zeta ( lDof );

            // Loop on data vector components
            for ( UInt icmp = 0; icmp < nbComp; ++icmp )
            {

                // Interpolating data at the nodal point
                sum = 0;
                for ( ID idof = 0; idof < nbDofPerElement; ++idof )  // Loop on local DOF on the adjacent element
                {
                    sum += vLoc[ icmp * nbDofPerElement + idof ] * M_refFE2->phi ( idof, x, y, z );
                }

                // Updating interpolating vector
                vI ( icmp * M_dof->numTotalDof() + M_dof->localToGlobalMap ( iElAd, lDof ) ) = sum;
            }
        }
    }
}

// ===================================================
// Private Methods
// ===================================================

template <typename MeshType>
void DOFInterface3Dto3D::updateFacetConnections ( const MeshType& mesh1, const markerID_Type& flag1,
                                                  const MeshType& mesh2, const markerID_Type& flag2, const Real& tol, const fct& coupled )
{

    UInt bdnF1 = mesh1.numBoundaryFacets(); // Number of boundary facets in mesh1
    UInt bdnF2 = mesh2.numBoundaryFacets(); // Number of boundary facets mesh2

    markerID_Type marker1;
    std::set<ID> facetsFlagged2;


    typedef typename MeshType::facetShape_Type GeoBShape; // Shape of the facets

    UInt nbVertexPerFacet = GeoBShape::S_numVertices; // Number of facet's vertices

    std::vector<Real>  v1 (nDimensions), v2 ( nDimensions );
    std::vector< std::vector<Real> > vertexVector (nbVertexPerFacet);

    LifeChrono  chrono;
    chrono.start();

    // select facets flagged with flag 2
    for ( ID ibF2 = 0; ibF2 < bdnF2; ++ibF2 )
        if ( flag2 == mesh2.boundaryFacet ( ibF2 ).markerID() )
        {
            facetsFlagged2.insert (ibF2);
        }

    // Loop on boundary facets on mesh1
    for ( ID ibF1 = 0; ibF1 < bdnF1; ++ibF1 )
    {

        // The facet marker
        marker1 = mesh1.boundaryFacet ( ibF1 ).markerID();

        // Is the facet on the interface?
        if ( marker1 == flag1 )
        {

            // Loop on facet vertices
            for ( ID iVeFa = 0; iVeFa < nbVertexPerFacet; ++iVeFa )
            {

                // Loop on vertex coordinates
                for ( ID j = 0; j < nDimensions; ++j )
                {
                    v1[ j ] = mesh1.boundaryFacet ( ibF1 ).point ( iVeFa ).coordinate ( j );
                }
                vertexVector[iVeFa] = v1;
            }
            // Loop on boundary facets on mesh2
            std::set<ID>::iterator it = facetsFlagged2.begin();
            for (; it != facetsFlagged2.end(); ++it )
            {
                ID ibF2 = *it;

                // Loop on facet vertices
                bool matched;
                ID iVeFa = 0;
                do
                {
                    // Loop on vertex coordinates
                    for ( ID j = 0; j < nDimensions; ++j )
                    {
                        v2[ j ] = mesh2.boundaryFacet ( ibF2 ).point ( iVeFa ).coordinate ( j );
                    }

                    // Loop on facet vertices on mesh1
                    ID ivefa = 0;
                    do
                    {
                        // Do the vertices match? if yes break the loop
                        matched = coupled ( vertexVector[ivefa], v2, tol );
                    }
                    while ( (++ivefa < nbVertexPerFacet) && !matched );
                }
                while ( (++iVeFa < nbVertexPerFacet) && matched);
                //! Do the facets match?
                if ( matched )
                {
                    std::pair<ID, ID> elc ( ibF1, ibF2 );
                    M_facetToFacetConnectionList.push_front ( elc );
                    facetsFlagged2.erase (it);
                    break; // Stop loop on boundary facets on mesh2
                }

            }
        }
    }
    chrono.stop();
}

//! This method builds the connections between DOF at the interface (M_dofToDofConnectionList container)
/*!
  \param mesh1 the mesh in which we want to make the computations
  \param dof1 the DOF object of the mesh in which we want to make the computations
  \param mesh2 the mesh which provides the data at the interface
  \param dof2 the DOF object of the mesh which provides the data at the interface
  \param tol tolerance for connecting points of both meshes at the interface
*/
template <typename Mesh>
void DOFInterface3Dto3D::updateDofConnections ( const Mesh& mesh1, const DOF& dof1,
                                                const Mesh& mesh2, const DOF& dof2, const Real& tol, const fct& coupled, Int const* const flag3)
{
    CurrentBoundaryFE feBd1 ( M_refFE1->boundaryFE(), getGeometricMap ( mesh1 ).boundaryMap() );
    CurrentBoundaryFE feBd2 ( M_refFE2->boundaryFE(), getGeometricMap ( mesh2 ).boundaryMap() );

    std::vector<Real> p1 ( nDimensions ), p2 ( nDimensions );

    // Loop on facets at the interface (matching facets)
    for ( Iterator i = M_facetToFacetConnectionList.begin(); i != M_facetToFacetConnectionList.end(); ++i )
    {
        feBd1.update ( mesh1.boundaryFacet ( i->first ), UPDATE_ONLY_CELL_NODES ); // Updating facet information on mesh1
        feBd2.update ( mesh2.boundaryFacet ( i->second ), UPDATE_ONLY_CELL_NODES ); // Updating facet information on mesh2

        std::vector<ID> localToGlobalMapOnBFacet1 = dof1.localToGlobalMapOnBdFacet (i->first);
        std::vector<ID> localToGlobalMapOnBFacet2 = dof2.localToGlobalMapOnBdFacet (i->second);

        for (ID lDof1 = 0; lDof1 < localToGlobalMapOnBFacet1.size(); lDof1++)
        {
            if ( flag3 != 0 && mesh1.boundaryFacet (i->first).point (lDof1).markerID() == *flag3)
            {
                continue;
            }
            ID gDof1 = localToGlobalMapOnBFacet1[lDof1];
            feBd1.coorMap ( p1[0], p1[1], p1[2], feBd1.refFE().xi ( lDof1 ), feBd1.refFE().eta ( lDof1 ) ); // Nodal coordinates on the current facet (mesh1)

            for (ID lDof2 = 0; lDof2 < localToGlobalMapOnBFacet2.size(); lDof2++)
            {
                ID gDof2 = localToGlobalMapOnBFacet2[lDof2];
                feBd2.coorMap ( p2[0], p2[1], p2[2], feBd2.refFE().xi ( lDof2 ), feBd2.refFE().eta ( lDof2 ) );

                if ( coupled ( p1, p2, tol ) )
                {
                    std::pair<ID, ID> locDof ( gDof1, gDof2 );
                    M_dofToDofConnectionList.push_front ( locDof ); // Updating the list of dof connections
                    break;
                }
            }
        }
    }

    // Updating the map containter with the connections
    for ( Iterator i = M_dofToDofConnectionList.begin(); i != M_dofToDofConnectionList.end(); ++i )
    {
        M_localDofMap[ i->first ] = i->second;
    }

    // Saving memory
    M_dofToDofConnectionList.clear();
}


}
#endif
