/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*!
  \file FESpace.hpp

  \version 1.0
  \author Gilles Fourestey
  \date 06/2007

  \version 1.2 - Added new constructors for "Standard FE Spaces" and some other minor changes.
  \author Cristiano Malossi<cristiano.malossi@epfl.ch>
  \date 06/2009

  \version 1.3 - Added interpolation features for a given solution vector
  \author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
  \date 08/2009

  \brief This file contains an abstract class for NavierStokes solvers.

*/

#ifndef _FESPACE_H_
#define _FESPACE_H_

#include <cmath>
#include <ext/slist>
#include <iomanip>
#include <sstream>
#include <utility>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <life/lifearray/pattern.hpp>
#include <life/lifearray/SimpleVect.hpp>

#include <life/lifecore/life.hpp>

#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/currentFE.hpp>
#include <life/lifefem/currentBdFE.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/sobolevNorms.hpp>

#include <life/lifemesh/partitionMesh.hpp>

using std::pair;

namespace LifeV
{

/*!
  \class FESpace

  Abstract class which defines the general structure of a NavierStokes solver.
  For each new NavierStokes solver  we have to implement the corresponding
  timeAdvance and an iterate methods

*/


inline Real
elemL22( boost::function<Real( Real, Real, Real, Real, UInt )> fct,
         const CurrentFE& fe, const Real t, const UInt nbcomp )
{
    Int ig;
    UInt ic;
    Real s = 0., f, x, y, z;
    for ( ig = 0;ig < fe.nbQuadPt;ig++ )
    {
        fe.coorQuadPt( x, y, z, ig );
        for ( ic = 0; ic < nbcomp; ic++ )
        {
            f = fct( t, x, y, z, ic + 1 );
            s += f * f * fe.weightDet( ig );
        }
    }
    return s;
}



template <typename Mesh, typename Map>
class FESpace
{

public:

    typedef boost::function<Real ( Real const&, Real const&, Real const&,
                                   Real const&, ID const& )> Function;


    typedef Mesh                         mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef Map                          map_type;

    /** @name Constructors, Destructor
     */
	//@{

    /*!
      \param data_file GetPot data file
      \param refFE_u reference FE for the velocity
      \param refFE_p reference FE for the pressure
      \param Qr_u volumic quadrature rule for the velocity
      \param bdQr_u surface quadrature rule for the velocity
      \param Qr_p volumic quadrature rule for the pressure
      \param bdQr_p surface quadrature rule for the pressure
      \param BCh_fluid boundary conditions for the fluid
      \param ord_bdf order of the bdf time advancing scheme and incremental pressure approach (default: Backward Euler)
    */

    FESpace(	partitionMesh<Mesh>&	mesh,
				const RefFE&			refFE,
				const QuadRule&			Qr,
				const QuadRule&			bdQr,
				const Int				fDim,
				Epetra_Comm&			comm
			);

    FESpace(	partitionMesh<Mesh>&	mesh,
				const std::string&		space,
				const Int				fDim,
				Epetra_Comm&			comm
			);

    FESpace(	mesh_ptrtype			mesh,
				const RefFE&			refFE,
				const QuadRule&			Qr,
				const QuadRule&			bdQr,
				const Int				fDim,
				Epetra_Comm&			comm
			);

    FESpace(	mesh_ptrtype			mesh,
				const std::string&		space,
				const Int				fDim,
				Epetra_Comm&			comm
			);

    //! Do nothing destructor
    virtual ~FESpace() {}

    //@}



    /** @name Get Functions
     */
	//@{

    //    const mesh_type& mesh()       {return M_mesh;}
    //    const mesh_type& mesh() const {return *M_mesh;}
    //    mesh_type& mesh() {return *M_mesh;}

    const mesh_ptrtype 	mesh()	const { return M_mesh; }
    mesh_ptrtype 		mesh()		  { return M_mesh; }

    //! Returns map
    const map_type&		map()	const { return M_map; }
    map_type&			map()		  { return M_map; }

	//! Returns the velocity dof
	const Dof&			dof()	const { return *M_dof; }
	Dof&				dof() 		  { return *M_dof; }

	//! Returns the current FE
	const CurrentFE&	fe()	const { return *M_fe; }
	CurrentFE& 			fe()		  { return *M_fe; }

	//! Returns the current boundary FE
	CurrentBdFE&		feBd()		  { return *M_feBd; }

	//! Returns the res FE
	const RefFE&		refFE()	const { return *M_refFE; }

	//! Returns the volumic quadratic rule
	const QuadRule&		qr() 	const { return *M_Qr; }

	//! Returns the surfasic quadratic rule
	const QuadRule&		bdQr()	const { return *M_bdQr; }

	//! Returns FE space dimension
    /*const*/ UInt		dim()      const { return M_dim; }
    /*const*/ UInt		fieldDim() const { return M_fieldDim; }

    //@}


	void initMap( const Int fDim, Epetra_Comm& comm );

    //! Interpolate a given velocity function nodally onto a velocity vector
    template<typename vector_type>
    void interpolate( const Function& fct, vector_type& vect, const Real time = 0. );

    template<typename vector_type>
    void interpolateBC( BCHandler& BCh, vector_type&    vect, const Real      time);

        //! takes into account a possible offset by a constant
    //! \param pexact the exact pressure as a function
    //! \param time the time
    //! \param relError Real* to store the relative error in
//     Real pErrorL2( const Function&            pexact,
//                    const ScalUnknown<Vector>& p,
//                    Real                       time,
//                    Real*                      relError = 0 );

    //! calculate L2 velocity error for given exact velocity function
    //! \param pexact the exact velocity as a function
    //! \param time the time
    //! \param relError Real* to store the relative error in

    // this computes vec = \int fct phi_i
    template<typename vector_type>
    void L2ScalarProduct(	const Function&		fct,
							vector_type&		vec,
							const Real			t
						);

    template<typename vector_type>
    Real L20Error(			const Function& 	fexact,
							const vector_type& 	vec,
							const Real 			time,
							Real* 				relError = 0 );

    template<typename function, typename vector_type>
    Real L2Error(			const function&		fexact,
							const vector_type&	vec,
							const Real			time,
							Real*				relError=0 );

    template<typename function, typename vector_type>
    Real H1Error(			const function&		fexact,
							const vector_type&	vec,
							const Real			time,
							Real*				relError=0 );


    template<typename vector_type>
    Real L2Norm( const vector_type& vec );

    template<typename vector_type>
    Real H1Norm( const vector_type& vec );
  
  
    // Added by S. Quinodoz 21.08.2009 
    // Interpolate the function given by a FE Vector "solutionVector" (not a real function)
    // and a cell "elementID" in ANY given point "pt" in the space
    //
    // !!! Until now, only for scalar FE (not checked)

    template<typename point_type, typename vector_type>
    Real FEinterpolateValue(const ID& elementID, const vector_type& solutionVector, const point_type& pt ) const;

    template<typename point_type, typename vector_type>
    Real FEinterpolateValueLocal(const ID& elementID, const vector_type& solutionVector, const point_type& pt ) const;

    template<typename point_type, typename vector_type>
    Real FEinterpolateGradient(const ID& elementID, const vector_type& solutionVector, const point_type& pt,
			       const unsigned int& component ) const;

    template<typename point_type, typename vector_type>
    Real FEinterpolateGradientLocal(const ID& elementID, const vector_type& solutionVector, const point_type& pt,
				    const unsigned int& component ) const;


    BasePattern::PatternType patternType();

private:

    //! copy constructor
    FESpace( const FESpace& fespace );

    //! Set space
    inline void setSpace( const std::string& space );

    //! Set space Map (useful for switch syntax with strings)
    enum spaceData{P1, P1_HIGH, P1Bubble, P2, P2_HIGH};
    std::map<std::string, spaceData>	M_spaceMap;

    //! reference to the mesh
    mesh_ptrtype						M_mesh;

    //! Reference FE for the velocity
    const RefFE*    					M_refFE;

    //! Quadrature rule for volumic elementary computations
    const QuadRule* 					M_Qr;

    //! Quadrature rule for surface elementary computations
    const QuadRule* 					M_bdQr;

    //! dimention of the field variable ( scalar/vector field)
    UInt								M_fieldDim;

    //! A shared pointer to the Dof object
    boost::shared_ptr<Dof>				M_dof;

    //! The number of total dofs
    UInt								M_dim;

    //! Current FE
    boost::shared_ptr<CurrentFE>		M_fe;
    boost::shared_ptr<CurrentBdFE>		M_feBd;

    //! Map
    map_type							M_map;

};

//
// IMPLEMENTATION
//


// Constructors
template <typename Mesh, typename Map>
FESpace<Mesh, Map>::
FESpace(	partitionMesh<Mesh>& 	mesh,
			const RefFE&         	refFE,
			const QuadRule&      	Qr,
			const QuadRule&      	bdQr,
			const Int            	fDim,
			Epetra_Comm&         	comm
		) :
        M_mesh			( mesh.mesh() ),
        M_refFE			( &refFE ),
        M_Qr			( &Qr ),
        M_bdQr			( &bdQr ),
        M_fieldDim		( fDim ),
        M_dof			( new Dof( *M_mesh, *M_refFE ) ),
        M_dim			( M_dof->numTotalDof() ),
        M_fe			( new CurrentFE  ( *M_refFE,              getGeoMap( *M_mesh ),               *M_Qr ) ),
        M_feBd			( new CurrentBdFE( M_refFE->boundaryFE(), getGeoMap( *M_mesh ).boundaryMap(), *M_bdQr ) ),
        M_map			( )
{
	Map map( *M_refFE, *M_mesh, comm );
	for ( UInt ii = 0; ii < M_fieldDim; ++ii )
		M_map += map;
}

template <typename Mesh, typename Map>
FESpace<Mesh, Map>::
FESpace(	partitionMesh<Mesh>&	mesh,
			const std::string&		space,
			const Int				fDim,
			Epetra_Comm&			comm
		) :
	M_mesh			( mesh.mesh() ),
    M_fieldDim		( fDim ),
    M_dof			( ),
    M_fe			( ),
    M_feBd			( ),
    M_map			( )
{
	// Set spaceMap
	M_spaceMap["P1"]		= P1;
	M_spaceMap["P1_HIGH"]	= P1_HIGH;
	M_spaceMap["P1Bubble"]	= P1Bubble;
	M_spaceMap["P2"]		= P2;
	M_spaceMap["P2_HIGH"]	= P2_HIGH;

	// Set space
	setSpace( space );

	// Set other quantities
	M_dof.reset( new Dof( *M_mesh, *M_refFE ) );
	M_dim = M_dof->numTotalDof();
	M_fe.reset( new CurrentFE( *M_refFE, getGeoMap( *M_mesh ), *M_Qr ) );
	M_feBd.reset( new CurrentBdFE( M_refFE->boundaryFE(), getGeoMap( *M_mesh ).boundaryMap(), *M_bdQr ) );

	// Build Map
	Map map( *M_refFE, *M_mesh, comm );
	for ( UInt ii = 0; ii < M_fieldDim; ++ii )
		M_map += map;
}

template <typename Mesh, typename Map>
FESpace<Mesh, Map>::
FESpace(	mesh_ptrtype			mesh,
			const RefFE&			refFE,
			const QuadRule&			Qr,
			const QuadRule&			bdQr,
			const Int				fDim,
			Epetra_Comm&			comm
		) :
    M_mesh			( mesh ),
    M_refFE			( &refFE ),
    M_Qr			( &Qr ),
    M_bdQr			( &bdQr ),
    M_fieldDim		( fDim ),
    M_dof			( new Dof( *M_mesh, *M_refFE ) ),
    M_dim			( M_dof->numTotalDof() ),
    M_fe			( new CurrentFE( *M_refFE, getGeoMap( *M_mesh ), *M_Qr ) ),
    M_feBd			( new CurrentBdFE( M_refFE->boundaryFE(), getGeoMap( *M_mesh ).boundaryMap(), *M_bdQr ) ),
    M_map			( )
{
	Map map( *M_refFE, *M_mesh, comm );
	for ( UInt ii = 0; ii < M_fieldDim; ++ii )
		M_map += map;
}

template <typename Mesh, typename Map>
FESpace<Mesh, Map>::
FESpace(	mesh_ptrtype			mesh,
			const std::string&		space,
			const Int				fDim,
			Epetra_Comm&			comm
		) :
	M_mesh			( mesh ),
    M_fieldDim		( fDim ),
    M_dof			( ),
    M_fe			( ),
    M_feBd			( ),
    M_map			( )
{
	// Set spaceMap
	M_spaceMap["P1"]		= P1;
	M_spaceMap["P1_HIGH"]	= P1_HIGH;
	M_spaceMap["P1Bubble"]	= P1Bubble;
	M_spaceMap["P2"]		= P2;
	M_spaceMap["P2_HIGH"]	= P2_HIGH;

	// Set space
	setSpace( space );

	// Set other quantities
	M_dof.reset( new Dof( *M_mesh, *M_refFE ) );
	M_dim = M_dof->numTotalDof();
	M_fe.reset( new CurrentFE( *M_refFE, getGeoMap( *M_mesh ), *M_Qr ) );
	M_feBd.reset( new CurrentBdFE( M_refFE->boundaryFE(), getGeoMap( *M_mesh ).boundaryMap(), *M_bdQr ) );

	// Build Map
	Map map( *M_refFE, *M_mesh, comm );
	for ( UInt ii = 0; ii < M_fieldDim; ++ii )
        M_map += map;
}

// Set FE space (default standard parameters)
template <typename Mesh, typename Map>
inline void
FESpace<Mesh, Map>::setSpace( const std::string& space )
{
	switch ( M_spaceMap[space] )
	{
		case P1 :

			#ifdef TWODIM
				M_refFE = &feTriaP1;
				M_Qr    = &quadRuleTria6pt;
				M_bdQr  = &quadRuleSeg3pt;
			#elif defined THREEDIM
				M_refFE = &feTetraP1;
				M_Qr    = &quadRuleTetra4pt;
				M_bdQr	= &quadRuleTria3pt;
			#endif

			break;

		case P1_HIGH :

			#ifdef TWODIM
				M_refFE = &feTriaP1;
				M_Qr    = &quadRuleTria6pt;
				M_bdQr  = &quadRuleSeg3pt;
			#elif defined THREEDIM
				M_refFE = &feTetraP1;
				M_Qr    = &quadRuleTetra15pt;
				M_bdQr	= &quadRuleTria4pt;
			#endif

			break;

		case P1Bubble :

			#ifdef TWODIM

				break;

			#elif defined THREEDIM
				M_refFE = &feTetraP1bubble;
				M_Qr    = &quadRuleTetra64pt;
				M_bdQr	= &quadRuleTria3pt;
			#endif

			break;

		case P2 :

			#ifdef TWODIM
				M_refFE = &feTriaP2;
				M_Qr    = &quadRuleTria6pt;
				M_bdQr  = &quadRuleSeg3pt;
			#elif defined THREEDIM
				M_refFE = &feTetraP2;
				M_Qr    = &quadRuleTetra15pt;
				M_bdQr	= &quadRuleTria3pt;
			#endif

		case P2_HIGH :

			#ifdef TWODIM
				M_refFE = &feTriaP2;
				M_Qr    = &quadRuleTria6pt;
				M_bdQr  = &quadRuleSeg3pt;
			#elif defined THREEDIM
				M_refFE = &feTetraP2;
				M_Qr    = &quadRuleTetra64pt;
				M_bdQr	= &quadRuleTria4pt;
			#endif

			break;

		default :

			std::cout << "!!! WARNING: Space " << space << "not implemented in FESpace class !!!" << std::endl;
	}
}


template <typename Mesh, typename Map>
template<typename vector_type>
void
FESpace<Mesh, Map>::interpolate( const Function& fct,
                                 vector_type&    vect,
                                 const Real      time)
{
    typedef typename mesh_type::ElementShape GeoShape ; // Element shape

    UInt nDofpV    = refFE().nbDofPerVertex; // number of Dof per vertex
    UInt nDofpE    = refFE().nbDofPerEdge;   // number of Dof per edge
    UInt nDofpF    = refFE().nbDofPerFace;   // number of Dof per face
//    UInt nDofpEl   = refFE().nbDofPerVolume; // number of Dof per Volume

    UInt nElemV    = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE    = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF    = GeoShape::numFaces;    // Number of element's faces

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element
//    UInt nDofElemF = nElemF * nDofpF; // number of face's Dof on a Element

    ID nbComp = M_fieldDim; // Number of components of the mesh velocity

    Real x, y, z;

    ID lDof;

    // Loop on elements of the mesh
    for ( ID iElem = 1; iElem <= mesh()->numElements(); ++iElem )
    {

        fe().updateJac( mesh()->element( iElem ) );
        ID elemId = mesh()->element( iElem ).localId();

        // Vertex based Dof
        if ( nDofpV )
        {
            // loop on element vertices
            for ( ID iVe = 1; iVe <= nElemV; ++iVe )
            {
                // Loop number of Dof per vertex
                for ( ID l = 1; l <= nDofpV; ++l )
                {
                    lDof = ( iVe - 1 ) * nDofpV + l; // Local dof in this element

                    // Nodal coordinates
                    fe().coorMap( x, y, z, refFE().xi( lDof - 1 ), refFE().eta( lDof - 1 ), refFE().zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        UInt iDof = icmp * dim() + dof().localToGlobal( elemId, lDof  );
                        vect.checkAndSet( iDof, fct( time, x, y, z, icmp + 1 ));
                    }
                }
            }
        }
        // Edge based Dof
        if ( nDofpE )
        {

            // loop on element edges
            for ( ID iEd = 1; iEd <= nElemE; ++iEd )
            {

                // Loop number of Dof per edge
                for ( ID l = 1; l <= nDofpE; ++l )
                {
                    lDof = nDofElemV + ( iEd - 1 ) * nDofpE + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    fe().coorMap( x, y, z, refFE().xi( lDof - 1 ), refFE().eta( lDof - 1 ), refFE().zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        UInt iDof = icmp * dim() + dof().localToGlobal( elemId, lDof  );
                        vect.checkAndSet( iDof, fct( time, x, y, z, icmp + 1 ));
                    }
                }
            }
        }

        // Face based Dof
        if ( nDofpF )
        {

            // loop on element faces
            for ( ID iFa = 1; iFa <= nElemF; ++iFa )
            {

                // Loop on number of Dof per face
                for ( ID l = 1; l <= nDofpF; ++l )
                {

                    lDof = nDofElemE + nDofElemV + ( iFa - 1 ) * nDofpF + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    fe().coorMap( x, y, z, refFE().xi( lDof - 1 ), refFE().eta( lDof - 1 ), refFE().zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        vect( icmp * dim() + dof().localToGlobal( elemId, lDof ) - 0 ) = fct( time, x, y, z, icmp + 1 );
                    }
                }
            }
        }
//         // Element based Dof
//         // Loop on number of Dof per Element
//         for ( ID l = 1; l <= nDofpEl; ++l )
//         {
//             lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

//             // Nodal coordinates
//             fe().coorMap( x, y, z, refFE().xi( lDof - 1 ), refFE().eta( lDof - 1 ), refFE().zeta( lDof - 1 ) );

//             // Loop on data vector components
//             for ( UInt icmp = 0; icmp < nbComp; ++icmp )
//             {
// //                vect( icmp * dim() + dof().localToGlobal( elemId, lDof ) - 0 ) = fct( time, x, y, z, icmp + 1 );
//             }
//         }
    }
}


template <typename Mesh, typename Map>
template<typename vector_type>
void
FESpace<Mesh, Map>::interpolateBC( BCHandler& BCh,
                                   vector_type&    vect,
                                   const Real      time)
{
    //
    // ID   nbComp   = M_fieldDim; // Number of components of the mesh velocity
    UInt totalDof = dof().numTotalDof();

    //
    ID idDof;

    if ( !BCh.bdUpdateDone() )
    {
        BCh.bdUpdate( *mesh(), feBd(), dof() );
    }


    for ( UInt ibc = 0; ibc < BCh.size(); ++ibc )
    {
        if (BCh[ ibc ].type() == Essential)
        {
            // Number of components involved in this boundary condition
            UInt nComp = BCh[ibc].numberOfComponents();

            //


            for ( ID i = 1; i <= BCh[ibc].list_size(); ++i )
            {
                // Coordinates of the node where we impose the value
                Real x = static_cast< const IdentifierEssential* >( BCh[ibc]( i ) ) ->x();
                Real y = static_cast< const IdentifierEssential* >( BCh[ibc]( i ) ) ->y();
                Real z = static_cast< const IdentifierEssential* >( BCh[ibc]( i ) ) ->z();

                for ( ID j = 1; j <= nComp; ++j )
                {
                    // Global Dof
                    idDof = BCh[ibc]( i ) ->id() + ( BCh[ibc].component( j ) - 1 ) * totalDof;
                    Real val = BCh[ibc]( time, x, y, z, BCh[ibc].component( j ) );

                    vect.checkAndSet(idDof,val);
                }
            }
        }
    }



}



// this computes vec = \int fct phi_i
template <typename Mesh, typename Map>
template<typename vector_type>
void
FESpace<Mesh, Map>::L2ScalarProduct( const Function& fct, vector_type& vec, const Real t)
{

    for ( UInt iVol = 1; iVol <= this->mesh()->numElements(); iVol++ )
    {
        this->fe().updateFirstDeriv( this->mesh()->element( iVol ) );

        Real f, x, y, z;

        Int i, inod, ig;
        UInt eleID = this->fe().currentLocalId();
        UInt ic;
        Real u_ig;

        for ( ic = 0; ic < M_fieldDim; ic++ )
        {
            for ( ig = 0; ig < this->fe().nbQuadPt; ig++ )
            {
                this->fe().coorQuadPt( x, y, z, ig );
                f = fct( t, x, y, z, ic + 1 );
                u_ig = 0.;
                for ( i = 0;i < this->fe().nbNode; i++ )
                {
                    inod = this->dof().localToGlobal( eleID, i + 1 ) + ic * dim();
                    u_ig = f*this->fe().phi( i, ig );
                    vec.sumIntoGlobalValues(inod,u_ig*this->fe().weightDet( ig ));
                }

            }
        }
    }

    vec.GlobalAssemble();
}



template <typename Mesh, typename Map>
template<typename vector_type>
Real
FESpace<Mesh, Map>::L20Error( const Function& fexact,
                              const vector_type& vec,
                              const Real time,
                              Real* relError )
{
    Real sum2      = 0.;
    Real sum1      = 0.;
    Real sum0      = 0.;
    Real sumExact2 = 0.;
    Real sumExact1 = 0.;

    for ( UInt iVol = 1; iVol <= this->mesh()->numElements(); iVol++ )
    {
        this->fe().updateFirstDeriv( this->mesh()->element( iVol ) );
        sum2 += elem_L2_diff_2( vec, fexact, this->fe(), this->dof(), time, M_fieldDim );
        sum1 += elem_integral_diff( vec, fexact, this->fe(), this->dof(), time, M_fieldDim );
        sum0 += this->fe().measure();
        if (relError)
        {
            sumExact2 += elem_f_L2_2( fexact, this->fe(), time, 1 );
            sumExact1 += elem_integral( fexact, this->fe(), time, 1 );
        }
    }


    Real sendbuff[5] = {sum0, sum1, sum2, sumExact1, sumExact2};
    Real recvbuff[5];

    this->map().Comm().SumAll(&sendbuff[0], &recvbuff[0], 5);


    //int me = this->map().Comm().MyPID();

//     if (me == 0)
//     {
    sum0      = recvbuff[0];
    sum1      = recvbuff[1];
    sum2      = recvbuff[2];
    sumExact1 = recvbuff[3];
    sumExact2 = recvbuff[4];

    Real absError = sqrt( sum2 - sum1*sum1/sum0 );

    if (relError)
    {
        Real normExact = sqrt( sumExact2 - sumExact1*sumExact1/sum0 );
        *relError = absError / normExact;
    }
    return absError;
}



    template <typename Mesh, typename Map>
    template<typename function, typename vector_type>
    Real
    FESpace<Mesh, Map>::L2Error( const function&    fexact,
                                 const vector_type& vec,
                                 const Real       time,
                                 Real*            relError )
    {
        Real normU       = 0.;
        Real sumExact    = 0.;

        for ( UInt iVol  = 1; iVol <= this->mesh()->numElements(); iVol++ )
        {
            this->fe().updateFirstDeriv( this->mesh()->element( iVol ) );

            normU += elem_L2_diff_2( vec, fexact,
                                     this->fe(),
                                     //p1,
                                     this->dof(),
                                     time,
                                     M_fieldDim );
            if (relError)
            {
                sumExact += elemL22( fexact,
                                     this->fe(),
                                     //p1,
                                     time,
                                     M_fieldDim );
            }
        }


    Real sendbuff[2] = {normU, sumExact};
    Real recvbuff[2];

    this->map().Comm().SumAll(&sendbuff[0], &recvbuff[0], 2);


//    int me = this->map().Comm().MyPID();

//     if (me == 0)
//     {
    normU    = recvbuff[0];
    sumExact = recvbuff[1];

    if (relError)
    {
        *relError = sqrt( normU / sumExact );
    }
//     }

    return sqrt( normU );
}


    template <typename Mesh, typename Map>
     template<typename function, typename vector_type>
     Real
     FESpace<Mesh, Map>::H1Error( const function&    fexact,
                                  const vector_type& vec,
                                  const Real       time,
                                  Real*            relError )
     {
         Real normU       = 0.;
         Real sumExact    = 0.;

         for ( UInt iVol  = 1; iVol <= this->mesh()->numElements(); iVol++ )
         {
             this->fe().updateFirstDeriv( this->mesh()->element( iVol ) );

             normU += elem_H1_diff_2( vec, fexact,
                                      this->fe(),
                                      //p1,
                                      this->dof(),
                                      time,
                                      M_fieldDim );
             if (relError)
             {
                 sumExact += elem_H1_2( fexact,
                                      this->fe(),
                                      //p1,
                                      time,
                                      M_fieldDim );
             }
         }


     Real sendbuff[2] = {normU, sumExact};
     Real recvbuff[2];

     this->map().Comm().SumAll(&sendbuff[0], &recvbuff[0], 2);


 //    int me = this->map().Comm().MyPID();

 //     if (me == 0)
 //     {
     normU    = recvbuff[0];
     sumExact = recvbuff[1];

     if (relError)
     {
         *relError = sqrt( normU / sumExact );
     }
 //     }

     return sqrt( normU );
 }




// template <typename Mesh, typename Map>
// template<typename vector_type>
// Real
// FESpace<Mesh, Map>::L2Norm(vector_type vec)
// {

//     Real norm = 0.;
//     Real tmp  = 0.;

//    for ( UInt ielem = 1; ielem <= this->mesh()->numElements(); ielem++ )
//     {
//         //UInt elem = M_FESpace.mesh()->element( ielem ).id();
//         this->fe().updateFirstDerivQuadPt( M_FESpace.mesh()->element( ielem ) );

//         //int    marker = M_FESpace.mesh()->element( ielem ).marker();
//         Real s      = 0;
//         Real volume = M_FESpace.fe().detJac(0);

//         for ( int ig = 0; ig < M_FESpace.fe().nbQuadPt; ++ig )
//         {
//             for ( int k = 0; k < M_FESpace.fe().nbNode; ++k )
//             {
//                 for ( int j = 0; j < M_fDim; ++j)
//                 {
//                     int i    = M_FESpace.fe().patternFirst(k);
//                     int idof = M_FESpace.dof().localToGlobal(M_FESpace.fe().currentLocalId(), i + 1);

//                     tmp = this->fe().phi( k, 0 , ig )*
//                         vec[idof + j*this->dim()];

//                     norm += this->fe().weigh(ig)*tmp*tmp;
//                 }
//             }
//         }


//         return norm;

// }

template <typename Mesh, typename Map>
template<typename vector_type>
Real
FESpace<Mesh, Map>::L2Norm( const vector_type& vec)
{
	//
    ID nbComp = M_fieldDim; // Number of components of the mesh velocity
    //
	Real norm = 0.;
	//
	for ( UInt ielem = 1; ielem <= this->mesh()->numElements(); ielem++ )
	{
		//UInt elem = M_FESpace.mesh()->element( ielem ).id();
		this->fe().updateJacQuadPt( this->mesh()->element( ielem ) );
		//
		norm += elem_L2_2( vec, this->fe(), this->dof(), nbComp );
	}

    Real sendbuff[1] = {norm};
    Real recvbuff[1];

    this->map().Comm().SumAll(&sendbuff[0], &recvbuff[0], 1);

    norm    = recvbuff[0];

    return sqrt( norm );

}


template <typename Mesh, typename Map>
template<typename vector_type>
Real
FESpace<Mesh, Map>::H1Norm(const vector_type& vec)
{
	//
    ID nbComp = M_fieldDim; // Number of components of the mesh velocity
    //
	Real norm = 0.;
	//
	for ( UInt ielem = 1; ielem <= this->mesh()->numElements(); ielem++ )
	{
		//UInt elem = M_FESpace.mesh()->element( ielem ).id();
		this->fe().updateFirstDerivQuadPt( this->mesh()->element( ielem ) );
		//
//		for ( UInt j = 0; j < nbComp; ++j)
		norm += elem_H1_2( vec, this->fe(), this->dof(), nbComp );
	}

    Real sendbuff[1] = {norm};
    Real recvbuff[1];

    this->map().Comm().SumAll(&sendbuff[0], &recvbuff[0], 1);

    norm    = recvbuff[0];

    return sqrt( norm );

}


 

//This function interpolates the value of the function given by
// - an element
// - the solution vector
// in a given point P: value = \sum_{dofs} solution_i \phi_i(P)
//
// Another faster way would have been to use a custom quadrature
// but the fact that the quadrature is already written inside the 
// refEle destroys all the hope of a faster implementation (the
// geoMap has to be redefined, and so has to be updated for the
// mesh, the refFE has to be rebuilt,...)

template<typename Mesh, typename Map>
template<typename point_type,typename vector_type>
Real
FESpace<Mesh, Map>::
FEinterpolateValue(const ID& elementID, const vector_type& solutionVector, const point_type& pt ) const
{
  // Assumptions for the types:
  // point_type has the STL-type [] accessor and the size() method

  // The vector has to be repeated
  ASSERT( solutionVector.getMaptype() == Repeated, " Vector must be repeated ");

  // Make sur everything is up to date
  M_fe->updateFirstDeriv( M_mesh->element( elementID ));

  // Map the point back to the ref FE    
  Real hat_x(0);  Real hat_y(0);  Real hat_z(0);
  Real x(0);  Real y(0);  Real z(0);

  if (pt.size()>=1)    x=pt[0];
  if (pt.size()>=2)    y=pt[1];
  if (pt.size()>=3)    z=pt[2];

  
  M_fe->coorBackMap(x,y,z,hat_x,hat_y,hat_z);

  // Store the number of local DoF
  UInt nDof(dof().numLocalDof());

  // Initialization
  Real value(0);

  // Loop over the DoFs
  for (UInt iter_dof(0); iter_dof<nDof ; ++iter_dof)
    {
      // The global ID of the selected dof
      // This should be changed in case of a VECTORIAL FE
      
      ID dofID(fe().patternFirst(iter_dof));
      ID globalDofID( dof().localToGlobal(elementID, dofID +1) ); // iter_dof -> dofID

      // Make the accumulation
      value += solutionVector[globalDofID] * M_refFE->phi(iter_dof, hat_x, hat_y, hat_z);
    };

  return value;
};



// Same function as before, with the difference that
// the solution vector is LOCAL (before it was the GLOBAL) one.

template<typename Mesh, typename Map>
template<typename point_type,typename vector_type>
Real
FESpace<Mesh, Map>::
FEinterpolateValueLocal(const ID& elementID, const vector_type& solutionVector, const point_type& pt ) const
{
  // Assumptions for the types:
  // point_type has the STL-type [] accessor and the size() method
  // vector_type has the STL_type [] accessor


  // Make sur everything is up to date
  M_fe->updateFirstDeriv( M_mesh->element( elementID ));

  // Map the point back to the ref FE    
  Real hat_x(0);  Real hat_y(0);  Real hat_z(0);

  Real x(0);  Real y(0);  Real z(0);

  if (pt.size()>=1)    x=pt[0];
  if (pt.size()>=2)    y=pt[1];
  if (pt.size()>=3)    z=pt[2];
  
  M_fe->coorBackMap(x,y,z,hat_x,hat_y,hat_z);
     
  // Store the number of local DoF
  UInt nDof(dof().numLocalDof());

  // Initialization
  Real value(0);

  // Loop over the DoFs
  for (UInt iter_dof(1); iter_dof<=nDof ; ++iter_dof)
    {
      // Make the accumulation
      value += solutionVector[iter_dof-1] * M_refFE->phi(iter_dof-1, hat_x, hat_y, hat_z);
    };

  return value;
};


template<typename Mesh, typename Map>
template<typename point_type,typename vector_type>
Real
FESpace<Mesh, Map>::
FEinterpolateGradient(const ID& elementID, const vector_type& solutionVector, const point_type& pt,
		      const unsigned int& component ) const
{
  // Assumptions for the types:
  // point_type has the STL-type [] accessor and the size() method


  // Make sur everything is up to date
  M_fe->updateFirstDeriv( M_mesh->element( elementID ));

  // Map the point back to the ref FE    
  Real hat_x(0);  Real hat_y(0);  Real hat_z(0);

  Real x(0);  Real y(0);  Real z(0);

  if (pt.size()>=1)    x=pt[0];
  if (pt.size()>=2)    y=pt[1];
  if (pt.size()>=3)    z=pt[2];
  
  M_fe->coorBackMap(x,y,z,hat_x,hat_y,hat_z);
  
  // Store the number of local DoF
  UInt nDof(dof().numLocalDof());

  // Initialization
  Real grad(0);

  // Loop over the DoFs
  for (UInt iter_dof(1); iter_dof<=nDof ; ++iter_dof)
    {
      // The global ID of the selected dof
      // This should be changed in case of a VECTORIAL FE
      ID globalDofID( dof().localToGlobal(elementID, iter_dof) );

      for (UInt iter_dim(0); iter_dim<3; ++iter_dim)
      {
	grad+= solutionVector(globalDofID) * M_refFE->dPhi(iter_dof-1,iter_dim,hat_x,hat_y,hat_z) 
	                                   * M_fe->pointInverseJacobian(hat_x,hat_y,hat_z,component,iter_dim);
      };
    };

  return grad;
};

template<typename Mesh, typename Map>
template<typename point_type,typename vector_type>
Real
FESpace<Mesh, Map>::
FEinterpolateGradientLocal(const ID& elementID, const vector_type& solutionVector, const point_type& pt,
			   const unsigned int& component ) const
{
  // Assumptions for the types:
  // point_type has the STL-type [] accessor and the size() method
  // vector_type has the STL_type [] accessor


  // Make sur everything is up to date
  M_fe->updateFirstDeriv( M_mesh->element( elementID ));

  // Map the point back to the ref FE    
  Real hat_x(0);  Real hat_y(0);  Real hat_z(0);

  Real x(0);  Real y(0);  Real z(0);

  if (pt.size()>=1)    x=pt[0];
  if (pt.size()>=2)    y=pt[1];
  if (pt.size()>=3)    z=pt[2];
  
  M_fe->coorBackMap(x,y,z,hat_x,hat_y,hat_z);
     
  // Store the number of local DoF
  UInt nDof(dof().numLocalDof());

  // Initialization
  Real grad(0);


  // Loop over the DoFs
  for (UInt iter_dof(1); iter_dof<=nDof ; ++iter_dof)
    {
     
      for (UInt iter_dim(0); iter_dim<3; ++iter_dim)
      {
	grad+= solutionVector[iter_dof-1] * M_refFE->dPhi(iter_dof-1,iter_dim,hat_x,hat_y,hat_z) 
	                                   * M_fe->pointInverseJacobian(hat_x,hat_y,hat_z,component,iter_dim);
      };
    };


  return grad;
};


} // end of the namespace
#endif
