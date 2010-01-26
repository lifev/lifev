//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 * @file
 * @brief This files contains the description and the implementation of the FESpace class.
 *
 * @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 * @date 00-06-2007
 *
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


//! FESpace - Short description here please!
/*!
 *  @author Gilles Fourestey
 *
 *  Here write a long and detailed description of the class.
 *
 */

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


    //! @name Methods
    //@{

    //! This method computes the interpolate value of a given FE function in a given point.
    /*!
     * The user of this function has to provide the element and the vector of the DOF values.
     * Given these informations and a point P, the function compute the value:
     * value = sum_{dof i} v(i) phi_i(P)
     * 
     * Note that the point P can be outside of the element considered (this is NOT checked).
     *
     * Warning: this method has been only checked in 3D.
     *
     *  @param elementID  The ID of the element considered. The ID is the local ID
     * (not the global one, see in the MeshEntity class).
     *  @param solutionVector  The vector containing the values in all nodes (not only of
     * the considered element).
     *  @param pt  The point where to perform the interpolation. Note that pt must allow 
     * for STL-type accessor [] and must have the size() method (typically a std::vector<Real>).
     *  @param component  The component for which the interpolation has to be performed
     * (usefull only for the vectorial FE, set to 0 if the FE is scalar).
     *  @return  The interpolated value in the point.
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template<typename point_type, typename vector_type>
    Real FEinterpolateValue(const ID& elementID, const vector_type& solutionVector, 
                            const point_type& pt, const UInt& component = 0) const;


    //! This method computes the interpolate value of a given FE function in a given point.
    /*!
     * This method is the same as FEinterpolateValue, but for the definition of the solution
     * vector. Here, the solution is given only for the degrees of freedom of the given
     * element. The parameter solutionVector is typically a std::vector<Real> containing
     * the values in the dofs.
     *
     * It is not possible to specify the component, the user has to provide directly the
     * values for the wanted component in the solutionVector.
     *
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template<typename point_type, typename vector_type>
    Real FEinterpolateValueLocal(const ID& elementID, const vector_type& solutionVector,
                                 const point_type& pt ) const;


    //! This method computes the interpolated gradient of a given FE function in a given point.
    /*!
     * This method is the same as FEinterpolateValue, but it computes the gradient insted of
     * the value. Therefor, the desired element of the gradient has to be expressed using the
     * parameter gradientElement.
     *
     * For example, if gradientElement=i and component=j, then the results corresponds to the
     * value of partial u_j / partial x_i.
     *
     * Warning: This method has not been tested deeply, buggy behaviour is possible.
     *
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template<typename point_type, typename vector_type>
    Real FEinterpolateGradient(const ID& elementID, const vector_type& solutionVector, const point_type& pt,
			       const UInt& gradientElement, const UInt& component=0 ) const;

    
    //! This method computes the interpolated gradient of a given FE function in a given point.
    /*!
     * This method is the same as FEinterpolateGradient, but it works with local values only,
     * just as FEinterpolateValueLocal.
     *
     * Warning: This method has not been tested deeply, buggy behaviour is possible.
     *
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template<typename point_type, typename vector_type>
    Real FEinterpolateGradientLocal(const ID& elementID, const vector_type& solutionVector, const point_type& pt,
				    const UInt& gradientElement ) const;


    //! This method enables to pass a solution from one FESpace to the present one.
    /*!
     * This method interpolates the values of a FE solution (given by the vector and the FESpace)
     * in the dofs of the present FESpace: we note \f$ v_i \f$ the component of the vector given in
     * argument and \f$ phi_i \f$ the basis functions of the FESpace given in argument. This defines a 
     * function on the domain given by: \f$ f(x) = \sum_i v_i \phi_i(x) \f$.
     * 
     * We search then the function \f$ g(x) \f$ belonging to the present FESpace such that \f$ g(x_j) = f(x_j) \f$
     * for all \f$ x_j \f$ dofs of the present FESpace. This function returns the vector  \f$ w_j \f$ such that
     * \f$ g(x) = sum_j w_j psi_j(x) \f$ where \f$ psi_j \f$ are the basis functions of the present space.
     *
     * Warning: It is NOT true in general that \f$ w_j = f(x_j) \f$ (it is for lagrangian FEs). For 
     * example, if we want to pass the solution from the P_1 to the P_1 with bubbles, then
     * the value associated to the bubble function will be always 0 even if the value of f
     * is not 0 at the location of the dof of the bubble!
     * 
     *  @param originalSpace  The space were the solution is defined originally.
     *  @param originalVector  The vector of the solution in the original space.
     *  @return The vector in the current FESpace. The vector is Unique (unless it is actually
     * the same space so that no interpolation is performed).
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template <typename vector_type>
    vector_type FeToFeInterpolate(const FESpace<mesh_type,map_type>& originalSpace,
    const vector_type& originalVector) const;


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


    BasePattern::PatternType patternType();

private:

    //! copy constructor
    FESpace( const FESpace& fespace );

     //! @name Private Methods
    //@{
    
     //! Set space
    inline void setSpace( const std::string& space );
 

    //! This is a specialized function called by FeToFeInterpolate method for P_2 to P_1 interpolation
    /*!
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template<typename vector_type>
    vector_type P2ToP1Interpolate(const FESpace<mesh_type,map_type>& originalSpace,
				  const vector_type& originalVector) const;


    //! This is a specialized function called by FeToFeInterpolate method for P_1 Bubble to P_1 interpolation
    /*!
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template<typename vector_type>
    vector_type P1bToP1Interpolate(const FESpace<mesh_type,map_type>& original_space,
				  const vector_type& original_vector) const;
    
    //! This is a specialized function called by FeToFeInterpolate method for P_1 to P_2 interpolation
    /*!
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template<typename vector_type>
    vector_type P1ToP2Interpolate(const FESpace<mesh_type,map_type>& original_space,
				  const vector_type& original_vector) const;
    

    //! This is a specialized function called by FeToFeInterpolate method for P_1 to P_1 Bubble interpolation
    /*!
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template<typename vector_type>
    vector_type P1ToP1bInterpolate(const FESpace<mesh_type,map_type>& original_space,
				   const vector_type& original_vector) const;


    //! This is a specialized function called by FeToFeInterpolate method for P_1 Bubble to P_2 interpolation
    /*!
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */
    
    template<typename vector_type>
    vector_type P1bToP2Interpolate(const FESpace<mesh_type,map_type>& original_space,
				   const vector_type& original_vector) const;


    //! This is a specialized function called by FESpace::FeToFeInterpolate method for P_1 Bubble to P_2 interpolation
    /*!
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template<typename vector_type>
    vector_type P2ToP1bInterpolate(const FESpace<mesh_type,map_type>& original_space,
				   const vector_type& original_vector) const;
    //@}




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

    UInt nDofpV    = refFE().nbDofPerVertex(); // number of Dof per vertex
    UInt nDofpE    = refFE().nbDofPerEdge();   // number of Dof per edge
    UInt nDofpF    = refFE().nbDofPerFace();   // number of Dof per face
    UInt nDofpEl   = refFE().nbDofPerVolume(); // number of Dof per Volume

    UInt nElemV    = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE    = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF    = GeoShape::numFaces;    // Number of element's faces

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element
    UInt nDofElemF = nElemF * nDofpF; // number of face's Dof on a Element

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
/*         for ( ID l = 1; l <= nDofpEl; ++l )
         {
             lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

             // Nodal coordinates
             fe().coorMap( x, y, z, refFE().xi( lDof - 1 ), refFE().eta( lDof - 1 ), refFE().zeta( lDof - 1 ) );

             // Loop on data vector components
             for ( UInt icmp = 0; icmp < nbComp; ++icmp )
             {
                vect( icmp * dim() + dof().localToGlobal( elemId, lDof ) - 0 ) = fct( time, x, y, z, icmp + 1 );
             }
	     }*/
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


// Warning! When using this function in parallel, the vector vec HAS TO BE REPEATED!
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



template<typename Mesh, typename Map>
template<typename point_type,typename vector_type>
Real
FESpace<Mesh, Map>::
FEinterpolateValue(const ID& elementID, const vector_type& solutionVector, const point_type& pt, const UInt& component ) const
{
    // The vector has to be repeated, so if it is not, we make is repeated and call this function again.
    if (solutionVector.getMaptype() != Repeated )
    {
        vector_type repeatedSolutionVector(solutionVector,Repeated);
        return FEinterpolateValue(elementID, repeatedSolutionVector, pt, component);
    };
    
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
    UInt totalDof(dof().numTotalDof());
    
    // Initialization
    Real value(0);
    
    // Loop over the DoFs
    for (UInt iter_dof(0); iter_dof<nDof ; ++iter_dof)
    {
        // The global ID of the selected dof
        // This should be changed in case of a VECTORIAL FE
        
        ID globalDofID(component*totalDof + dof().localToGlobal(elementID, iter_dof +1) ); // iter_dof -> dofID
        
        // Make the accumulation
        value += solutionVector[globalDofID] * M_refFE->phi(iter_dof, hat_x, hat_y, hat_z);
    };
    
    return value;
};


template<typename Mesh, typename Map>
template<typename point_type,typename vector_type>
Real
FESpace<Mesh, Map>::
FEinterpolateValueLocal(const ID& elementID, const vector_type& solutionVector, const point_type& pt ) const
{
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
        // Make the accumulation
        value += solutionVector[iter_dof] * M_refFE->phi(iter_dof, hat_x, hat_y, hat_z);
    };
    
    return value;
};


template<typename Mesh, typename Map>
template<typename point_type,typename vector_type>
Real
FESpace<Mesh, Map>::
FEinterpolateGradient(const ID& elementID, const vector_type& solutionVector, const point_type& pt,
                      const UInt& gradientElement, const UInt& component ) const
{
    // The vector has to be repeated, so if it is not, we make is repeated and call this function again.
    if (solutionVector.getMaptype() != Repeated )
    {
        vector_type repeatedSolutionVector(solutionVector,Repeated);
        return FEinterpolateGradient(elementID, repeatedSolutionVector, pt, gradientElement, component);
    };
    
    
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
    UInt totalDof(dof().numTotalDof());
    
    // Initialization
    Real grad(0);
    
    // Loop over the DoFs
    for (UInt iter_dof(0); iter_dof<nDof ; ++iter_dof)
    {
        // The global ID of the selected dof
        ID globalDofID(totalDof * component + dof().localToGlobal(elementID, iter_dof+1) );
        
        for (UInt iter_dim(0); iter_dim<3; ++iter_dim)
        {
            grad+= solutionVector(globalDofID) * M_refFE->dPhi(iter_dof,iter_dim,hat_x,hat_y,hat_z)
                * M_fe->pointInverseJacobian(hat_x,hat_y,hat_z,gradientElement,iter_dim);
        };
    };
    
    return grad;
};


template<typename Mesh, typename Map>
template<typename point_type,typename vector_type>
Real
FESpace<Mesh, Map>::
FEinterpolateGradientLocal(const ID& elementID, const vector_type& solutionVector, const point_type& pt,
			   const UInt& gradientElement ) const
{
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
    for (UInt iter_dof(0); iter_dof<nDof ; ++iter_dof)
    {
        
        for (UInt iter_dim(0); iter_dim<3; ++iter_dim)
        {
            grad+= solutionVector[iter_dof] * M_refFE->dPhi(iter_dof,iter_dim,hat_x,hat_y,hat_z)
                * M_fe->pointInverseJacobian(hat_x,hat_y,hat_z,gradientElement,iter_dim);
        };
    };
    
    return grad;
};


template<typename Mesh, typename Map>
template<typename vector_type>
vector_type
FESpace<Mesh,Map>::
FeToFeInterpolate(const FESpace<mesh_type,map_type>& OriginalSpace,
		  const vector_type& OriginalVector) const
{
  // This method just check that everything is alright and then call the
  // appropriate method __To__Interpolate if needed.

  // First, check that the interpolation is possible
  ASSERT(fieldDim() == OriginalSpace.fieldDim(),"Incompatible field dimension for interpolation");
  ASSERT(refFE().shape == OriginalSpace.refFE().shape , "Incompatible element shape for interpolation");

  // If the spaces are the same, just return the original vector
  if (refFE().type == OriginalSpace.refFE().type)
  {
    return OriginalVector;
  };

  // Then, check that the original vector is repeated.
  // If not, make it repeated and recall this method.
  if ( OriginalVector.getMaptype() == Unique )
  {
    const vector_type OriginalRepeated(OriginalVector,Repeated); 
    return FeToFeInterpolate(OriginalSpace,OriginalRepeated);
  };

  // Distinguish the other cases

  if (refFE().type == FE_P1_3D)
  {
    if (OriginalSpace.refFE().type == FE_P1bubble_3D)
    {
      return P1bToP1Interpolate(OriginalSpace,OriginalVector);
    } 
    else if (OriginalSpace.refFE().type == FE_P2_3D)
    {
      return P2ToP1Interpolate(OriginalSpace,OriginalVector);
    }
    else
    {
      ERROR_MSG(" The interpolation with this host space has not been yet implemented. Please, add it!");
    };
  }
  else if (refFE().type == FE_P1bubble_3D)
  {
    if (OriginalSpace.refFE().type == FE_P1_3D)
    {
      return P1ToP1bInterpolate(OriginalSpace,OriginalVector);
    } 
    else if (OriginalSpace.refFE().type == FE_P2_3D)
    {
      return P2ToP1bInterpolate(OriginalSpace,OriginalVector);
    }
    else
    {
      ERROR_MSG(" The interpolation with this host space has not been yet implemented. Please, add it!");
    };
  }
  else if (refFE().type == FE_P2_3D)
  {
    if (OriginalSpace.refFE().type == FE_P1bubble_3D)
    {
      return P1bToP2Interpolate(OriginalSpace,OriginalVector);
    } 
    else if (OriginalSpace.refFE().type == FE_P1_3D)
    {
      return P1ToP2Interpolate(OriginalSpace,OriginalVector);
    }
    else
    {
      ERROR_MSG(" The interpolation with this host space has not been yet implemented. Please, add it!");
    };
  }
  else
  {
    ERROR_MSG(" The interpolation with this space has not been yet implemented. Please, add it!");
  };


};


// This method returns a vector with Unique map.
template<typename Mesh, typename Map>
template<typename vector_type>
vector_type
FESpace<Mesh,Map>::
P2ToP1Interpolate(const FESpace<mesh_type,map_type>& OriginalSpace,
		  const vector_type& OriginalVector) const
{
  // Create a zero vector.
  vector_type Interpolated(map(),Repeated);
  Interpolated *= 0.0;
 
  // Some constants (avoid recomputing them each time used)
  UInt FieldDim(fieldDim());
  
  UInt numVolumes(mesh()->numVolumes());
  
  UInt totalDofsOriginal(OriginalSpace.dof().numTotalDof());
  UInt totalDofsPresent(dof().numTotalDof());

  // Loop over the elements to get the values
  for ( ID iElem = 1; iElem <= numVolumes ; ++iElem )
  {
    UInt elemId (mesh()->volume(iElem).localId());
   
    // In the file /lifefem/refEle.hpp, we can see that the dofs pour P1
    // are the first ones of the P2 dofs. We have 4 P1 dofs to report per
    // field dimension.
    
    for (UInt iComponent(0); iComponent< FieldDim; ++iComponent)
    {
      for (UInt iP1dof(0); iP1dof<4; ++iP1dof)
      {
	ID globalDofID_original(iComponent*totalDofsOriginal + OriginalSpace.dof().localToGlobal(elemId, iP1dof +1));
	ID globalDofID_present(iComponent*totalDofsPresent + dof().localToGlobal(elemId, iP1dof +1) );

	Real value =  OriginalVector[globalDofID_original];
	Interpolated[globalDofID_present]= value;
      };
    };
  };

  // Here we do need to use the combine mode "Insert": the default combine mode
  // is "Add" and then, as we pass several times on the same DoF, the values would be added, what would
  // be wrong. Instead, we keep only one value (which anyway should be always the same for continuous
  // finite element).
  vector_type return_vector(Interpolated,Unique,Insert);
  return return_vector;

};


template<typename Mesh, typename Map>
template<typename vector_type>
vector_type
FESpace<Mesh,Map>::
P1bToP1Interpolate(const FESpace<mesh_type,map_type>& OriginalSpace,
		  const vector_type& OriginalVector) const
{
// Create a zero vector.
  vector_type Interpolated(map(),Repeated);
  Interpolated *= 0.0;
 
  // Some constants (avoid recomputing them each time used)
  UInt FieldDim(fieldDim());
  
  UInt numVolumes(mesh()->numVolumes());
  
  UInt totalDofsOriginal(OriginalSpace.dof().numTotalDof());
  UInt totalDofsPresent(dof().numTotalDof());

  // Loop over the elements to get the values
  for ( ID iElem = 1; iElem <= numVolumes ; ++iElem )
  {
    UInt elemId (mesh()->volume(iElem).localId());
   
    // In the file /lifefem/refEle.hpp, we can see that the dofs pour P1
    // are the first ones of the P1Bubble dofs. We have 4 P1 dofs to report per
    // field dimension.
    
    for (UInt iComponent(0); iComponent< FieldDim; ++iComponent)
    {
      for (UInt iP1dof(0); iP1dof<4; ++iP1dof)
      {
	ID globalDofID_original(iComponent*totalDofsOriginal + OriginalSpace.dof().localToGlobal(elemId, iP1dof +1));
	ID globalDofID_present(iComponent*totalDofsPresent + dof().localToGlobal(elemId, iP1dof +1) );

	Real value =  OriginalVector[globalDofID_original];
	Interpolated[globalDofID_present]= value;
      };
    };
  };

  // Here we do need to use the combine mode "Insert": the default combine mode
  // is "Add" and then, as we pass several times on the same DoF, the values would be added, what would
  // be wrong. Instead, we keep only one value (which anyway should be always the same for continuous
  // finite element).
  vector_type return_vector(Interpolated,Unique,Insert);
  return return_vector;

};


template<typename Mesh, typename Map>
template<typename vector_type>
vector_type
FESpace<Mesh,Map>::
P1ToP2Interpolate(const FESpace<mesh_type,map_type>& OriginalSpace,
		  const vector_type& OriginalVector) const
{
 // Create a zero vector.
  vector_type Interpolated(map(),Repeated);
  Interpolated *= 0.0;
 
  // Some constants (avoid recomputing them each time used)
  UInt FieldDim(fieldDim());
  
  UInt numVolumes(mesh()->numVolumes());
  
  UInt totalDofsOriginal(OriginalSpace.dof().numTotalDof());
  UInt totalDofsPresent(dof().numTotalDof());

  // A vector to store the values to set in the dofs
  std::vector<Real> DofValues(10,0.0);

  // Loop over the elements to get the values
  for ( ID iElem = 1; iElem <= numVolumes ; ++iElem )
  {
    UInt elemId (mesh()->volume(iElem).localId());
   
    // In the file /lifefem/refEle.hpp, we can see that the dofs pour P1
    // are the first ones of the P2 dofs. We have then to recompute the values
    // in the "face" dofs of the P2 element.
    
    for (UInt iComponent(0); iComponent< FieldDim; ++iComponent)
    {
      // Get the values in the vertices
      for (UInt iP1dof(0); iP1dof<4; ++iP1dof)
      {
	ID globalDofID_original(iComponent*totalDofsOriginal+OriginalSpace.dof().localToGlobal(elemId,iP1dof +1));

	DofValues[iP1dof]  =  OriginalVector[globalDofID_original];
      };
      
      // Compute the values in the faces (!this vector starts from 0, not from 1 as the dofs numeration)
      DofValues[4] = 0.5*(DofValues[0]+DofValues[1]);
      DofValues[5] = 0.5*(DofValues[1]+DofValues[2]);
      DofValues[6] = 0.5*(DofValues[0]+DofValues[2]);
      DofValues[7] = 0.5*(DofValues[0]+DofValues[3]);
      DofValues[8] = 0.5*(DofValues[1]+DofValues[3]);
      DofValues[9] = 0.5*(DofValues[2]+DofValues[3]);

      // Now set them
      for (UInt iP2dof(0); iP2dof<10; ++iP2dof)
      {
	ID globalDofID_present(iComponent*totalDofsPresent+dof().localToGlobal(elemId,iP2dof +1));

	Interpolated[globalDofID_present] = DofValues[iP2dof];
      };
     
    };
  };

  // Here we do need to use the combine mode "Insert": the default combine mode
  // is "Add" and then, as we pass several times on the same DoF, the values would be added, what would
  // be wrong. Instead, we keep only one value (which anyway should be always the same for continuous
  // finite element).
  vector_type return_vector(Interpolated,Unique,Insert);
  return return_vector;

};

template<typename Mesh, typename Map>
template<typename vector_type>
vector_type
FESpace<Mesh,Map>::
P1ToP1bInterpolate(const FESpace<mesh_type,map_type>& OriginalSpace,
		  const vector_type& OriginalVector) const
{
 // Create a zero vector.
  vector_type Interpolated(map(),Repeated);
  Interpolated *= 0.0;
 
  // Some constants (avoid recomputing them each time used)
  UInt FieldDim(fieldDim());
  
  UInt numVolumes(mesh()->numVolumes());
  
  UInt totalDofsOriginal(OriginalSpace.dof().numTotalDof());
  UInt totalDofsPresent(dof().numTotalDof());

  // A vector to store the values to set in the dofs
  std::vector<Real> DofValues(5,0.0);

  // Loop over the elements to get the values
  for ( ID iElem = 1; iElem <= numVolumes ; ++iElem )
  {
    UInt elemId (mesh()->volume(iElem).localId());
   
    // In the file /lifefem/refEle.hpp, we can see that the dofs pour P1
    // are the first ones of the P1Bubble dofs. The value of the Bubble is
    // zero when interpolating a P1 functions.
    
    for (UInt iComponent(0); iComponent< FieldDim; ++iComponent)
    {
      // Get the values in the vertices
      for (UInt iP1dof(0); iP1dof<4; ++iP1dof)
      {
	ID globalDofID_original(iComponent*totalDofsOriginal+OriginalSpace.dof().localToGlobal(elemId,iP1dof +1));

	DofValues[iP1dof]  =  OriginalVector[globalDofID_original];
      };
      
      // This is the value for the bubble function.
      DofValues[4] = 0.0;
      
      // Now set them
      for (UInt iP1bdof(0); iP1bdof<5; ++iP1bdof)
      {
	ID globalDofID_present(iComponent*totalDofsPresent+dof().localToGlobal(elemId,iP1bdof +1));

	Interpolated[globalDofID_present] = DofValues[iP1bdof];
      };
     
    };
  };

  // Here we do need to use the combine mode "Insert": the default combine mode
  // is "Add" and then, as we pass several times on the same DoF, the values would be added, what would
  // be wrong. Instead, we keep only one value (which anyway should be always the same for continuous
  // finite element).
  vector_type return_vector(Interpolated,Unique,Insert);
  return return_vector;

};

template<typename Mesh, typename Map>
template<typename vector_type>
vector_type
FESpace<Mesh,Map>::
P1bToP2Interpolate(const FESpace<mesh_type,map_type>& OriginalSpace,
		  const vector_type& OriginalVector) const
{
  // Create a zero vector.
  vector_type Interpolated(map(),Repeated);
  Interpolated *= 0.0;
 
  // Some constants (avoid recomputing them each time used)
  UInt FieldDim(fieldDim());
  
  UInt numVolumes(mesh()->numVolumes());
  
  UInt totalDofsOriginal(OriginalSpace.dof().numTotalDof());
  UInt totalDofsPresent(dof().numTotalDof());

  // A vector to store the values to set in the dofs
  std::vector<Real> DofValues(10,0.0);

  // Loop over the elements to get the values
  for ( ID iElem = 1; iElem <= numVolumes ; ++iElem )
  {
    UInt elemId (mesh()->volume(iElem).localId());
   
    // In the file /lifefem/refEle.hpp, we can see that the dofs pour P1
    // are the first ones of the P2 dofs. We have then to recompute the values
    // in the "face" dofs of the P2 element. Moreover, the bubble is zero on the
    // faces of the element. It gives then no contribution.
    
    for (UInt iComponent(0); iComponent< FieldDim; ++iComponent)
    {
      // Get the values in the vertices
      for (UInt iP1dof(0); iP1dof<4; ++iP1dof)
      {
	ID globalDofID_original(iComponent*totalDofsOriginal+OriginalSpace.dof().localToGlobal(elemId,iP1dof +1));

	DofValues[iP1dof]  =  OriginalVector[globalDofID_original];
      };
      
      // Compute the values in the faces (!this vector starts from 0, not from 1 as the dofs numeration)
      DofValues[4] = 0.5*(DofValues[0]+DofValues[1]);
      DofValues[5] = 0.5*(DofValues[1]+DofValues[2]);
      DofValues[6] = 0.5*(DofValues[0]+DofValues[2]);
      DofValues[7] = 0.5*(DofValues[0]+DofValues[3]);
      DofValues[8] = 0.5*(DofValues[1]+DofValues[3]);
      DofValues[9] = 0.5*(DofValues[2]+DofValues[3]);

      // Now set them
      for (UInt iP2dof(0); iP2dof<10; ++iP2dof)
      {
	ID globalDofID_present(iComponent*totalDofsPresent+dof().localToGlobal(elemId,iP2dof +1));

	Interpolated[globalDofID_present] = DofValues[iP2dof];
      };
     
    };
  };

  // Here we do need to use the combine mode "Insert": the default combine mode
  // is "Add" and then, as we pass several times on the same DoF, the values would be added, what would
  // be wrong. Instead, we keep only one value (which anyway should be always the same for continuous
  // finite element).
  vector_type return_vector(Interpolated,Unique,Insert);
  return return_vector;

};

template<typename Mesh, typename Map>
template<typename vector_type>
vector_type
FESpace<Mesh,Map>::
P2ToP1bInterpolate(const FESpace<mesh_type,map_type>& OriginalSpace,
		  const vector_type& OriginalVector) const
{
  // Create a zero vector.
  vector_type Interpolated(map(),Repeated);
  Interpolated *= 0.0;
 
  // Some constants (avoid recomputing them each time used)
  UInt FieldDim(fieldDim());
  
  UInt numVolumes(mesh()->numVolumes());
  
  UInt totalDofsOriginal(OriginalSpace.dof().numTotalDof());
  UInt totalDofsPresent(dof().numTotalDof());

  // A vector to store the values to set in the dofs
  std::vector<Real> DofValues(5,0.0);

  // Loop over the elements to get the values
  for ( ID iElem = 1; iElem <= numVolumes ; ++iElem )
  {
    UInt elemId (mesh()->volume(iElem).localId());
   
    // This is a tricky case. The problem is that the P2 functions
    // have an influence on the gravity center. Moreover, the P1Bubble
    // FE is NOT a Lagrangian FE, so we need to make a combinaison of
    // the values in the different points. But first of all, we extract
    // the P1 values in the vertices.
    
    for (UInt iComponent(0); iComponent< FieldDim; ++iComponent)
    {
      // Get the values in the vertices
      for (UInt iP1dof(0); iP1dof<4; ++iP1dof)
      {
	ID globalDofID_original(iComponent*totalDofsOriginal+OriginalSpace.dof().localToGlobal(elemId,iP1dof +1));

	DofValues[iP1dof]  =  OriginalVector[globalDofID_original];
      };
      
      // Here we have to compute the value of the interpolated function
      // in the gravity center of the cell.
      std::vector<Real> gravityCenter(3,0);
      for (UInt iterVertices(0); iterVertices<4; ++iterVertices)
      {
          gravityCenter[0] += mesh()->volume(iElem).point(iterVertices+1).coordinate(0+1)/4.0;
          gravityCenter[1] += mesh()->volume(iElem).point(iterVertices+1).coordinate(1+1)/4.0;
          gravityCenter[2] += mesh()->volume(iElem).point(iterVertices+1).coordinate(2+1)/4.0;
      };

      Real gravityCenterValue(0);
      gravityCenterValue = FEinterpolateValue(elemId,OriginalVector,gravityCenter,iComponent);
      
      // Here we do the combinaison because it is not Lagrangian.
      DofValues[4] = 256*(gravityCenterValue - (DofValues[0]+DofValues[1]+DofValues[2]+DofValues[3])/4.0);
     
      // Now set them
      for (UInt iP1bdof(0); iP1bdof<5; ++iP1bdof)
      {
	ID globalDofID_present(iComponent*totalDofsPresent+dof().localToGlobal(elemId,iP1bdof +1));

        Interpolated[globalDofID_present] = DofValues[iP1bdof];
      };
     
    };
  };

  // Here we do need to use the combine mode "Insert": the default combine mode
  // is "Add" and then, as we pass several times on the same DoF, the values would be added, what would
  // be wrong. Instead, we keep only one value (which anyway should be always the same for continuous
  // finite element).
  vector_type return_vector(Interpolated,Unique,Insert);
  return return_vector;
};

} // end of the namespace
#endif
