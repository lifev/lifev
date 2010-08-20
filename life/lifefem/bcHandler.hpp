//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2096 EPFL, Politecnico di Milano, INRIA
               2006-2010 EPFL, Politecnico di Milano

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
/**
   @file
   @brief File containing a class for BC handling

   @author Miguel Fernandez
   @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>

   @date 2004-10-11
*/

#ifndef BCHandler_H
#define BCHandler_H 1

#include <life/lifefem/bcCond.hpp>

namespace LifeV {

/*! @page bc_page Boundary Conditions

  @section bc_identifier Conventions

  Here we explain what is the convention for the boundary conditions is LifeV.

  The main issue is to understand what is the flag that are attributed to the degrees of freedom not located on a vertex. In that case, the flag of the smallest entity including the degree of freedom is used.

  This means that:
  <ol>
  <li> A degree of freedom based on a vertex will use the flag of that vertex.
  <li> A degree of freedom based on an edge (but not on a vertex) will use the flag of that edge.
  <li> A degree of freedom based on a face (but not on a vertex or an edge) will use the flag of that face.
  <li> ...
  </ol>

  The user of the library has then the responsability to define the good flags for each mesh entity. There is no priority for any of the boundary conditions. This is an important issue has it can create buggy behaviours that are difficult to trace back. We provide here an example.

  Suppose that we want to keep a fluid inside a cube.

  @image html cube.png "The domain"

  Then, no penetration boundary conditions have to be imposed on all the sides, otherwise the fluid will escape from the container. No penetration boundary conditions are different depending on the face/edge/vertex that we are considering, as shown on the next picture (for a corner of the cube).

   @image html bc_cube_types.png "Boundary conditions to be imposed"

   To impose these conditions, we might have defined in the mesh file flags for all the vertices. Here, we might have defined them as in the next figure.

   @image html bc_cube_flags_I_P1.png "A choice for the flags"

   However, if P2 finite elements are used, the degrees of freedom in the middle of each edge are concerned by the boundary conditions. By default, when the mesh is built, the flag attributed to the edge is the minimum of the flags of the vertices defining it. The flag of the edge is then used as the flag of the P2 degree of freedom, as shown in the next picture.

   @image html bc_cube_flags_I_P2.png "The resulting flags for the P2 degrees of freedom"

   We can see that this will result in a wrong behaviour, as the boundary conditions of the edges are spread on the faces of the cube, what is not what we wanted. If we are aware of the default behaviour, we can reorder the flags to get boundary conditions that are closer to what we wanted, as shown in the next two figures.

   @image html bc_cube_flags_II_P1.png "A better choice for the flags ..."

   @image html bc_cube_flags_II_P2.png "... and the resulting flags for the P2 degrees of freedom"

   However, we see on the last picture that there are still two degrees of freedom with the wrong boundary conditions, near the corner. On the front face, we would like the upper "3" to be "1" as the boundary condition for the face should apply. The library cannot guess that it has to take the flag "1"! This is why in this situation, the user has to provide a flag for this edge in the mesh file.

   @image html bc_cube_flags_III_P2.png "The correct flag repartition"

   <b>Remark</b>: The mesh shown in the previous examples should be avoided, not only because it makes "complicated" boundary conditions for the P2 elements, but also because the tetrahedron in the corner is "locked" if the P1 elements are choosen.
 */


//! BCHandler - The BC handler class
/*!
   @class BCHandler

   Container for @c BCBase classes

   @c BCHandler is a container for the boundary condition (BC) classes. It just
   uses an stl vector to store @c BCBase objets. The usage is simple: the
   user creates the data functors from the user defined functions:

   @verbatim
   BCFunctionBase gv( g );
   BCFunctionMixte gp( h, q );
   @endverbatim
   Then he/she specifies the number of BC and uses the @c BCHandler for
   creating the actual @c BCBase objects and storing them.

   @verbatim
   BCHandler bcV( 3 );
   bcV.addBC( "inlet", 10, Essential, Full, gv, 3 );
   bcV.addBC( "outflow", 11, Mixte, Scalar );
   bcV.addBC( "wall", 10, Essential, Full, gv, 3 );
   @endverbatim
*/
class BCHandler
{
public:

    //! @name Type & enum definitions
    //@{

    typedef std::vector<BCBase>::iterator       BCBase_Iterator;
    typedef std::vector<BCBase>::const_iterator BCBase_ConstIterator;

    /*!
      @typedef enum
      Hints about the nature of the BCs
    */
    enum BCHints
    {
        HINT_BC_NONE,          //!< no hint
        HINT_BC_ONLY_ESSENTIAL //!< BCs are only essential
    };

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
      @param nbc the number of BC to be stored
      @param hint hint about their nature
    */
    BCHandler( const ID& nbc = 0 , const BCHints& hint = HINT_BC_NONE );

    //! Copy constructor
    /*!
     * @param BCh BCHandler
     */
    BCHandler( const BCHandler & BCh );

    //@}


    //! @name Operators
    //@{

    //! assignment operator
    BCHandler& operator= ( const BCHandler& BCh );

    //! extracting a BC in the list
          BCBase& operator[] ( const Index_t& );

    //! extracting a BC in the list
    const BCBase& operator[] ( const Index_t& ) const;

    //@}


    //! @name Methods
    //@{

    //! Add new BC to the list (user defined function)
    /*!
      @param name the name of the boundary condition
      @param flag the mesh flag identifying the part of the mesh where the boundary condtion applies
      @param type the boundary condition type: Natural, Essential, Mixte
      @param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
      @param bcf the function holding the user defined function involved in this boundary condition
      @param std::vector<ID> storing the list of components involved in this boundary condition
    */
    void addBC( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCFunctionBase& bcf,
                const std::vector<ID>& comp );


    //! Add new BC to the list without specified components for Scalar, Tangential or Normal mode problems (user defined function).
    /*!
       @param name the name of the boundary condition
       @param flag the mesh flag identifying the part of the mesh where the boundary condtion applies
       @param type the boundary condition type: Natural, Essential, Mixte
       @param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
       @param bcf the function holding the user defined function involved in this boundary condition
    */
    void addBC( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCFunctionBase& bcf );


    //! Add new BC to the list without list of components for Full mode problems  (user defined function)
    /*!
      @param name the name of the boundary condition
      @param flag the mesh flag identifying the part of the mesh where the boundary condiion applies
      @param type the boundary condition type: Natural, Essential, Mixte
      @param mode the boundary condition mode: Scalar, Full, Normal, Tangential
      @param bcf the function holding the user defined function involved in this boundary condition
      @param nComp the number of componets involved in this boundary condition
    */
    void addBC( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCFunctionBase& bcf,
                const UInt& nComp );


    //! Add new BC to the list (data vector)
    /*!
      @param name the name of the boundary condition
      @param flag the mesh flag identifying the part of the mesh where the boundary condtion applies
      @param type the boundary condition type: Natural, Essential, Mixte
      @param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
      @param bcv data vector
      @param std::vector<ID> storing the list of components involved in this boundary condition
    */
    void addBC( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCVectorBase& bcv,
                const std::vector<ID>& comp );


    //! Add new BC to the list  without specified components for Scalar, Tangential or Normal  mode problemst (data vector)
    /*!
      @param name the name of the boundary condition
      @param flag the mesh flag identifying the part of the mesh where the boundary condtion applies
      @param type the boundary condition type: Natural, Essential, Mixte
      @param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
      @param bcv data vector
    */
    void addBC( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCVectorBase& bcv );


    //! Add new BC to the list without list of components for Full mode problems t (data vector)
    /*!
      @param name the name of the boundary condition
      @param flag the mesh flag identifying the part of the mesh where the boundary condiion applies
      @param type the boundary condition type: Natural, Essential, Mixte
      @param mode the boundary condition mode: Scalar, Full, Normal, Tangential
      @param bcv data vector
      @param nComp the number of componets involved in this boundary condition
    */
    void addBC( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCVectorBase& bcv,
                const UInt& nComp );

    void addBC( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCFunctionUDepBase& bcf );


    //! modify the boundary condition @c name with condition @c bcv
    void modifyBC( std::string const& name, BCFunctionBase& bcv );

    //! modify the boundary condition @c name with condition @c bcv
    void modifyBC( std::string const& name, BCVectorBase& bcv );

    void modifyBC( std::string const& __name, BCFunctionUDepBase& __bcf );

    //! modify the boundary condition @c lab with condition @c bcv
    void modifyBC( Int lab, BCFunctionBase& bcv );

    //! modify the boundary condition @c lab with condition @c bcv
    void modifyBC( Int lab, BCVectorBase& bcv );

    void modifyBC( Int lab, BCFunctionUDepBase& __bcf );

    //! Build the boundary stuff
    /*!
      This method udapte the BC classes checking the markers
      on the boundary. It builds the list of identifiers depending
      on the BC type.
      It updates @c BCHandler objects not @c Dof objects
      @param mesh the mesh
      @param feBd the current finite element on the boundary
      @param dof give the local to global table
    */
    template <typename Mesh>
    void bdUpdate( Mesh& mesh, CurrentBdFE& feBd, const Dof& dof );

    //! merges two different BCHandler (with their offsets) into one
    void merge( BCHandler& bch );

    //! Display the content of the variables
    void showMe( bool verbose = false, std::ostream & out = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! the offset is used on all BCs
    void setOffset( const UInt& offset );

    //! specific BC offset
    void setOffset( const std::string& name, Int offset );

    //@}


    //! @name Get Methods
    //@{

    //! Extract the type of the boundary condition with flag aFlag.
    /*!
     * Useful for weak imposition of boundary conditions (e.g. in DG-FEM).
     */
    BCType boundaryType(const EntityFlag& aFlag) const;

    //! Extract a BC in the list according to its flag
    BCBase& GetBCWithFlag(const EntityFlag& );
    const BCBase& GetBCWithFlag(const EntityFlag&) const;

    //! Get a vector list of BC with specific type. The list contains the BCName of the BC.
    /*!
      @param type - the BC type to be inserted in the list
    */
    std::vector<BCName> getBCWithType( const BCType& type );

    //! Get the number of BC with specific type.
    /*!
      @param type - the BC type to be counted
    */
    UInt getNumberBCWithType( const BCType& type );

    // get the BC number from the list using its name
    UInt getBCbyName(const std::string _BCName) const;

    UInt offset() const;

    //! Iterators for the begining of the BC list
    BCBase_Iterator begin();

    //! Iterators for the end of the BC list
    BCBase_Iterator end();

    //! How many BC stored?
    Index_t size() const;

    //! Is there no BC stored?
    bool empty() const;

    //! returns true if the bdUpdate has been done before
    bool bdUpdateDone() const;

    /*!
      @brief returns true if all the stored BC are of Essential type

      It throws a @c logic_error exception if state is not consistent with
      hint.
    */
    bool hasOnlyEssential() const;

    //@}

private:

    /**
       look for the BC named @c name.

       It throws an @c invalid_argument exception if @c name is not
       found.

       @return a pointer to @c BCBase
    */
    BCBase* findBC( const std::string& __name );

    BCBase* findBC( const Int& lab );

    //! determine whether BCs in M_bcList are only essential
    bool listHasOnlyEssential() const;

    void sumOffsets();

    //! stores status (essential or not) of the boundary associated to a marker
    class EssentialStatus
    {
    public:
        EssentialStatus() : M_normal( false ), M_tangential( false )
        {
            for(UInt iComp=0; iComp<nDimensions; ++iComp)
                M_components[iComp] = false;
        }

        void setComponent(UInt comp) { M_components[comp-1] = true; }

        void setAllComponents()
        {
            for(UInt iComp=0; iComp<nDimensions; ++iComp)
                M_components[iComp] = true;
        }

        void setNormal() { M_normal = true; }
        void setTangential() { M_tangential = true; }
        void setDirectional() {M_directional = true;}

        bool isEssential() const
        {
            bool result = true;
            for (UInt iComp=0; iComp<nDimensions; ++iComp)
                result &= M_components[iComp];
            result = result || ( M_normal && M_tangential );
            return result;
        }

    private:

        bool M_components[NDIM];
        bool M_normal;
        bool M_tangential;
        bool M_directional;
    };

    //! Number of BC to be stored: true if the bdUpdate has been done
    bool M_bdUpdateDone;

    //! hints about the nature of the BCs
    BCHints M_hint;

    //! vector list holding the stored BC
    std::vector<BCBase> M_bcList;

    UInt M_offset;

    //! set of markers which are in the mesh but not in the list
    std::set<EntityFlag> M_notFoundMarkers;

};

// ===================================================
// Template methods implementations
// ===================================================
template <typename Mesh>
void
BCHandler::bdUpdate( Mesh& mesh, CurrentBdFE& feBd, const Dof& dof )
{
	typedef typename Mesh::ElementShape GeoShape;

    // Some useful local variables, to save some typing
        UInt nDofPerVert = dof.fe.nbDofPerVertex(); // number of Dof per vertices
        UInt nDofPerEdge = dof.fe.nbDofPerEdge();   // number of Dof per edges
        UInt nDofPerFace = dof.fe.nbDofPerFace();   // number of Dof per faces

    UInt numBElements = mesh.numBElements();    // number of boundary elements

    EntityFlag marker; //will store the marker of each geometric entity
    EntityFlag elementMarker; //will store the marker of the element

    typedef typename GeoShape::GeoBShape GeoBShape;

    UInt nBElementV = GeoBShape::numVertices; // Number of boundary element's vertices
    UInt nBElementE = GeoBShape::numEdges;    // Number of boundary element's edges

    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges

    UInt nDofBElV = nDofPerVert * nBElementV; // number of vertex's Dof on a boundary element
    UInt nDofBElE = nDofPerEdge * nBElementE; // number of edge's Dof on a boundary element
#ifdef TWODIM
    UInt nDofBEl = nDofBElV + nDofBElE; // number of total Dof on a boundary element
#elif defined THREEDIM
    UInt nDofBEl = nDofBElV + nDofBElE + nDofPerFace; // number of total Dof on a boundary element
#endif
    UInt nDofElemV = nElemV * nDofPerVert; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofPerEdge; // number of edge's Dof on a Element

    SimpleVect<ID> bdltg( nDofBEl );
    typedef std::vector<BCBase>::iterator Iterator;
    Iterator where;
    std::vector<Iterator> whereList;

    UInt iElAd, iVeEl, iBElEl, iEdEl;
    ID lDof, gDof;
    Real x, y, z;

    std::set<EntityFlag> notFoundMarkersCurrent;

    // ===================================================
    // Loop on boundary faces
    // ===================================================
    for ( ID iBoundaryElement = 1 ; iBoundaryElement <= numBElements; ++iBoundaryElement )
        {
            iElAd = mesh.bElement( iBoundaryElement ).ad_first();  // id of the element adjacent to the face
            iBElEl = mesh.bElement( iBoundaryElement ).pos_first(); // local id of the face in its adjacent element

            feBd.updateMeas( mesh.bElement( iBoundaryElement ) );  // updating finite element information
            elementMarker = mesh.bElement( iBoundaryElement ).marker(); // We keep the element marker (added by Gwenol Grandperrin)

            // ===================================================
            // Vertex based Dof
            // ===================================================
            if ( nDofPerVert )
                {

                    // loop on boundary elements vertices
                    for ( ID iVeBEl = 1; iVeBEl <= nBElementV; ++iVeBEl )
                        {

                            marker = mesh.bElement( iBoundaryElement ).point( iVeBEl ).marker(); // vertex marker
#ifdef TWODIM
                            iVeEl = GeoShape::eToP( iBElEl, iVeBEl ); // local vertex number (in element)
#elif defined THREEDIM
                            iVeEl = GeoShape::fToP( iBElEl, iVeBEl ); // local vertex number (in element)
#endif

                            // Finding this marker on the BC list
                            whereList.clear();
                            where = M_bcList.begin();
                            while ( ( where = find( where, M_bcList.end(), marker ) ) != M_bcList.end() )
                                {
                                    whereList.push_back( where );
                                    ++where;
                                }
                            if ( whereList.size() == 0 )
                                {
                                    notFoundMarkersCurrent.insert(marker);
                                }

                            // Loop number of Dof per vertex
                            for ( ID l = 1; l <= nDofPerVert; ++l )
                                {

                                    lDof = ( iVeBEl - 1 ) * nDofPerVert + l ; // local Dof

                                    // global Dof
                                    gDof = dof.localToGlobal( iElAd, ( iVeEl - 1 ) * nDofPerVert + l );
                                    bdltg( lDof ) = gDof; // local to global on this face

                                    // Adding identifier
                                    for ( UInt i = 0 ; i < whereList.size(); ++i )
                                        {
                                            where = whereList[ i ];
                                            //                         std::cout << "type = " << where->type() << std::endl;
                                            switch ( where->type() )
                                                {
                                                case Essential:
                                                    // Which kind of data ?
                                                    if ( where->dataVector() )
                                                        { // With data vector
                                                            where->addIdentifier( new IdentifierBase( gDof ) ); // We only need the dof number
                                                        }
                                                    else
                                                        { // With user defined functions
                                                            feBd.coorMap( x, y, z, feBd.refFE.xi( lDof - 1 ), feBd.refFE.eta( lDof - 1 ) );
                                                            where->addIdentifier( new IdentifierEssential( gDof, x, y, z ) );
                                                        }
                                                    break;
                                                case Natural:
                                                    if ( where->dataVector() )
                                                        { // With data
                                                            switch ( where->pointerToBCVector() ->type() )
                                                                {
                                                                case 0:
                                                                    // if the BC is a function or a vector which values
                                                                    // don't need to be integrated
                                                                    //std::cout << "adding natural identifier" << std::endl;
                                                                    where->addIdentifier( new IdentifierNatural( gDof ) );
                                                                    break;
                                                                case 1:  // if the BC is a vector of values to be integrated
                                                                    break;
                                                                case 2:  // if the BC is a vector of values to be integrated
                                                                    break;
                                                                default:
                                                                    ERROR_MSG( "This boundary condition type is not yet implemented" );
                                                                }
                                                        }
                                                    else
                                                        {}
                                                    break;
                                                case Mixte:
                                                    // Why kind of data ?
                                                    // vincent please check again for your Mixte-FE it doesn't work for Q1
                                                    //       if ( where->dataVector()  ) { // With data vector
                                                    //        where->addIdentifier( new IdentifierNatural(gDof) );
                                                    //       }
                                                    break;
                                                case Flux:
                                                    if ( where->dataVector() )
                                                        { // With data
                                                            switch ( where->pointerToBCVector() ->type() )
                                                                {
                                                                case 0:
                                                                    // if the BC is a function or a vector which values
                                                                    // don't need to be integrated
                                                                    where->addIdentifier( new IdentifierNatural( gDof ) );
                                                                    break;
                                                                case 1:  // if the BC is a vector of values to be integrated
                                                                    where->addIdentifier( new IdentifierNatural( gDof ) );
                                                                    break;
                                                                case 2:  // if the BC is a vector of values to be integrated
                                                                    break;
                                                                default:
                                                                    ERROR_MSG( "This boundary condition type is not yet implemented" );
                                                                }
                                                        }
                                                    else
                                                        {
                                                            //                                        feBd.coorMap( x, y, z, feBd.refFE.xi( lDof - 1 ), feBd.refFE.eta( lDof - 1 ) );
                                                            //                                        where->addIdentifier( new IdentifierEssential( gDof, x, y, z ) );
//                                                             std::cout << "adding flux identifier " << gDof << std::endl;
//                                                             where->addIdentifier( new IdentifierNatural( gDof ) );
                                                        }
                                                  break;
						case Resistance:

						  if ( where->dataVector() )

						  where->addIdentifierIdGlobal( new IdentifierNatural( gDof ) );

						  else
						    { // With user defined functions
						      // feBd.coorMap( x, y, z, feBd.refFE.xi( lDof - 1 ), feBd.refFE.eta( lDof - 1 ) );
						      //where->addIdentifierIdGlobal( new IdentifierEssential( gDof, x, y, z ) );
						    }
                                                    break;
						default:
					            ERROR_MSG( "This boundary condition type is not yet implemented" );
						}
					}
                                }
                        }
                }


            // ===================================================
            // Edge based Dof
            // ===================================================
            if ( nDofPerEdge )
                {
                    // loop on boundary element's edges
                    for ( ID iLocalBElement = 1; iLocalBElement <= nBElementE; ++iLocalBElement )
                        {
#ifdef TWODIM
                            iEdEl = iBElEl; // local edge number (in element)
#elif defined THREEDIM
                            iEdEl = GeoShape::fToE( iBElEl, iLocalBElement ).first; // local edge number (in element)
#endif
                            marker = mesh.boundaryEdge( mesh.localEdgeId( iElAd, iEdEl ) ).marker(); // edge marker
                            //if(marker!= elementMarker){continue;}
                            // Finding this marker on the BC list
                            whereList.clear();
                            where = M_bcList.begin();
                            while ( ( where = find( where, M_bcList.end(), marker ) ) != M_bcList.end() )
                                {
                                    whereList.push_back( where );
                                    ++where;
                                }
                            if ( whereList.size() == 0 )
                                {
                                    notFoundMarkersCurrent.insert(marker);
                                }

                            // Loop number of Dof per edge
                            for ( ID l = 1; l <= nDofPerEdge; ++l )
                                {

                                    lDof = nDofBElV + ( iLocalBElement - 1 ) * nDofPerEdge + l ; // local Dof
                                    gDof = dof.localToGlobal( iElAd, nDofElemV + ( iEdEl - 1 ) * nDofPerEdge + l ); // global Dof
                                    bdltg( lDof ) = gDof; // local to global on this face

                                    // Adding identifier
                                    for ( UInt i = 0 ; i < whereList.size(); ++i )
                                        {
                                            where = whereList[ i ];
                                            switch ( where->type() )
                                                {
                                                case Essential:
                                                    // Which kind of data ?
                                                    if ( where->dataVector() )
                                                        { // With data vector
                                                            where->addIdentifier( new IdentifierBase( gDof ) );
                                                        }
                                                    else
                                                        { // With user defined functions
                                                            feBd.coorMap( x, y, z, feBd.refFE.xi( lDof - 1 ), feBd.refFE.eta( lDof - 1 ) );
                                                            where->addIdentifier( new IdentifierEssential( gDof, x, y, z ) );
                                                        }
                                                    break;
                                                case Natural:
                                                    // Which kind of data ?
                                                    if ( where->dataVector() )
                                                        { // With data vector
                                                            switch ( where->pointerToBCVector() ->type() )
                                                                {
                                                                case 0:
                                                                    // if the BC is a function or a vector which values
                                                                    // don't need to be integrated
                                                                    where->addIdentifier( new IdentifierNatural( gDof ) );
                                                                    break;
                                                                case 1:  // if the BC is a vector of values to be integrated
                                                                    break;
                                                                case 2:  // if the BC is a vector of values to be integrated
                                                                    break;
                                                                default:
                                                                    ERROR_MSG( "This boundary condition type is not yet implemented" );
                                                                }
                                                        }
                                                    break;
                                                case Mixte:
                                                    // Which kind of data ?
                                                    if ( where->dataVector() )
                                                        { // With data vector
                                                            where->addIdentifier( new IdentifierNatural( gDof ) );
                                                        }
                                                    break;
                                                case Flux:
						  break;
						case Resistance:
						  // Which kind of data ?
						  if ( where->dataVector() )
						    { // With data vector
						      where->addIdentifier( new IdentifierNatural( gDof ) );

						      where->addIdentifierIdGlobal(new IdentifierBase( gDof ) );
						    }
						  break;
                                                default:
						  ERROR_MSG( "This boundary condition type is not yet implemented" );
                                                }
                                        }
                                }
                        }
                }


#ifndef TWODIM
            // ===================================================
            // Face based Dof
            // ===================================================
            marker = mesh.bElement( iBoundaryElement ).marker(); // edge marker

            // Finding this marker on the BC list
            whereList.clear();
            where = M_bcList.begin();

            while ( ( where = find( where, M_bcList.end(), marker ) ) != M_bcList.end() )
                {
                    whereList.push_back( where );
                    ++where;
                }
            if ( whereList.size() == 0 )
                {
                    notFoundMarkersCurrent.insert(marker);
                }

            // Adding identifier
            for ( UInt i = 0 ; i < whereList.size(); ++i )
                {
                    where = whereList[ i ];
                    switch ( where->type() )
                        {
                        case Essential:
                            // Loop on number of Dof per face
                            for ( ID l = 1; l <= nDofPerFace; ++l )
                                {
                                    lDof = nDofBElE + nDofBElV + l; // local Dof
                                    gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + ( iBElEl - 1 ) * nDofPerFace + l ); // global Dof
                                    // Why kind of data ?
                                    if ( where->dataVector() )
                                        { // With data vector
                                            where->addIdentifier( new IdentifierBase( gDof ) );
                                        }
                                    else
                                        { // With user defined functions
                                            feBd.coorMap( x, y, z, feBd.refFE.xi( lDof - 1 ), feBd.refFE.eta( lDof - 1 ) );
                                            where->addIdentifier( new IdentifierEssential( gDof, x, y, z ) );
                                        }
                                }
                            break;
                        case Natural:

                            // Why kind of data ?
                            // vincent please check again for your Mixte-FE it doesn't work for Q1
                            if ( where->dataVector() )
                                { // With data vector
                                    UInt type = where->pointerToBCVector()->type() ;
                                    if ( type == 0 )
                                        {
                                            // if the BC is a vector which values don't need to be integrated
                                            for ( ID l = 1; l <= nDofPerFace; ++l )
                                                {
                                                    lDof = nDofBElE + nDofBElV + l; // local Dof
                                                    gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + ( iBElEl - 1 ) * nDofPerFace + l ); // global Dof
                                                    where->addIdentifier( new IdentifierNatural( gDof ) );
                                                }
                                        }
                                    else if ( (type == 1) || (type == 2) )
                                        {
                                            // Loop on number of Dof per face
                                            for ( ID l = 1; l <= nDofPerFace; ++l )
                                                {
                                                    lDof = nDofBElE + nDofBElV + l; // local Dof
                                                    gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + ( iBElEl - 1 ) * nDofPerFace + l ); // global Dof
                                                    std::cout << lDof << " to " << gDof << std::endl;
                                                    bdltg( lDof ) = gDof; // local to global on this face
                                                }
                                            where->addIdentifier( new IdentifierNatural( iBoundaryElement, bdltg ) );
                                        }

                                    else
                                        ERROR_MSG( "This BCVector type is not yet implemented" );

                                }
                            else
                                {
                                    // Loop on number of Dof per face
                                    for ( ID l = 1; l <= nDofPerFace; ++l )
                                        {
                                            lDof = nDofBElE + nDofBElV + l; // local Dof
                                            gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + ( iBElEl - 1 ) * nDofPerFace + l ); // global Dof
                                            bdltg( lDof ) = gDof; // local to global on this face
                                        }
                                    where->addIdentifier( new IdentifierNatural( iBoundaryElement, bdltg ) );
                                }
                            break;
                        case Mixte:
                            // Why kind of data ?
                            // vincent please check again for your Mixte-FE it doesn't work for Q1
                            // if ( where->dataVector()  ) { // With data vector
                            //   for (ID l=1; l<=nDofPerFace; ++l) {
                            //     lDof = nDofBElE + nDofBElV + l; // local Dof
                            //     gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + (iBElEl-1)*nDofPerFace + l); // global Dof
                            //     where->addIdentifier( new IdentifierNatural(gDof) );
                            //   }
                            // }
                            // else {
                            // Loop on number of Dof per face
                            for ( ID l = 1; l <= nDofPerFace; ++l )
                                {
                                    lDof = nDofBElE + nDofBElV + l; // local Dof
                                    gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + ( iBElEl - 1 ) * nDofPerFace + l ); // global Dof
                                    bdltg( lDof ) = gDof; // local to global on this face
                                }
                            where->addIdentifier( new IdentifierNatural( iBoundaryElement, bdltg ) );
                            // }
                            break;
                        case Flux:
//                             std::cout << iBoundaryElement << " = ";
//                             for (Int ii = 0; ii < bdltg.size(); ++ii)
//                                 std::cout << bdltg[ii] << " ";
//                             std::cout << std::endl;
                            where->addIdentifier( new IdentifierNatural( iBoundaryElement, bdltg ) );
                            break;

			case Resistance:
                          if ( where->dataVector()  )
			  {

			   for ( ID l = 1; l <= nDofPerFace; ++l )
			     {
			       lDof = nDofBElE + nDofBElV + l; // local Dof
			       gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + ( iBElEl - 1 ) * nDofPerFace + l ); // global Dof
			       where->addIdentifierIdGlobal(new IdentifierBase( gDof ) );
			       bdltg( lDof ) = gDof; // local to global on this face
			     }
			   where->addIdentifier( new IdentifierNatural( iBoundaryElement, bdltg ) );
			  }
			  break;
			default:
			  ERROR_MSG( "This boundary condition type is not yet implemented" );
                        }
                }

#endif
        }

    std::set<EntityFlag> notFoundMarkersNew;

    for( std::set<EntityFlag>::iterator it = notFoundMarkersCurrent.begin(); it != notFoundMarkersCurrent.end(); ++it )
        {
	  if ( M_notFoundMarkers.find( *it ) == M_notFoundMarkers.end() )
	    {
	      notFoundMarkersNew.insert( *it );
	    }
        }


    if( notFoundMarkersNew.size() > 0 )
      {

#ifdef DEBUG
          Debug(5010) <<
              "WARNING -- BCHandler::bdUpdate()\n" <<
              "  boundary degrees of freedom with the following markers\n" <<
              "  have no boundary condition set: ";
	for( std::set<EntityFlag>::iterator it = notFoundMarkersNew.begin();
	     it != notFoundMarkersNew.end(); ++it )
                {
                    Debug(5010) << *it << " ";
                }
            Debug(5010) << "\n";
#endif
        }

    M_notFoundMarkers = notFoundMarkersCurrent;

    whereList.clear();
    // ============================================================================
    // There is no more identifiers to add to the boundary conditions
    // We finalise de set of identifiers by transfering it elements to a std::vector
    // ============================================================================
    for ( Iterator it = M_bcList.begin(); it != M_bcList.end(); ++it )
        {
            it->finalise();
	    it->finaliseIdGlobal();
        }

    M_bdUpdateDone = true;
} // bdUpdate

} // namespace LifeV

#endif /* BCHandler_H */
