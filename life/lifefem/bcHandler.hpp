/* -*- mode: c++ -*-

This file is part of the LifeV library

Author(s): Miguel Fernandez
Christophe Prud'homme <christophe.prudhomme@epfl.ch>
Date: 2004-10-11

Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano

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
/**
   \file bcHandler.hpp
   \author Miguel Fernandez
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>

   \date 2004-10-11
*/
#ifndef __BCHandler_H
#define __BCHandler_H 1

#include <life/lifefem/bcCond.hpp>

namespace LifeV
{
/**
   \class BCHandler

   Container for \c BCBase classes

   \c BCHandler is a container for the boundary condition (BC) classes. It just
   uses an stl vector to store \c BCBase objets. The usage is simple: the
   user creates the data functors from the user defined functions:

   \verbatim
   BCFunctionBase gv( g );
   BCFunctionMixte gp( h, q );
   \endverbatim
   Then he/she specifies the number of BC and uses the \c BCHandler for
   creating the actual \c BCBase objects and storing them.

   \verbatim
   BCHandler bcV( 3 );
   bcV.addBC( "inlet", 10, Essential, Full, gv, 3 );
   bcV.addBC( "outflow", 11, Mixte, Scalar );
   bcV.addBC( "wall", 10, Essential, Full, gv, 3 );
   \endverbatim
*/
class BCHandler
{
public:

    /** @name Typedefs
     */
    //@{

    typedef std::vector<BCBase>::iterator Iterator;
    typedef std::vector<BCBase>::const_iterator ConstIterator;

    /*!
      \typedef enum
      Hints about the nature of the BCs
    */
    typedef enum BCHints
    {
        HINT_BC_NONE,          //!< no hint
        HINT_BC_ONLY_ESSENTIAL //!< BCs are only essential
    };

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /*!
      Constructor.
      \param nbc the number of BC to be stored
      \param hint hint about their nature
    */
    BCHandler( const ID& nbc = 0 , BCHints hint = HINT_BC_NONE );

    //@}


    BCHandler( const BCHandler & BCh );

    /** @name Operator overloads
     */
    //@{

    //! assignment operator
    BCHandler& operator = ( const BCHandler &BCh );
    //! extracting a BC in the list
    BCBase& operator[] ( const Index_t& );

    //! extracting a BC in the list
    const BCBase& operator[] ( const Index_t& ) const;

    //@}

    /** @name Accessors
     */
    //@{

    //! How many BC stored?
    Index_t size() const;

    //! Is there no BC stored?
    bool empty() const;

    //! Iterators for the begining and end of the BC list
    Iterator begin()
        {
            return M_bcList.begin();
        }
    Iterator end()
        {
            return M_bcList.end();
        }

    //! returns true if the bdUpdate has been done before
    bool bdUpdateDone() const;

    /*!
      \brief returns true if all the stored BC are of Essential type

      It throws a \c logic_error exception if state is not consistent with
      hint.
    */

    bool hasOnlyEssential() const;

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    //! add new BC to the list (user defined function)
    /*!
      \param name the name of the boundary condition
      \param flag the mesh flag identifying the part of the mesh where the boundary condtion applies
      \param type the boundary condition type: Natural, Essential, Mixte
      \param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
      \param bcf the function holding the user defined function involved in this boundary condition
      \param std::vector<ID> storing the list of components involved in this boundary condition
    */
    void addBC( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCFunctionBase& bcf,
                const std::vector<ID>& comp );

    /**
       Add new BC to the list without specified components for Scalar,
       Tangential or Normal mode problems (user defined function).

       \param name the name of the boundary condition
       \param flag the mesh flag identifying the part of the mesh where the boundary condtion applies
       \param type the boundary condition type: Natural, Essential, Mixte
       \param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
       \param bcf the function holding the user defined function involved in this boundary condition
    */
    void addBC( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCFunctionBase& bcf );


    //! add new BC to the list without list of components for Full mode problems  (user defined function)
    /*!
      \param name the name of the boundary condition
      \param flag the mesh flag identifying the part of the mesh where the boundary condiion applies
      \param type the boundary condition type: Natural, Essential, Mixte
      \param mode the boundary condition mode: Scalar, Full, Normal, Tangential
      \param bcf the function holding the user defined function involved in this boundary condition
      \param nComp the number of componets involved in this boundary condition
    */
    void addBC( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCFunctionBase& bcf,
                const UInt& nComp );


    //! add new BC to the list (data vector)
    /*!
      \param name the name of the boundary condition
      \param flag the mesh flag identifying the part of the mesh where the boundary condtion applies
      \param type the boundary condition type: Natural, Essential, Mixte
      \param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
      \param bcv data vector
      \param std::vector<ID> storing the list of components involved in this boundary condition
    */
    void addBC( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCVectorBase& bcv,
                const std::vector<ID>& comp );


    //! add new BC to the list  without specified components for Scalar, Tangential or Normal  mode problemst (data vector)
    /*!
      \param name the name of the boundary condition
      \param flag the mesh flag identifying the part of the mesh where the boundary condtion applies
      \param type the boundary condition type: Natural, Essential, Mixte
      \param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
      \param bcv data vector
    */
    void addBC( const std::string& name,
                const EntityFlag& flag,
                const BCType& type,
                const BCMode& mode,
                BCVectorBase& bcv );


    //! add new BC to the list without list of components for Full mode problems t (data vector)
    /*!
      \param name the name of the boundary condition
      \param flag the mesh flag identifying the part of the mesh where the boundary condiion applies
      \param type the boundary condition type: Natural, Essential, Mixte
      \param mode the boundary condition mode: Scalar, Full, Normal, Tangential
      \param bcv data vector
      \param nComp the number of componets involved in this boundary condition
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


    //! modify the boundary condition \c name with condition \c bcv
    void modifyBC( std::string const& name, BCFunctionBase& bcv );

    //! modify the boundary condition \c name with condition \c bcv
    void modifyBC( std::string const& name, BCVectorBase& bcv );

    void modifyBC( std::string const& __name, BCFunctionUDepBase& __bcf );

    //! modify the boundary condition \c lab with condition \c bcv
    void modifyBC( int lab, BCFunctionBase& bcv );

    //! modify the boundary condition \c lab with condition \c bcv
    void modifyBC( int lab, BCVectorBase& bcv );

    void modifyBC( int lab, BCFunctionUDepBase& __bcf );


    //! Build the boundary stuff
    /*!
      \param mesh the mesh
      \param feBd the current finite element on the boundary
      \param dof give the local to global table
      This method udapte the BC classes checking the markers
      on the boundary. It builds the list of identifiers depending
      on the BC type
    */
    template <typename Mesh>
    void bdUpdate( Mesh& mesh, CurrentBdFE& feBd, const Dof& dof );

    /*! Extract the type of the boundary condition with flag aFlag. Useful for weak
      imposition of boundary conditions (e.g. in DG-FEM).
    */
    BCType boundaryType(const EntityFlag& aFlag) const;

    /*! Extract a BC in the list according to its flag
     */
    BCBase& GetBCWithFlag(const EntityFlag& );
    const BCBase& GetBCWithFlag(const EntityFlag&) const;

    // get the BC number from the list using its name

    UInt BCHandler::getBCbyName(const std::string _BCName) const;

    //! output
    std::ostream& showMe( bool verbose = false, std::ostream & out = std::cout ) const;

    //@}

private:
    //! number of BC to be stored

    //! true if the bdUpdate has been done
    bool M_bdUpdateDone;

    //! hints about the nature of the BCs
    BCHints M_hint;

    //! vector list holding the stored BC
    std::vector<BCBase> M_bcList;

    /**
       look for the BC named \c name.

       It throws an \c invalid_argument exception if \c name is not
       found.

       \return a pointer to \c BCBase
    */
    BCBase* findBC( std::string const& __name );

    BCBase* findBC( int lab );

    //! determine whether BCs in M_bcList are only essential
    bool listHasOnlyEssential() const;

    //! stores status (essential or not) of the boundary associated to a marker
    class EssentialStatus
    {
    public:
        EssentialStatus() : M_normal( false ), M_tangential( false )
            {
                for(UInt iComp=0; iComp<nDimensions; ++iComp)
                {
                    M_components[iComp] = false;
                }
            }
        void setComponent(UInt comp) { M_components[comp-1] = true; }
        void setAllComponents()
            {
                for(UInt iComp=0; iComp<nDimensions; ++iComp)
                {
                    M_components[iComp] = true;
                }
            }
        void setNormal() { M_normal = true; }
        void setTangential() { M_tangential = true; }
        bool isEssential() const
            {
                bool result = true;
                for (UInt iComp=0; iComp<nDimensions; ++iComp)
                {
                    result &= M_components[iComp];
                }
                result = result || ( M_normal && M_tangential );
                return result;
            }
    private:
        bool M_components[nDimensions];
        bool M_normal;
        bool M_tangential;
    };

    //! set of markers which are in the mesh but not in the list
    std::set<EntityFlag> M_notFoundMarkers;

};

/********************************************************************
                      IMPLEMENTATIONS
********************************************************************/

/**
   The update method must be a method of \c BCHandler.
   
   It updates \c BCHandler objects not \c Dof objects
   Build the boundary stuff
*/
template <typename Mesh>
void
BCHandler::bdUpdate( Mesh& mesh, CurrentBdFE& feBd, const Dof& dof )
{
    
    typedef typename Mesh::VolumeShape GeoShape;
    
    // Some useful local variables, to save some typing
    UInt nDofpV = dof.fe.nbDofPerVertex; // number of Dof per vertices
    UInt nDofpE = dof.fe.nbDofPerEdge;   // number of Dof per edges
    UInt nDofpF = dof.fe.nbDofPerFace;   // number of Dof per faces
    
    UInt bdnF = mesh.numBFaces();    // number of faces on boundary
    
    EntityFlag marker; //will store the marker of each geometric entity
    
    typedef typename GeoShape::GeoBShape GeoBShape;
    
    UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
    UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges
    
    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges

    UInt nDofFV = nDofpV * nFaceV; // number of vertex's Dof on a face
    UInt nDofFE = nDofpE * nFaceE; // number of edge's Dof on a face

    UInt nDofF = nDofFV + nDofFE + nDofpF; // number of total Dof on a face

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element

    SimpleVect<ID> bdltg( nDofF );
    typedef std::vector<BCBase>::iterator Iterator;
    Iterator where;
    std::vector<Iterator> whereList;

    UInt iElAd, iVeEl, iFaEl, iEdEl;
    ID lDof, gDof;
    Real x, y, z;

    std::set<EntityFlag> notFoundMarkersCurrent;

    // ===================================================
    // Loop on boundary faces
    // ===================================================
    for ( ID ibF = 1 ; ibF <= bdnF; ++ibF )
    {
    
        iElAd = mesh.boundaryFace( ibF ).ad_first();  // id of the element adjacent to the face
        iFaEl = mesh.boundaryFace( ibF ).pos_first(); // local id of the face in its adjacent element
        feBd.updateMeas( mesh.boundaryFace( ibF ) );  // updating finite element information

        // ===================================================
        // Vertex based Dof
        // ===================================================
        if ( nDofpV )
        {
        
            // loop on face vertices
            for ( ID iVeFa = 1; iVeFa <= nFaceV; ++iVeFa )
            {

                marker = mesh.boundaryFace( ibF ).point( iVeFa ).marker(); // vertex marker
                iVeEl = GeoShape::fToP( iFaEl, iVeFa ); // local vertex number (in element)
        
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
                for ( ID l = 1; l <= nDofpV; ++l )
                {
            
                    lDof = ( iVeFa - 1 ) * nDofpV + l ; // local Dof
                    gDof = dof.localToGlobal( iElAd, ( iVeEl - 1 ) * nDofpV + l ); // global Dof
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
                                // Why kind of data ?
                                // vincent please check again for your Mixte-FE it doesn't work for Q1
                                //       if ( where->dataVector()  ) { // With data vector
                                //        where->addIdentifier( new IdentifierNatural(gDof) );
                                //       }
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
        if ( nDofpE )
        {
        
            // loop on face edges
            for ( ID iEdFa = 1; iEdFa <= nFaceE; ++iEdFa )
            {
        
                iEdEl = GeoShape::fToE( iFaEl, iEdFa ).first; // local edge number (in element)
                marker = mesh.boundaryEdge( mesh.localEdgeId( iElAd, iEdEl ) ).marker(); // edge marker
        
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
                for ( ID l = 1; l <= nDofpE; ++l )
                {
            
                    lDof = nDofFV + ( iEdFa - 1 ) * nDofpE + l ; // local Dof
                    gDof = dof.localToGlobal( iElAd, nDofElemV + ( iEdEl - 1 ) * nDofpE + l ); // global Dof
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
                            default:
                                ERROR_MSG( "This boundary condition type is not yet implemented" );
                        }
                    }
                }
            }
        }
    
    
    
        // ===================================================
        // Face based Dof
        // ===================================================
        marker = mesh.boundaryFace( ibF ).marker(); // edge marker
    
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
                    for ( ID l = 1; l <= nDofpF; ++l )
                    {
                        lDof = nDofFE + nDofFV + l; // local Dof
                        gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + ( iFaEl - 1 ) * nDofpF + l ); // global Dof
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
                            for ( ID l = 1; l <= nDofpF; ++l )
                            {
                                lDof = nDofFE + nDofFV + l; // local Dof
                                gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + ( iFaEl - 1 ) * nDofpF + l ); // global Dof
                                where->addIdentifier( new IdentifierNatural( gDof ) );
                            }
                        }
                        else if ( (type == 1) || (type == 2) ) 
                        {
                            // Loop on number of Dof per face
                            for ( ID l = 1; l <= nDofpF; ++l )
                            {
                                lDof = nDofFE + nDofFV + l; // local Dof
                                gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + ( iFaEl - 1 ) * nDofpF + l ); // global Dof
                                bdltg( lDof ) = gDof; // local to global on this face
                            }
                            where->addIdentifier( new IdentifierNatural( ibF, bdltg ) );
                        }
            
                        else
                            ERROR_MSG( "This BCVector type is not yet implemented" );
            
                    }
                    else
                    {
                        // Loop on number of Dof per face
                        for ( ID l = 1; l <= nDofpF; ++l )
                        {
                            lDof = nDofFE + nDofFV + l; // local Dof
                            gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + ( iFaEl - 1 ) * nDofpF + l ); // global Dof
                            bdltg( lDof ) = gDof; // local to global on this face
                        }
                        where->addIdentifier( new IdentifierNatural( ibF, bdltg ) );
                    }
                    break;
                case Mixte:
                    // Why kind of data ?
                    // vincent please check again for your Mixte-FE it doesn't work for Q1
                    // if ( where->dataVector()  ) { // With data vector
                    //   for (ID l=1; l<=nDofpF; ++l) {
                    //     lDof = nDofFE + nDofFV + l; // local Dof
                    //     gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + (iFaEl-1)*nDofpF + l); // global Dof
                    //     where->addIdentifier( new IdentifierNatural(gDof) );
                    //   }
                    // }
                    // else {
                    // Loop on number of Dof per face
                    for ( ID l = 1; l <= nDofpF; ++l )
                    {
                        lDof = nDofFE + nDofFV + l; // local Dof
                        gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + ( iFaEl - 1 ) * nDofpF + l ); // global Dof
                        bdltg( lDof ) = gDof; // local to global on this face
                    }
                    where->addIdentifier( new IdentifierNatural( ibF, bdltg ) );
                    // }
                    break;
                default:
                    ERROR_MSG( "This boundary condition type is not yet implemented" );
            }
        }
    }
    
    std::set<EntityFlag> notFoundMarkersNew;
    for( std::set<EntityFlag>::iterator it = notFoundMarkersCurrent.begin();
         it != notFoundMarkersCurrent.end(); ++it )
    {
        if ( M_notFoundMarkers.find( *it ) == M_notFoundMarkers.end() ) {
            notFoundMarkersNew.insert( *it );
        }
    }

    
    if( notFoundMarkersNew.size() > 0 )
    {
        std::cerr <<
            "WARNING -- BCHandler::bdUpdate()\n" <<
            "  boundary degrees of freedom with the following markers\n" <<
            "  have no boundary condition set: ";
        for( std::set<EntityFlag>::iterator it = notFoundMarkersNew.begin();
             it != notFoundMarkersNew.end(); ++it )
        {
            std::cerr << *it << " ";
        }
        std::cerr << std::endl;
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
    }
    
    M_bdUpdateDone = true;
} // bdUpdate
  
} // namespace LifeV

#endif /* __BCHandler_H */

