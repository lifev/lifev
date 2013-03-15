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
    @brief File containing BCHandler class for handling boundary conditions

    @author Miguel Fernandez <miguel.fernandez@inria.fr>
    @contributor Christophe Prud'homme <christophe.prudhomme@epfl.ch>
    @contributor Mauro Perego <perego.mauro@gmail.com>
    @maintainer Mauro Perego <perego.mauro@gmail.com>

    @date 10-11-2004
 *///@HEADER

#ifndef BCHANDLER_H
#define BCHANDLER_H 1

#include <lifev/core/fem/BCBase.hpp>

namespace LifeV
{

//! BCHandler - class for handling boundary conditions
/*!
   @author Miguel Fernandez

   Container for @c BCBase classes

   @c BCHandler is a container of @c BCBase objects.
   It instantiates the needed @c BCBase objects and provide them with the data needed.
   As an example, if a @c BCBase object is related with an Essential boundary condition, then the global DOFs IDs associated with that essential boundary conditions will be passed to the object.

   Boundary conditions are added using the method @c addBC(@c name, @c flag, @c type, ...).
   @c type is a @c bcType_Type and describes the type of boundary condition (@c Essential, @c Natural, etc.).
   @c flag is used to determine on which elements the boundary condition should be prescribed.

   In general @c flag refers to the boundary elements' marker.
   The only exception are the @c bcType_Type @c EssentialEdges and @c EssentialVertices.
   In these two cases the flag refers respectively to the Edges' markers and to the Vertices' markers.
   In general @c EssentialEdges and @c EssentialVertices types are not recommended and they should be used only by advanced users.

   For this reasons boundary element markers are needed (and they better be specified in the mesh files).


   The typical way to use the class is to create function handlers:

   @verbatim
   BCFunctionBase bcFuncEss( solOnBoundary ), bcFuncNat( ZeroFunction);

   @endverbatim

   Then, instantiate a BCHandler object and add boundary conditions (this will create @c BCBase objects).

   @verbatim
   BCHandler bcHandler;
   bcHandler.addBC( "inlet", 10, Essential, Full, bcFuncEss, 3 );
   bcHandler.addBC( "outlet", 20, Natural, Full, bcFuncNat, 3 );
   bcHandler.addBC( "wall", 30, Essential, Full, bcFuncEss, 3 );
   @endverbatim

   Number 10, 20, 30 refers to the markers on the mesh boundary elements (facets: faces in 3D, edges in 2D).
   More information about boundary conditions will be found in BCBase.hpp

   Then, update the boundary conditions, using the function bcUpdate
   @verbatim
   bcHandler.bcUpdate(mesh, boundaryFe, dof);
   @endverbatim

   This method will looks for the markers in the mesh and will fill the boundary conditions container BCBase, with the appropriate data structures (global DOF id.. coordinates.. ecc)

   Now all is set to prescribe boundary conditions using one of the functions in bcManager.hpp
 */



class BCHandler
{
public:

    //! @name Public Types
    //@{

    typedef std::vector<BCBase>::iterator       bcBaseIterator_Type;
    typedef std::vector<BCBase>::const_iterator bcBaseConstIterator_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    BCHandler();

    //! Copy constructor
    /*!
       @warning This is not a copy constructor, since bcUpdate is set to false,
               as a consequence of the fact that also BCBase "copy" constructor is not
                actually a copy constructor
       @param bcHandler BCHandler
     */
    BCHandler ( const BCHandler& bcHandler );

    //! Destructor
    ~BCHandler();

    //@}


    //! @name Operators
    //@{

    //! Assignment operator
    /*!
      @param bcHandler BCHandler
      @return Reference to a new BCHandler with the same
                content of BCHandler bcHandler
     */
    BCHandler& operator= ( const BCHandler& bcHandler );

    //! Extract a BC in the list
    /*!
      @param i Position index
      @return Reference to the BCBase object in M_bcList at position i
     */
    BCBase& operator[] ( const ID& );

    //! Extract a BC in the list, const
    /*!
      @param i Position index
      @return Const reference to the BCBase object in M_bcList at position i
     */
    const BCBase& operator[] ( const ID& ) const;

    //@}


    //! @name Methods
    //@{

    //! Add new BC to the list for Component or Directional mode problems (user defined function case)
    /*!
      @param name The name of the boundary condition
      @param flag The mesh flag identifying the part of the mesh where the boundary condition applies
      @param type The boundary condition type: Natural, Robin, Flux, Resistance, Periodic, Essential, EssentialEdges, EssentialVertices
      @param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential, Directional
      @param bcFunction  The container holding the user defined function involved in this boundary condition
      @param components storing the list of components involved in this boundary condition
    */
    void addBC ( const std::string& name,
                 const bcFlag_Type& flag,
                 const bcType_Type& type,
                 const bcMode_Type& mode,
                 BCFunctionBase& bcFunction,
                 const bcComponentsVec_Type& components );


    //! Add new BC to the list for Scalar, Tangential or Normal mode problems (user defined function case)
    /*!
      @param name The name of the boundary condition
      @param flag The mesh flag identifying the part of the mesh where the boundary condition applies
      @param type The boundary condition type: Natural, Robin, Flux, Resistance, Periodic, Essential, EssentialEdges, EssentialVertices
      @param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential, Directional
      @param bcFunction  The container holding the user defined function involved in this boundary condition
    */
    void addBC ( const std::string& name,
                 const bcFlag_Type& flag,
                 const bcType_Type& type,
                 const bcMode_Type& mode,
                 BCFunctionBase& bcFunction );


    //! Add new BC to the list for Full mode problems  (user defined function case)
    /*!
     @param name The name of the boundary condition
     @param flag The mesh flag identifying the part of the mesh where the boundary condition applies
     @param type The boundary condition type: Natural, Robin, Flux, Resistance, Periodic, Essential, EssentialEdges, EssentialVertices
     @param mode The boundary condition mode: Scalar, Full, Component, Normal, Tangential, Directional
     @param bcFunction The container holding the user defined function involved in this boundary condition
     @param numberOfComponents The number of components involved in this boundary condition
    */
    void addBC ( const std::string& name,
                 const bcFlag_Type& flag,
                 const bcType_Type& type,
                 const bcMode_Type& mode,
                 BCFunctionBase& bcFunction,
                 const UInt& numberOfComponents );


    //! Add new BC to the list for Component or Directional mode problems (data vector case)
    /*!
      @param name The name of the boundary condition
      @param flag The mesh flag identifying the part of the mesh where the boundary condition applies
      @param type The boundary condition type: Natural, Robin, Flux, Resistance, Periodic, Essential, EssentialEdges, EssentialVertices
      @param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential, Directional
      @param bcVector The container holding the finite element vector involved in this boundary condition
      @param components storing the list of components involved in this boundary condition
    */
    void addBC ( const std::string& name,
                 const bcFlag_Type& flag,
                 const bcType_Type& type,
                 const bcMode_Type& mode,
                 BCVectorBase& bcVector,
                 const bcComponentsVec_Type& components );


    //! Add new BC to the list for Scalar, Tangential or Normal  mode problems (data vector case)
    /*!
      @param name The name of the boundary condition
      @param flag The mesh flag identifying the part of the mesh where the boundary condition applies
      @param type The boundary condition type: Natural, Robin, Flux, Resistance, Periodic, Essential, EssentialEdges, EssentialVertices
      @param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential, Directional
      @param bcVector The container holding the finite element vector involved in this boundary condition
     */
    void addBC ( const std::string& name,
                 const bcFlag_Type& flag,
                 const bcType_Type& type,
                 const bcMode_Type& mode,
                 BCVectorBase& bcVector );


    //! Add new BC to the list for Full mode problems (data vector case)
    /*!
      @param name The name of the boundary condition
      @param flag The mesh flag identifying the part of the mesh where the boundary condition applies
      @param type The boundary condition type: Natural, Robin, Flux, Resistance, Periodic, Essential, EssentialEdges, EssentialVertices
      @param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential, Directional
      @param bcVector The container holding the finite element vector involved in this boundary condition
      @param numberOfComponents The number of components involved in this boundary condition
    */
    void addBC ( const std::string& name,
                 const bcFlag_Type& flag,
                 const bcType_Type& type,
                 const bcMode_Type& mode,
                 BCVectorBase& bcVector,
                 const UInt& nComp );


    //! Add new BC to the list for Scalar, Tangential or Normal mode problems (user defined function case, with function depending on U)
    /*!
      @param name The name of the boundary condition
      @param flag The mesh flag identifying the part of the mesh bcBaseIterator the boundary condition applies
      @param type The boundary condition type: Natural, Robin, Flux, Resistance, Periodic, Essential, EssentialEdges, EssentialVertices
      @param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential, Directional
      @param  bcFunctionFEVectorDependent  The container holding the user defined function, depending on a FE vector, involved in this boundary condition
    */
    void addBC ( const std::string& name,
                 const bcFlag_Type& flag,
                 const bcType_Type& type,
                 const bcMode_Type& mode,
                 BCFunctionUDepBase&  bcFunctionFEVectorDependent );

    //! Add new BC to the list using a BCBase object
    /*!
      @param bcBase a BCBase object
    */
    void addBC ( const BCBase& bcBase );


    //! Modify the boundary condition @c name, assigning the function @c bcFunction
    /*!
      @param name The name of the boundary condition to be modified
      @param bcFunction The container holding the user defined function which will replace the existing one

     */
    void modifyBC ( std::string const& name, BCFunctionBase const& bcFunction );

    //! Modify the boundary condition @c assigning the FE vector in @c bcVector
    /*!
      @param name The name of the boundary condition to be modified
      @param bcVector The container holding the user FE vector which will replace the existing one
     */
    void modifyBC ( std::string const& name, BCVectorBase const& bcVector );

    //! Modify the boundary condition @c name, assigning the function in @c  bcFunctionFEVectorDependent
    /*!
      @param name The name of the boundary condition to be modified
      @param  bcFunctionFEVectorDependent The container holding the user defined function, depending on an FE vector, which will replace the existing one
     */
    void modifyBC ( std::string const& name, BCFunctionUDepBase const&  bcFunctionFEVectorDependent );

    //! Modify the boundary condition associated with flag @c aFlag, assigning the function in @c bcFunction
    /*!
      @param aFlag The flag associated with the boundary condition to be modified
      @param bcFunction The container holding the user defined function which will replace the existing one
     */
    void modifyBC ( bcFlag_Type const& aFlag, BCFunctionBase const& bcFunction );

    //! Modify the boundary condition associated with flag @c aFlag, assigning the FE vector in @c bcVector
    /*!
      @param aFlag The flag associated with the boundary condition to be modified
      @param bcVector The container holding the user FE vector which will replace the existing one
     */
    void modifyBC ( bcFlag_Type const& aFlag, BCVectorBase const& bcVector );

    //! Modify the boundary condition associated with flag @c aFlag, assigning the function in @c  bcFunctionFEVectorDependent
    /*!
      @param aFlag The flag associated with the boundary condition to be modified
      @param  bcFunctionFEVectorDependent The container holding the user defined function, depending on an FE vector, which will replace the existing one
     */
    void modifyBC ( bcFlag_Type const& aFlag, BCFunctionUDepBase const&  bcFunctionFEVectorDependent );

    //! Modify the boundary condition associated with flag @c aFlag, assigning the type in @c bcType
    /*!
      @param aFlag The flag associated with the boundary condition to be modified
      @param bcVector The type of b.c. which will replace the existing one
     */
    void modifyBC ( bcFlag_Type const& aFlag, bcType_Type const& bcType );

    //! Update all the boundary conditions
    /*!
      This method update the BC classes checking the markers
      on the boundary.
      Except for EssentialEdges and EssentialVertices type of boundary conditions,
      the flags associated with the boundary conditions will be searched among the markers of the boundary element (facets).

      Then, all the DOFs belonging to a matching boundary element will be associated with the appropriate boundary conditions.
      It is possible then the same DOF is shared by different boundary conditions.
      In the case of essential boundary conditions, the largest condition will be prescribed on the shared DOF (largest in the ordering given by the operator< in BCBase.hpp)
      In particular, if two Essential boundary conditions share the same DOF, it will be prescribed the condition with the largest flag.
      This behavior is due to the fact that the largest boundary condition is the last to be prescribed.

      Finally M_bcUpdateDone is set to true, and it is possible to prescribed boundary conditions using functions in BCManage.hpp.

      @param mesh The mesh
      @param boundaryFE Current finite element on the boundary
      @param dof Container of the local to global map of DOF id
    */
    template <typename Mesh>
    void bcUpdate ( Mesh& mesh, CurrentBoundaryFE& boundaryFE, const DOF& dof );

    //! Merges the boundary condition bcHandler (with its offset) with the stored one
    /*!
      @param bch BCHandler
     */
    void merge ( BCHandler& bcHandler );

    //! Display the content of the variables
    /*!
      @param verbose The verbosity (default: true)
      @param out The ostream output (default: std::cout)
     */
    void showMe ( bool verbose = false, std::ostream& out = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set @c offset in all boundary conditions
    /*!
      @param offset The boundary condition offset
     */
    void setOffset ( const UInt& offset );

    //! Set @c offset in boundary conditions @c name
    /*!
      @param name The name of the boundary condition
      @param offset The boundary condition offset
     */
    void setOffset ( const std::string& name, Int offset );

    //@}


    //! @name Get Methods
    //@{

    //! Extract a BC in the list according to its flag
    /*!
      @param aFlag The flag associated with the boundary condition
      @return constant Reference to the  boundary condition associated with flag aFlag
     */
    BCBase& findBCWithFlag (const bcFlag_Type& aFlag);


    //! Extract a BC in the list according to its flag (non const)
    /*!
      @param aFlag The flag associated with the boundary condition
      @return constant Reference to the  boundary condition associated with flag aFlag
     */
    const BCBase& findBCWithFlag (const bcFlag_Type& aFlag) const;

    //! Get a vector list of BC with specific type. The list contains the bcName_Type of the BC.
    /*!
      @param aType The BC type to be inserted in the list
    */
    std::vector<bcName_Type> findAllBCWithType ( const bcType_Type& aType ) const;


    //! Get the number of boundary conditions with type @c aType
    /*!
      @param aType The BC type to be counted
      @return Number of boundary conditions with type aType
    */
    UInt numberOfBCWithType ( const bcType_Type& aType ) const;


    // Get the BC index from the list using its name
    /*!
      @param name The name of the boundary condition
      @return The index of the boundary condition
     */
    ID findBCIndexWithName (const bcName_Type& name) const;

    // Get the BC from the list using its name
    /*!
      @param name The name of the boundary condition
      @return The index of the boundary condition
     */
    const BCBase* findBCWithName (const bcName_Type& name) const;


    // Get the offset of the boundary conditions
    /*!
      @return The offset of the boundary conditions
     */
    inline UInt offset() const
    {
        return M_offset;
    }


    //! Iterator of the beginning of the boundary elements list
    /*!
      @return The Iterator of the beginning of M_bcList
     */
    inline bcBaseIterator_Type begin()
    {
        return M_bcList.begin();
    }
    inline bcBaseConstIterator_Type begin() const
    {
        return M_bcList.begin();
    }


    //! Iterator of the end of the boundary elements list
    /*!
      @return The Iterator of the end of M_bcList
     */
    inline bcBaseIterator_Type end()
    {
        return M_bcList.end();
    }
    inline bcBaseConstIterator_Type end() const
    {
        return M_bcList.end();
    }


    //! Number of the stored boundary conditions
    /*!
     @return the number of boundary conditions stored
     */
    inline UInt size() const
    {
        return M_bcList.size();
    }


    //! Determine whether no boundary conditions are stored
    /*!
     @return true only if no boundary conditions are stored
     */
    inline bool empty() const
    {
        return M_bcList.empty();
    }


    //! Determine whether bcUpdate has been done before
    /*!
     @return true only if bcUpdate has been done before
     */
    inline bool bcUpdateDone() const
    {
        return M_bcUpdateDone;
    }


    //!Determine whether all the stored boundary conditions have EssentialXXX type
    /*!  It throws a @c logic_error exception if state is not consistent with hint.
      @return true only if the stored boundary conditions are of EssentialXXX type
    */
    bool hasOnlyEssential() const;

    //@}

private:

    //! Find the BC named @c name
    /*!
       It throws an @c invalid_argument exception if @c name is not
       found.

       @return A pointer to @c BCBase
    */
    BCBase* findBC ( const std::string& name );


    //! Find the BC named @c aFlag
    /*!
       It throws an @c invalid_argument exception if @c aFlag is not
       found.

       @return A pointer to @c BCBase
    */
    BCBase* findBC ( const bcFlag_Type& aFlag);

    //! Sum the M_offset to boundary conditions offsets
    void M_sumOffsets();

    //! true only if the bcUpdate has been done
    bool M_bcUpdateDone;

    //! vector list holding the stored BCs
    std::vector<BCBase> M_bcList;

    //! offset
    UInt M_offset;

    //! set of markers which are in the mesh but not in the list
    std::set<bcFlag_Type> M_notFoundMarkers;

};


// ===================================================
// Template methods implementations
// ===================================================
template <typename Mesh>
void
BCHandler::bcUpdate ( Mesh& mesh, CurrentBoundaryFE& boundaryFE, const DOF& dof )
{
    typedef typename Mesh::elementShape_Type geoShape_Type;
    typedef typename geoShape_Type::GeoBShape geoBShape_Type;

    // Some useful local variables
    UInt nDofPerVert = dof.localDofPattern().nbDofPerVertex(); // number of DOF per vertices
    UInt nBElemVertices = geoBShape_Type::S_numVertices; // Number of boundary element's vertices

    UInt nDofPerEdge = dof.localDofPattern().nbDofPerEdge(); // number of DOF per vertices
    UInt nBElemEdges = geoBShape_Type::S_numEdges; // Number of boundary element's vertices

    bcFlag_Type marker; //will store the marker of each geometric entity
    bcFlag_Type elementMarker; //will store the marker of the element

    // BCbase Iterator
    bcBaseIterator_Type bcBaseIterator;

    // Iterators which point to the beginning of EssentialEdges and EssentialVertices conditions
    bcBaseIterator_Type beginEssential, beginEssVertices, beginEssEdges;

    // coordinates of DOFs points
    Real x, y, z;

    //sets containing markers of boundary elements which have not been found in the boundary conditions' flags
    std::set<bcFlag_Type> notFoundMarkers;

    bool marker_found (false);

    //We set the iterators which points to the beginning of essential types bc.
    //(Notice that essential bc are positioned at the end of M_bcList, in the order
    // Essential, EssentialEdges, EssentialVertices)

    beginEssential = M_bcList.begin();
    while ( (beginEssential != M_bcList.end() ) && (beginEssential->type() < Essential) )
    {
        beginEssential++;
    }

    beginEssEdges = beginEssential;
    while ( (beginEssEdges != M_bcList.end() ) && (beginEssEdges->type() < EssentialEdges) )
    {
        beginEssEdges++;
    }

    beginEssVertices = beginEssEdges;
    while ( (beginEssVertices != M_bcList.end() ) && (beginEssVertices->type() < EssentialVertices) )
    {
        beginEssVertices++;
    }

    // ===================================================
    // Loop on boundary faces
    // ===================================================
    for ( ID iBoundaryElement = 0 ; iBoundaryElement < mesh.numBoundaryFacets(); ++iBoundaryElement )
    {
        // ===================================================================================
        // construction of localToGlobalMapOnBElem (this part should be moved in DOF.hpp)
        // ===================================================================================

        boundaryFE.updateMeas ( mesh.boundaryFacet ( iBoundaryElement ) ); // updating finite element information
        elementMarker = mesh.boundaryFacet ( iBoundaryElement ).markerID(); // We keep the element marker


        //vector containing the local to global map on each element
        std::vector<ID> localToGlobalMapOnBElem = dof.localToGlobalMapOnBdFacet (iBoundaryElement);

        // =============================================================
        // Insertion of boundary conditions defined on boundary Elements
        // =============================================================

        //looking for elementMarker in boundary conditions' flags, and updating the boundary conditions data accordingly

        marker_found = false;
        bcBaseIterator = M_bcList.begin();
        while ( ( bcBaseIterator = std::find ( bcBaseIterator, beginEssEdges, elementMarker ) ) != beginEssEdges )
        {
            marker_found = true;
            switch ( bcBaseIterator->type() )
            {
                case EssentialVertices:
                case EssentialEdges:
                    ERROR_MSG ("EssentialVertices and EssentialEdges types should not be prescribed on face labels");
                    break;

                case Essential:
                    for (ID lDof = 0; lDof < localToGlobalMapOnBElem.size(); lDof++)
                    {
                        ID gDof = localToGlobalMapOnBElem[ lDof]; // global DOF

                        //providing Essential boundary conditions with needed data (global DOF id and their coordinates)
                        if ( bcBaseIterator->isDataAVector() )
                        {
                            bcBaseIterator->addBCIdentifier ( new BCIdentifierBase ( gDof) );
                        }
                        else
                        {
                            // With user defined functions
                            boundaryFE.coorMap ( x, y, z, boundaryFE.refFE.xi ( lDof ), boundaryFE.refFE.eta ( lDof ) );
                            bcBaseIterator->addBCIdentifier ( new BCIdentifierEssential ( gDof, x, y, z ) );
                        }
                    }
                    break;

                case Natural:
                    //providing Natural boundary conditions with global DOFs on element
                    if ( bcBaseIterator->isDataAVector() )
                    {
                        // With data vector
                        UInt type = bcBaseIterator->pointerToBCVector()->type() ;
                        if ( type == 0 )  // if the BC is a vector which values don't need to be integrated
                        {
                            for (ID lDof = 0; lDof < localToGlobalMapOnBElem.size(); lDof++)
                            {
                                bcBaseIterator->addBCIdentifier ( new BCIdentifierNatural ( localToGlobalMapOnBElem[lDof] ) );
                            }
                        }
                        else if ( type <= 3 )  // FE vector which needed to be integrated on the boundary
                        {
                            bcBaseIterator->addBCIdentifier ( new BCIdentifierNatural ( iBoundaryElement, localToGlobalMapOnBElem ) );
                        }
                        else
                        {
                            ERROR_MSG ( "This BCVector type is not yet implemented" );
                        }
                    }
                    else
                    {
                        bcBaseIterator->addBCIdentifier ( new BCIdentifierNatural ( iBoundaryElement, localToGlobalMapOnBElem) );
                    }
                    break;

                case Robin:
                    //providing Robin boundary conditions with global DOFs on element
                    bcBaseIterator->addBCIdentifier ( new BCIdentifierNatural ( iBoundaryElement, localToGlobalMapOnBElem ) );
                    break;
                case Resistance:
                    //providing Resistance boundary conditions with global DOFs on element
                    if ( bcBaseIterator->isDataAVector() )
                    {
                        bcBaseIterator->addBCIdentifier ( new BCIdentifierNatural ( iBoundaryElement, localToGlobalMapOnBElem ) );
                    }
                    break;
                case Flux:
                    //providing Flux boundary conditions with global DOFs on element
                    bcBaseIterator->addBCIdentifier ( new BCIdentifierNatural ( iBoundaryElement, localToGlobalMapOnBElem ) );
                    break;
                default:
                    ERROR_MSG ("BC not yet implemented");
                    break;
            }
            bcBaseIterator++;
        }

        //inserting found and not found markers in respective containers.
        if (!marker_found)
        {
            notFoundMarkers.insert (elementMarker);
        }


        // ===================================================
        // Insertion of EssentialEdges bc
        // ===================================================

        //looking for which EssentialEdges are involved in this element.
        if ( (beginEssEdges != beginEssVertices) && (nDofPerVert || nDofPerEdge) )
        {
            ID iAdjacentElem = mesh.boundaryFacet ( iBoundaryElement ).firstAdjacentElementIdentity(); // id of the element adjacent to the face
            ID iElemBElement = mesh.boundaryFacet ( iBoundaryElement ).firstAdjacentElementPosition(); // local id of the face in its adjacent element

            //loop on boundary element edges
            for ( ID iBElemEdge = 0; iBElemEdge < nBElemEdges; ++iBElemEdge )
            {
                //index of edge in element
                ID iElemEdge = geoShape_Type::faceToEdge ( iElemBElement, iBElemEdge ).first;

                //marker on boundary edge
                marker = mesh.boundaryEdge ( mesh.localEdgeId ( iAdjacentElem, iElemEdge ) ).markerID(); // edge marker

                //indices of edge's vertices
                UInt iEdgeFirstVert =  geoBShape_Type::edgeToPoint ( iBElemEdge, 0);
                UInt iEdgeSecondVert = geoBShape_Type::edgeToPoint ( iBElemEdge, 1);

                // Finding this marker on the BC list
                bcBaseIterator = beginEssEdges;
                while ( ( bcBaseIterator = std::find ( bcBaseIterator, beginEssVertices, marker ) ) != beginEssVertices )
                {
                    //vector of boundary element DOFs on current edge
                    std::vector<ID> vecEdgeDofs (nDofPerEdge + 2 * nDofPerVert);

                    ID iTotalEdgeDof = 0; //local index of total DOF on current edge

                    //loop on Vertices DOFs
                    for ( ID iVertexDof = 0; iVertexDof < nDofPerVert; ++iVertexDof )
                    {
                        vecEdgeDofs[iTotalEdgeDof++] = iEdgeFirstVert * nDofPerVert + iVertexDof ; // local Dof
                        vecEdgeDofs[iTotalEdgeDof++] = iEdgeSecondVert * nDofPerVert + iVertexDof ;
                    }

                    //loop on Edge DOFs
                    for ( ID iEdgeDof = 0; iEdgeDof < nDofPerEdge; ++iEdgeDof )
                    {
                        vecEdgeDofs[iTotalEdgeDof++] = nDofPerVert * nBElemVertices + iBElemEdge * nDofPerEdge + iEdgeDof ;    // local Dof
                    }


                    //loop on total DOFs of current edge
                    for (iTotalEdgeDof = 0; iTotalEdgeDof < vecEdgeDofs.size(); iTotalEdgeDof++)
                    {
                        UInt lDof = vecEdgeDofs[iTotalEdgeDof]; //local DOF on boundary element
                        UInt gDof = localToGlobalMapOnBElem[ lDof ]; // global DOF

                        //providing the boundary conditions with needed data
                        if ( bcBaseIterator->isDataAVector() )
                        {
                            bcBaseIterator->addBCIdentifier ( new BCIdentifierBase ( gDof ) );
                        }
                        else
                        {
                            // With user defined functions
                            boundaryFE.coorMap ( x, y, z, boundaryFE.refFE.xi ( lDof ), boundaryFE.refFE.eta ( lDof ) );
                            bcBaseIterator->addBCIdentifier ( new BCIdentifierEssential ( gDof, x, y, z ) );
                        }
                    }
                    bcBaseIterator++;
                }
            }
        }


        // ===================================================
        // Insertion of EssentialVertices bc
        // ===================================================
        if ( (beginEssVertices != M_bcList.end() ) && nDofPerVert)
        {
            //loop on boundary element vertices
            for ( ID iBElemVert = 0; iBElemVert < nBElemVertices; ++iBElemVert )
            {
                marker = mesh.boundaryFacet ( iBoundaryElement ).point ( iBElemVert ).markerID(); // vertex marker

                // Finding this marker on the BC list
                bcBaseIterator = beginEssVertices;
                while ( ( bcBaseIterator = std::find ( bcBaseIterator, M_bcList.end(), marker ) ) != M_bcList.end() )
                {
                    for ( ID iVertexDof = 0; iVertexDof < nDofPerVert; ++iVertexDof )
                    {
                        UInt lDof = iBElemVert * nDofPerVert + iVertexDof ; // local Dof
                        UInt gDof = localToGlobalMapOnBElem[lDof]; // global Dof

                        //providing the boundary conditions with needed data
                        if ( bcBaseIterator->isDataAVector() )
                        {
                            bcBaseIterator->addBCIdentifier ( new BCIdentifierBase ( gDof ) );
                        }
                        else
                        {
                            // With user defined functions
                            boundaryFE.coorMap ( x, y, z, boundaryFE.refFE.xi ( lDof ), boundaryFE.refFE.eta ( lDof ) );
                            bcBaseIterator->addBCIdentifier ( new BCIdentifierEssential ( gDof, x, y, z ) );
                        }
                    }
                    bcBaseIterator++;
                }
            }
        }
    }


#ifdef DEBUG

    if ( notFoundMarkers.size() > 0 )
    {


        debugStream (5010) <<
                           "WARNING -- BCHandler::bcUpdate()\n" <<
                           "  boundary degrees of freedom with the following markers\n" <<
                           "  have no boundary condition set: ";
        for ( std::set<bcFlag_Type>::iterator it = notFoundMarkers.begin();
                it != notFoundMarkers.end(); ++it )
        {
            debugStream (5010) << *it << " ";
        }
        debugStream (5010) << "\n";

    }
#endif

    // ============================================================================
    // There is no more identifiers to add to the boundary conditions
    // We finalize the set of identifiers by transferring it elements to a std::vector
    // ============================================================================
    for ( bcBaseIterator = M_bcList.begin(); bcBaseIterator != M_bcList.end(); ++bcBaseIterator )
    {
        bcBaseIterator->copyIdSetIntoIdVector();
    }

    M_bcUpdateDone = true;
} // bcUpdate

} // namespace LifeV

#endif /* BCHANDLER_H */
