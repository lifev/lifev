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
    @brief File contains BCManageNormal class for handling normal essential boundary conditions

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @contributor Mauro Perego <perego.mauro@gmail.com>
    @maintainer Mauro Perego <perego.mauro@gmail.com>

    @date 12-02-2009
 *///@HEADER

#ifndef BCMANAGENORMAL_H
#define BCMANAGENORMAL_H

#include <lifev/core/fem/BCHandler.hpp>

namespace LifeV
{

//! BCManageNormal - class for handling normal essential boundary conditions
/*!
   @author Gwenol Grandperrin

  The purpose of this class is to store the data and provide the methods needed to
  create the rotation matrix that help to impose normal essential boundary conditions.
  */


template<typename MatrixType>
class BCManageNormal
{
public:

    //! @name Public Types
    //@{
    typedef std::map<ID, ID> flagsMap_Type;
    typedef std::map<ID, Vector > versorsMap_Type;
    typedef MatrixType matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
    typedef boost::shared_ptr<MapEpetra> epetraMapPtr_Type;
    typedef boost::shared_ptr<VectorEpetra> epetraVectorPtr_type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    BCManageNormal();

    //! Copy constructor
    /*!
     * All the stored data are copied so that the pointers of bcManageNormal and this do not share the same objects
       @param bcManageNormal BCManageNormal
     */
    BCManageNormal ( BCManageNormal const& bcManageNormal );


    //! Destructor
    ~BCManageNormal();

    //@}


    //! @name Operators
    //@{

    //! Assignment operator
    /*!
      Stored data are copied so that the pointers of bcManageNormal and of the returned class do not share the same objects
      @param bcManageNormal BCManageNormal
      @return Reference to a new BCManageNormal with the same content of bcManageNormal
     */
    BCManageNormal& operator= ( const BCManageNormal& bcManageNormal );

    //@}


    //! @name Methods
    //@{

    //! Store boundary DOFs id and coordinates related to the boundary conditions
    /*!
        This method store the coordinates and the flags of the DOFs on the boundary.
        In the directional case, it will also store the normals.
        @param boundaryCondition A BCBase object
        @param time The actual time
     */
    void init (const BCBase& boundaryCondition, const Real& time);

    //! Build the rotation matrix
    /*!
     *  This method calculates the normal and tangential vectors and builds the rotation matrix
        @param dof The DOF object
        @param currentBdFE the current boundary finite element
        @param systemMatrix The system matrix
        @param offset The boundary condition offset
        @param commPtr pointer to Epetra_Comm object
     */
    void build (const DOF& dof, CurrentFEManifold& currentBdFE, matrix_Type& systemMatrix, UInt offset, MapEpetra::comm_ptrtype& commPtr);

    //! Build the rotation matrix
    /*!
     *  This method calculates the normal and tangential vectors and builds the rotation matrix
        @param mesh The mesh
        @param dof The DOF object
        @param currentBdFE the current boundary finite element
        @param systemMatrix The system matrix
        @param offset The boundary condition offset
        @param commPtr pointer to Epetra_Comm object
     */
    template<typename MeshType>
    void build (const MeshType& mesh, const DOF& dof, CurrentFEManifold& currentBdFE, const MapEpetra& map, UInt offset);

    //! This function modify the system matrix to apply a change of basis from the Cartesian coordinate system to the local coordinate system given by tangents and normals to the boundary
    /*!
        This method is used to change the basis of the DOF. Then it will be possible to prescribe
        an Essential boundary condition in the (t1, t2, n) direction by simply prescribing it along the x,y and z axis respectively.
        In fact it will perform the operations A := R*A*Rt and R*b, where A and be are the system matrix and right hand side,
        while R is the rotation matrix which transform the Cartesian system into the system of normals and tangents to the boundary.
        @param systemMatrix the matrix of the problem
        @param rightHandSide the right hand side
     */
    template <typename VectorType>
    void bcShiftToNormalTangentialCoordSystem (matrix_Type& systemMatrix, VectorType& rightHandSide) const;

    //! This function modify the system matrix to apply a change of basis from the local coordinate system given by tangents and normals to the boundary to the Cartesian coordinate system
    /*!
        This method is used to do the opposite operation of the method bcShiftToNormalTangentialCoordSystem.
        In fact it will perform the operations A := Rt*A*R and Rt*b, where A and be are the system matrix and right hand side,
        and R is the rotation matrix which transform the Cartesian system into the system of normals and tangents to the boundary.
        @param systemMatrix the matrix of the problem
        @param rightHandSide the right hand side
     */
    template <typename VectorType>
    void bcShiftToCartesianCoordSystem (matrix_Type& systemMatrix, VectorType& rightHandSide) const;

    //! This method computes the normals to the domain boundary
    /*!
        @warning the computed normals can be not compatible with the incompressibility condition
        @param dof The DOF object containing the local-to-global map
        @param currentBdFE The current finite element on the boundary
        @param normals The output vector
        @param mesh The mesh
     */
    template <typename VectorType, typename MeshType>
    void computeIntegratedNormals (const DOF& dof, CurrentFEManifold& currentBdFE, VectorType& normals,  const MeshType& mesh);

    //! Export in the vtk format the tangential and normal vectors (t1, t2, n) for each DOF on the domain boundary
    /*!
        @param fileName The name of the file where the data are exported
     */
    void exportToParaview (std::string fileName) const;


    //@}

private:

    //! @name Private Methods
    //@{
    //! stores the coordinates of the mesh vertices involved in the normal boundary conditions
    /*!
       @param mesh The mesh
     */
    template<typename MeshType>
    void M_calculateCoordinates (const MeshType& mesh);

    //! Add point to the map of flags
    /*!
         @param ID boundaryCondition A BCBase object
         @param time The actual time
      */
    void M_addBoundaryPoint (const ID& dofId, const ID& flag);

    //! Add point to the map of flags
    /*!
        @param ID dofId The point global id
        @param nx The x component of the unit vector
        @param ny The y component of the unit vector
        @param nz The z component of the unit vector
     */
    void M_addVersor (const ID& dofId, const Real& vx, const Real& vy, const Real& vz);

    //! Calculate the normal vectors
    /*!
        @param dof the dof class
        @param currentBdFE the current boundary finite element
     */
    template<typename MeshType>
    void M_calculateNormals (const MeshType& mesh, const DOF& dof, CurrentFEManifold& currentBdFE);

    //! Store in *M_normalPtr the versors given by the user
    /*!
        This function normalizes the versors in case they are not already normalized
     */
    void M_storeGivenVersors();

    //! Compute the tangential vectors
    void M_calculateTangentVectors();

    //! This function builds the rotation matrix for the change the basis.
    /*!
        @param systemMatrix the matrix of the problem
        @param offset that will be used if there is more than one unknown to recover the global ID
     */
    void M_buildRotationMatrix (const MapEpetra& map, UInt offset);
    //@}

    //! true when there are stored normals
    bool          M_dataBuilt;

    //! Shared pointer to the rotation matrix and its transpose
    matrixPtr_Type M_rotationMatrixPtr;
    matrixPtr_Type M_rotationMatrixTransposePtr;

    //! Shared pointer to the local Map
    epetraMapPtr_Type M_localMapEpetraPtr;

    //! Shared pointer to the vector of the first tangents to the domain boundary
    epetraVectorPtr_type M_firstTangentPtr;

    //! Shared pointer to the vector of the second tangents to the domain boundary
    epetraVectorPtr_type M_secondTangentPtr;

    //! Shared pointer to the vector of the normals to the domain boundary
    epetraVectorPtr_type M_normalPtr;

    //! Shared coordinates of the point
    epetraVectorPtr_type M_coordPtr;

    //! Number of degree of freedom (of one component)
    UInt          M_numDof;

    //! Number of degree of freedom (of one component) involved in boundary conditions
    UInt          M_numInvoledDof;

    //! flag of the boundary elements
    flagsMap_Type      M_flags;

    //! Versors that are given by the user
    versorsMap_Type    M_givenVersors;
};

//==============================================
//               IMPLEMENTATION
//==============================================

//==============================================
// Constructors and Destructor
//==============================================

//Empty Constructor
template<typename MatrixType>
BCManageNormal<MatrixType>::BCManageNormal() :
    M_dataBuilt (false),
    M_numDof (0),
    M_numInvoledDof (0)
{
    // Nothing to be done here
}

//Copy Constructor
template<typename MatrixType>
BCManageNormal<MatrixType>::BCManageNormal ( const BCManageNormal& bcManageNormal ) :
    M_dataBuilt (bcManageNormal.M_dataBuilt),
    M_rotationMatrixPtr (new matrix_Type (*bcManageNormal.M_rotationMatrixPtr) ),
    M_localMapEpetraPtr (new MapEpetra (*bcManageNormal.M_localMapEpetraPtr) ),
    M_firstTangentPtr (new VectorEpetra (*bcManageNormal.M_firstTangentPtr) ),
    M_secondTangentPtr (new VectorEpetra (*bcManageNormal.M_secondTangentPtr) ),
    M_normalPtr (new VectorEpetra (*bcManageNormal.M_normalPtr) ),
    M_coordPtr (new VectorEpetra (*bcManageNormal.M_coordPtr) ),
    M_numDof (bcManageNormal.M_numDof),
    M_numInvoledDof (bcManageNormal.M_numInvoledDof),
    M_flags (bcManageNormal.M_flags),
    M_givenVersors (bcManageNormal.M_givenVersors)
{
    // Nothing to be done here
}

// Destructor.
template<typename MatrixType>
BCManageNormal<MatrixType>::~BCManageNormal()
{
    // Nothing to be done here
}

//==============================================
// Operators
//==============================================

// Assignment operator
template<typename MatrixType>
BCManageNormal<MatrixType>&
BCManageNormal<MatrixType>::operator= ( const BCManageNormal& bcManageNormal )
{
    if (this != &bcManageNormal)
    {
        M_dataBuilt (bcManageNormal.M_dataBuilt);
        M_rotationMatrixPtr.reset (new matrix_Type (*bcManageNormal.M_rotationMatrixPtr) );
        M_localMapEpetraPtr.reset (new MapEpetra (*bcManageNormal.M_localMapEpetraPtr) );
        M_firstTangentPtr.reset (new VectorEpetra (*bcManageNormal.M_firstTangentPtr) );
        M_secondTangentPtr.reset (new VectorEpetra (*bcManageNormal.M_secondTangentPtr) );
        M_normalPtr.reset (new VectorEpetra (*bcManageNormal.M_normalPtr) );
        M_coordPtr.reset (new VectorEpetra (*bcManageNormal.M_coordPtr) );
        M_numDof = bcManageNormal.M_numDof;
        M_numInvoledDof = bcManageNormal.M_numInvoledDof;
        M_flags = bcManageNormal.M_flags;
        M_givenVersors = bcManageNormal.M_givenVersors;
    }
    return *this;
}

//==============================================
// Methods
//==============================================

template<typename MatrixType>
void BCManageNormal<MatrixType>::init (const BCBase& boundaryCondition, const Real& time)
{
    // Loop on BC identifiers
    for ( ID i = 0; i < boundaryCondition.list_size(); ++i )
    {
        const BCIdentifierEssential* pId = static_cast< const BCIdentifierEssential* > ( boundaryCondition[ i ] );

        if (boundaryCondition.mode() == Directional)
        {
            const BCFunctionDirectional* pBcF = static_cast<const BCFunctionDirectional*> ( boundaryCondition.pointerToFunctor() );
            Real nx (pBcF->vectFct (time, pId->x(), pId->y(), pId->z(), 0) );
            Real ny (pBcF->vectFct (time, pId->x(), pId->y(), pId->z(), 1) );
            Real nz (pBcF->vectFct (time, pId->x(), pId->y(), pId->z(), 2) );

            M_addVersor (boundaryCondition[i]->id(), nx, ny, nz);
        }
        else
        {
            M_addBoundaryPoint (boundaryCondition[i]->id(), boundaryCondition.flag() );
        }
    }
    M_dataBuilt = true; //Since vectors has been given we must apply the basis change.
}

template<typename MatrixType>
template<typename MeshType>
void BCManageNormal<MatrixType>::build (const MeshType& mesh, const DOF& dof, CurrentFEManifold& currentBdFE, const MapEpetra& map, UInt offset)
{
    if (M_dataBuilt)
    {
        //-----------------------------------------------------
        // STEP 1: Building the map
        //-----------------------------------------------------

        //First we build the map
        UInt nbPoints (M_flags.size() + M_givenVersors.size() );
        UInt i (0);
        std::vector<Int> idList (3 * nbPoints); //3 times because we want 3 coordinates

        //We store the number of Degrees of Freedom
        M_numDof = dof.numTotalDof();

        //Creating the iterators to explore the maps
        flagsMap_Type::iterator mapIt;
        versorsMap_Type::iterator mapIt2;

        //Building the list
        for ( mapIt = M_flags.begin() ; mapIt != M_flags.end(); mapIt++ )
        {
            idList[i] = (*mapIt).first;
            idList[i + nbPoints] = (*mapIt).first + M_numDof;
            idList[i + 2 * nbPoints] = (*mapIt).first + 2 * M_numDof;
            ++i;
        }
        for ( mapIt2 = M_givenVersors.begin() ; mapIt2 != M_givenVersors.end(); mapIt2++ )
        {
            idList[i] = (*mapIt2).first;
            idList[i + nbPoints] = (*mapIt2).first + M_numDof;
            idList[i + 2 * nbPoints] = (*mapIt2).first + 2 * M_numDof;
            ++i;
        }

        M_localMapEpetraPtr.reset ( new MapEpetra (-1, 3 * nbPoints, &idList[0], map.commPtr() ) );

        //-----------------------------------------------------
        // STEP 2: Compute normals and tangents
        //-----------------------------------------------------
        M_calculateNormals (mesh, dof, currentBdFE);
        M_storeGivenVersors();

        //this is used only for exporting the normals in vtk format.. should be put elsewhere?
        M_calculateCoordinates (mesh);
        M_calculateTangentVectors();
        M_buildRotationMatrix (map, offset);
    }
}

template<typename MatrixType>
template <typename VectorType>
void BCManageNormal<MatrixType>::bcShiftToNormalTangentialCoordSystem (matrix_Type& systemMatrix, VectorType& rightHandSide) const
{
    if (M_dataBuilt)
    {
        //Shift to tangential system

        //C = R*A
        matrix_Type C (systemMatrix.map(), systemMatrix.meanNumEntries() );
        M_rotationMatrixPtr->multiply (false, systemMatrix, false, C);

        //A = C*Rt"
        matrix_Type D (systemMatrix.map(), systemMatrix.meanNumEntries() );
        C.multiply (false, *M_rotationMatrixTransposePtr, false, D);
        systemMatrix.swapCrsMatrix (D);

        //b = R*b
        VectorType c (rightHandSide);
        M_rotationMatrixPtr->multiply (false, c, rightHandSide);
    }
}

template<typename MatrixType>
template <typename VectorType>
void BCManageNormal<MatrixType>::bcShiftToCartesianCoordSystem (matrix_Type& systemMatrix, VectorType& rightHandSide) const
{
    if (M_dataBuilt)
    {
        // C = Rt*A;
        matrix_Type C (systemMatrix.map(), systemMatrix.meanNumEntries() );
        M_rotationMatrixTransposePtr->multiply (false, systemMatrix, false, C);

        // A = C*R";
        matrix_Type D (systemMatrix.map(), systemMatrix.meanNumEntries() );
        C.multiply (false, *M_rotationMatrixPtr, false, D);
        systemMatrix.swapCrsMatrix (D);

        // b = Rt*b;
        VectorType c (rightHandSide);
        M_rotationMatrixTransposePtr->multiply (false, c, rightHandSide);
    }
}

template<typename MatrixType>
template<typename VectorType, typename MeshType>
void BCManageNormal<MatrixType>::computeIntegratedNormals (const DOF& dof, CurrentFEManifold& currentBdFE, VectorType& normals,  const MeshType& mesh)
{

    //-----------------------------------------------------
    // STEP 1: Calculating the normals
    //-----------------------------------------------------

    VectorType repNormals (normals.map(), Repeated);
    //Loop on the Faces
    for ( UInt iFace = 0; iFace < mesh.numBoundaryFacets(); ++iFace )
    {
        //Update the currentBdFE with the face data
        currentBdFE.update ( mesh.boundaryFacet ( iFace ), UPDATE_NORMALS | UPDATE_W_ROOT_DET_METRIC );
        UInt nDofF = currentBdFE.nbFEDof();

        //For each node on the face
        for (UInt icheck = 0; icheck < nDofF; ++icheck)
        {
            ID idf = dof.localToGlobalMapByBdFacet (iFace, icheck);

            //If the face exists and the point is on this processor
            if (M_flags.find (idf) != M_flags.end() )
            {
                ID flag = M_flags[idf];

                //if the normal is not already calculated
                //and the marker correspond to the flag of the point
                if ( (flag == mesh.boundaryFacet (iFace).markerID() ) || (flag == 0) )
                {
                    //Warning: the normal is taken in the first Gauss point
                    //since the normal is the same over the triangle
                    //(not true in the case of quadratic and bilinear maps)
                    Real nx (currentBdFE.normal (0, 0) );
                    Real ny (currentBdFE.normal (1, 0) );
                    Real nz (currentBdFE.normal (2, 0) );

                    //We get the area
                    Real area (currentBdFE.measure() );

                    //We update the normal component of the boundary point
                    (repNormals) [idf] += nx * area;
                    (repNormals) [idf + dof.numTotalDof()] += ny * area;
                    (repNormals) [idf + 2 * dof.numTotalDof()] += nz * area;
                }
            }
        }
    }

    //-----------------------------------------------------
    // STEP 2: Gathering the data from others processors
    //-----------------------------------------------------

    normals = VectorType (repNormals, Unique);

    //-----------------------------------------------------
    // STEP 3: Normalizing the vectors
    //-----------------------------------------------------

    //We obtain the ID of the element
    Int NumMyElements = normals.map().map (Unique)->NumMyElements();
    std::vector<Int> MyGlobalElements (NumMyElements);
    normals.map().map (Unique)->MyGlobalElements (& (MyGlobalElements[0]) );

    //We normalize the normal
    Real norm;
    UInt id;

    //Need to run only over the first third of MyGlobalElements
    //(the larger values are the y and z components)
    for ( Int i (0); i < NumMyElements / static_cast<Int> ( nDimensions ); ++i )
    {
        id = MyGlobalElements[i];
        Real nx ( (normals) [id] );
        Real ny ( (normals) [id + dof.numTotalDof()] );
        Real nz ( (normals) [id + 2 * dof.numTotalDof()] );
        norm = std::sqrt ( nx * nx + ny * ny + nz * nz );
        (normals) [id]            /= norm;
        (normals) [id + dof.numTotalDof()]   /= norm;
        (normals) [id + 2 * dof.numTotalDof()] /= norm;
    }
}

//TODO this function can be improved using the vtk writers
template<typename MatrixType>
void BCManageNormal<MatrixType>::exportToParaview (std::string fileName) const
{
    if (M_dataBuilt)
    {
        fileName.append ("_proc");
        std::ostringstream ossMyPid;
        ossMyPid << M_localMapEpetraPtr->comm().MyPID();
        fileName.append ( ossMyPid.str() );
        fileName.append (".vtk");
        std::ofstream file (fileName.c_str() );

        //Is the file open?
        if (file.fail() )
        {
            std::cerr << "Error: The file " << fileName << " is not opened " << std::endl;
        }
        else
        {
            //We obtain the ID of the element
            Int NumMyElements = M_localMapEpetraPtr->map (Unique)->NumMyElements();
            std::vector<Int> MyGlobalElements (NumMyElements);
            M_localMapEpetraPtr->map (Unique)->MyGlobalElements (&MyGlobalElements[0]);
            ID idof (0);

            //Writing the header
            file << "# vtk DataFile Version 2.0" << std::endl;
            file << "Normal directions" << std::endl;
            file << "ASCII" << std::endl;

            //Writing the points
            file << "DATASET POLYDATA" << std::endl;
            file << "POINTS " << M_numInvoledDof << " float" << std::endl;

            //Need to run only over the first third of MyGlobalElements
            //(the larger values are the y and z components)
            for ( Int i (0); i < NumMyElements / static_cast<Int> ( nDimensions ); ++i )
            {
                idof = MyGlobalElements[i];

                file << (*M_coordPtr) [idof] << "\t";
                file << (*M_coordPtr) [idof + M_numDof] << "\t";
                file << (*M_coordPtr) [idof + 2 * M_numDof] << std::endl;
            }

            //Starting the data part of the file
            file << "POINT_DATA " << M_numInvoledDof << std::endl;

            //Writing t1
            file << "VECTORS cell_tangent_1 float" << std::endl;

            //Need to run only over the first third of MyGlobalElements
            //(the larger values are the y and z components)
            for ( Int i (0); i < NumMyElements / static_cast<Int> ( nDimensions ); ++i )
            {
                idof = MyGlobalElements[i];

                file << (*M_firstTangentPtr) [idof] << "\t";
                file << (*M_firstTangentPtr) [idof + M_numDof] << "\t";
                file << (*M_firstTangentPtr) [idof + 2 * M_numDof] << std::endl;
            }

            //Writing t2
            file << "VECTORS cell_tangent_2 float" << std::endl;

            //Need to run only over the first third of MyGlobalElements
            //(the larger values are the y and z components)
            for ( Int i (0); i < NumMyElements / static_cast<Int> ( nDimensions ); ++i )
            {
                idof = MyGlobalElements[i];

                file << (*M_secondTangentPtr) [idof] << "\t";
                file << (*M_secondTangentPtr) [idof + M_numDof] << "\t";
                file << (*M_secondTangentPtr) [idof + 2 * M_numDof] << std::endl;
            }

            //Writing n
            file << "VECTORS cell_normals float" << std::endl;

            //Need to run only over the first third of MyGlobalElements
            //(the larger values are the y and z components)
            for ( Int i (0); i < NumMyElements / static_cast<Int> ( nDimensions ); ++i )
            {
                idof = MyGlobalElements[i];

                file << (*M_normalPtr) [idof] << "\t";
                file << (*M_normalPtr) [idof + M_numDof] << "\t";
                file << (*M_normalPtr) [idof + 2 * M_numDof] << std::endl;
            }

            //Closing the file
            file.close();
        }
    }
}

//==============================================
// Private Methods
//==============================================

template< typename MatrixType>
template< typename MeshType >
void BCManageNormal<MatrixType>::M_calculateCoordinates (MeshType const& mesh)
{
    M_coordPtr.reset ( new VectorEpetra (*M_localMapEpetraPtr, Unique) );

    //We obtain the ID of the element
    Int NumMyElements = M_localMapEpetraPtr->map (Unique)->NumMyElements();
    std::vector<Int> MyGlobalElements (NumMyElements);
    M_localMapEpetraPtr->map (Unique)->MyGlobalElements (&MyGlobalElements[0]);

    UInt id;

    //Need to run only over the first third of MyGlobalElements
    //(the larger values are the y and z components)
    for ( Int i (0); i < NumMyElements / static_cast<Int> ( nDimensions ); ++i )
    {
        id = MyGlobalElements[i];

        for (UInt j (0); j < mesh.pointList.size(); ++j)
        {
            if (id == mesh.pointList[j].id() )
            {
                (*M_coordPtr) [id]            = mesh.pointList[j].x();
                (*M_coordPtr) [id + M_numDof]   = mesh.pointList[j].y();
                (*M_coordPtr) [id + 2 * M_numDof] = mesh.pointList[j].z();
                j = mesh.pointList.size();
            }
        }
    }
}

template<typename MatrixType>
void BCManageNormal<MatrixType>::M_addBoundaryPoint (const ID& idof, const ID& flag)
{
    M_flags.insert (std::pair<ID, ID> (idof, flag) );
}


template<typename MatrixType>
void BCManageNormal<MatrixType>::M_addVersor (const ID& idof, const Real& vx, const Real& vy, const Real& vz)
{
    Vector n (3);
    n[0] = vx;
    n[1] = vy;
    n[2] = vz;
    M_givenVersors.insert (std::pair<ID, Vector> (idof, n) );
}

template< typename MatrixType>
template< typename MeshType>
void BCManageNormal< MatrixType>::M_calculateNormals (const MeshType& mesh, const DOF& dof, CurrentFEManifold& currentBdFE)
{
    //-----------------------------------------------------
    // STEP 1: Calculating the normals
    //-----------------------------------------------------

    M_normalPtr.reset ( new VectorEpetra (*M_localMapEpetraPtr, Repeated) );
    //(*M_normalPtr)*=0;
    computeIntegratedNormals (dof, currentBdFE, *M_normalPtr, mesh);

    //-----------------------------------------------------
    // STEP 4: Cleaning the memory
    //-----------------------------------------------------

    M_flags.clear();
}

//! Save the imposed normal vectors.
template<typename MatrixType>
void BCManageNormal<MatrixType>::M_storeGivenVersors()
{
    //-----------------------------------------------------
    // STEP 1: Retrieve the normals
    //-----------------------------------------------------

    //We obtain the ID of the element
    Int NumMyElements = M_localMapEpetraPtr->map (Unique)->NumMyElements();
    std::vector<Int> MyGlobalElements (NumMyElements);
    M_localMapEpetraPtr->map (Unique)->MyGlobalElements (&MyGlobalElements[0]);

    //We normalize the normal
    Real norm;
    UInt id;

    //Need to run only over the first third of MyGlobalElements
    //(the larger values are the y and z components)
    for ( Int i (0); i < NumMyElements / static_cast<Int> ( nDimensions ); ++i )
    {
        id = MyGlobalElements[i];

        if ( M_givenVersors.find (id) != M_givenVersors.end() )
        {
            Real nx ( (M_givenVersors) [id][0] );
            Real ny ( (M_givenVersors) [id][1] );
            Real nz ( (M_givenVersors) [id][2] );
            norm = std::sqrt ( nx * nx + ny * ny + nz * nz );
            (*M_normalPtr) [id]            = nx / norm;
            (*M_normalPtr) [id + M_numDof]   = ny / norm;
            (*M_normalPtr) [id + 2 * M_numDof] = nz / norm;
        }
    }

    //-----------------------------------------------------
    // STEP 2: Cleaning the memory
    //-----------------------------------------------------

    M_givenVersors.clear();
}

//! Calculate the tangential vectors for each triad in the vector triad.
template<typename MatrixType>
void BCManageNormal<MatrixType>::M_calculateTangentVectors()
{
    //-----------------------------------------------------
    // STEP 1: Initialization
    //-----------------------------------------------------

    //We obtain the ID of the element
    Int NumMyElements = M_localMapEpetraPtr->map (Unique)->NumMyElements();
    std::vector<Int> MyGlobalElements (NumMyElements);
    M_localMapEpetraPtr->map (Unique)->MyGlobalElements (&MyGlobalElements[0]);

    //Building the tangential vectors
    M_firstTangentPtr.reset ( new VectorEpetra (*M_localMapEpetraPtr, Unique) );
    (*M_firstTangentPtr) *= 0;
    M_secondTangentPtr.reset ( new VectorEpetra (*M_localMapEpetraPtr, Unique) );
    (*M_secondTangentPtr) *= 0;

    //We are going to use the loop to count the number
    //of imposed DOF because since the VectorEpetra is
    //unique now it is possible that we have less DOF
    M_numInvoledDof = 0;

    //-----------------------------------------------------
    // STEP 2: Calculation of the tangential vectors
    //-----------------------------------------------------

    // Real norm;
    UInt id;

    //Need to run only over the first third of MyGlobalElements
    //(the larger values are the y and z components)
    for ( Int i (0); i < NumMyElements / static_cast<Int> ( nDimensions ); ++i )
    {
        id = MyGlobalElements[i];

        //Counting the number of DOF
        M_numInvoledDof++;

        //We take max{|n x i|,|n x j|,|n x k|}
        //            =max{sqrt(ny^2+nz^2),sqrt(nx^2+nz^2),sqrt(nx^2+ny^2)}
        //            =max{r1,r2,r3}
        Real nx ( (*M_normalPtr) [id] );
        Real ny ( (*M_normalPtr) [id + M_numDof] );
        Real nz ( (*M_normalPtr) [id + 2 * M_numDof] );
        Real nxi = std::sqrt (ny * ny + nz * nz);
        Real nxj = std::sqrt (nx * nx + nz * nz);
        Real nxk = std::sqrt (nx * nx + ny * ny);

        if ( (nxi >= nxj) && (nxi >= nxk) ) //max = |n x i|
        {
            //We create t1
            (*M_firstTangentPtr) [id]            = 0;
            (*M_firstTangentPtr) [id + M_numDof]   = nz / nxi;
            (*M_firstTangentPtr) [id + 2 * M_numDof] = -ny / nxi;

            //We create t2
            (*M_secondTangentPtr) [id]            = -nxi;
            (*M_secondTangentPtr) [id + M_numDof]   = nx * ny / nxi;
            (*M_secondTangentPtr) [id + 2 * M_numDof] = nx * nz / nxi;
        }
        else if ( (nxj >= nxi) && (nxj >= nxk) ) //max = |n x j|
        {
            //We create t1
            (*M_firstTangentPtr) [id]            = -nz / nxj;
            (*M_firstTangentPtr) [id + M_numDof]   = 0;
            (*M_firstTangentPtr) [id + 2 * M_numDof] = nx / nxj;

            //We create t2
            (*M_secondTangentPtr) [id]            = nx * ny / nxj;
            (*M_secondTangentPtr) [id + M_numDof]   = -nxj;
            (*M_secondTangentPtr) [id + 2 * M_numDof] = ny * nz / nxj;
        }
        else //max = |n x k|
        {
            //We create t1
            (*M_firstTangentPtr) [id]            = ny / nxk;
            (*M_firstTangentPtr) [id + M_numDof]   = -nx / nxk;
            (*M_firstTangentPtr) [id + 2 * M_numDof] = 0;

            //We create t2
            (*M_secondTangentPtr) [id]            = nx * nz / nxk;
            (*M_secondTangentPtr) [id + M_numDof]   = ny * nz / nxk;
            (*M_secondTangentPtr) [id + 2 * M_numDof] = -nxk;
        }
    }
}

template<typename MatrixType>
void BCManageNormal<MatrixType>::M_buildRotationMatrix (const MapEpetra& map, UInt offset)
{
    //Initialization of the map to store the normal vectors
    std::map< ID, std::vector< Real > >::iterator mapIt;

    //Creating the matrix
    M_rotationMatrixPtr.reset ( new matrix_Type(map, 3) );

    //Adding one to the diagonal
    M_rotationMatrixPtr->insertOneDiagonal();

    static const Int nbRows (3);
    static const Int nbCols (3);
    std::vector<Real*> values (nbCols);
    Int Indices[3];
    for ( Int n = 0; n < nbCols; ++n )
    {
        values[n] = new Real[nbRows];
        for ( Int m = 0; m < nbRows; ++m )
        {
            values[n][m] = 0.0;
        }
    }

    std::vector<Int> rows;
    std::vector<Int> cols;

    //We obtain the ID of the element
    Int NumMyElements = M_localMapEpetraPtr->map (Unique)->NumMyElements();
    std::vector<Int> MyGlobalElements (NumMyElements);
    M_localMapEpetraPtr->map (Unique)->MyGlobalElements (&MyGlobalElements[0]);

    UInt id;

    //Need to run only over the first third of MyGlobalElements
    //(the larger values are the y and z components)
    for ( Int i (0); i < NumMyElements / static_cast<Int> ( nDimensions ); ++i )
    {
        id = MyGlobalElements[i];

        //...Except for the nodes where we make the rotation
        //Global Dof
        //idDof = boundaryCondition( i ) ->id() + boundaryCondition.component( j ) * totalDof;
        //i,j and k take values in [0,totalDof)
        Indices[0] = id + offset;
        Indices[1] = id + M_numDof + offset;
        Indices[2] = id + 2 * M_numDof + offset;

        cols.clear();
        cols.push_back (Indices[0]);
        cols.push_back (Indices[1]);
        cols.push_back (Indices[2]);

        rows.clear();
        rows.push_back (Indices[0]);
        rows.push_back (Indices[1]);
        rows.push_back (Indices[2]);

        //Line i (first tangential vector)
        values[0][0] = (*M_firstTangentPtr) [id];
        values[1][0] = (*M_firstTangentPtr) [id + M_numDof];
        values[2][0] = (*M_firstTangentPtr) [id + 2 * M_numDof];

        //-1 because we added one to the diagonal
        M_rotationMatrixPtr->addToCoefficient (Indices[0], Indices[0], -1.0);

        //Line j (second tangential vector)
        values[0][1] = (*M_secondTangentPtr) [id];
        values[1][1] = (*M_secondTangentPtr) [id + M_numDof];
        values[2][1] = (*M_secondTangentPtr) [id + 2 * M_numDof];

        //-1 because we added one to the diagonal
        M_rotationMatrixPtr->addToCoefficient (Indices[1], Indices[1], -1.0);

        //Line k (normal vector)
        values[0][2] = (*M_normalPtr) [id];
        values[1][2] = (*M_normalPtr) [id + M_numDof];
        values[2][2] = (*M_normalPtr) [id + 2 * M_numDof];

        //-1 because we added one to the diagonal
        M_rotationMatrixPtr->addToCoefficient (Indices[2], Indices[2], -1.0);

        M_rotationMatrixPtr->addToCoefficients (nbCols, nbRows, cols, rows, &values[0]);
    }

    for ( Int n = 0; n < nbRows; ++n )
    {
        delete[] values[n];
    }

    M_rotationMatrixPtr->globalAssemble();
    M_rotationMatrixTransposePtr = M_rotationMatrixPtr->transpose();
}

} //end of namespace LifeV

#endif // BCMANAGENORMAL_H
