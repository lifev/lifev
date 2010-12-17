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
    @brief File contains BCNormalManager class for handling normal essential boundary conditions

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @contributor Mauro Perego <perego.mauro@gmail.com>
    @maintainer Mauro Perego <perego.mauro@gmail.com>

    @date 12-02-2009
 *///@HEADER

#ifndef _BCNORMALMANAGER_
#define _BCNORMALMANAGER_

#include <life/lifefem/bcHandler.hpp>

namespace LifeV
{

//! BCNormalManager - class for handling normal essential boundary conditions
/*!
   @author Gwenol Grandperrin

  The purpose of this class is to store the data and provide the methods needed to
  create the rotation matrix that help to impose normal essential boundary conditions.
  */


//TODO remove M_meshPtr and make the class template only on MatrixType
template<typename MeshType,typename MatrixType>
class BCNormalManager
{
public:

    //! @name Public Types
    //@{
    typedef std::map<ID,ID> FlagsMap;  //deprecated
    typedef std::map<ID, Vector > NormalsMap;  //deprecated

    typedef std::map<ID,ID> flagsMap_Type;
    typedef std::map<ID, Vector > versorsMap_Type;
    typedef MatrixType matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
    typedef boost::shared_ptr<EpetraMap> epetraMapPtr_Type;
    typedef boost::shared_ptr<EpetraVector> epetraVectorPtr_type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    BCNormalManager();

    //!Constructor
    /*!
     * @warning This constructor is deprecated and will be removed soon
     * @param mesh The mesh where the normal should be calculated
     */
    BCNormalManager(const MeshType& mesh);


    //! Copy constructor
    /*!
     * All the stored data are copied so that the pointers of bcNormalManager and this do not share the same objects
       @param bcNormalManager BCNormalManager
     */
    BCNormalManager( BCNormalManager const& bcNormalManager );


    //! Destructor
    ~BCNormalManager();

    //@}


    //! @name Operators
    //@{

    //! Assignment operator
    /*!
      Stored data are copied so that the pointers of bcNormalManager and of the returned class do not share the same objects
      @param bcNormalManager BCNormalManager
      @return Reference to a new BCNormalManager with the same content of bcNormalManager
     */
    BCNormalManager& operator= ( const BCNormalManager& bcNormalManager );

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
    void init(const BCBase& boundaryCondition, const Real& time);

    //! Build the rotation matrix
    /*!
     *  This method calculates the normal and tangential vectors and builds the rotation matrix
        @param dof The Dof object
        @param currentBdFE the current boundary finite element
        @param systemMatrix The system matrix
        @param offset The boundary condition offset
        @param commPtr pointer to Epetra_Comm object
     */
    void build(const Dof& dof, CurrentBdFE& currentBdFE, matrix_Type& systemMatrix, UInt offset, EpetraMap::comm_ptrtype& commPtr);


    //! Build the rotation matrix
    /*!
     *  This method calculates the normal and tangential vectors and builds the rotation matrix
        @param mesh The mesh
    	@param dof The Dof object
    	@param currentBdFE the current boundary finite element
    	@param systemMatrix The system matrix
    	@param offset The boundary condition offset
    	@param commPtr pointer to Epetra_Comm object
     */
    void build(const MeshType& mesh, const Dof& dof, CurrentBdFE& currentBdFE, matrix_Type& systemMatrix, UInt offset, EpetraMap::comm_ptrtype& commPtr);


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
    void bcShiftToNormalTangentialCoordSystem(matrix_Type& systemMatrix, VectorType& rightHandSide) const;


    //! This function modify the system matrix to apply a change of basis from the local coordinate system given by tangents and normals to the boundary to the Cartesian coordinate system
    /*!
        This method is used to do the opposite operation of the method bcShiftToNormalTangentialCoordSystem.
        In fact it will perform the operations A := Rt*A*R and Rt*b, where A and be are the system matrix and right hand side,
        and R is the rotation matrix which transform the Cartesian system into the system of normals and tangents to the boundary.
        @param systemMatrix the matrix of the problem
        @param rightHandSide the right hand side
     */
    template <typename VectorType>
    void bcShiftToCartesianCoordSystem(matrix_Type& systemMatrix, VectorType& rightHandSide) const;


    //! This method computes the normals to the domain boundary
    /*!
    	@warning the computed normals can be not compatible with the incompressibility condition
        @param dof The Dof object containing the local-to-global map
        @param currentBdFE The current finite element on the boundary
        @param normals The output vector
        @param mesh The mesh
     */
    template <typename VectorType>
    void computeIntegratedNormals(const Dof& dof, CurrentBdFE& currentBdFE, VectorType& normals,  const MeshType& mesh);



    //! Export in the vtk format the tangential and normal vectors (t1, t2, n) for each DOF on the domain boundary
    /*!
        @param fileName The name of the file where the data are exported
     */
    void exportToParaview(std::string fileName) const;


    //@}

private:

    //! @name Private Methods
    //@{
    //! stores the coordinates of the mesh vertices involved in the normal boundary conditions
    /*!
       @param mesh The mesh
     */
    void M_calculateCoordinates(const MeshType& mesh);


    //! Add point to the map of flags
    /*!
    	 @param ID boundaryCondition A BCBase object
    	 @param time The actual time
      */
    void M_addBoundaryPoint(const ID& dofId,const ID& flag);


    //! Add point to the map of flags
    /*!
        @param ID dofId The point global id
        @param nx The x component of the unit vector
        @param ny The y component of the unit vector
        @param nz The z component of the unit vector
     */
    void M_addVersor(const ID& dofId,const Real& vx,const Real& vy, const Real& vz);


    //! Calculate the normal vectors
    /*!
        @param dof the dof class
        @param currentBdFE the current boundary finite element
     */
    void M_calculateNormals(const MeshType& mesh, const Dof& dof,CurrentBdFE& currentBdFE);


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
    void M_buildRotationMatrix(matrix_Type& systemMatrix, UInt offset=0);
    //@}


    //! true when there are stored normals
    bool          M_dataBuilt;

    //! The mesh, this member will be removed soon
    const MeshType* const M_meshPtr;

    //! Shared pointer to the rotation matrix
    matrixPtr_Type M_rotationMatrixPtr;

    //! Shared pointer to the local Map
    epetraMapPtr_Type M_localEpetraMapPtr;

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
template<typename MeshType,typename MatrixType>
BCNormalManager<MeshType,MatrixType>::BCNormalManager():
        M_dataBuilt(false),
        M_meshPtr(0),
        M_numDof(0),
        M_numInvoledDof(0)
{

}

//Constructor
template<typename MeshType,typename MatrixType>
BCNormalManager<MeshType,MatrixType>::BCNormalManager(const MeshType& mesh):
        M_dataBuilt(false),
        M_meshPtr(&mesh),
        M_numDof(0),
        M_numInvoledDof(0)
{

}

//Copy Constructor
template<typename MeshType,typename MatrixType>
BCNormalManager<MeshType,MatrixType>::BCNormalManager( const BCNormalManager & bcNormalManager ):
        M_dataBuilt(bcNormalManager.M_dataBuilt),
        M_rotationMatrixPtr(new matrix_Type(*bcNormalManager.M_rotationMatrixPtr) ),
        M_localEpetraMapPtr(new EpetraMap(*bcNormalManager.M_localEpetraMapPtr) ),
        M_firstTangentPtr(new EpetraVector(*bcNormalManager.M_firstTangentPtr) ),
        M_secondTangentPtr(new EpetraVector(*bcNormalManager.M_secondTangentPtr) ),
        M_normalPtr(new EpetraVector(*bcNormalManager.M_normalPtr) ),
        M_coordPtr(new EpetraVector(*bcNormalManager.M_coordPtr) ),
        M_numDof(bcNormalManager.M_numDof),
        M_numInvoledDof(bcNormalManager.M_numInvoledDof),
        M_flags(bcNormalManager.M_flags),
        M_givenVersors(bcNormalManager.M_givenVersors)
{

}

// Destructor.
template<typename MeshType,typename MatrixType>
BCNormalManager<MeshType,MatrixType>::~BCNormalManager()
{

}

//==============================================
// Operators
//==============================================



// Assignment operator
template<typename MeshType,typename MatrixType>
BCNormalManager<MeshType,MatrixType>&
BCNormalManager<MeshType,MatrixType>::operator= ( const BCNormalManager & bcNormalManager )
{
    if (this != &bcNormalManager)
    {
        M_dataBuilt(bcNormalManager.M_dataBuilt);
        M_rotationMatrixPtr.reset(new matrix_Type(*bcNormalManager.M_rotationMatrixPtr) );
        M_localEpetraMapPtr.reset(new EpetraMap(*bcNormalManager.M_localEpetraMapPtr) );
        M_firstTangentPtr.reset(new EpetraVector(*bcNormalManager.M_firstTangentPtr) );
        M_secondTangentPtr.reset(new EpetraVector(*bcNormalManager.M_secondTangentPtr) );
        M_normalPtr.reset(new EpetraVector(*bcNormalManager.M_normalPtr) );
        M_coordPtr.reset(new EpetraVector(*bcNormalManager.M_coordPtr) );
        M_numDof = bcNormalManager.M_numDof;
        M_numInvoledDof = bcNormalManager.M_numInvoledDof;
        M_flags = bcNormalManager.M_flags;
        M_givenVersors = bcNormalManager.M_givenVersors;
    }
    return *this;
}


//==============================================
// Methods
//==============================================



template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::init(const BCBase& boundaryCondition,const Real& time)
{
    // Loop on BC identifiers
    for ( ID i = 1; i <= boundaryCondition.list_size(); ++i )
    {
        const IdentifierEssential* pId = static_cast< const IdentifierEssential* >( boundaryCondition( i ) );

        if (boundaryCondition.mode()==Directional)
        {
            const BCFunctionDirectional* pBcF = static_cast<const BCFunctionDirectional*>( boundaryCondition.pointerToFunctor() );
            Real nx(pBcF->vectFct(time,pId->x(),pId->y(),pId->z(),1));
            Real ny(pBcF->vectFct(time,pId->x(),pId->y(),pId->z(),2));
            Real nz(pBcF->vectFct(time,pId->x(),pId->y(),pId->z(),3));

            M_addVersor(boundaryCondition(i)->id(),nx,ny,nz);
        }
        else
        {
            M_addBoundaryPoint(boundaryCondition(i)->id(),boundaryCondition.flag());
        }
    }
    M_dataBuilt = true; //Since vectors has been given we must apply the basis change.
}

template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::build(const Dof& dof,CurrentBdFE& currentBdFE,MatrixType& systemMatrix, UInt offset,EpetraMap::comm_ptrtype& commPtr)
{
    build(M_meshPtr, dof, currentBdFE, systemMatrix, offset, commPtr);
}


template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::build(const MeshType& mesh, const Dof& dof,CurrentBdFE& currentBdFE,MatrixType& systemMatrix, UInt offset,EpetraMap::comm_ptrtype& commPtr)
{
    if (M_dataBuilt)
    {
        //-----------------------------------------------------
        // STEP 1: Building the map
        //-----------------------------------------------------

        //First we build the map
        UInt nbPoints(M_flags.size()+M_givenVersors.size());
        UInt i(0);
        int idList[3*nbPoints]; //3 times because we want 3 coordinates

        //We store the number of Degrees of Freedom
        M_numDof = dof.numTotalDof();

        //Creating the iterators to explore the maps
        flagsMap_Type::iterator mapIt;
        versorsMap_Type::iterator mapIt2;

        //Building the list
        for ( mapIt=M_flags.begin() ; mapIt != M_flags.end(); mapIt++ )
        {
            idList[i] = (*mapIt).first;
            idList[i+nbPoints] = (*mapIt).first+M_numDof;
            idList[i+2*nbPoints] = (*mapIt).first+2*M_numDof;
            ++i;
        }
        for ( mapIt2=M_givenVersors.begin() ; mapIt2 != M_givenVersors.end(); mapIt2++ )
        {
            idList[i] = (*mapIt2).first;
            idList[i+nbPoints] = (*mapIt2).first+M_numDof;
            idList[i+2*nbPoints] = (*mapIt2).first+2*M_numDof;
            ++i;
        }

        M_localEpetraMapPtr.reset( new EpetraMap(-1,3*nbPoints,idList,1,commPtr) );

        //-----------------------------------------------------
        // STEP 2: Compute normals and tangents
        //-----------------------------------------------------
        M_calculateNormals(mesh, dof,currentBdFE);
        M_storeGivenVersors();

        //this is used only for exporting the normals in vtk format.. should be put elsewhere?
        M_calculateCoordinates(mesh);
        M_calculateTangentVectors();
        M_buildRotationMatrix(systemMatrix,offset);
    }
}



template<typename MeshType,typename MatrixType> template <typename VectorType>
void BCNormalManager<MeshType,MatrixType>::bcShiftToNormalTangentialCoordSystem(matrix_Type& systemMatrix, VectorType& rightHandSide) const
{
    if (M_dataBuilt)
    {
        //Shift to tangential system

        int errCode(0);

        //C = R*A
        matrix_Type C(systemMatrix.getMap(), systemMatrix.getMeanNumEntries());
        errCode = M_rotationMatrixPtr->Multiply(false,systemMatrix,false,C);

        //A = C*Rt"
        matrix_Type D(systemMatrix.getMap(), systemMatrix.getMeanNumEntries());
        errCode = C.Multiply(false,*M_rotationMatrixPtr,true,D);
        //std::cout<< errCode <<std::endl;
        systemMatrix.swapCrsMatrix(D);

        //b = R*b
        VectorType c(rightHandSide);
        errCode = M_rotationMatrixPtr->Multiply(false,c,rightHandSide);
    }
}


template<typename MeshType,typename MatrixType> template <typename VectorType>
void BCNormalManager<MeshType,MatrixType>::bcShiftToCartesianCoordSystem(matrix_Type& systemMatrix, VectorType& rightHandSide) const
{
    if (M_dataBuilt)
    {
        int errCode(0);

        // C = Rt*A;
        matrix_Type C(systemMatrix.getMap(), systemMatrix.getMeanNumEntries());
        errCode = M_rotationMatrixPtr->Multiply(true,systemMatrix,false,C);

        // A = C*R";
        matrix_Type D(systemMatrix.getMap(), systemMatrix.getMeanNumEntries());
        errCode = C.Multiply(false,*M_rotationMatrixPtr,false,D);
        systemMatrix.swapCrsMatrix(D);

        // b = Rt*b;
        VectorType c(rightHandSide);
        errCode = M_rotationMatrixPtr->Multiply(true,c,rightHandSide);
    }
}



template<typename MeshType, typename MatrixType> template< typename VectorType>
void BCNormalManager<MeshType, MatrixType>::computeIntegratedNormals(const Dof& dof,CurrentBdFE& currentBdFE, VectorType& normals,  const MeshType& mesh)
{

    //-----------------------------------------------------
    // STEP 1: Calculating the normals
    //-----------------------------------------------------


    VectorType repNormals(normals.getMap(), Repeated);
    //Loop on the Faces
    for ( UInt iFace = 1; iFace<= mesh.numBElements(); ++iFace )
    {
        //Update the currentBdFE with the face data
        currentBdFE.updateMeasNormalQuadPt( mesh.bElement( iFace ) );
        ID idFace = mesh.bElement( iFace ).id();
        UInt nDofF = currentBdFE.nbNode();

        //For each node on the face
        for (UInt icheck = 1; icheck<= nDofF; ++icheck)
        {
            bool idFaceExist(false); //Is the face in the array?
            ID idf = dof.localToGlobalByFace(idFace,icheck,idFaceExist);

            //If the face exists and the point is on this processor
            if (idFaceExist && (M_flags.find(idf) != M_flags.end()))
            {
                ID flag = M_flags[idf];

                //if the normal is not already calculated
                //and the marker correspond to the flag of the point
                if ((flag == mesh.bElement(iFace).marker())||(flag == 0))
                {
                    //Warning: the normal is taken in the first Gauss point
                    //since the normal is the same over the triangle
                    //(not true in the case of quadratic and bilinear maps)
                    Real nx(currentBdFE.normal(0,0));
                    Real ny(currentBdFE.normal(1,0));
                    Real nz(currentBdFE.normal(2,0));

                    //We get the area
                    Real area(currentBdFE.measure());

                    //We update the normal component of the boundary point
                    (repNormals)[idf] += nx * area;
                    (repNormals)[idf+dof.numTotalDof()] += ny * area;
                    (repNormals)[idf+2*dof.numTotalDof()] += nz * area;
                }
            }
        }
    }

    //-----------------------------------------------------
    // STEP 2: Gathering the data from others processors
    //-----------------------------------------------------

    normals = VectorType(repNormals,Unique);

    //-----------------------------------------------------
    // STEP 3: Normalizing the vectors
    //-----------------------------------------------------

    //We obtain the ID of the element
    int NumMyElements = normals.getMap().getMap(Unique)->NumMyElements();
    int MyGlobalElements[NumMyElements];
    normals.getMap().getMap(Unique)->MyGlobalElements(MyGlobalElements);

    //We normalize the normal
    Real norm;
    UInt id;

    for ( int i(0); i<NumMyElements; ++i )
    {
        id = MyGlobalElements[i];

        //The id must be smaller than M_numDof
        //(the larger values are the y and z components)
        if (id <= dof.numTotalDof())
        {
            Real nx( (normals)[id] );
            Real ny( (normals)[id+dof.numTotalDof()] );
            Real nz( (normals)[id+2*dof.numTotalDof()] );
            norm = sqrt( nx*nx + ny*ny + nz*nz );
            (normals)[id]            /= norm;
            (normals)[id+dof.numTotalDof()]   /= norm;
            (normals)[id+2*dof.numTotalDof()] /= norm;
        }
    }
}

//TODO this function can be improved using the vtk writers
template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::exportToParaview(std::string fileName) const
{
    if (M_dataBuilt)
    {
        fileName.append("_proc");
        std::ostringstream ossMyPid;
        ossMyPid << M_localEpetraMapPtr->Comm().MyPID();
        fileName.append( ossMyPid.str() );
        fileName.append(".vtk");
        std::ofstream file(fileName.c_str());

        //Is the file open?
        if (file.fail())
        {
            std::cerr << "Error: The file " << fileName << " is not opened " << std::endl;
        }
        else
        {
            //We obtain the ID of the element
            int NumMyElements = M_localEpetraMapPtr->getMap(Unique)->NumMyElements();
            int MyGlobalElements[NumMyElements];
            M_localEpetraMapPtr->getMap(Unique)->MyGlobalElements(MyGlobalElements);
            ID idof(0);

            //Writing the header
            file << "# vtk DataFile Version 2.0" << std::endl;
            file << "Normal directions" << std::endl;
            file << "ASCII" << std::endl;

            //Writing the points
            file << "DATASET POLYDATA" << std::endl;
            file << "POINTS " << M_numInvoledDof << " float" << std::endl;
            for ( int i(0); i<NumMyElements; ++i )
            {
                idof = MyGlobalElements[i];

                //The id must be smaller than M_numDof
                //(the larger values are the y and z components)
                if (idof<=M_numDof)
                {
                    file << (*M_coordPtr)[idof] << "\t";
                    file << (*M_coordPtr)[idof+M_numDof] << "\t";
                    file << (*M_coordPtr)[idof+2*M_numDof] << std::endl;
                }
            }

            //Starting the data part of the file
            file << "POINT_DATA " << M_numInvoledDof << std::endl;

            //Writing t1
            file << "VECTORS cell_tangent_1 float" << std::endl;
            for ( int i(0); i<NumMyElements; ++i )
            {
                idof = MyGlobalElements[i];

                //The id must be smaller than M_numDof
                //(the larger values are the y and z components)
                if (idof<=M_numDof)
                {
                    file << (*M_firstTangentPtr)[idof] << "\t";
                    file << (*M_firstTangentPtr)[idof+M_numDof] << "\t";
                    file << (*M_firstTangentPtr)[idof+2*M_numDof] << std::endl;
                }
            }

            //Writing t2
            file << "VECTORS cell_tangent_2 float" << std::endl;
            for ( int i(0); i<NumMyElements; ++i )
            {
                idof = MyGlobalElements[i];

                //The id must be smaller than M_numDof
                //(the larger values are the y and z components)
                if (idof<=M_numDof)
                {
                    file << (*M_secondTangentPtr)[idof] << "\t";
                    file << (*M_secondTangentPtr)[idof+M_numDof] << "\t";
                    file << (*M_secondTangentPtr)[idof+2*M_numDof] << std::endl;
                }
            }

            //Writing n
            file << "VECTORS cell_normals float" << std::endl;
            for ( int i(0); i<NumMyElements; ++i )
            {
                idof = MyGlobalElements[i];

                //The id must be smaller than M_numDof
                //(the larger values are the y and z components)
                if (idof<=M_numDof)
                {
                    file << (*M_normalPtr)[idof] << "\t";
                    file << (*M_normalPtr)[idof+M_numDof] << "\t";
                    file << (*M_normalPtr)[idof+2*M_numDof] << std::endl;
                }
            }

            //Closing the file
            file.close();
        }
    }
}



//==============================================
// Private Methods
//==============================================


template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::M_calculateCoordinates(MeshType const& mesh)
{
    M_coordPtr.reset( new EpetraVector(*M_localEpetraMapPtr,Unique) );

    //We obtain the ID of the element
    int NumMyElements = M_localEpetraMapPtr->getMap(Unique)->NumMyElements();
    int MyGlobalElements[NumMyElements];
    M_localEpetraMapPtr->getMap(Unique)->MyGlobalElements(MyGlobalElements);

    UInt id;
    for ( int i(0); i<NumMyElements; ++i )
    {
        id = MyGlobalElements[i];

        //The id must be smaller than M_numDof
        if (id<=M_numDof)
        {
            for (UInt j(0); j<mesh.pointList.size(); ++j)
            {
                if (id==mesh.pointList[j].id())
                {
                    (*M_coordPtr)[id]            = mesh.pointList[j].x();
                    (*M_coordPtr)[id+M_numDof]   = mesh.pointList[j].y();
                    (*M_coordPtr)[id+2*M_numDof] = mesh.pointList[j].z();
                    j=mesh.pointList.size();
                }
            }
        }

    }
}


template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::M_addBoundaryPoint(const ID& idof,const ID& flag)
{
    M_flags.insert(std::pair<ID,ID>(idof,flag));
}


template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::M_addVersor(const ID& idof,const Real& vx,const Real& vy, const Real& vz)
{
    Vector n(3);
    n[0] = vx;
    n[1] = vy;
    n[2] = vz;
    M_givenVersors.insert(std::pair<ID, Vector>(idof,n));
}


template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::M_calculateNormals(const MeshType& mesh, const Dof& dof,CurrentBdFE& currentBdFE)
{
    //-----------------------------------------------------
    // STEP 1: Calculating the normals
    //-----------------------------------------------------

    M_normalPtr.reset ( new EpetraVector(*M_localEpetraMapPtr,Repeated) );
    //(*M_normalPtr)*=0;
    computeIntegratedNormals(dof, currentBdFE, *M_normalPtr, mesh);


    //-----------------------------------------------------
    // STEP 4: Cleaning the memory
    //-----------------------------------------------------

    M_flags.clear();
}

//! Save the imposed normal vectors.
template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::M_storeGivenVersors()
{
    //-----------------------------------------------------
    // STEP 1: Retrieve the normals
    //-----------------------------------------------------

    //We obtain the ID of the element
    int NumMyElements = M_localEpetraMapPtr->getMap(Unique)->NumMyElements();
    int MyGlobalElements[NumMyElements];
    M_localEpetraMapPtr->getMap(Unique)->MyGlobalElements(MyGlobalElements);

    //We normalize the normal
    Real norm;
    UInt id;
    for ( int i(0); i<NumMyElements; ++i )
    {
        id = MyGlobalElements[i];

        //The id must be smaller than M_numDof
        if ( ( id<=M_numDof ) && ( M_givenVersors.find(id)!=M_givenVersors.end() ) )
        {
            Real nx( (M_givenVersors)[id][0] );
            Real ny( (M_givenVersors)[id][1] );
            Real nz( (M_givenVersors)[id][2] );
            norm = sqrt( nx*nx + ny*ny + nz*nz );
            (*M_normalPtr)[id]            = nx/norm;
            (*M_normalPtr)[id+M_numDof]   = ny/norm;
            (*M_normalPtr)[id+2*M_numDof] = nz/norm;
        }
    }

    //-----------------------------------------------------
    // STEP 2: Cleaning the memory
    //-----------------------------------------------------

    M_givenVersors.clear();
}

//! Calculate the tangential vectors for each triad in the vector triad.
template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::M_calculateTangentVectors()
{
    //-----------------------------------------------------
    // STEP 1: Initialization
    //-----------------------------------------------------

    //We obtain the ID of the element
    int NumMyElements = M_localEpetraMapPtr->getMap(Unique)->NumMyElements();
    int MyGlobalElements[NumMyElements];
    M_localEpetraMapPtr->getMap(Unique)->MyGlobalElements(MyGlobalElements);

    //Building the tangential vectors
    M_firstTangentPtr.reset ( new EpetraVector(*M_localEpetraMapPtr,Unique) );
    (*M_firstTangentPtr)*=0;
    M_secondTangentPtr.reset ( new EpetraVector(*M_localEpetraMapPtr,Unique) );
    (*M_secondTangentPtr)*=0;

    //We are going to use the loop to count the number
    //of imposed DOF because since the EpetraVector is
    //unique now it is possible that we have less DOF
    M_numInvoledDof = 0;

    //-----------------------------------------------------
    // STEP 2: Calculation of the tangential vectors
    //-----------------------------------------------------

    // Real norm;
    UInt id;
    for ( int i(0); i<NumMyElements; ++i )
    {
        id = MyGlobalElements[i];

        //The id must be smaller than M_numDof
        //(the larger values are the y and z components)
        if (id<=M_numDof)
        {
            //Counting the number of DOF
            M_numInvoledDof++;

            //We take max{|n x i|,|n x j|,|n x k|}
            //			=max{sqrt(ny^2+nz^2),sqrt(nx^2+nz^2),sqrt(nx^2+ny^2)}
            //			=max{r1,r2,r3}
            Real nx( (*M_normalPtr)[id] );
            Real ny( (*M_normalPtr)[id+M_numDof] );
            Real nz( (*M_normalPtr)[id+2*M_numDof] );
            Real nxi=sqrt(ny*ny+nz*nz);
            Real nxj=sqrt(nx*nx+nz*nz);
            Real nxk=sqrt(nx*nx+ny*ny);

            if ((nxi>=nxj)&&(nxi>=nxk)) //max = |n x i|
            {
                //We create t1
                (*M_firstTangentPtr)[id]            = 0;
                (*M_firstTangentPtr)[id+M_numDof]   = nz/nxi;
                (*M_firstTangentPtr)[id+2*M_numDof] = -ny/nxi;

                //We create t2
                (*M_secondTangentPtr)[id]            = -nxi;
                (*M_secondTangentPtr)[id+M_numDof]   = nx*ny/nxi;
                (*M_secondTangentPtr)[id+2*M_numDof] = nx*nz/nxi;
            }
            else if ((nxj>=nxi)&&(nxj>=nxk)) //max = |n x j|
            {
                //We create t1
                (*M_firstTangentPtr)[id]            = -nz/nxj;
                (*M_firstTangentPtr)[id+M_numDof]   = 0;
                (*M_firstTangentPtr)[id+2*M_numDof] = nx/nxj;

                //We create t2
                (*M_secondTangentPtr)[id]            = nx*ny/nxj;
                (*M_secondTangentPtr)[id+M_numDof]   = -nxj;
                (*M_secondTangentPtr)[id+2*M_numDof] = ny*nz/nxj;
            }
            else //max = |n x k|
            {
                //We create t1
                (*M_firstTangentPtr)[id]            = ny/nxk;
                (*M_firstTangentPtr)[id+M_numDof]   = -nx/nxk;
                (*M_firstTangentPtr)[id+2*M_numDof] = 0;

                //We create t2
                (*M_secondTangentPtr)[id]            = nx*nz/nxk;
                (*M_secondTangentPtr)[id+M_numDof]   = ny*nz/nxk;
                (*M_secondTangentPtr)[id+2*M_numDof] = -nxk;
            }
        }
    }
}


template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::M_buildRotationMatrix(matrix_Type& systemMatrix, UInt offset)
{
    //Initialization of the map to store the normal vectors
    std::map< ID,std::vector< Real > >::iterator mapIt;

    //Creating the matrix
    M_rotationMatrixPtr.reset( new matrix_Type(systemMatrix.getMap(), systemMatrix.getMeanNumEntries() ) );

    //Adding one to the diagonal
    M_rotationMatrixPtr->insertOneDiagonal();

    int nbRows(3);
    int nbCols(3);
    double* values[nbCols];
    int Indices[3];
    for ( int n = 0; n < nbCols; ++n )
    {
        values[n] = new Real[nbRows];
        for ( int m = 0; m < nbRows; ++m )
        {
            values[n][m] = 0.0;
        }
    }

    std::vector<int> rows;
    std::vector<int> cols;

    //We obtain the ID of the element
    int NumMyElements = M_localEpetraMapPtr->getMap(Unique)->NumMyElements();
    int MyGlobalElements[NumMyElements];
    M_localEpetraMapPtr->getMap(Unique)->MyGlobalElements(MyGlobalElements);

    UInt id;
    for ( int i(0); i<NumMyElements; ++i )
    {
        id = MyGlobalElements[i];

        //The id must be smaller than M_numDof
        if (id<=M_numDof)
        {
            //...Except for the nodes where we make the rotation
            //Global Dof
            //idDof = boundaryCondition( i ) ->id() + ( boundaryCondition.component( j ) - 1 ) * totalDof;
            //i,j and k take values in [1,totalDof]
            Indices[0] = id + offset - 1;
            Indices[1] = id + M_numDof + offset - 1;
            Indices[2] = id + 2 * M_numDof + offset - 1;

            cols.clear();
            cols.push_back(Indices[0]);
            cols.push_back(Indices[1]);
            cols.push_back(Indices[2]);

            rows.clear();
            rows.push_back(Indices[0]);
            rows.push_back(Indices[1]);
            rows.push_back(Indices[2]);


            //Line i (first tangential vector)
            values[0][0] = (*M_firstTangentPtr)[id];
            values[1][0] = (*M_firstTangentPtr)[id+M_numDof];
            values[2][0] = (*M_firstTangentPtr)[id+2*M_numDof];

            //-1 because we added one to the diagonal
            M_rotationMatrixPtr->set_mat_inc(Indices[0],Indices[0],-1.0);

            //Line j (second tangential vector)
            values[0][1] = (*M_secondTangentPtr)[id];
            values[1][1] = (*M_secondTangentPtr)[id+M_numDof];
            values[2][1]= (*M_secondTangentPtr)[id+2*M_numDof];

            //-1 because we added one to the diagonal
            M_rotationMatrixPtr->set_mat_inc(Indices[1],Indices[1],-1.0);

            //Line k (normal vector)
            values[0][2] = (*M_normalPtr)[id];
            values[1][2] = (*M_normalPtr)[id+M_numDof];
            values[2][2] = (*M_normalPtr)[id+2*M_numDof];

            //-1 because we added one to the diagonal
            M_rotationMatrixPtr->set_mat_inc(Indices[2],Indices[2],-1.0);

            M_rotationMatrixPtr->set_mat_inc(nbCols, nbRows, cols, rows, values);
        }

    }

    for ( int n = 0; n < nbRows; ++n )
    {
        delete values[n];
    }

    M_rotationMatrixPtr->GlobalAssemble();
}


} //end of namespace LifeV

#endif
