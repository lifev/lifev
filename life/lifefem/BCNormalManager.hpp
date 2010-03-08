/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2009-12-02

  Copyright (C) 2009

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
  @file BCNormalManager.h
  @brief Class for imposing boundary condtitions.
  @version 1.0
  @author Gwenol Grandperrin
  @date 12/2009

  This file contains a class that manages all the data needed to
  create the rotation matrix that help to impose normal Dirichlet
  boundary conditions.
*/

#ifndef _BCNORMALMANAGER_
#define _BCNORMALMANAGER_

#include <Epetra_Comm.h>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifearray/EpetraVector.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifefem/currentBdFE.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <string>
#include <map>
#include <sstream>

namespace LifeV
{

template<typename MeshType,typename MatrixType>
class BCNormalManager
{
public:
    //Some typedef
    typedef std::map<ID,ID> FlagsMap;
    typedef std::map<ID,std::vector<Real> > NormalsMap;

    //Constructor/destructor
    BCNormalManager(const MeshType& mesh);
    ~BCNormalManager();

    //Assembling the data
    template<typename DataType>
    void init(const BCBase& BCb,const DataType& t);
    void addBoundaryPoint(const ID& idof,const ID& flag);
    void addNormalPoint(const ID& idof,const Real& vx,const Real& vy, const Real& vz);
    void build(const Dof& dof,CurrentBdFE& bdfem,MatrixType& A, UInt offset,const Epetra_Comm& comm);

    //Exporting the data
    void exportToParaview(std::string fileName) const;

    //To change the basis
    template <typename VectorType>
    void bcShiftToNormalTangentialCoordSystem(MatrixType& A, VectorType& b);
    template <typename VectorType>
    void bcShiftToCartesianCoordSystem(MatrixType& A, VectorType& b);
    template <typename VectorType>
    void computeIntegratedNormals(const Dof& dof,CurrentBdFE& bdfem, VectorType& normals,  const MeshType& mesh);

private:
    bool          M_dataBuilt;      // true if some normal has been set
    const MeshType* const M_mesh;   // Mesh that are related to the normal
    MatrixType*   M_rotMat;         // Rotation matrix
    EpetraMap*    M_idMap;          // Map use on this processor
    EpetraVector* M_tangent1;       // components of the 1st tangential vector
    EpetraVector* M_tangent2;       // components of the 2nd tangential vector
    EpetraVector* M_normal;         // components of the normal
    EpetraVector* M_coordinates;    // coordinates of the point
    UInt          M_numDof;         // Number of DOF
    UInt          M_numImposedDof;  // Number of DOF that are imposed
    FlagsMap      M_flags;          // flag of the boundary
    NormalsMap    M_imposedNormals; // Normal that are imposed by the user

    //Internal methods to build the triad
    void M_calculateCoordinates();
    void M_calculateNormals(const Dof& dof,CurrentBdFE& bdfem);
    void M_saveImposedNormals();
    void M_calculateTangentVectors();
    void M_createRotationMatrix(MatrixType& A, UInt offset=0);

    //Just a small utility function to make
    //the conversion between int and std::string
    std::string M_toString(const int& n) const;

};

//==============================================
//               IMPLEMENTATION
//==============================================

//! Constructor.
/*!
    @param mesh The mesh where the normal should be calculated
 */
template<typename MeshType,typename MatrixType>
BCNormalManager<MeshType,MatrixType>::BCNormalManager(const MeshType& mesh)
    :M_dataBuilt(false),M_mesh(&mesh),M_rotMat(0),M_idMap(0),M_tangent1(0),M_tangent2(0)
    ,M_normal(0),M_coordinates(0),M_numDof(0),M_numImposedDof(0),M_flags()
{

}

//! Destructor.
template<typename MeshType,typename MatrixType>
BCNormalManager<MeshType,MatrixType>::~BCNormalManager()
{
    if(M_tangent1!=NULL) //check if the pointers have been initialized
    {
        delete M_tangent1;
        delete M_tangent2;
        delete M_normal;
        delete M_coordinates;
        delete M_idMap;
    }
    if(M_rotMat!=NULL)
    {
        delete M_rotMat;
    }
}

//! Retrieve the DOF of the boundary
/*!
    This method will get the coordinates and the flag of the DOF on the boundary.
    In the directional case, it will also set the normal.
    @param BCb a BCBase object
    @param t the actual time
 */
template <typename MeshType,typename MatrixType> template<typename DataType>
void BCNormalManager<MeshType,MatrixType>::init(const BCBase& BCb,const DataType& t)
{
    // Loop on BC identifiers
    for ( ID i = 1; i <= BCb.list_size(); ++i )
    {
        const IdentifierEssential* pId = static_cast< const IdentifierEssential* >( BCb( i ) );

        if(BCb.mode()==Directional){
            const BCFunctionDirectional* pBcF = static_cast<const BCFunctionDirectional*>( BCb.pointerToFunctor() );
            Real nx(pBcF->vectFct(t,pId->x(),pId->y(),pId->z(),1));
            Real ny(pBcF->vectFct(t,pId->x(),pId->y(),pId->z(),2));
            Real nz(pBcF->vectFct(t,pId->x(),pId->y(),pId->z(),3));

            addNormalPoint(BCb(i)->id(),nx,ny,nz);
        }
        else
        {
            addBoundaryPoint(BCb(i)->id(),BCb.flag());
        }
    }
    M_dataBuilt = true; //Since vectors has been given we must apply the basis change.
}

//! This method add a dof in which we want to impose a Normal Dirichlet boundary condition.
/*!
    @param idof id of the dof
    @param flag id of the boundary flag
 */
template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::addBoundaryPoint(const ID& idof,const ID& flag)
{
    M_flags.insert(pair<ID,ID>(idof,flag));
}

//! This method add a dof in which we want to impose a Normal Dirichlet boundary condition with a given vector.
/*!
    @param idof id of the dof
    @param vx x component of the vector
    @param vy y component of the vector
    @param vz z component of the vector
 */
template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::addNormalPoint(const ID& idof,const Real& vx,const Real& vy, const Real& vz)
{
    std::vector<Real> n(3,0.0);
    n[0] = vx;
    n[1] = vy;
    n[2] = vz;
    M_imposedNormals.insert(pair<ID,std::vector<Real> >(idof,n));
}

//! This method calculate all the normal and tangential vectors.
/*!
    This method will also export the normal and tangential vectors for paraview into the file normalAndTangentialDirections.vtk.
    @param dof the dof class
    @param bdfem the current boundary finite element
    @param comm Epetra_Comm object
 */
template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::build(const Dof& dof,CurrentBdFE& bdfem,MatrixType& A, UInt offset,const Epetra_Comm& comm)
{
    if(M_dataBuilt)
    {
        //-----------------------------------------------------
        // STEP 1: Building the map
        //-----------------------------------------------------

        //First we build the map
        UInt nbPoints(M_flags.size()+M_imposedNormals.size());
        UInt i(0);
        int idList[3*nbPoints]; //3 times because we want 3 coordinates

        //We store the number of Degrees of Freedom
        M_numDof = dof.numTotalDof();

        //Creating the iterators to explore the maps
        std::map< ID,ID >::iterator mapIt;
        std::map< ID,std::vector<Real> >::iterator mapIt2;

        //Building the list
        for ( mapIt=M_flags.begin() ; mapIt != M_flags.end(); mapIt++ )
        {
            idList[i] = (*mapIt).first;
            idList[i+nbPoints] = (*mapIt).first+M_numDof;
            idList[i+2*nbPoints] = (*mapIt).first+2*M_numDof;
            ++i;
        }
        for ( mapIt2=M_imposedNormals.begin() ; mapIt2 != M_imposedNormals.end(); mapIt2++ )
        {
            idList[i] = (*mapIt2).first;
            idList[i+nbPoints] = (*mapIt2).first+M_numDof;
            idList[i+2*nbPoints] = (*mapIt2).first+2*M_numDof;
            ++i;
        }

        M_idMap = new EpetraMap(-1,3*nbPoints,idList,1,comm);

        //-----------------------------------------------------
        // STEP 2: Compute normals and tangents
        //-----------------------------------------------------
        M_calculateNormals(dof,bdfem);
        M_saveImposedNormals();
        M_calculateCoordinates();
        M_calculateTangentVectors();
        M_createRotationMatrix(A,offset);
    }
}

//! This function creates the rotation matrix R which change the basis.
/*!
    @param A the matrix of the problem
    @param offset that will be used if there is more than one unknown to recover the global ID
 */
template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::M_createRotationMatrix(MatrixType& A, UInt offset)
{
    //std::cout << "Creating the R matrix... ";
    //Initialization of the map to store the normal vectors
    std::map< ID,std::vector< Real > >::iterator mapIt;

    //Creating the matrix
    M_rotMat = new MatrixType(A.getMap(), A.getMeanNumEntries());

    //Adding one to the diagonal
    M_rotMat->insertOneDiagonal();

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
    int NumMyElements = M_idMap->getMap(Unique)->NumMyElements();
    int MyGlobalElements[NumMyElements];
    M_idMap->getMap(Unique)->MyGlobalElements(MyGlobalElements);

    UInt id;
	for ( UInt i(0);i<NumMyElements;++i )
    {
        id = MyGlobalElements[i];

        //The id must be smaller than M_numDof
        //(the larger values are the y and z components)
        if(id<=M_numDof)
        {
            //...Except for the nodes where we make the rotation
            //Global Dof
            //idDof = BCb( i ) ->id() + ( BCb.component( j ) - 1 ) * totalDof;
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
            values[0][0] = (*M_tangent1)[id];
            values[1][0] = (*M_tangent1)[id+M_numDof];
            values[2][0] = (*M_tangent1)[id+2*M_numDof];

            //-1 because we added one to the diagonal
            M_rotMat->set_mat_inc(Indices[0],Indices[0],-1.0);

            //Line j (second tangential vector)
            values[0][1] = (*M_tangent2)[id];
            values[1][1] = (*M_tangent2)[id+M_numDof];
            values[2][1]= (*M_tangent2)[id+2*M_numDof];

            //-1 because we added one to the diagonal
            M_rotMat->set_mat_inc(Indices[1],Indices[1],-1.0);

            //Line k (normal vector)
            values[0][2] = (*M_normal)[id];
            values[1][2] = (*M_normal)[id+M_numDof];
            values[2][2] = (*M_normal)[id+2*M_numDof];

            //-1 because we added one to the diagonal
            M_rotMat->set_mat_inc(Indices[2],Indices[2],-1.0);

            M_rotMat->set_mat_inc(nbCols, nbRows, cols, rows, values);
        }

    }

    for ( int n = 0; n < nbRows; ++n )
    {
        delete values[n];
    }

    M_rotMat->GlobalAssemble();
    //M_rotMat->removeZeros();

    //M_rotMat->spy("R");
}

//! Export the vectors (\tau_1,\tau_2,n) for each DOF on the boundary
/*!
    @param filename name of the file to store the triad informations
 */
template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::exportToParaview(std::string fileName) const
{
    if(M_dataBuilt)
    {
        fileName.append("_proc");
        fileName.append( M_toString( M_idMap->Comm().MyPID() ) );
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
            int NumMyElements = M_idMap->getMap(Unique)->NumMyElements();
            int MyGlobalElements[NumMyElements];
            M_idMap->getMap(Unique)->MyGlobalElements(MyGlobalElements);
            ID idof(0);

            //Writing the header
            file << "# vtk DataFile Version 2.0" << std::endl;
            file << "Normal directions" << std::endl;
            file << "ASCII" << std::endl;

            //Writing the points
            file << "DATASET POLYDATA" << std::endl;
            file << "POINTS " << M_numImposedDof << " float" << std::endl;
            for ( UInt i(0);i<NumMyElements;++i )
            {
                idof = MyGlobalElements[i];

                //The id must be smaller than M_numDof
                //(the larger values are the y and z components)
                if(idof<=M_numDof)
                {
                    file << (*M_coordinates)[idof] << "\t";
                    file << (*M_coordinates)[idof+M_numDof] << "\t";
                    file << (*M_coordinates)[idof+2*M_numDof] << std::endl;
                }
            }

            //Starting the data part of the file
            file << "POINT_DATA " << M_numImposedDof << std::endl;

            //Writing t1
            file << "VECTORS cell_tangent_1 float" << std::endl;
            for ( UInt i(0);i<NumMyElements;++i )
            {
                idof = MyGlobalElements[i];

                //The id must be smaller than M_numDof
                //(the larger values are the y and z components)
                if(idof<=M_numDof)
                {
                    file << (*M_tangent1)[idof] << "\t";
                    file << (*M_tangent1)[idof+M_numDof] << "\t";
                    file << (*M_tangent1)[idof+2*M_numDof] << std::endl;
                }
            }

            //Writing t2
            file << "VECTORS cell_tangent_2 float" << std::endl;
            for ( UInt i(0);i<NumMyElements;++i )
            {
                idof = MyGlobalElements[i];

                //The id must be smaller than M_numDof
                //(the larger values are the y and z components)
                if(idof<=M_numDof)
                {
                    file << (*M_tangent2)[idof] << "\t";
                    file << (*M_tangent2)[idof+M_numDof] << "\t";
                    file << (*M_tangent2)[idof+2*M_numDof] << std::endl;
                }
            }

            //Writing n
            file << "VECTORS cell_normals float" << std::endl;
            for ( UInt i(0);i<NumMyElements;++i )
            {
                idof = MyGlobalElements[i];

                //The id must be smaller than M_numDof
                //(the larger values are the y and z components)
                if(idof<=M_numDof)
                {
                    file << (*M_normal)[idof] << "\t";
                    file << (*M_normal)[idof+M_numDof] << "\t";
                    file << (*M_normal)[idof+2*M_numDof] << std::endl;
                }
            }

            //Closing the file
            file.close();
        }

    }
}

//! Retrieve the coordinates of the points of each dof
template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::M_calculateCoordinates()
{
    M_coordinates = new EpetraVector(*M_idMap,Unique);

    //We obtain the ID of the element
    int NumMyElements = M_idMap->getMap(Unique)->NumMyElements();
    int MyGlobalElements[NumMyElements];
    M_idMap->getMap(Unique)->MyGlobalElements(MyGlobalElements);

    UInt id;
	for ( UInt i(0);i<NumMyElements;++i )
    {
        id = MyGlobalElements[i];

        //The id must be smaller than M_numDof
        //(the larger values are the y and z components)
        if(id<=M_numDof)
        {
            for(UInt j(0);j<M_mesh->pointList.size();++j)
            {
                if(id==M_mesh->pointList[j].id())
                {
                    (*M_coordinates)[id]            = M_mesh->pointList[j].x();
                    (*M_coordinates)[id+M_numDof]   = M_mesh->pointList[j].y();
                    (*M_coordinates)[id+2*M_numDof] = M_mesh->pointList[j].z();
                    j=M_mesh->pointList.size();
                }
            }
        }

    }
}

//! Calculate the normal vectors
/*!
    @param dof the dof class
    @param bdfem the current boundary finite element
 */
template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::M_calculateNormals(const Dof& dof,CurrentBdFE& bdfem)
{
    // Author:	Gwenol Grandperrin
    // Date:	23.09.09
    //
    // This function calculate the normal vectors
    // and store the component in the triad vector

    //-----------------------------------------------------
    // STEP 1: Calculating the normals
    //-----------------------------------------------------

    M_normal = new EpetraVector(*M_idMap,Repeated);
    //(*M_normal)*=0;
    computeIntegratedNormals(dof, bdfem, *M_normal, *M_mesh);


    //-----------------------------------------------------
    // STEP 4: Cleaning the memory
    //-----------------------------------------------------

    M_flags.clear();
}

//! Save the imposed normal vectors.
template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::M_saveImposedNormals()
{
    //-----------------------------------------------------
    // STEP 1: Retrieve the normals
    //-----------------------------------------------------

    //We obtain the ID of the element
    int NumMyElements = M_idMap->getMap(Unique)->NumMyElements();
    int MyGlobalElements[NumMyElements];
    M_idMap->getMap(Unique)->MyGlobalElements(MyGlobalElements);

    //We normalize the normal
    Real norm;
    UInt id;
	for ( UInt i(0);i<NumMyElements;++i )
    {
        id = MyGlobalElements[i];

        //The id must be smaller than M_numDof
        //(the larger values are the y and z components)
        if(id<=M_numDof && M_imposedNormals.find(id)!=M_imposedNormals.end())
        {
            Real nx( (M_imposedNormals)[id][0] );
            Real ny( (M_imposedNormals)[id][1] );
            Real nz( (M_imposedNormals)[id][2] );
            norm = sqrt( nx*nx + ny*ny + nz*nz );
            (*M_normal)[id]            = nx/norm;
            (*M_normal)[id+M_numDof]   = ny/norm;
            (*M_normal)[id+2*M_numDof] = nz/norm;
        }
    }

    //-----------------------------------------------------
    // STEP 2: Cleaning the memory
    //-----------------------------------------------------

    M_imposedNormals.clear();
}

//! Calculate the tangential vectors for each triad in the vector triad.
template<typename MeshType,typename MatrixType>
void BCNormalManager<MeshType,MatrixType>::M_calculateTangentVectors()
{
    //-----------------------------------------------------
    // STEP 1: Initialization
    //-----------------------------------------------------

    //We obtain the ID of the element
    int NumMyElements = M_idMap->getMap(Unique)->NumMyElements();
    int MyGlobalElements[NumMyElements];
    M_idMap->getMap(Unique)->MyGlobalElements(MyGlobalElements);

    //Building the tangential vectors
    M_tangent1 = new EpetraVector(*M_idMap,Unique);
    (*M_tangent1)*=0;
    M_tangent2 = new EpetraVector(*M_idMap,Unique);
    (*M_tangent2)*=0;

    //We are going to use the loop to count the number
    //of imposed DOF because since the EpetraVector is
    //unique now it is possible that we have less DOF
    M_numImposedDof = 0;

    //-----------------------------------------------------
    // STEP 2: Calculation of the tangential vectors
    //-----------------------------------------------------

    Real norm;
    UInt id;
	for ( UInt i(0);i<NumMyElements;++i )
    {
        id = MyGlobalElements[i];

        //The id must be smaller than M_numDof
        //(the larger values are the y and z components)
        if(id<=M_numDof)
        {
            //Counting the number of DOF
            M_numImposedDof++;

            //We take max{|n x i|,|n x j|,|n x k|}
            //			=max{sqrt(ny^2+nz^2),sqrt(nx^2+nz^2),sqrt(nx^2+ny^2)}
            //			=max{r1,r2,r3}
            Real nx( (*M_normal)[id] );
            Real ny( (*M_normal)[id+M_numDof] );
            Real nz( (*M_normal)[id+2*M_numDof] );
            Real nxi=sqrt(ny*ny+nz*nz);
            Real nxj=sqrt(nx*nx+nz*nz);
            Real nxk=sqrt(nx*nx+ny*ny);

            if((nxi>=nxj)&&(nxi>=nxk)) //max = |n x i|
            {
                //We create t1
                (*M_tangent1)[id]            = 0;
                (*M_tangent1)[id+M_numDof]   = nz/nxi;
                (*M_tangent1)[id+2*M_numDof] = -ny/nxi;

                //We create t2
                (*M_tangent2)[id]            = -nxi;
                (*M_tangent2)[id+M_numDof]   = nx*ny/nxi;
                (*M_tangent2)[id+2*M_numDof] = nx*nz/nxi;
            }
            else if((nxj>=nxi)&&(nxj>=nxk)) //max = |n x j|
            {
                //We create t1
                (*M_tangent1)[id]            = -nz/nxj;
                (*M_tangent1)[id+M_numDof]   = 0;
                (*M_tangent1)[id+2*M_numDof] = nx/nxj;

                //We create t2
                (*M_tangent2)[id]            = nx*ny/nxj;
                (*M_tangent2)[id+M_numDof]   = -nxj;
                (*M_tangent2)[id+2*M_numDof] = ny*nz/nxj;
            }
            else //max = |n x k|
            {
                //We create t1
                (*M_tangent1)[id]            = ny/nxk;
                (*M_tangent1)[id+M_numDof]   = -nx/nxk;
                (*M_tangent1)[id+2*M_numDof] = 0;

                //We create t2
                (*M_tangent2)[id]            = nx*nz/nxk;
                (*M_tangent2)[id+M_numDof]   = ny*nz/nxk;
                (*M_tangent2)[id+2*M_numDof] = -nxk;
            }
        }
    }
}

template<typename MeshType,typename MatrixType>
std::string BCNormalManager<MeshType,MatrixType>::M_toString(const int& n) const
{
    std::ostringstream oss;
    oss << n;
    return oss.str();
}

//! This function will change A to fit the (t1,t2,n) basis
/*!
    This method is used to change the basis of the DOF locally. Then to apply
    A Dirichlet boundary condition in the t1,t2,n direction it is enough
    to apply a Dirichlet boundary condition to the x,y and z component respectively
    In fact it will perform the operations A:=R*A*Rt and R*b
    @param A the matrix of the problem
    @param b the right hand side
 */
template<typename MeshType,typename MatrixType> template <typename VectorType>
void BCNormalManager<MeshType,MatrixType>::bcShiftToNormalTangentialCoordSystem(MatrixType& A, VectorType& b)
{
    if(M_dataBuilt)
    {
        //std::cout << "Shift to tangential system" << std::endl;

        int errCode(0);

        //std::cout << "C = R*A" << std::endl;
        MatrixType C(A.getMap(), A.getMeanNumEntries());
        errCode = M_rotMat->Multiply(false,A,false,C);
        //std::cout<< errCode <<std::endl;

        //C.spy("RA");

        //std::cout << "A = C*Rt" << std::endl;
        MatrixType D(A.getMap(), A.getMeanNumEntries());
        errCode = C.Multiply(false,*M_rotMat,true,D);
        //std::cout<< errCode <<std::endl;
        A.swapCrsMatrix(D);

        //Here we have an error of 9.7e-04
        //After a test we have D.spy()==A.spy()
        //A.spy("RAR");

        //std::cout << "b = R*b" << std::endl;
        VectorType c(b);
        errCode = M_rotMat->Multiply(false,c,b);
        //std::cout<< errCode <<std::endl;
    }
}

//! This function will change A to fit the cartesian basis
/*!
    This method is used to do the opposite operation of the method bcShiftToNormalTangentialCoordSystem.
    In fact it will perform the operations A:=Rt*A*R and Rt*b
    @param A the matrix of the problem
    @param b the right hand side
 */
template<typename MeshType,typename MatrixType> template <typename VectorType>
void BCNormalManager<MeshType,MatrixType>::bcShiftToCartesianCoordSystem(MatrixType& A, VectorType& b)
{
    if(M_dataBuilt)
    {
        //std::cout << "Shift to cartesian system" << std::endl;

        int errCode(0);

        //std::cout << "C = Rt*A" << std::endl;
        MatrixType C(A.getMap(), A.getMeanNumEntries());
        errCode = M_rotMat->Multiply(true,A,false,C);
        //std::cout<< errCode <<std::endl;

        //std::cout << "A = C*R" << std::endl;
        MatrixType D(A.getMap(), A.getMeanNumEntries());
        errCode = C.Multiply(false,*M_rotMat,false,D);
        //std::cout<< errCode <<std::endl;
        A.swapCrsMatrix(D);

        //std::cout << "b = Rt*b" << std::endl;
        VectorType c(b);
        errCode = M_rotMat->Multiply(true,c,b);
        //std::cout<< errCode <<std::endl;
    }
}


//! This method computes the normals to a surface
/*!
    @param Dof the dof containing the local-to-global map and the number of dofs of the volume whose boundary
    is the surface of interest.
    @param bdfem the finite element of the boundary
    @param normals the templated output vector
    @param mesh the templated mesh
 */
template<typename MeshType, typename MatrixType> template< typename VectorType>
void BCNormalManager<MeshType, MatrixType>::computeIntegratedNormals(const Dof& dof,CurrentBdFE& bdfem, VectorType& normals,  const MeshType& mesh)
{

    //-----------------------------------------------------
    // STEP 1: Calculating the normals
    //-----------------------------------------------------


    VectorType repNormals(normals.getMap(), Repeated);
    //Loop on the Faces
    for ( UInt iFace = 1; iFace<= mesh.numBElements(); ++iFace )
    {
        //Update the bdfem with the face data
        bdfem.updateMeasNormalQuadPt( mesh.bElement( iFace ) );
        ID idFace = mesh.bElement( iFace ).id();
        UInt nDofF = bdfem.nbNode;

        //For each node on the face
        for (UInt icheck = 1; icheck<= nDofF; ++icheck)
        {
            bool idFaceExist(false); //Is the face in the array?
            ID idf = dof.localToGlobalByFace(idFace,icheck,idFaceExist);
            //std::cout << "idf = " << idf << std::endl;

            //If the face exists and the point is on this processor
            if(idFaceExist && (M_flags.find(idf) != M_flags.end()))
            {
                ID flag = M_flags[idf];

                //std::cout << "Flag: " << flag << std::endl;

                //if the normal is not already calculated
                //and the marker correspond to the flag of the point
                if((flag == mesh.bElement(iFace).marker())||(flag == 0))
                {
                    //Warning we take the normal in the first gauss point
                    //since the normal is the same over the triangle
                    Real nx(bdfem.normal(0,0));
                    Real ny(bdfem.normal(1,0));
                    Real nz(bdfem.normal(2,0));

                    //We get the area
                    Real area(bdfem.measure());

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

	for ( UInt i(0); i<NumMyElements; ++i )
    {
        id = MyGlobalElements[i];

        //The id must be smaller than M_numDof
        //(the larger values are the y and z components)
        if(id <= dof.numTotalDof())
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


} //end of namespace LifeV

#endif
