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

    @contributor Simone Rossi <simone.rossi@epfl.ch>
    @maintainer Simone Rossi <simone.rossi@epfl.ch>

    This file contains a set of base utilities used to read vectorial data (mainly fiber
    and sheet directions) from different formats to VectorEpetra objects.
    Also other useful functions are collected here.
 */

#ifndef HEARTUTILITY_H
#define HEARTUTILITY_H 1

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

namespace LifeV
{

// Predeclaration

namespace ElectrophysiologyUtility
{

//! HeartUtility
/*!
 *  @author(s) Simone Rossi
 *
 *  \c HeartUtility contains methods for reading vectorial data  from different formats to VectorEpetra objects.
 *
 */

//! @name Methods
//@{

//! Read fiber vector field from HDF5 file
/*!
 * @param fiberVector VectorEpetra object for storing the vector field
 * @param fileName    Name of the HDF5 file to read from
 * @param localMesh   Pointer to the mesh
 */
template<typename Mesh> inline void importFibers (  boost::shared_ptr<VectorEpetra> fiberVector, const std::string& fileName, boost::shared_ptr< Mesh > localMesh, const std::string& filePath = "./" )
{
    typedef Mesh                                                                          mesh_Type;
    typedef ExporterData<mesh_Type>                                                       exporterData_Type;
    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh<LifeV::LinearTetra> > >  filterPtr_Type;
    typedef LifeV::ExporterHDF5< RegionMesh<LinearTetra> >                                hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                                            hdf5FilterPtr_Type;


    boost::shared_ptr<Epetra_Comm> comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > fiberSpace ( new FESpace< mesh_Type, MapEpetra > ( localMesh, "P1", 3, comm ) );

    exporterData_Type impData (exporterData_Type::VectorField, "fibers.00000", fiberSpace,
                               fiberVector, UInt (0), exporterData_Type::UnsteadyRegime);

    filterPtr_Type importer ( new hdf5Filter_Type() );
    importer -> setMeshProcId ( localMesh, comm -> MyPID() );
    importer -> setPrefix (fileName);
    importer -> setPostDir (filePath);
    importer -> readVariable (impData);
    importer -> closeFile();

}

//! Read scalar field from HDF5 file
/*!
 * @param vector      VectorEpetra object for storing the scalar field
 * @param fileName    Name of the HDF5 file to read from
 * @param fieldName   Name of the scalar field in the HDF5 file
 * @param localMesh   Pointer to the mesh
 */
template<typename Mesh> inline void importScalarField (  boost::shared_ptr<VectorEpetra> vector, const std::string& fileName, const std::string& fieldName, boost::shared_ptr< Mesh > localMesh  )
{
    typedef Mesh                                                                         mesh_Type;
    typedef ExporterData<mesh_Type>                                                      exporterData_Type;
    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh<LifeV::LinearTetra> > > filterPtr_Type;
    typedef LifeV::ExporterHDF5< RegionMesh<LinearTetra> >                               hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                                           hdf5FilterPtr_Type;


    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > feSpace ( new FESpace< mesh_Type, MapEpetra > ( localMesh, "P1", 1, comm ) );

    exporterData_Type impData (exporterData_Type::ScalarField, fieldName + ".00000", feSpace,
                               vector, UInt (0), exporterData_Type::UnsteadyRegime);

    filterPtr_Type importer ( new hdf5Filter_Type() );
    importer -> setMeshProcId ( localMesh, comm -> MyPID() );
    importer -> setPrefix (fileName);
    importer -> readVariable (impData);
    importer -> closeFile();

}

//! Read vector field from HDF5 file
/*!
 * @param vector      VectorEpetra object for storing the scalar field
 * @param fileName    Name of the HDF5 file to read from
 * @param fieldName   Name of the vector field in the HDF5 file
 * @param localMesh   Pointer to the mesh
 */
template<typename Mesh> inline void importVectorField (  boost::shared_ptr<VectorEpetra> vector, const std::string& fileName, const std::string& fieldName, boost::shared_ptr< Mesh > localMesh  , const std::string& postDir = "./")
{
    typedef Mesh                                                                         mesh_Type;
    typedef ExporterData<mesh_Type>                                                      exporterData_Type;
    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh<LifeV::LinearTetra> > > filterPtr_Type;
    typedef LifeV::ExporterHDF5< RegionMesh<LinearTetra> >                               hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                                           hdf5FilterPtr_Type;


    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > feSpace ( new FESpace< mesh_Type, MapEpetra > ( localMesh, "P1", 3, comm ) );

    exporterData_Type impData (exporterData_Type::VectorField, fieldName + ".00000", feSpace,
                               vector, UInt (0), exporterData_Type::UnsteadyRegime);

    filterPtr_Type importer ( new hdf5Filter_Type() );
    importer -> setMeshProcId ( localMesh, comm -> MyPID() );
    importer -> setPrefix (fileName);
    importer -> setPostDir (postDir);
    importer -> readVariable (impData);
    importer -> closeFile();

}

//! Read fiber field from text file
/*!
 * @param fiberVector VectorEpetra object for storing the vector field
 * @param fileName    Name of the HDF5 file to read from
 * @param filePath    Path of the HDF5 file to read from
 * @param format      The fibers can be in the text file in two different formats:
 *
 * format 0 = fibers saved as (fx, fy, fz) in each row
 *
 * format 1 = fibers saved as fx in each row for all the mesh
 *                            fy in each row for all the mesh
 *                            fz in each row for all the mesh
 */
inline void importFibersFromTextFile ( boost::shared_ptr<VectorEpetra> fiberVector, std::string fileName, std::string filePath, int format = 0 )
{
    typedef VectorEpetra                                    vector_Type;
    typedef boost::shared_ptr<vector_Type>                  vectorPtr_Type;

    std::ifstream fibers ( (filePath + fileName).c_str() );

    UInt NumGlobalElements =  fiberVector -> size();
    std::vector<Real> fiber_global_vector (NumGlobalElements);

    // Importing fibers
    for ( UInt i = 0; i < NumGlobalElements; ++i)
    {
        fibers >> fiber_global_vector[i];
        if ( fiber_global_vector[i] == 0 )
        {
            std::cout << "\nzero component!!!! \t";
            std::cout << "in: " << filePath  + fileName << "\n";

        }
    }
    int n = (*fiberVector).epetraVector().MyLength();
    int d = n / 3;
    int i (0);
    int j (0);
    int k (0);
    int offset = (*fiberVector).size() / 3;

    for (int l = 0; l < d; ++l)
    {
        i = (*fiberVector).blockMap().GID (l);
        j = (*fiberVector).blockMap().GID (l + d);
        k = (*fiberVector).blockMap().GID (l + 2 * d);
        if ( format == 0 )
        {
            (*fiberVector) [i] = fiber_global_vector[3 * i];
            (*fiberVector) [j] = fiber_global_vector[3 * i + 1];
            (*fiberVector) [k] = fiber_global_vector[3 * i + 2];
        }
        else
        {

            (*fiberVector) [i] = fiber_global_vector[ i ];
            (*fiberVector) [j] = fiber_global_vector[ i + offset ];
            (*fiberVector) [k] = fiber_global_vector[ i + 2 * offset ];
        }

        //normalizing
        Real norm = std::sqrt ( (*fiberVector) [i] * (*fiberVector) [i] + (*fiberVector) [j] * (*fiberVector) [j] + (*fiberVector) [k] * (*fiberVector) [k] );
        if ( norm != 0 )
        {
            (*fiberVector) [i] = (*fiberVector) [i] / norm;
            (*fiberVector) [j] = (*fiberVector) [j] / norm;
            (*fiberVector) [k] = (*fiberVector) [k] / norm;
        }
        else
        {
            std::cout << "\n\nThe fiber vector in the node: " << i << " has component:";
            std::cout << "\nx: " << fiber_global_vector [i];
            std::cout << "\ny: " << fiber_global_vector [i + offset];
            std::cout << "\nz: " << fiber_global_vector [i + 2 * offset];
            std::cout << "\nI will put it to: (f_x, f_y, f_z) = (1, 0, 0)\n\n";

            (*fiberVector) [i] = 1.;
            (*fiberVector) [j] = 0.;
            (*fiberVector) [k] = 0.;
        }



    }

    fiber_global_vector.clear();


}

//! Setup fiber field from unidirectional VectorSmall object
/*!
 * @param fiberVector    VectorEpetra object for storing the vector field
 * @param fiberDirection Direction of fiber vectors as a VectorSmall
 */
inline void setupFibers ( VectorEpetra& fiberVector, VectorSmall<3>& fiberDirection)
{
    int n1 = fiberVector.epetraVector().MyLength();
    int d1 = n1 / 3;
    fiberVector *= 0;
    int i (0);
    int j (0);
    int k (0);

    for ( int l (0); l < d1; l++)
    {
        i = fiberVector.blockMap().GID (l);
        j = fiberVector.blockMap().GID (l + d1);
        k = fiberVector.blockMap().GID (l + 2 * d1);
        fiberVector [i] = fiberDirection[0];
        fiberVector [j] = fiberDirection[1];
        fiberVector [k] = fiberDirection[2];
    }

}

//! Setup fiber field from unidirectional std::vector object
/*!
 * @param fiberVector    VectorEpetra object for storing the vector field
 * @param fiberDirection Direction of fiber vectors as a VectorSmall
 */
inline void setupFibers ( VectorEpetra& fiberVector, std::vector<Real>& fiberDirection)
{
    VectorSmall<3> fiberVectorSmall;
    fiberVectorSmall[0] = fiberDirection.at (0);
    fiberVectorSmall[1] = fiberDirection.at (1);
    fiberVectorSmall[2] = fiberDirection.at (2);
    setupFibers (fiberVector, fiberVectorSmall);
}

//! Setup fiber field from three real components
/*!
 * @param fiberVector VectorEpetra object for storing the vector field
 * @param fx          First component of vector
 * @param fy          Second component of vector
 * @param fz          Third component of vector
 */
inline void setupFibers ( VectorEpetra& fiberVector, Real fx, Real fy, Real fz)
{
    VectorSmall<3> fiberVectorSmall;
    fiberVectorSmall[0] = fx;
    fiberVectorSmall[1] = fy;
    fiberVectorSmall[2] = fz;
    setupFibers (fiberVector, fiberVectorSmall);
}


//! On the boundary with the specified flag, set the wanted value
/*!
 * @param vec         VectorEpetra object where we want to impose the boundary value
 * @param fullMesh    Pointer to the non partitioned mesh
 * @param value       value to set on that boundary
 * @param flag        flag of the boundary
 */
inline void setValueOnBoundary ( VectorEpetra& vec, boost::shared_ptr<  RegionMesh<LinearTetra> > fullMesh, Real value, UInt flag)
{

    for ( Int j (0); j < vec.epetraVector().MyLength() ; ++j )
    {
        if ( fullMesh -> point ( vec.blockMap().GID (j) ).markerID() == flag )
        {
            if ( vec.blockMap().LID ( vec.blockMap().GID (j) ) != -1 )
            {
                (vec) ( vec.blockMap().GID (j) ) = value;
            }
        }
    }
}


//! On the boundary with the specified flags, set the wanted value
/*!
 * @param vec         VectorEpetra object where we want to impose the boundary value
 * @param fullMesh    Pointer to the non partitioned mesh
 * @param value       value to set on that boundary
 * @param flags        flags of the boundary
 */
inline void setValueOnBoundary ( VectorEpetra& vec, boost::shared_ptr<  RegionMesh<LinearTetra> > fullMesh, Real value, std::vector<UInt> flags)
{
    for ( int j (0); j < vec.epetraVector().MyLength() ; ++j )
    {
        for ( UInt k (0); k < flags.size(); k++ )
        {
            if ( fullMesh -> point ( vec.blockMap().GID (j) ).markerID() == flags.at (k) )
            {
                if ( vec.blockMap().LID ( vec.blockMap().GID (j) ) != -1 )
                {
                    (vec) ( vec.blockMap().GID (j) ) = value;
                }
            }
        }
    }
}

//! Rescale a scalar field to be between requested bounds
/*!
 * @param vector      VectorEpetra object that contains the scalar field
 * @param minValue    Minimum value
 * @param maxValue    Maximum value
 * @param scaleFactor Additional scaling factor (defaults to 1)
 */
inline void rescaleVector ( VectorEpetra& vector, Real minValue, Real maxValue, Real scaleFactor = 1.0 )
{
    vector -= minValue;
    vector *= ( scaleFactor / ( maxValue - minValue ) );
}

//! Rescale a scalar field by a constant factor
/*!
 * @param vector      VectorEpetra object that contains the scalar field
 * @param scaleFactor Additional scaling factor (defaults to 1)
 */
inline void rescaleVector ( VectorEpetra& vector, Real scaleFactor = 1.0 )
{
    Real max = vector.maxValue();
    Real min = vector.minValue();
    rescaleVector ( vector, min, max, scaleFactor);
}

//! Rescale a scalar field by a constant factor on a boundary
/*!
 * @param vector      VectorEpetra object that contains the scalar field
 * @param fullMesh    pointer to the non partitioned mesh
 * @param flag        flag of the boundary
 * @param scaleFactor Additional scaling factor (defaults to 1)
 */
inline void rescaleVectorOnBoundary ( VectorEpetra& vector, boost::shared_ptr<  RegionMesh<LinearTetra> > fullMesh, UInt flag, Real scaleFactor = 1.0 )
{
    for ( Int j (0); j < vector.epetraVector().MyLength() ; ++j )
    {
        if ( fullMesh -> point ( vector.blockMap().GID (j) ).markerID() == flag )
        {
            if ( vector.blockMap().LID ( vector.blockMap().GID (j) ) != -1 )
            {
                (vector) ( vector.blockMap().GID (j) ) *= scaleFactor;
            }
        }
    }
}

//! Normalizes a vector field to unit length
/*!
 * @param vector      VectorEpetra object that contains the vector field
 */
inline void normalize ( VectorEpetra& vector )
{
    int n1 = vector.epetraVector().MyLength();
    int d1 = n1 / 3;
    int i (0);
    int j (0);
    int k (0);

    for ( int l (0); l < d1; l++)
    {
        i = vector.blockMap().GID (l);
        j = vector.blockMap().GID (l + d1);
        k = vector.blockMap().GID (l + 2 * d1);
        Real norm = std::sqrt ( vector[i] * vector[i] + vector[j] * vector[j] + vector[k] * vector[k] );
        if ( norm != 0 )
        {
            (vector) [i] = (vector) [i] / norm;
            (vector) [j] = (vector) [j] / norm;
            (vector) [k] = (vector) [k] / norm;
        }
        else
        {
            std::cout << "\n\nThe vector in the node: " << i << " has component:";
            std::cout << "\nx: " <<  vector[i];
            std::cout << "\ny: " <<  vector[j];
            std::cout << "\nz: " <<  vector[k];
            std::cout << "\nI will put it to: (v_x, v_y, v_z) = (1, 0, 0)\n\n";
        }

    }
}


//! Add random component to the fibers
/*!
 * @param fiberVector    VectorEpetra object for storing the vector field
 * @param fiberDirection Direction of fiber vectors as a VectorSmall
 */
inline void addNoiseToFibers ( VectorEpetra& fiberVector, Real magnitude = 0.01,  std::vector<bool> component =  std::vector<bool> (3, true) )
{
    int n1 = fiberVector.epetraVector().MyLength();
    int d1 = n1 / 3;
    int i (0);
    int j (0);
    int k (0);

    for ( int l (0); l < d1; l++)
    {
        if (component[0])
        {
            i = fiberVector.blockMap().GID (l);
            fiberVector [i] += magnitude * (0.01 * ( (std::rand() % 100 ) - 50.0 ) );
        }
        if (component[1])
        {
            j = fiberVector.blockMap().GID (l + d1);
            fiberVector [j] += magnitude * (0.01 * ( (std::rand() % 100 ) - 50.0 ) );
        }
        if (component[2])
        {
            k = fiberVector.blockMap().GID (l + 2 * d1);
            fiberVector [k] += magnitude * (0.01 * ( (std::rand() % 100 ) - 50.0 ) );
        }
    }

    ElectrophysiologyUtility::normalize (fiberVector);
}



typedef boost::function < Real (const Real&  t,
                                const Real&  x,
                                const Real&  y,
                                const Real&  z,
                                const ID&    i ) >   function_Type;
template<typename Mesh> inline void setFibersFromFunction (  boost::shared_ptr<VectorEpetra> vector, boost::shared_ptr< Mesh > localMesh, function_Type f)
{
    typedef Mesh                                                                         mesh_Type;
    typedef ExporterData<mesh_Type>                                                      exporterData_Type;
    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh<LifeV::LinearTetra> > > filterPtr_Type;
    typedef LifeV::ExporterHDF5< RegionMesh<LinearTetra> >                               hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                                           hdf5FilterPtr_Type;


    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > feSpace ( new FESpace< mesh_Type, MapEpetra > ( localMesh, "P1", 3, comm ) );

    feSpace->interpolate (
        static_cast<FESpace<RegionMesh<LinearTetra>, MapEpetra>::function_Type> (f),
        *vector, 0.0);
}
//@}

} // namespace HeartUtility

} // namespace LifeV

#endif /* HEARTUTILITY_H */

