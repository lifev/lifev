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
 @brief   File containing a class for coloring meshes.
 
 @author D. Forti
 @date 06-2016
 */

#ifndef MESHCOLORING_H
#define MESHCOLORING_H

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{
    
//! MeshData - class for coloring mesh.
/*!
 @author D. Forti
 
 The class color the mesh
 */

class MeshColoring
{
public:
    
    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
    
    typedef Epetra_Comm comm_Type;
    typedef boost::shared_ptr< comm_Type > commPtr_Type;
    
    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr< vector_Type > vectorPtr_Type;

    //! @name Constructors & Destructor
    //@{
    
    //! Empty Constructor
    MeshColoring( const commPtr_Type& communicator );
    
    //! Destructor
    ~MeshColoring();
    
    //@}
    
    
    //! @name Methods
    //@{
    
    //! Print node to volume connections
    void printNodesToVolumes ( );
    
    //! Print volume to node connections
    void printVolumesToNodes ( );
    
    //! Print number of element neighbors per node
    void printNumNeighbors ( );
    
    //! Print volume to volume connections
    void printVolumesToVolumes ( );
    
    //! Print vector of colors
    void printColors ( );
    
    //! Check if colors have been assigned correctly
    /*!
     @return integer to check if there have been errors
     */
    int performCheck ( );
    
    //! Initial setup of the class
    void setup ( );
    
    //! Perform the coloring of the mesh
    void colorMesh ( );
    
    //! Vector of colors for the assembly
    void colorsForAssembly ( );
    
    //! Vector of colors used in the assembly
    void printColorsForAssembly ( );
    
    //! Info about coloring
    void printInfo ( );

    //@}
    
    
    //! @name Set Methods
    //@{
    
    //! Set the mesh
    /*!
     @param mesh mesh to be colored
     */
    
    void setMesh ( const meshPtr_Type& local_mesh );
    
    //@}
    
    
    //! @name Get Methods
    //@{
    
    //! Get the VectorEpetra with colors
    /*!
     @return Vector of colors
     */
    vectorPtr_Type getColors ( ) { return M_vectorColors; }
    
    std::vector<std::vector<UInt>> getColorsForAssembly ( ) { return M_colorsForAssembly; }
    
    //@}
    
private:
    
    // mesh
    meshPtr_Type M_mesh;

    // communicator
    commPtr_Type M_comm;
    
    // Element to Vertices connectivity
    std::vector<std::vector<UInt>> M_volumesToNodes;
    
    // Vertices to Element connectivity
    std::vector<std::vector<UInt>> M_nodesToVolumes;
    
    std::vector<UInt> M_numNeighbors;

    UInt M_maxNumNeighbors;
    
    UInt M_indexNodeMaxNumNeighbors;
    
    // Vector with colors
    std::vector<int> M_colors;
    
    // Number of colors
    UInt M_numColors;
    
    // finding neighbors
    std::vector<std::vector<UInt>> M_volumesToVolumes;
    
    // epetra vector with colors
    vectorPtr_Type M_vectorColors;
    
    std::vector<std::vector<UInt>> M_colorsForAssembly;

};


} // namespace LifeV

#endif  /* MESHCOLORING_H */
