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
  @brief This file provides an interface for post-processing with Paraview
     by writing files in VTK XML format
  @date 11-2010

  VTK_XML is now inherited from Exporter.

  Usage: two steps
  - first: add the variables using addVariable
  - second: call postProcess(  );

  @author Tiziano Passerini <tiziano@mathcs.emory.edu>
  @maintainer Tiziano Passerini <tiziano@mathcs.emory.edu>
*/

#ifndef EXPORTER_VTK_H
#define EXPORTER_VTK_H 1

#include <life/lifefilters/Exporter.hpp>
#include <life/lifefem/DOF.hpp>
#include <life/lifefem/CurrentFE.hpp>

namespace LifeV
{

template<typename MeshType>
class VTK_XML : public Exporter<MeshType> {

public:

    enum VTK_CELL{
        // Linear cells
        VTK_VERTEX = 1,
        VTK_POLY_VERTEX = 2,
        VTK_LINE = 3,
        VTK_POLY_LINE = 4,
        VTK_TRIANGLE = 5,
        VTK_TRIANGLE_STRIP = 6,
        VTK_POLYGON = 7,
        VTK_PIXEL = 8,
        VTK_QUAD = 9,
        VTK_TETRA = 10,
        VTK_VOXEL = 11,
        VTK_HEXAHEDRON = 12,
        VTK_WEDGE = 13,
        VTK_PYRAMID = 14,

        // Quadratic, isoparametric cells
        VTK_QUADRATIC_EDGE = 21,
        VTK_QUADRATIC_TRIANGLE = 22,
        VTK_QUADRATIC_QUAD = 23,
        VTK_QUADRATIC_TETRA = 24,
        VTK_QUADRATIC_HEXAHEDRON = 25

    };

    typedef MeshType                       mesh_Type;
    typedef Exporter<MeshType>             super;
    typedef typename super::meshPtr_Type   meshPtr_Type;
    typedef DOF*                           dofPtr_Type;
    typedef CurrentFE*                     fePtr_Type;
    typedef std::multimap<ExporterData::Where,ExporterData>::iterator
            whereToDataMapIterator_Type;

    //! Empty constructor
    VTK_XML();

    /* *
       Constructor for VTK_XML

       \param dfile the GetPot data file where you must provide an [exporter] section with:
       "start" (start index for filenames 0 for 000, 1 for 001 etc.),
       "save" (how many time steps per postprocessing)
       "multimesh" (=true if the mesh has to be saved at each post-processing step)

       \param mesh the mesh

       \param the prefix for the case file (ex. "test" for test.case)
     */
    VTK_XML(const GetPot& data_file, meshPtr_Type meshPtr, dofPtr_Type dofPtr, fePtr_Type M_fePtr,
            const std::string prefix);


    /* *
       Post-process the variables added to the list

       \param time the solver time
     */
    virtual void postProcess(const Real& time);

    //! Import data from previous simulations at a certain time
    /*!
       @param Time the time of the data to be imported
     */
    virtual UInt importFromTime( const Real& /*Time*/ ) { return -1; }

    //! Import data from previous simulations
    /*!
       @param time the solver time
       @param dt time step used to rebuild the history up to now
    */
    virtual void import(const Real& /*Tstart*/, const Real& /*dt*/) {}

    //! Read  only last timestep
    virtual void import(const Real& /*Tstart*/) {}


private:

    UInt whichCellType();

    void M_wr_ascii();
    void composeHeader();
    void composeFooter();
    void composeTypeDataHeader(const ExporterData::Where& where,
                               std::stringstream& typeDataHeaderStringStream);
    void composeTypeDataFooter(const ExporterData::Where& where,
                               std::stringstream& typeDataFooterStringStream);
    void composeDataArrays(const ExporterData::Where& where,
                           std::stringstream& dataArraysStringStream);

    virtual void M_rd_scalar( ExporterData& /*dvar*/ ) {}
    virtual void M_rd_vector( ExporterData& /*dvar*/ ) {}

    meshPtr_Type M_meshPtr;
    dofPtr_Type  M_dofPtr;
    fePtr_Type   M_fePtr;

    std::stringstream M_headerStringStream, M_footerStringStream,
    M_pointDataHeaderStringStream, M_cellDataHeaderStringStream,
    M_pointDataFooterStringStream, M_cellDataFooterStringStream;
};


//
// Implementation
//

template<typename Mesh>
VTK_XML<Mesh>::VTK_XML():
    super()
{
}


template<typename Mesh> VTK_XML<Mesh >::VTK_XML(
        const GetPot& data_file, meshPtr_Type meshPtr, dofPtr_Type dofPtr,
        fePtr_Type M_fePtr, const std::string prefix)
    :
    super(data_file, prefix),
    M_meshPtr(meshPtr), // pointer to mesh
    M_dofPtr(dofPtr),
    M_fePtr(M_fePtr)
{
    composeHeader();
    composeFooter();
}


template <typename Mesh>
UInt VTK_XML<Mesh>::whichCellType()
{
    UInt vtkCellType;

    switch ( M_fePtr->refFE.type )
    {
    case FE_P1_2D:
    case FE_P1bubble_2D:
        vtkCellType = VTK_TRIANGLE;
        break;
    case FE_P2_2D:
        vtkCellType = VTK_QUADRATIC_TRIANGLE;
        break;
    case FE_P1_3D:
        vtkCellType = VTK_TETRA;
        break;
    case FE_P2_3D:
    case FE_P2tilde_3D:
        vtkCellType = VTK_QUADRATIC_TETRA; // maybe to be modified (P2 on linear tetra...) see also the geomap...
        break;
    case FE_Q1_2D:
        vtkCellType = VTK_QUAD;
        break;
    case FE_Q2_2D:
        vtkCellType = VTK_QUADRATIC_QUAD;
        break;
    case FE_Q0_3D:
    case FE_Q1_3D:
        vtkCellType = VTK_HEXAHEDRON;
        break;
    case FE_Q2_3D:
        vtkCellType = VTK_QUADRATIC_HEXAHEDRON;
        break;
    default:
        std::cout << "WARNING: the element is not yet implemented in vtk_wrtrs.h\n";
    }

    return vtkCellType;
}


template<typename Mesh>
void VTK_XML<Mesh>::postProcess(const Real& /*time*/)
{
    typedef std::list< ExporterData >::iterator Iterator;

    this->computePostfix();

    std::size_t found( this->M_postfix.find( "*" ) );
    if ( found == string::npos )
    {
        std::cout << "  x-  VTK post-processing ...        " << std::flush;

        // redundant. should be done just the first time
        composeTypeDataHeader(ExporterData::Node, M_pointDataHeaderStringStream);
        composeTypeDataFooter(ExporterData::Node, M_pointDataFooterStringStream);
        composeTypeDataHeader(ExporterData::Cell, M_cellDataHeaderStringStream);
        composeTypeDataFooter(ExporterData::Cell, M_cellDataFooterStringStream);

        std::ofstream vtkFile( ( this->M_post_dir+this->M_prefix+this->M_postfix+".vtu").c_str() );
        std::stringstream buffer("");

        Chrono chrono;
        chrono.start();

        vtkFile << M_headerStringStream.str();

        vtkFile << M_pointDataHeaderStringStream.str();
        composeDataArrays(ExporterData::Node, buffer);
        vtkFile << buffer.str();
        vtkFile << M_pointDataFooterStringStream.str();

        buffer.str("");

        vtkFile << M_cellDataHeaderStringStream.str();
        composeDataArrays(ExporterData::Cell, buffer);
        vtkFile << buffer.str();
        vtkFile << M_cellDataFooterStringStream.str();

        vtkFile << M_footerStringStream.str();

        vtkFile.close();

        chrono.stop();

        M_pointDataHeaderStringStream.str("");
        M_pointDataFooterStringStream.str("");
        M_cellDataHeaderStringStream.str("");
        M_cellDataFooterStringStream.str("");
        std::cout << "      done in " << chrono.diff() << " s." << std::endl;
    }
}


template <typename Mesh>
void
VTK_XML<Mesh>::composeDataArrays(const ExporterData::Where& where,
                                 std::stringstream& dataArraysStringStream)
{
    Debug(8000) << "\n[VTK_XML::composeDataArrays] where = " << where << "\n";
    // Debug(8000) << "\n[VTK_XML::composeDataArrays] dataArraysStringStream = " << dataArraysStringStream;
    // std::stringstream dataArraysStringStream;

    whereToDataMapIterator_Type it;
    std::pair<whereToDataMapIterator_Type, whereToDataMapIterator_Type> rangeFound;

    rangeFound = this->M_whereToDataMap.equal_range(where);

    if( rangeFound.first != rangeFound.second )
    {
        Debug(8000) << "\n[VTK_XML::composeDataArrays] found data\n";

        dataArraysStringStream.setf(std::ios_base::fixed);
        dataArraysStringStream.precision(5);
        dataArraysStringStream.width(12);

        for (it = rangeFound.first; it != rangeFound.second; ++it)
        {
            switch( it->second.type() )
            {
            case ExporterData::ScalarData:

                dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
                        << it->second.variableName() << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";

                for (UInt i=0; i<it->second.size(); ++i){
                    dataArraysStringStream << it->second( i ) << " ";
                }
                break;
            case ExporterData::VectorData:

                dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
                        << it->second.variableName() << "\" NumberOfComponents=\""
                        << nDimensions << "\" format=\"ascii\">\n";

                for (UInt i=0; i<it->second.size(); ++i){
                    for (UInt icoor=0; icoor< nDimensions;++icoor){
                        dataArraysStringStream << it->second( i + icoor * it->second.size() ) << " ";
                    }
                }
                break;
            case ExporterData::TensorData:

                dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
                        << it->second.variableName() << "\" NumberOfComponents=\""
                        << nDimensions*nDimensions << "\" format=\"ascii\">\n";

                for (UInt i=0; i<it->second.size(); ++i){
                    for (UInt icoor=0; icoor< nDimensions;++icoor){
                        for (UInt jcoor=0; jcoor< nDimensions;++jcoor){
                            dataArraysStringStream << it->second( i * nDimensions * nDimensions + icoor * nDimensions + jcoor ) << " ";
                        }
                    }
                }
                break;
            }
            dataArraysStringStream << "\n\t\t\t\t</DataArray>\n";
        }
        // return dataArraysStringStream;
    }
}


template <typename Mesh>
void VTK_XML<Mesh>::composeHeader()
{
    Debug(8000) << "\n[VTK_XML::composeHeader]\n";

    ID numVertices = M_meshPtr->numVertices();
    ID numElements = M_meshPtr->numElements();
    UInt numLocalDof = M_dofPtr->numLocalDof();

    UInt numPoints = M_dofPtr->numTotalDof();
    UInt numNonVertexPoints = M_dofPtr->numTotalDof() - numVertices;

    std::vector<Vector> coordinatesOfNonVertexPoints( nDimensions, ZeroVector(numPoints) );

    Real x, y, z;

    // ASSERT( vtkFile, "Error: Output file cannot be open" );

    //header part of the file
    M_headerStringStream << "<?xml version=\"1.0\"?>\n";
    M_headerStringStream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    M_headerStringStream << "\t<UnstructuredGrid>\n";

    M_headerStringStream << "\t\t<Piece NumberOfPoints=\"" << numPoints << "\""
            << " NumberOfCells=\"" << numElements << "\">\n";


    M_headerStringStream << "\t\t\t<Points>\n";
    M_headerStringStream << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"" << nDimensions
            << "\" format=\"ascii\">\n";
    // Vertex based Dof: the coordinates are available from the Point List
    for(ID iVertex=1; iVertex <= numVertices; ++iVertex)
    {
        for (ID jCoor=1; jCoor <= nDimensions; ++jCoor)
        {
            M_headerStringStream << M_meshPtr->point(iVertex).coordinate(jCoor) << " ";
        }
    }

    // Now I store the coordinates of the supplementary nodes in a temporary vector
    for ( UInt iElement = 1; iElement <= numElements; ++iElement )
    {
        M_fePtr->updateJac( M_meshPtr->element( iElement ) );
        for ( UInt iPoint = M_dofPtr->numLocalVertices(); iPoint < M_dofPtr->numLocalDof(); ++iPoint )
        {
            M_fePtr->coorMap( x, y, z, M_fePtr->refFE.xi( iPoint ), M_fePtr->refFE.eta( iPoint ), M_fePtr->refFE.zeta( iPoint ) );
            UInt index = M_dofPtr->localToGlobalMap( iElement, iPoint ) - numVertices;
            coordinatesOfNonVertexPoints[0][index] = x;
            coordinatesOfNonVertexPoints[1][index] = y;
            coordinatesOfNonVertexPoints[2][index] = z;
        }
    }
    for ( UInt iPoint = 0; iPoint < numNonVertexPoints; ++iPoint )
    {
        for( UInt iCoor = 0; iCoor < nDimensions; ++iCoor )
            M_headerStringStream << coordinatesOfNonVertexPoints[iCoor][ iPoint ] << " ";
    }
    M_headerStringStream << "\n\t\t\t\t</DataArray>\n";
    M_headerStringStream << "\t\t\t</Points>\n";

    // connectivity
    // cells_size = nldof*(nV+1);
    M_headerStringStream << "\t\t\t<Cells>\n";
    M_headerStringStream << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    for (UInt iElement=1; iElement <= numElements; ++iElement)
    {
        for( UInt jPoint=0; jPoint < numLocalDof; ++jPoint)
        {
            // vtkFile.width(8);
            M_headerStringStream << M_dofPtr->localToGlobalMap( iElement, jPoint ) << " ";
        }
    }
    M_headerStringStream << "\n\t\t\t\t</DataArray>\n";

    M_headerStringStream << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";

    for (UInt iElement=1; iElement <= numElements; ++iElement)
    {
        M_headerStringStream << iElement * numLocalDof << " ";
    }

    M_headerStringStream << "\n\t\t\t\t</DataArray>\n";

    M_headerStringStream << "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";

    // definition of cells types
    UInt vtkCellType = whichCellType();
    for (UInt iElement=1; iElement <= numElements; ++iElement)
    {
        M_headerStringStream << vtkCellType << " ";
    }

    M_headerStringStream << "\n\t\t\t\t</DataArray>\n";
    M_headerStringStream << "\t\t\t</Cells>\n";

//    vtkFile << "\t\t</Piece>\n";
//    vtkFile << "\t</UnstructuredGrid>\n";
//    vtkFile << "</VTKFile>\n";

//    vtkFile.close();
}


template <typename Mesh>
void VTK_XML<Mesh>::composeFooter()
{
    Debug(8000) << "\n[VTK_XML::composeFooter]\n";

    //footer part of the file
    M_footerStringStream << "\t\t</Piece>\n";
    M_footerStringStream << "\t</UnstructuredGrid>\n";
    M_footerStringStream << "</VTKFile>\n";

}


template <typename Mesh>
void
VTK_XML<Mesh>::composeTypeDataHeader(const ExporterData::Where& where,
                                     std::stringstream& dataHeaderStringStream)
{
    Debug(8000) << "\n[VTK_XML::composeTypeDataHeader] where = " << where << "\n";
    // std::stringstream M_dataHeaderStringStream;

    whereToDataMapIterator_Type it;
    std::pair<whereToDataMapIterator_Type, whereToDataMapIterator_Type> rangeFound;

    // every dvar will produce a different data set
    // I'm just separating CELL data from POINT data
    rangeFound = this->M_whereToDataMap.equal_range(where);

    if( rangeFound.first != rangeFound.second )
    {
        Debug(8000) << "\n[VTK_XML::composeTypeDataHeader] found data\n";

        std::string whereString;
        switch( where )
        {
        case ExporterData::Node:
            whereString = "PointData";
            break;
        case ExporterData::Cell:
            whereString = "CellData";
            break;
        }
        dataHeaderStringStream << "\t\t\t<" << whereString << " ";
/*
        for (it = rangeFound.first; it != rangeFound.second; ++it)
        {
            switch( it->second.type() )
            {
            case ExporterData::ScalarData:
                dataHeaderStringStream << "Scalars=\"" << it->second.variableName() << "\" ";
                break;
            case ExporterData::VectorData:
                dataHeaderStringStream << "Vectors=\"" << it->second.variableName() << "\" ";
                break;
            case ExporterData::TensorData:
                M_pointDataHeaderStringStream << "Tensors=\"" << it->second.variableName() << "\" ";
                break;
            }
        }
*/
        dataHeaderStringStream << ">\n";
    }

    // return M_dataHeaderStringStream;
}


template <typename Mesh>
void
VTK_XML<Mesh>::composeTypeDataFooter(const ExporterData::Where& where,
                                     std::stringstream& dataFooterStringStream)
{
    // std::stringstream dataFooterStringStream;

    whereToDataMapIterator_Type it;
    std::pair<whereToDataMapIterator_Type, whereToDataMapIterator_Type> rangeFound;

    // every dvar will produce a different data set
    // I'm just separating CELL data from POINT data
    rangeFound = this->M_whereToDataMap.equal_range(where);

    if( rangeFound.first != rangeFound.second )
    {
        Debug(8000) << "\n[VTK_XML::composeTypeDataFooter] where = " << where << "\n";

        std::string whereString;
        switch( where )
        {
        case ExporterData::Node:
            whereString = "PointData";
            break;
        case ExporterData::Cell:
            whereString = "CellData";
            break;
        }
        dataFooterStringStream << "\t\t\t</" << whereString << ">\n";

    }
    // return dataFooterStringStream;
}

}
#endif
