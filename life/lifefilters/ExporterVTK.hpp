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

  ExporterVTK is now inherited from Exporter.

  Usage: two steps
  - first: add the variables using addVariable
  - second: call postProcess(  );

  @author Tiziano Passerini <tiziano@mathcs.emory.edu>
  @maintainer Tiziano Passerini <tiziano@mathcs.emory.edu>
*/

#ifndef EXPORTER_VTK_H
#define EXPORTER_VTK_H 1

#include <life/lifefilters/Exporter.hpp>
#include <life/lifefem/FESpace.hpp>

namespace LifeV
{

template<typename MeshType>
class ExporterVTK : public Exporter<MeshType> {

public:

    //! @name Public Types
    //@{

    /*! @enum VTK_CELL
        A list of the cell types admitted in VTK
     */
    enum VTK_CELL {
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
    typedef boost::shared_ptr<FESpace>     feSpacePtr_Type;
    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    ExporterVTK();

    /* *
       Constructor for ExporterVTK

       \param data_file the GetPot data file where you must provide an [exporter] section with:
       "start" (start index for filenames 0 for 000, 1 for 001 etc.),
       "save" (how many time steps per postprocessing)
       "multimesh" (=true if the mesh has to be saved at each post-processing step)

       \param the prefix for the output file (ex. "test" for test.vtu)
     */
//    ExporterVTK(const GetPot& data_file, meshPtr_Type meshPtr, dofPtr_Type dofPtr, fePtr_Type M_fePtr,
//            const std::string prefix);
    ExporterVTK(const GetPot& data_file, const std::string prefix);

    //! Copy constructor
    // ExporterVTK( const ExporterVTK& example );

    //! Destructor
    ~ExporterVTK() {}

    //@}

    //! @name Methods
    //@{
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

    //@}

    //! @name Get methods
    //@{

    //! returns the type of the map to use for the VectorEpetra
    MapEpetraType mapType() const;

    //@}

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

    virtual void readScalar( ExporterData& /*dvar*/ ) {}
    virtual void readVector( ExporterData& /*dvar*/ ) {}

    meshPtr_Type    M_meshPtr;
    // dofPtr_Type     M_dofPtr;
    feSpacePtr_Type M_feSpacePtr;

    std::stringstream M_headerStringStream, M_footerStringStream,
    M_pointDataHeaderStringStream, M_cellDataHeaderStringStream,
    M_pointDataFooterStringStream, M_cellDataFooterStringStream;
};


//
// Implementation
//

template<typename Mesh>
ExporterVTK<Mesh>::ExporterVTK():
        super()
{
}


template<typename Mesh> ExporterVTK<Mesh >::ExporterVTK(
    const GetPot& data_file,
    // meshPtr_Type meshPtr, dofPtr_Type dofPtr, fePtr_Type M_fePtr,
    const std::string prefix)
        :
        super(data_file, prefix)
//    M_meshPtr(meshPtr), // pointer to mesh
//    M_dofPtr(dofPtr),
//    M_fePtr(M_fePtr)
{
    // composeHeader();
    // composeFooter();
}


template <typename Mesh>
UInt ExporterVTK<Mesh>::whichCellType()
{
    ASSERT( M_feSpacePtr.get(), "\nA pointer to a valid FE object is required!");

    UInt vtkCellType;

    switch ( M_feSpacePtr->fe().refFE().type() )
    {
    case FE_P1_2D:
        // case FE_P1bubble_2D:
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
void ExporterVTK<Mesh>::postProcess(const Real& /*time*/)
{
    // typedef std::list< ExporterData >::iterator Iterator;

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

        std::ofstream vtkFile( ( this->M_postDir+this->M_prefix+this->M_postfix+".vtu").c_str() );
        std::stringstream buffer("");

        LifeChrono chrono;
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
ExporterVTK<Mesh>::composeDataArrays(const ExporterData::Where& where,
                                     std::stringstream& dataArraysStringStream)
{
    Debug(8000) << "\n[ExporterVTK::composeDataArrays] where = " << where << "\n";
    // Debug(8000) << "\n[ExporterVTK::composeDataArrays] dataArraysStringStream = " << dataArraysStringStream;
    // std::stringstream dataArraysStringStream;

    typename super::iterator_Type it;
    std::pair<typename super::iterator_Type, typename super::iterator_Type> rangeFound;

    rangeFound = this->M_whereToDataMap.equal_range(where);

    if ( rangeFound.first != rangeFound.second )
    {
        Debug(8000) << "\n[ExporterVTK::composeDataArrays] found data\n";

        dataArraysStringStream.setf(std::ios_base::fixed);
        dataArraysStringStream.precision(5);
        dataArraysStringStream.width(12);

        for (it = rangeFound.first; it != rangeFound.second; ++it)
        {
            switch ( it->second.type() )
            {
            case ExporterData::ScalarField:

                dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
                        << it->second.variableName() << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";

                for (UInt i=0; i<it->second.size(); ++i) {
                    dataArraysStringStream << it->second( i ) << " ";
                }
                break;
            case ExporterData::VectorField:

                dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
                        << it->second.variableName() << "\" NumberOfComponents=\""
                        << nDimensions << "\" format=\"ascii\">\n";

                for (UInt i=0; i<it->second.size(); ++i) {
                    for (UInt icoor=0; icoor< nDimensions; ++icoor) {
                        dataArraysStringStream << it->second( i + icoor * it->second.size() ) << " ";
                    }
                }
                break;
                /*case ExporterData::TensorField:

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
                    break;*/
            }
            dataArraysStringStream << "\n\t\t\t\t</DataArray>\n";
        }
        // return dataArraysStringStream;
    }
}


template <typename Mesh>
void ExporterVTK<Mesh>::composeHeader()
{
    ASSERT( M_meshPtr.get(), "\nA pointer to a valid mesh object is required!");
    ASSERT( M_feSpacePtr.get(), "\nA pointer to a valid FESpace object is required!");

    Debug(8000) << "\n[ExporterVTK::composeHeader]\n";

    ID numVertices = M_meshPtr->numVertices();
    ID numElements = M_meshPtr->numElements();
    UInt numLocalDof = M_feSpacePtr->dof().numLocalDof();

    UInt numPoints = M_feSpacePtr->dof().numTotalDof();
    UInt numNonVertexPoints = M_feSpacePtr->dof().numTotalDof() - numVertices;

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
    for (ID iVertex=1; iVertex <= numVertices; ++iVertex)
    {
        for (ID jCoor=1; jCoor <= nDimensions; ++jCoor)
        {
            M_headerStringStream << M_meshPtr->point(iVertex).coordinate(jCoor) << " ";
        }
    }

    // Now I store the coordinates of the supplementary nodes in a temporary vector
    for ( UInt iElement = 1; iElement <= numElements; ++iElement )
    {
        M_feSpacePtr->fe().updateJac( M_meshPtr->element( iElement ) );
        for ( UInt iPoint = M_feSpacePtr->dof().numLocalVertices()+1; iPoint <= M_feSpacePtr->dof().numLocalDof(); ++iPoint )
        {
            M_feSpacePtr->fe().coorMap( x, y, z,
                                        M_feSpacePtr->fe().refFE().xi( iPoint ), M_feSpacePtr->fe().refFE().eta( iPoint ), M_feSpacePtr->fe().refFE().zeta( iPoint ) );
            UInt index = M_feSpacePtr->dof().localToGlobal( iElement, iPoint ) - numVertices;
            coordinatesOfNonVertexPoints[0][index] = x;
            coordinatesOfNonVertexPoints[1][index] = y;
            coordinatesOfNonVertexPoints[2][index] = z;
        }
    }
    for ( UInt iPoint = 0; iPoint < numNonVertexPoints; ++iPoint )
    {
        for ( UInt iCoor = 0; iCoor < nDimensions; ++iCoor )
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
        for ( UInt jPoint=1; jPoint <= numLocalDof; ++jPoint)
        {
            // vtkFile.width(8);
            M_headerStringStream << M_feSpacePtr->dof().localToGlobal( iElement, jPoint ) << " ";
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
void ExporterVTK<Mesh>::composeFooter()
{
    Debug(8000) << "\n[ExporterVTK::composeFooter]\n";

    //footer part of the file
    M_footerStringStream << "\t\t</Piece>\n";
    M_footerStringStream << "\t</UnstructuredGrid>\n";
    M_footerStringStream << "</VTKFile>\n";

}


template <typename Mesh>
void
ExporterVTK<Mesh>::composeTypeDataHeader(const ExporterData::Where& where,
                                         std::stringstream& dataHeaderStringStream)
{
    Debug(8000) << "\n[ExporterVTK::composeTypeDataHeader] where = " << where << "\n";
    // std::stringstream M_dataHeaderStringStream;

    typename super::iterator_Type it;
    std::pair<typename super::iterator_Type, typename super::iterator_Type> rangeFound;

    // every dvar will produce a different data set
    // I'm just separating CELL data from POINT data
    rangeFound = this->M_whereToDataMap.equal_range(where);

    if ( rangeFound.first != rangeFound.second )
    {
        Debug(8000) << "\n[ExporterVTK::composeTypeDataHeader] found data\n";

        std::string whereString;
        switch ( where )
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
                    case ExporterData::ScalarField:
                        dataHeaderStringStream << "Scalars=\"" << it->second.variableName() << "\" ";
                        break;
                    case ExporterData::VectorField:
                        dataHeaderStringStream << "Vectors=\"" << it->second.variableName() << "\" ";
                        break;
                    case ExporterData::TensorField:
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
ExporterVTK<Mesh>::composeTypeDataFooter(const ExporterData::Where& where,
                                         std::stringstream& dataFooterStringStream)
{
    // std::stringstream dataFooterStringStream;

    typename super::iterator_Type it;
    std::pair<typename super::iterator_Type, typename super::iterator_Type> rangeFound;

    // every dvar will produce a different data set
    // I'm just separating CELL data from POINT data
    rangeFound = this->M_whereToDataMap.equal_range(where);

    if ( rangeFound.first != rangeFound.second )
    {
        Debug(8000) << "\n[ExporterVTK::composeTypeDataFooter] where = " << where << "\n";

        std::string whereString;
        switch ( where )
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
