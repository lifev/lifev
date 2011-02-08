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
  @brief This file provides an interface for post-processing with VTK/Paraview
     by writing files in VTK XML format
  @date 11-2010

  ExporterVTK is inherited from Exporter.

  Usage: two steps
  - first: add the variables using addVariable
  - second: call postProcess(  );

  @author Tiziano Passerini <tiziano@mathcs.emory.edu>
  @maintainer Tiziano Passerini <tiziano@mathcs.emory.edu>
 */

#ifndef EXPORTERVTK_H
#define EXPORTERVTK_H 1

#include <life/lifefilters/Exporter.hpp>

namespace LifeV
{

/**
 * @class ExporterVTK
 * @brief ExporterVTK data exporter
 */
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

    typedef MeshType                          mesh_Type;
    typedef Exporter<MeshType>                super;
    typedef typename super::meshPtr_Type      meshPtr_Type;
    typedef typename super::feSpacePtr_Type   feSpacePtr_Type;
    typedef typename super::exporterData_Type exporterData_Type;
    typedef typename super::feTypeToDataIdMap_Type::iterator
                                              feTypeToDataIdMap_Type_Type;
    //@}

    //! @name Constructors & Destructor
    //@{

    //! Default constructor
    ExporterVTK();

    /*!
       Constructor for ExporterVTK

       \param data_file the GetPot data file where you must provide an [exporter] section with:
       "start" (start index for filenames 0 for 000, 1 for 001 etc.),
       "save" (how many time steps per postprocessing)
       "multimesh" (=true if the mesh has to be saved at each post-processing step)

       \param the prefix for the output file (ex. "test" for test.vtu)
     */
    ExporterVTK(const GetPot& data_file, const std::string prefix);

private:
    //! Copy constructor
    ExporterVTK( const ExporterVTK& example ) {}

public:
    //! Destructor
    ~ExporterVTK() {}

    //@}

    //! @name Public Methods
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
    virtual MapEpetraType mapType() const;
    //@}

private:

    //! @name Private methods
    //@{
    /*!
       \param _feSpacePtr a point to the FE space descriptor
       \return the identifier of the VTK elementary cell
     */
    UInt whichCellType( const feSpacePtr_Type & _feSpacePtr );

    /*!
       This method fills a buffer for the *.pvtu file. We will have a PVTU file for each
       ExporterData object, basically listing the names of the VTU files containing the
       "partitioned" data.

       \param dvar the ExporterData object
       \param[out] pVTUStringStream the stringstream object (a file buffer)
     */
    void composePVTUStream( const exporterData_Type& dvar, std::stringstream& pVTUStringStream );

    /*!
       This method creates the data structures needed by each processor to retrieve the
       point global IDs and coordinates.

       \param _feSpacePtr a pointer to the FE Space descriptor
       \param[out] globalToLocalPointsMap the key of this map is the global ID of the point,
              the value is the position in the local data structure
       \param[out] localToGlobalPointsMap the key of this map is the position in the local data structure,
              the value is the global ID of the point
       \param[out] coordinatesOfPoints a nDimensionsxnumPoints matrix, storing the nDimensions coordinates
              of the numPoints points
     */
    void createMaps( const feSpacePtr_Type & _feSpacePtr,
                     std::map<UInt, UInt>& globalToLocalPointsMap,
                     std::map<UInt, UInt>& localToGlobalPointsMap,
                     std::vector<Vector>& coordinatesOfPoints );

    /*!
       This method fills a buffer for the header part of a *.vtu file.

       \param numPoints we need to specify the dimension of the data set written in the VTU file
       \param[out] vtuHeaderStringStream the stringstream object (a file buffer)
     */
    void composeVTUHeaderStream( UInt numPoints,
                                 std::stringstream& vtuHeaderStringStream );

    /*!
       This method fills a buffer for the footer part of a *.vtu file.

       \param[out] vtuHeaderStringStream the stringstream object (a file buffer)
     */
    void composeVTUFooterStream( std::stringstream& vtuFooterStringStream );

    /*!
       This method writes in a buffer (according to the VTK XML format)
       1- the list of points' coordinates
       2- the connectivity of the cells
       3- the offsets of the connectivity list
       4- the type of each cell

       \param _feSpacePtr a pointer to the FE Space descriptor
       \param globalToLocalPointsMap the key of this map is the global ID of the point,
              the value is the position in the local data structure
       \param coordinatesOfPoints a nDimensionsxnumPoints matrix, storing the nDimensions coordinates
              of the numPoints points
     */
    void composeVTUGeoStream( const feSpacePtr_Type & _feSpacePtr,
                              const std::map<UInt, UInt>& globalToLocalPointsMap,
                              const std::vector<Vector>& coordinatesOfPoints,
                              std::stringstream& vtuGeoStringStream );

    /*!
       This method writes in a buffer (according to the VTK XML format) the header
       of the Data section of the file. This section may be named PointData or
       CellData, and to allow for several DataArrays to be grouped inside a single
       PointData or CellData environment, we need a separate method taking care
       of the header (footer)

       \param where Node or Cell
       \param typeDataHeaderStringStream  the stringstream object (a file buffer)
     */
    void composeTypeDataHeaderStream(const typename exporterData_Type::WhereEnum& where,
                                     std::stringstream& typeDataHeaderStringStream);

    /*!
       \sa composeTypeDataHeaderStream

       \param where Node or Cell
       \param typeDataHeaderStringStream  the stringstream object (a file buffer)
     */
    void composeTypeDataFooterStream(const typename exporterData_Type::WhereEnum& where,
                                     std::stringstream& typeDataFooterStringStream);
    /*!
       This method writes in a buffer (according to the VTK XML format) the content
       of the Data section of the file.

       \param dvar the ExporterData object
       \param localToGlobalPointsMap a map to query the dvar object for global IDs
       \param dataArraysStringStream the stringstream object (a file buffer)
     */
    void composeDataArrayStream(const exporterData_Type& dvar,
                                const std::map<UInt,UInt>& localToGlobalPointsMap,
                                std::stringstream& dataArraysStringStream);

    // to be checked - do not use for now
    void composeDataArrayStream(const typename exporterData_Type::WhereEnum& where,
                           std::stringstream& dataArraysStringStream);

    virtual void readScalar( ExporterData<mesh_Type>& /*dvar*/ ) {}
    virtual void readVector( ExporterData<mesh_Type>& /*dvar*/ ) {}
    //@}

    //! @name Private members
    //@{
    //@}

};


// ==============
// Implementation
// ==============

// ==============
// Constructors
// ==============

template<typename Mesh>
ExporterVTK<Mesh>::ExporterVTK():
super()
{
}


template<typename Mesh> ExporterVTK<Mesh >::ExporterVTK(
                const GetPot& data_file,
                const std::string prefix)
                :
                super(data_file, prefix)
{
}


// =====================
// Public methods
// =====================

template<typename Mesh>
void ExporterVTK<Mesh>::postProcess(const Real& /*time*/)
{
    // typedef std::list< ExporterData >::iterator Iterator;

    this->computePostfix();

    std::size_t found( this->M_postfix.find( "*" ) );
    if ( found == string::npos )
    {
        if (!this->M_procId) std::cout << "\t[VTK post-processing] ...        " << std::endl;

        LifeChrono chrono;
        chrono.start();
        for (typename super::dataVectorIterator_Type i=this->M_dataVector.begin(); i != this->M_dataVector.end(); ++i)
        {
            std::ofstream vtkFile;
            std::stringstream buffer("");

            // a unique PVTU file is produced by the leader process
            if(this->M_procId==0)
            {
                composePVTUStream(*i, buffer);

                std::ofstream vtkPFile;
                vtkPFile.open( ( this->M_postDir+this->M_prefix+"_" + i->variableName()+this->M_postfix+".pvtu").c_str() );
                vtkPFile << buffer.str();
                vtkPFile.close();

                buffer.str("");
            }


            // redundant. should be done just the first time
            std::map<UInt,UInt> ltgPointsMap, gtlPointsMap;
            std::vector<Vector> pointsCoords;
            createMaps( i->feSpacePtr(), gtlPointsMap, ltgPointsMap, pointsCoords );

            composeVTUHeaderStream( gtlPointsMap.size(), buffer );
            composeVTUGeoStream( i->feSpacePtr(), gtlPointsMap, pointsCoords, buffer );

            composeTypeDataHeaderStream(i->where(), buffer);
            composeDataArrayStream(*i, ltgPointsMap, buffer);
            composeTypeDataFooterStream(i->where(), buffer);

            composeVTUFooterStream( buffer );

            // each process writes its own file
            vtkFile.open( ( this->M_postDir+this->M_prefix+"_" + i->variableName()+this->M_postfix+"."+this->M_procId+".vtu").c_str() );
            vtkFile << buffer.str();
            vtkFile.close();

            buffer.str("");

            /*
        composeTypeDataHeaderStream(exporterData_Type::Node, M_pointDataHeaderStringStream);
        composeTypeDataFooterStream(exporterData_Type::Node, M_pointDataFooterStringStream);
        composeTypeDataHeaderStream(exporterData_Type::Cell, M_cellDataHeaderStringStream);
        composeTypeDataFooterStream(exporterData_Type::Cell, M_cellDataFooterStringStream);

        vtkFile << M_pointDataHeaderStringStream.str();
        composeDataArrayStream(exporterData_Type::Node, buffer);
        vtkFile << buffer.str();
        vtkFile << M_pointDataFooterStringStream.str();

        buffer.str("");

        vtkFile << M_cellDataHeaderStringStream.str();
        composeDataArrayStream(exporterData_Type::Cell, buffer);
        vtkFile << buffer.str();
        vtkFile << M_cellDataFooterStringStream.str();

            M_headerStringStream.str("");
            M_dataHeaderStringStream.str("");
            M_dataFooterStringStream.str("");
            M_footerStringStream.str("");

        M_pointDataHeaderStringStream.str("");
        M_pointDataFooterStringStream.str("");
        M_cellDataHeaderStringStream.str("");
        M_cellDataFooterStringStream.str("");
             */
        }
        chrono.stop();
        std::cout << "      done in " << chrono.diff() << " s." << std::endl;
    }
}


// ===================
// Get methods
// ===================

template<typename MeshType>
MapEpetraType ExporterVTK<MeshType>::mapType() const
{
    return Repeated;
}


// ===================
// Private methods
// ===================

template <typename Mesh>
UInt ExporterVTK<Mesh>::whichCellType( const feSpacePtr_Type & _feSpacePtr )
{
    ASSERT( _feSpacePtr.get(), "\nA pointer to a valid FE object is required!");

    UInt vtkCellType;

    switch ( _feSpacePtr->fe().refFE().type() )
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


template <typename Mesh>
void
ExporterVTK<Mesh>::composeDataArrayStream(const exporterData_Type& dvar,
                                          const std::map<UInt,UInt>& localToGlobalPointsMap,
                                          std::stringstream& dataArraysStringStream)
{
    dataArraysStringStream.setf(std::ios_base::fixed);
    dataArraysStringStream.precision(5);
    dataArraysStringStream.width(12);

    UInt start   = dvar.start();
    UInt numDOF  = dvar.numDOF();
    UInt numPoints( localToGlobalPointsMap.size() );

    switch ( dvar.fieldType() )
    {
        case exporterData_Type::ScalarField:

            dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
            << dvar.variableName() << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";

            for (UInt i=0; i<numPoints; ++i) {
                Int id = localToGlobalPointsMap.find(i)->second;
                dataArraysStringStream << dvar( start + id ) << " ";
            }
            break;
        case exporterData_Type::VectorField:

            dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
            << dvar.variableName() << "\" NumberOfComponents=\""
            << nDimensions << "\" format=\"ascii\">\n";

            for (UInt i=0; i<numPoints; ++i) {
                for (UInt icoor=0; icoor< dvar.typeDim(); ++icoor) {
                    Int id = localToGlobalPointsMap.find(i)->second;
                    dataArraysStringStream << dvar( start + id + icoor * numDOF ) << " ";
                }
            }
            break;
            /*case typename exporterData_Type::TensorField:

                    dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
                            << dvar.variableName() << "\" NumberOfComponents=\""
                            << nDimensions*nDimensions << "\" format=\"ascii\">\n";

                    for (UInt i=0; i<dvar.size(); ++i){
                        for (UInt icoor=0; icoor< nDimensions;++icoor){
                            for (UInt jcoor=0; jcoor< nDimensions;++jcoor){
                                dataArraysStringStream << it->second( i * nDimensions * nDimensions + icoor * nDimensions + jcoor ) << " ";
                            }
                        }
                    }
                    break;*/
    }
    dataArraysStringStream << "\n\t\t\t\t</DataArray>\n";
    // return dataArraysStringStream;
}


/*
    preliminary attempt at managing simultaneously all data associated to Nodes
    as opposed to data associated to Cells.

    Untested for now.
 */
template <typename Mesh>
void
ExporterVTK<Mesh>::composeDataArrayStream(const typename exporterData_Type::WhereEnum& where,
                                     std::stringstream& dataArraysStringStream)
{


    typename super::iterator_Type it;
    std::pair<typename super::iterator_Type, typename super::iterator_Type> rangeFound;

    rangeFound = this->M_whereToDataMap.equal_range(where);

    if ( rangeFound.first != rangeFound.second )
    {
        Debug(8000) << "\n[ExporterVTK::composeDataArrayStream] found data\n";

        dataArraysStringStream.setf(std::ios_base::fixed);
        dataArraysStringStream.precision(5);
        dataArraysStringStream.width(12);

        for (it = rangeFound.first; it != rangeFound.second; ++it)
        {
            switch ( it->second.type() )
            {
                case exporterData_Type::ScalarField:

                    dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
                    << it->second.variableName() << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";

                    for (UInt i=1; i<=it->second.size(); ++i) {
                        dataArraysStringStream << it->second( i ) << " ";
                    }
                    break;
                case exporterData_Type::VectorField:

                    dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
                    << it->second.variableName() << "\" NumberOfComponents=\""
                    << nDimensions << "\" format=\"ascii\">\n";

                    for (UInt i=1; i<=it->second.size(); ++i) {
                        for (UInt icoor=0; icoor< nDimensions; ++icoor) {
                            dataArraysStringStream << it->second( i + icoor * it->second.size() ) << " ";
                        }
                    }
                    break;
                    /*case typename exporterData_Type::TensorField:

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
void ExporterVTK<Mesh>::composePVTUStream(const exporterData_Type& dvar,
                                          std::stringstream& pVTUStringStream)
{
    // ASSERT( vtkFile, "Error: Output file cannot be open" );

    //header part of the file
    pVTUStringStream << "<?xml version=\"1.0\"?>\n";
    pVTUStringStream << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pVTUStringStream << "\t<PUnstructuredGrid GhostLevel=\"0\">\n";

    pVTUStringStream << "\t\t<PPoints>\n";
    pVTUStringStream << "\t\t\t<PDataArray type=\"Float32\" NumberOfComponents=\"" << nDimensions
                    << "\" format=\"ascii\"/>\n";
    pVTUStringStream << "\t\t</PPoints>\n";

    // connectivity
    // cells_size = nldof*(nV+1);
    pVTUStringStream << "\t\t<PCells>\n";
    pVTUStringStream << "\t\t\t<PDataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"/>\n";

    pVTUStringStream << "\t\t\t<PDataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"/>\n";

    pVTUStringStream << "\t\t\t<PDataArray type=\"Int32\" Name=\"types\" format=\"ascii\"/>\n";

    pVTUStringStream << "\t\t</PCells>\n";

    std::string whereString;
    switch ( dvar.where() )
    {
        case exporterData_Type::Node:
            whereString = "PPointData";
            break;
        case exporterData_Type::Cell:
            whereString = "PCellData";
            break;
    }
    pVTUStringStream << "\t\t<" << whereString << " ";

    pVTUStringStream << ">\n";

    switch ( dvar.fieldType() )
    {
        case exporterData_Type::ScalarField:

            pVTUStringStream << "\t\t\t<PDataArray type=\"Float32\" Name=\""
            << dvar.variableName() << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";

            break;
        case exporterData_Type::VectorField:

            pVTUStringStream << "\t\t\t<PDataArray type=\"Float32\" Name=\""
            << dvar.variableName() << "\" NumberOfComponents=\""
            << nDimensions << "\" format=\"ascii\">\n";

            break;
    }
    pVTUStringStream << "\t\t\t</PDataArray>\n";

    switch ( dvar.where() )
    {
        case exporterData_Type::Node:
            whereString = "PPointData";
            break;
        case exporterData_Type::Cell:
            whereString = "PCellData";
            break;
    }
    pVTUStringStream << "\t\t</" << whereString << ">\n";

    for( int iProc = 0; iProc < dvar.feSpacePtr()->map().comm().NumProc(); ++iProc )
    {
        std::stringstream fileName( ( this->M_postDir+this->M_prefix+"_" + dvar.variableName()+this->M_postfix+"."+iProc+".vtu").c_str() );

        //footer part of the file
        pVTUStringStream << "\t\t<Piece Source=\"" << fileName.str() << "\"/>\n";
    }

    pVTUStringStream << "\t</PUnstructuredGrid>\n";
    pVTUStringStream << "</VTKFile>\n";

}


template <typename Mesh>
void ExporterVTK<Mesh>::createMaps( const feSpacePtr_Type & _feSpacePtr,
                                    std::map<UInt, UInt>& globalToLocalPointsMap,
                                    std::map<UInt, UInt>& localToGlobalPointsMap,
                                    std::vector<Vector>& coordinatesOfPoints )
{
    ASSERT( this->M_mesh.get(), "\nA pointer to a valid mesh object is required!");
    ASSERT( _feSpacePtr.get(), "\nA pointer to a valid FESpace object is required!");

    ID numVertices = this->M_mesh->numVertices();
    ID numElements = this->M_mesh->numElements();

    // careful: the vertex map in the mesh is repeated. to know how many non vertex dofs I have
    // in the partitioned mesh I need to look at repeated maps
    UInt numPoints = _feSpacePtr->map().map(Repeated)->NumMyElements() / _feSpacePtr->fieldDim();

    coordinatesOfPoints.resize( nDimensions, ZeroVector(numPoints) );

    Real x, y, z;

    // The global ID of the considered DOF
    UInt globalPointId(0);
    // The local ID of the considered DOF
    UInt positionInPartitionedMesh(0);
    // Helper iterator
    std::pair< std::map<UInt,UInt>::iterator, bool > returnType;

    // Vertex based Dof: the coordinates are available from the Point List
    for (ID iVertex=1; iVertex <= numVertices; ++iVertex)
    {
        globalPointId = this->M_mesh->point(iVertex).id();
        localToGlobalPointsMap.insert( std::pair<UInt,UInt>( positionInPartitionedMesh, globalPointId ) );
        globalToLocalPointsMap.insert( std::pair<UInt,UInt>( globalPointId, positionInPartitionedMesh ) );

        for (ID jCoor=0; jCoor < nDimensions; ++jCoor)
        {
            coordinatesOfPoints[jCoor][positionInPartitionedMesh] =
                            this->M_mesh->point(iVertex).coordinate(jCoor+1);
        }
        ++positionInPartitionedMesh;
    }
    ASSERT( positionInPartitionedMesh == numVertices, "didn't store all vertices in the maps");
    ASSERT( localToGlobalPointsMap.size() == globalToLocalPointsMap.size(),
            "problem in storing the local to global and global to local maps" );

    // Now I store the coordinates of the supplementary nodes in a temporary vector
    for ( UInt iElement = 1; iElement <= numElements; ++iElement )
    {
        _feSpacePtr->fe().updateJac( this->M_mesh->element( iElement ) );
        for ( UInt iPoint = _feSpacePtr->dof().numLocalVertices(); iPoint < _feSpacePtr->dof().numLocalDof(); ++iPoint )
        {
            _feSpacePtr->fe().coorMap( x, y, z,
                                       _feSpacePtr->fe().refFE().xi( iPoint ),
                                       _feSpacePtr->fe().refFE().eta( iPoint ),
                                       _feSpacePtr->fe().refFE().zeta( iPoint ) );

            globalPointId = _feSpacePtr->dof().localToGlobal( iElement, iPoint+1 );

            returnType = globalToLocalPointsMap.insert( std::pair<UInt,UInt>( globalPointId,
                                                                             positionInPartitionedMesh ) );

            // by looping over mesh elements I get repetition in the list of points
            // update the maps only when finding a new point
            if( returnType.second )
            {
                localToGlobalPointsMap.insert( std::pair<UInt,UInt>( positionInPartitionedMesh, globalPointId ) );
                coordinatesOfPoints[0][positionInPartitionedMesh] = x;
                coordinatesOfPoints[1][positionInPartitionedMesh] = y;
                coordinatesOfPoints[2][positionInPartitionedMesh] = z;
                ++positionInPartitionedMesh;
            }

        }
    }
    ASSERT( positionInPartitionedMesh == numPoints, "didn't store all points in the maps" );
    ASSERT( localToGlobalPointsMap.size() == globalToLocalPointsMap.size(),
            "problem in storing the local to global and global to local maps" );
    ASSERT( localToGlobalPointsMap.size() == numPoints,
            "problem in storing the local to global and global to local maps" );

}


template <typename Mesh>
void ExporterVTK<Mesh>::composeVTUGeoStream( const feSpacePtr_Type & _feSpacePtr,
                                             const std::map<UInt, UInt>& globalToLocalPointsMap,
                                             const std::vector<Vector>& coordinatesOfPoints,
                                             std::stringstream& vtuGeoStringStream )
{
    ASSERT( this->M_mesh.get(), "\nA pointer to a valid mesh object is required!");
    ASSERT( _feSpacePtr.get(), "\nA pointer to a valid FESpace object is required!");

    Debug(8000) << "\n[ExporterVTK::composeVTUHeaderStream]\n";

    UInt numPoints = globalToLocalPointsMap.size();
    UInt numElements = this->M_mesh->numElements();
    UInt numLocalDof = _feSpacePtr->dof().numLocalDof();

    vtuGeoStringStream << "\t\t\t<Points>\n";
    vtuGeoStringStream << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"" << nDimensions
                    << "\" format=\"ascii\">\n";

    for ( UInt iPoint = 0; iPoint < numPoints; ++iPoint )
    {
        for ( UInt iCoor = 0; iCoor < nDimensions; ++iCoor )
            vtuGeoStringStream << coordinatesOfPoints[iCoor][ iPoint ] << " ";
    }
    vtuGeoStringStream << "\n\t\t\t\t</DataArray>\n";
    vtuGeoStringStream << "\t\t\t</Points>\n";

    // connectivity
    // cells_size = nldof*(nV+1);
    vtuGeoStringStream << "\t\t\t<Cells>\n";
    vtuGeoStringStream << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    for (UInt iElement=1; iElement <= numElements; ++iElement)
    {
        for ( UInt jPoint=1; jPoint <= numLocalDof; ++jPoint)
        {
            // vtkFile.width(8);
            // UInt globalElementId( this->M_mesh->element(iElement).id() );
            UInt globalPointId( _feSpacePtr->dof().localToGlobal( iElement, jPoint ) );
            ASSERT( globalToLocalPointsMap.find( globalPointId )!=globalToLocalPointsMap.end(),
                    "didn't find a local ID for global point" );
            UInt localId( globalToLocalPointsMap.find( globalPointId )->second );
            vtuGeoStringStream << localId << " ";
        }
    }
    vtuGeoStringStream << "\n\t\t\t\t</DataArray>\n";

    vtuGeoStringStream << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";

    for (UInt iElement=1; iElement <= numElements; ++iElement)
    {
        vtuGeoStringStream << iElement * numLocalDof << " ";
    }

    vtuGeoStringStream << "\n\t\t\t\t</DataArray>\n";

    vtuGeoStringStream << "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";

    // definition of cells types
    UInt vtkCellType = whichCellType(_feSpacePtr);
    for (UInt iElement=1; iElement <= numElements; ++iElement)
    {
        vtuGeoStringStream << vtkCellType << " ";
    }

    vtuGeoStringStream << "\n\t\t\t\t</DataArray>\n";
    vtuGeoStringStream << "\t\t\t</Cells>\n";

}


template <typename Mesh>
void ExporterVTK<Mesh>::composeVTUHeaderStream( UInt numPoints,
                                                std::stringstream& vtuHeaderStringStream )
{

    //header part of the file
    vtuHeaderStringStream << "<?xml version=\"1.0\"?>\n";
    vtuHeaderStringStream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    vtuHeaderStringStream << "\t<UnstructuredGrid>\n";
    vtuHeaderStringStream << "\t\t<Piece NumberOfPoints=\"" << numPoints << "\""
                          << " NumberOfCells=\"" << this->M_mesh->numElements() << "\">\n";


}


template <typename Mesh>
void ExporterVTK<Mesh>::composeVTUFooterStream( std::stringstream& vtuFooterStringStream )
{
    Debug(8000) << "\n[ExporterVTK::composeFooter]\n";

    //footer part of the file
    vtuFooterStringStream << "\t\t</Piece>\n";
    vtuFooterStringStream << "\t</UnstructuredGrid>\n";
    vtuFooterStringStream << "</VTKFile>\n";


}


template <typename Mesh>
void
ExporterVTK<Mesh>::composeTypeDataHeaderStream(const typename exporterData_Type::WhereEnum& where,
                                         std::stringstream& dataHeaderStringStream)
                                         {
    Debug(8000) << "\n[ExporterVTK::composeTypeDataHeaderStream] where = " << where << "\n";
    // std::stringstream M_dataHeaderStringStream;
    /*
    typename super::iterator_Type it;
    std::pair<typename super::iterator_Type, typename super::iterator_Type> rangeFound;

    // every dvar will produce a different data set
    // I'm just separating CELL data from POINT data
    rangeFound = this->M_whereToDataMap.equal_range(where);

    if ( rangeFound.first != rangeFound.second )
    {
        Debug(8000) << "\n[ExporterVTK::composeTypeDataHeaderStream] found data\n";
     */
    std::string whereString;
    switch ( where )
    {
        case exporterData_Type::Node:
            whereString = "PointData";
            break;
        case exporterData_Type::Cell:
            whereString = "CellData";
            break;
    }
    dataHeaderStringStream << "\t\t\t<" << whereString << " ";
    /*
                for (it = rangeFound.first; it != rangeFound.second; ++it)
                {
                    switch( it->second.type() )
                    {
                    case typename exporterData_Type::ScalarField:
                        dataHeaderStringStream << "Scalars=\"" << it->second.variableName() << "\" ";
                        break;
                    case typename exporterData_Type::VectorField:
                        dataHeaderStringStream << "Vectors=\"" << it->second.variableName() << "\" ";
                        break;
                    case typename exporterData_Type::TensorField:
                        M_pointDataHeaderStringStream << "Tensors=\"" << it->second.variableName() << "\" ";
                        break;
                    }
                }
     */
    dataHeaderStringStream << ">\n";
    // }

    // return M_dataHeaderStringStream;
                                         }


template <typename Mesh>
void
ExporterVTK<Mesh>::composeTypeDataFooterStream(const typename exporterData_Type::WhereEnum& where,
                                         std::stringstream& dataFooterStringStream)
                                         {
    // std::stringstream dataFooterStringStream;
    /*
    typename super::iterator_Type it;
    std::pair<typename super::iterator_Type, typename super::iterator_Type> rangeFound;

    // every dvar will produce a different data set
    // I'm just separating CELL data from POINT data
    rangeFound = this->M_whereToDataMap.equal_range(where);

    if ( rangeFound.first != rangeFound.second )
    {
        Debug(8000) << "\n[ExporterVTK::composeTypeDataFooterStream] where = " << where << "\n";
     */
    std::string whereString;
    switch ( where )
    {
        case exporterData_Type::Node:
            whereString = "PointData";
            break;
        case exporterData_Type::Cell:
            whereString = "CellData";
            break;
    }
    dataFooterStringStream << "\t\t\t</" << whereString << ">\n";

    //    }
    // return dataFooterStringStream;
}

}
#endif // define EXPORTERVTK_H
