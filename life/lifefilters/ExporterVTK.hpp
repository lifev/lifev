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
#include <life/lifecore/EncoderBase64.hpp>

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

    /*! @enum EXPORT_MODE
        The export modes currently supported are ascii and binary
     */
    enum EXPORT_MODE {
        ASCII_EXPORT = 1,
        BINARY_EXPORT = 2
    };

    /*! @enum FLOAT_PRECISION
        Currently supported are single and double precision
     */
    enum FLOAT_PRECISION {
        SINGLE_PRECISION = 1,
        DOUBLE_PRECISION = 2
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
    ~ExporterVTK();
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
    virtual void import(const Real& Tstart);

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
    void createPointsMaps( const feSpacePtr_Type & _feSpacePtr,
                           std::map<UInt, UInt>& globalToLocalPointsMap,
                           std::map<UInt, UInt>& localToGlobalPointsMap,
                           std::vector<Vector>& coordinatesOfPoints );

    /*!
       This method fills a buffer for the VTK collection (a *.pvd file).

       \param time the current time
       \param[out] vtkCollectionStringStream the stringstream object (a file buffer)
     */
    void composeVTKCollection( const std::string& variableName,
                               std::stringstream& vtkCollectionStringStream);

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
                              //const std::map<UInt, UInt>& localToGlobalPointsMap,
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
       \param localToGlobalMap a map to query the dvar object for global IDs
       \param dataArraysStringStream the stringstream object (a file buffer)
     */
    void composeDataArrayStream(const exporterData_Type& dvar,
                                const std::map<UInt,UInt>& localToGlobalMap,
                                std::stringstream& dataArraysStringStream);

    // to be checked - do not use for now
    void composeDataArrayStream(const typename exporterData_Type::WhereEnum& where,
                           std::stringstream& dataArraysStringStream);

    virtual void readScalar( ExporterData<mesh_Type>& /*dvar*/ ) {}
    virtual void readVector( ExporterData<mesh_Type>& /*dvar*/ ) {}
    void readVTUFiles( exporterData_Type& dvar );
    void readBinaryData( const std::string& line, std::vector<Real>& values, const UInt& numBits );
    void readASCIIData( const std::string& line, std::vector<Real>& values );
    //@}

    //! @name Private members
    //@{
    EXPORT_MODE M_exportMode;

    FLOAT_PRECISION M_floatPrecision;

    std::map< std::string, std::list<std::string> > M_pvtuFiles;

    UInt M_numImportProc;
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
super(),
M_exportMode(ASCII_EXPORT),
M_floatPrecision( DOUBLE_PRECISION ),
M_numImportProc( 0 )
{
}


template<typename Mesh> ExporterVTK<Mesh >::ExporterVTK(
                const GetPot& data_file,
                const std::string prefix)
                :
                super(data_file, prefix)
{
    switch( data_file("exporter/exportMode",1) )
    {
        case 1:
            M_exportMode = ASCII_EXPORT;
            break;
        case 2:
            M_exportMode = BINARY_EXPORT;
            break;

    }
    switch( data_file("exporter/floatPrecision",2) )
    {
        case 1:
            M_floatPrecision = SINGLE_PRECISION;
            break;
        case 2:
            M_floatPrecision = DOUBLE_PRECISION;
            break;

    }
    M_numImportProc = data_file("exporter/numImportProc",1);

}


// ==============
// Destructor
// ==============
template<typename Mesh>
ExporterVTK<Mesh>::~ExporterVTK()
{
    // Right before dying, the exporter will write a VTK collection in a pvd file
    // This allows ParaView to load a single file for the entire time series,
    // and each frame will be associated to the correct time value.
    // This will happen only if at least one time step was actually post-processed
    if( this->M_timeSteps.size() )
    {
        std::stringstream buffer("");

        // a unique time collection is produced by the leader process
        if(this->M_procId==0)
        {
            for (typename super::dataVectorIterator_Type iData=this->M_dataVector.begin();
                            iData != this->M_dataVector.end(); ++iData)
            {
                composeVTKCollection( iData->variableName(), buffer );

                std::ofstream vtkCollectionFile;
                vtkCollectionFile.open( ( this->M_postDir+this->M_prefix+"_" + iData->variableName() +
                                +".pvd").c_str() );
                vtkCollectionFile << buffer.str();
                vtkCollectionFile.close();

                buffer.str("");
            }
        }
    }

}

// =====================
// Public methods
// =====================

template<typename Mesh>
void ExporterVTK<Mesh>::postProcess(const Real& time)
{
    // typedef std::list< ExporterData >::iterator Iterator;

    this->computePostfix();

    std::size_t found( this->M_postfix.find( "*" ) );
    if ( found == string::npos )
    {
        if (!this->M_procId) std::cout << "  X-  ExporterVTK post-processing" << std::flush;

        LifeChrono chrono;
        chrono.start();

        this->M_timeSteps.push_back(time);

        for (typename super::dataVectorIterator_Type iData=this->M_dataVector.begin();
             iData != this->M_dataVector.end(); ++iData)
        {
            std::ofstream vtkFile;
            std::stringstream buffer("");

            // a unique PVTU file + a time collection is produced by the leader process
            if(this->M_procId==0)
            {
                composePVTUStream(*iData, buffer);

                std::string vtkPFileName( this->M_postDir+this->M_prefix+"_" + iData->variableName()+
                                          this->M_postfix+".pvtu" );
                std::ofstream vtkPFile;
                vtkPFile.open( vtkPFileName.c_str() );
                vtkPFile << buffer.str();
                vtkPFile.close();

                buffer.str("");

                this->M_pvtuFiles[iData->variableName()].push_back(vtkPFileName);

            }


            // redundant. should be done just the first time
            std::map<UInt,UInt> ltgPointsMap, gtlPointsMap;
            std::vector<Vector> pointsCoords;
            createPointsMaps( iData->feSpacePtr(), gtlPointsMap, ltgPointsMap, pointsCoords );

            composeVTUHeaderStream( gtlPointsMap.size(), buffer );
            composeVTUGeoStream( iData->feSpacePtr(), gtlPointsMap, pointsCoords, buffer );

            composeTypeDataHeaderStream(iData->where(), buffer);
            composeDataArrayStream(*iData, ltgPointsMap, buffer);
            composeTypeDataFooterStream(iData->where(), buffer);

            composeVTUFooterStream( buffer );

            // each process writes its own file
            vtkFile.open( ( this->M_postDir+this->M_prefix+"_" + iData->variableName()+
                            this->M_postfix+"."+this->M_procId+".vtu").c_str() );
            vtkFile << buffer.str();
            vtkFile.close();

            buffer.str("");

        }
        chrono.stop();
        if (!this->M_procId) std::cout << "...done in " << chrono.diff() << " s." << std::endl;
    }
}


template<typename Mesh>
void ExporterVTK<Mesh>::import(const Real& /*time*/)
{
    //this->M_timeSteps.push_back(time);

    this->computePostfix();

    assert( this->M_postfix != "*****" );

    if (!this->M_procId) std::cout << "  X-  ExporterVTK importing ..."<< std::endl;

    LifeChrono chrono;
    chrono.start();
    for (typename super::dataVectorIterator_Type iData=this->M_dataVector.begin();
                    iData != this->M_dataVector.end(); ++iData)
    {
        this->readVTUFiles(*iData);
    }
    chrono.stop();
    if (!this->M_procId) std::cout << "      done in " << chrono.diff() << " s." << std::endl;

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

    UInt vtkCellType(0);

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
            vtkCellType = VTK_QUADRATIC_TETRA;
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
            if (!this->M_procId)
                std::cout << "WARNING: the element is not yet implemented in ExporterVTK\n";
            abort();

    }

    return vtkCellType;
}


template <typename Mesh>
void
ExporterVTK<Mesh>::composeDataArrayStream(const exporterData_Type& dvar,
                                          const std::map<UInt,UInt>& localToGlobalMap,
                                          std::stringstream& dataArraysStringStream)
{
    UInt start        ( dvar.start() );
    UInt numGlobalDOF ( dvar.numDOF() );
    UInt numMyDOF     ( localToGlobalMap.size() );

    std::stringstream dataToBeEncoded; dataToBeEncoded.str("");
    std::string encodedDataString;

    std::string formatString;
    std::stringstream nComponents;
    nComponents << dvar.fieldDim();

    int32_type sizeOfFloat;
    std::string floatTypeString;
    switch( M_floatPrecision )
    {
        case SINGLE_PRECISION:
            ASSERT( sizeof(float) == 4, "\nThis piece of code assumes sizeof(float) == 4" )
            sizeOfFloat = sizeof(float);
            floatTypeString = "Float32";
            break;
        case DOUBLE_PRECISION:
            ASSERT( sizeof(Real) == 8, "\nThis piece of code assumes sizeof(float) == 4" )
            sizeOfFloat = sizeof(Real);
            floatTypeString = "Float64";
            break;
        default:
            abort();
    }

    int32_type lengthOfRawData( dvar.fieldDim()*numMyDOF*sizeOfFloat );

    switch( M_exportMode )
    {
        case ASCII_EXPORT:
            formatString = "ascii";
            break;
        case BINARY_EXPORT:
            formatString = "binary";
            dataToBeEncoded.write( reinterpret_cast<char *>( &lengthOfRawData ),
                                   sizeof(int32_type) );
            lengthOfRawData += sizeof(int32_type);
            break;
    }

    dataArraysStringStream << "\t\t\t\t<DataArray type=\"" << floatTypeString << "\" Name=\""
                    << dvar.variableName() << "\" NumberOfComponents=\""
                    << nComponents.str() << "\" format=\"" << formatString << "\">\n";

    switch( M_exportMode )
    {
        case ASCII_EXPORT:
            dataArraysStringStream.setf(std::ios_base::fixed);
            dataArraysStringStream.precision(5);
            dataArraysStringStream.width(12);

            for (UInt iDOF=0; iDOF<numMyDOF; ++iDOF)
            {
                Int id = localToGlobalMap.find(iDOF)->second;
                for (UInt iCoor=0; iCoor< dvar.fieldDim(); ++iCoor)
                {
                    // the tensor case is not ready yet
                    //                    for (UInt jcoor=0; jcoor< dvar.fieldDim() / nDimensions; ++jcoor)
                    //                    {
                    dataArraysStringStream << dvar( start + id +
                                                    iCoor * numGlobalDOF ) << " ";
                    //                                                        + jcoor * numGlobalDOF * dvar.fieldDim() ) << " ";
                    //                    }
                }
            }
            break;
        case BINARY_EXPORT:
            for (UInt iDOF=0; iDOF<numMyDOF; ++iDOF)
            {
                Int id = localToGlobalMap.find(iDOF)->second;
                for (UInt iCoor=0; iCoor< dvar.fieldDim(); ++iCoor)
                {
                    // the tensor case is not ready yet
                    //                    for (UInt jcoor=0; jcoor< dvar.fieldDim() / nDimensions; ++jcoor)
                    //                    {
                    if( M_floatPrecision == SINGLE_PRECISION )
                    {
                        float value( dvar( start + id + iCoor * numGlobalDOF ) );
                        //                                                        + jcoor * numGlobalDOF * dvar.fieldDim() ) << " ";
                        dataToBeEncoded.write( reinterpret_cast<const char *>(&value), sizeof(float) );
                        //                    }
                    }
                    else
                    {
                        Real value( dvar( start + id + iCoor * numGlobalDOF ) );
                        //                                                        + jcoor * numGlobalDOF * dvar.fieldDim() ) << " ";
                        dataToBeEncoded.write( reinterpret_cast<const char *>(&value), sizeof(Real) );
                        //                    }
                    }
                }
            }

            encodedDataString = base64_encode(reinterpret_cast<const unsigned char *>( dataToBeEncoded.str().c_str() ),
                                              lengthOfRawData );
            dataArraysStringStream << encodedDataString;

            break;
    }

    dataArraysStringStream << "\n\t\t\t\t</DataArray>\n";

    dataArraysStringStream << "\t\t\t\t<DataArray type=\"Int32\" NumberOfComponents=\"1\" "
                    << "Name=\"GlobalId\" format=\"ascii\">\n";
    for ( UInt iDOF = 0; iDOF < numMyDOF; ++iDOF )
    {
        dataArraysStringStream << localToGlobalMap.find(  iDOF )->second << " ";
    }
    dataArraysStringStream << "\n\t\t\t\t</DataArray>\n";
}


template <typename Mesh>
void
ExporterVTK<Mesh>::readVTUFiles( exporterData_Type& dvar )
{
    ASSERT( M_numImportProc, "The number of pieces to be loaded was not specified." );

    UInt numPoints, numCells;
    std::vector<Real> inputValues;
    std::vector<Real> localDOF;

    UInt start        ( dvar.start() );
    UInt numGlobalDOF ( dvar.numDOF() );

    dvar.feSpacePtr()->map().map(Repeated)->Print( std::cout );
    dvar.feSpacePtr()->map().map(Unique)->Print( std::cout );

    // Each processor will read all the files, and fill just its own component of the vectors
    for( UInt iProc = 0; iProc < M_numImportProc; ++iProc )
    {
        std::string filename( this->M_postDir + this->M_prefix + "_" + dvar.variableName() +
                              this->M_postfix + "." + iProc + ".vtu" );
        std::ifstream inputFile( filename.c_str() );

        if (!this->M_procId) std::cout << "\tfile "<< filename << std::endl;

        ASSERT(inputFile.good(), std::stringstream("There is an error while reading " +
                                                   filename).str().c_str() );

        // file parsing: line by line
        std::string line;
        size_t found;
        std::stringstream parseLine;

        while ( getline( inputFile, line ) )
        {
            // this is essentially a consistency check: the number of DOF is explicitly
            // written in the VTK files. We will check that the number of values read
            // from file matches this number
            found = line.find( "NumberOfPoints" );
            if ( found != std::string::npos )
            {
                // place the get pointer at the first " after NumberOfPoints
                found = line.find( "\"", found, 1 );
                // after the " we'll find the number: parse the substring
                parseLine.str( line.substr(found+1) );
                parseLine >> numPoints;

                // do the same for NumberOfCells
                found = line.find( "NumberOfCells" );
                found = line.find( "\"", found, 1 );
                parseLine.str( line.substr(found+1) );
                parseLine >> numCells;
            }

            // load all PointData arrays
            if ( line.find( "<PointData" ) != std::string::npos )
            {
                while ( getline( inputFile, line ) )
                {
                    if ( line.find( dvar.variableName() ) != std::string::npos )
                    {
                        inputValues.resize( dvar.fieldDim()*numPoints );

                        if ( line.find( "binary" ) != std::string::npos )
                        {
                            UInt numBitsFloat;
                            found = line.find( "Float" );
                            parseLine.str( line.substr(found+5) );
                            parseLine >> numBitsFloat;
                            getline( inputFile, line );
                            readBinaryData( line, inputValues, numBitsFloat );
                        }
                    }
                    if ( line.find( "GlobalId" ) != std::string::npos )
                    {
                        localDOF.resize( numPoints );
                        getline( inputFile, line );
                        readASCIIData( line, localDOF );
                    }
                }

                for (UInt iPoint=0; iPoint<numPoints; ++iPoint)
                {
                    Int id = localDOF[iPoint];
                    if( dvar.feSpacePtr()->map().map(Repeated)->MyGID( id ) )
                    {
                        std::cout << "\nProcessor " << this->M_procId
                                        << " will take care of (" << dvar.variableName()
                                        << ") Global ID " << id << std::endl;
                        for (UInt iCoor=0; iCoor< dvar.fieldDim(); ++iCoor)
                        {
                            dvar( start + id + iCoor * numGlobalDOF ) =
                                            inputValues[ iPoint * dvar.fieldDim() + iCoor ];
                        }
                    }
                }
            }
        }
    }
}


template <typename Mesh>
void
ExporterVTK<Mesh>::readBinaryData( const std::string& line, std::vector<Real>& values, const UInt& numBits )
{
    std::stringstream decodedData, dataToBeDecoded; dataToBeDecoded.str("");
    std::string decodedDataString;
    UInt sizeOfFloat, sizeOfVector( values.size() ), lengthOfRawData;

    switch( numBits )
    {
        case 32:
            ASSERT( sizeof(float) == 4, "\nThis piece of code assumes sizeof(float) == 4" );
            sizeOfFloat = sizeof(float);
            break;
        case 64:
            ASSERT( sizeof(Real) == 8, "\nThis piece of code assumes sizeof(Real) == 8" );
            sizeOfFloat = sizeof(Real);
            break;
        default:
            ASSERT( false, "unmanaged float type" );
            break;
    }

    lengthOfRawData = sizeOfVector*sizeOfFloat + sizeof(int32_type);

    // read from file a block of char
    //char* inputRawData = new char[ lengthOfRawData ];
    //iFile.read( inputRawData, lengthOfRawData );

    // assign the block of char to a stringstream (to convert it into a string)
    dataToBeDecoded.write( line.c_str(), line.size() );

    // perform the string decoding and store the string in a stringstream
    // (to access it "bitwise")
    decodedDataString = base64_decode( dataToBeDecoded.str() );
    decodedData.str( decodedDataString );

    //std::cout << "\nlengthOfRawData = " << lengthOfRawData << ", line.size() = " << line.size()
    //            << ", decodedDataString.size() = " << decodedDataString.size() << std::endl;
    ASSERT( lengthOfRawData == decodedDataString.size(), "unexpected line length" );


    // the first value in the string is the size of the subsequent bunch of data
    int32_type* inputInt = new int32_type;
    decodedData.read( reinterpret_cast< char *>( inputInt ),
                      sizeof(int32_type) );

    ASSERT( *inputInt - sizeOfVector, "Inconsistent size of data!" );

    switch( numBits )
    {
        case 32:
            {
                float* inputFloat = new float[sizeOfVector];
                decodedData.read( reinterpret_cast< char *>( inputFloat ),
                                  sizeOfVector*sizeOfFloat );
                std::vector<float> inputValuesTemp(sizeOfVector);
                inputValuesTemp.assign( inputFloat, inputFloat+sizeOfVector );
                values.assign( inputValuesTemp.begin(), inputValuesTemp.end() );
                //for (UInt i = 0; i < inputValues.size(); i++)
                //  inputValues[i] = inputValuesTemp[i];
                break;
            }
        case 64:
            {
                Real* inputReal = new Real[sizeOfVector];
                decodedData.read( reinterpret_cast< char *>( inputReal ),
                                  sizeOfVector*sizeOfFloat );
                values.assign( inputReal, inputReal+sizeOfVector );
                break;
            }
    }
    // output data to screen to verify/compare the results
    //for (UInt i = 0; i < values.size(); i++)
    //  std::cout << values[i] << " " << std::flush;
    //std::cout << std::endl;
}


template <typename Mesh>
void
ExporterVTK<Mesh>::readASCIIData( const std::string& line, std::vector<Real>& values )
{
    std::stringstream readData( line );

    // simply parse the line to fill the vector of values
    for (UInt i = 0; i < values.size(); i++)
      readData >> values[i];

    // output data to screen to verify/compare the results
    //for (UInt i = 0; i < values.size(); i++)
    //  std::cout << values[i] << " " << std::flush;
    //std::cout << std::endl;

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
            switch ( it->second.fieldType() )
            {
                case exporterData_Type::ScalarField:

                    dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
                    << it->second.variableName() << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";

                    for (UInt iValue=1; iValue<=it->second.size(); ++iValue) {
                        dataArraysStringStream << it->second( iValue ) << " ";
                    }
                    break;
                case exporterData_Type::VectorField:

                    dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
                    << it->second.variableName() << "\" NumberOfComponents=\""
                    << nDimensions << "\" format=\"ascii\">\n";

                    for (UInt iValue=1; iValue<=it->second.size(); ++iValue) {
                        for (UInt iCoor=0; iCoor< nDimensions; ++iCoor) {
                            dataArraysStringStream << it->second( iValue + iCoor * it->second.size() ) << " ";
                        }
                    }
                    break;
                    /*case typename exporterData_Type::TensorField:

                    dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
                            << it->second.variableName() << "\" NumberOfComponents=\""
                            << nDimensions*nDimensions << "\" format=\"ascii\">\n";

                    for (UInt i=0; i<it->second.size(); ++i){
                        for (UInt iCoor=0; iCoor< nDimensions;++iCoor){
                            for (UInt jcoor=0; jcoor< nDimensions;++jcoor){
                                dataArraysStringStream << it->second( i * nDimensions * nDimensions +
                                 iCoor * nDimensions + jcoor ) << " ";
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

    std::string floatTypeString;
    switch( M_floatPrecision )
    {
        case SINGLE_PRECISION:
            floatTypeString = "Float32";
            break;
        case DOUBLE_PRECISION:
            floatTypeString = "Float64";
            break;
        default:
            abort();
    }

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

    std::string formatString;
    std::stringstream nComponents;
    nComponents << dvar.fieldDim();
    switch( M_exportMode )
    {
        case ASCII_EXPORT:
            formatString = "ascii";
            break;
        case BINARY_EXPORT:
            formatString = "binary";
            break;
    }
    pVTUStringStream << "\t\t\t<PDataArray type=\"" << floatTypeString << "\" Name=\""
                    << dvar.variableName() << "\" NumberOfComponents=\""
                    << nComponents.str() << "\" format=\"" << formatString << "\">\n";

    pVTUStringStream << "\t\t\t</PDataArray>\n";

    pVTUStringStream << "\t\t\t<PDataArray type=\"Int32\" Name=\"GlobalId\" NumberOfComponents=\"1\" "
                    << "format=\"ascii\">\n";
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
        std::stringstream fileName( ( this->M_postDir+this->M_prefix+"_" + dvar.variableName()+
                                      this->M_postfix+"."+iProc+".vtu").c_str() );

        //footer part of the file
        pVTUStringStream << "\t\t<Piece Source=\"" << fileName.str() << "\"/>\n";
    }

    pVTUStringStream << "\t</PUnstructuredGrid>\n";
    pVTUStringStream << "</VTKFile>\n";

}


template <typename Mesh>
void ExporterVTK<Mesh>::createPointsMaps( const feSpacePtr_Type & _feSpacePtr,
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
    for (ID iVertex=0; iVertex < numVertices; ++iVertex)
    {
        globalPointId = this->M_mesh->point(iVertex).id();
        localToGlobalPointsMap.insert( std::pair<UInt,UInt>( positionInPartitionedMesh, globalPointId ) );
        globalToLocalPointsMap.insert( std::pair<UInt,UInt>( globalPointId, positionInPartitionedMesh ) );

        for (ID jCoor=0; jCoor < nDimensions; ++jCoor)
        {
            coordinatesOfPoints[jCoor][positionInPartitionedMesh] =
                            this->M_mesh->point(iVertex).coordinate(jCoor);
        }
        ++positionInPartitionedMesh;
    }
    ASSERT( positionInPartitionedMesh == numVertices, "didn't store all vertices in the maps");
    ASSERT( localToGlobalPointsMap.size() == globalToLocalPointsMap.size(),
            "problem in storing the local to global and global to local maps" );

    // Now I store the coordinates of the supplementary nodes in a temporary vector
    for ( UInt iElement = 0; iElement < numElements; ++iElement )
    {
        _feSpacePtr->fe().updateJac( this->M_mesh->element( iElement ) );
        for ( UInt iPoint = _feSpacePtr->dof().numLocalVertices();
              iPoint < _feSpacePtr->dof().numLocalDof(); ++iPoint )
        {
            _feSpacePtr->fe().coorMap( x, y, z,
                                       _feSpacePtr->fe().refFE().xi( iPoint ),
                                       _feSpacePtr->fe().refFE().eta( iPoint ),
                                       _feSpacePtr->fe().refFE().zeta( iPoint ) );

            globalPointId = _feSpacePtr->dof().localToGlobalMap( iElement, iPoint );

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
void ExporterVTK<Mesh>::composeVTKCollection( const std::string& variableName,
                                              std::stringstream& vtkCollectionStringStream )
{
    ASSERT( this->M_timeSteps.size(), "No time values to be saved in the VTK collection!");
    ASSERT( this->M_pvtuFiles[variableName].size(), "No file names to be saved in the VTK collection!");
    ASSERT( this->M_pvtuFiles[variableName].size() == this->M_timeSteps.size(),
            "The number of post-processed files does not match the number of time steps in the list!" );

    //header part of the file
    vtkCollectionStringStream << "<?xml version=\"1.0\"?>\n";
    vtkCollectionStringStream << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    vtkCollectionStringStream << "\t<Collection>\n";

    typedef std::list<Real>::const_iterator realIterator;
    realIterator iTime = this->M_timeSteps.begin();
    typedef std::list<std::string>::const_iterator stringIterator;
    stringIterator iFileName = this->M_pvtuFiles[variableName].begin();
    for ( ; iTime != this->M_timeSteps.end(); ++iTime, ++iFileName)
    {
        vtkCollectionStringStream << "\t\t<DataSet timestep=\"" << *iTime << "\" group=\"\" part=\"0\" "
                        << "file=\"" << *iFileName << "\" />\n";
    }

    vtkCollectionStringStream << "\t</Collection>\n";
    vtkCollectionStringStream << "</VTKFile>\n";

}


template <typename Mesh>
void ExporterVTK<Mesh>::composeVTUGeoStream( const feSpacePtr_Type & _feSpacePtr,
                                             const std::map<UInt, UInt>& globalToLocalPointsMap,
                                             //const std::map<UInt, UInt>& localToGlobalPointsMap,
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

    for (UInt iElement=0; iElement < numElements; ++iElement)
    {
        for ( UInt jPoint=0; jPoint < numLocalDof; ++jPoint)
        {
            // vtkFile.width(8);
            // UInt globalElementId( this->M_mesh->element(iElement).id() );
            UInt globalPointId( _feSpacePtr->dof().localToGlobalMap( iElement, jPoint ) );
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

    dataHeaderStringStream << ">\n";
}


template <typename Mesh>
void
ExporterVTK<Mesh>::composeTypeDataFooterStream(const typename exporterData_Type::WhereEnum& where,
                                         std::stringstream& dataFooterStringStream)
{
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
}

}
#endif // define EXPORTERVTK_H
