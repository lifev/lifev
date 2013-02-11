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

#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/util/EncoderBase64.hpp>

namespace LifeV
{

/**
 * @class ExporterVTK
 * @brief ExporterVTK data exporter
 */
template<typename MeshType>
class ExporterVTK : public Exporter<MeshType>
{

public:

    //! @name Public Types
    //@{

    /*! @enum VTK_CELL
        A list of the cell types admitted in VTK
     */
    enum VTK_CELL
    {
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
    enum EXPORT_MODE
    {
        ASCII_EXPORT = 1,
        BINARY_EXPORT = 2
    };

    /*! @enum FLOAT_PRECISION
        Currently supported are single and double precision
     */
    enum FLOAT_PRECISION
    {
        SINGLE_PRECISION = 1,
        DOUBLE_PRECISION = 2
    };

    typedef MeshType                          mesh_Type;
    typedef Exporter<mesh_Type>               super;
    typedef typename super::meshPtr_Type      meshPtr_Type;
    typedef typename super::feSpacePtr_Type   feSpacePtr_Type;
    typedef typename super::exporterData_Type exporterData_Type;
    typedef typename super::feTypeToDataIdMap_Type::iterator
    feTypeToDataIdMap_Type_Type;
    typedef const typename exporterData_Type::WhereEnum& where_Type;
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
    ExporterVTK (const GetPot& data_file, const std::string prefix);

private:
    //! Copy constructor
    ExporterVTK ( const ExporterVTK& example );

public:
    //! Destructor
    virtual ~ExporterVTK();
    //@}

    //! @name Public Methods
    //@{

    /* *
       Post-process the variables added to the list

       \param time the solver time
     */
    virtual void postProcess (const Real& time);

    //! Import data from previous simulations at a certain time
    /*!
      @param Time the time of the data to be imported

      Not yet implemented for ExporterVTK
    */
    virtual UInt importFromTime ( const Real& /*Time*/ )
    {
        ERROR_MSG ( "ExporterVTK::importFromTime has not yet been implemented.");
        return -1;
    }

    //! Import data from previous simulations and rebuild the internal time counters
    /*!
      @param importTime the time of the snapshot to be imported
      @param dt the time step, is used to rebuild the history up to now

      Not yet implemented for ExporterVTK
    */
    virtual void import (const Real& /*Tstart*/, const Real& /*dt*/)
    {
        ERROR_MSG ( "ExporterVTK::importFromTime has not yet been implemented.");
    }

    //! Import data from previous simulations
    /*!
      @param importTime the time of the snapshot to be imported
    */
    virtual void import (const Real& Tstart);

    //! temporary: the method should work form the Exporter class
    void exportPID ( boost::shared_ptr<MeshType> /*mesh*/, boost::shared_ptr<Epetra_Comm> /*comm*/ )
    {
        std::cerr << "  X-  exportPID is not working with VTK (missing P0 element support)" << std::endl;
    }

    //! Set data from file.
    /*!
       @param dataFile data file.
       @param section section in the data file.
     */
    virtual void setDataFromGetPot ( const GetPot& dataFile, const std::string& section = "exporter" );

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
    UInt whichCellType ( const feSpacePtr_Type& _feSpacePtr );

    /*!
       This method fills a buffer for the *.pvtu file. We will have a PVTU file for each
       ExporterData object, basically listing the names of the VTU files containing the
       "partitioned" data.

       \param dvar the ExporterData object
       \param[out] pVTUStringStream the stringstream object (a file buffer)
     */
    void composePVTUStream ( const exporterData_Type& dvar, std::stringstream& pVTUStringStream );

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
    void createPointsMaps ( const feSpacePtr_Type& _feSpacePtr,
                            std::map<UInt, UInt>& globalToLocalPointsMap,
                            std::map<UInt, UInt>& localToGlobalPointsMap,
                            std::vector<Vector>& coordinatesOfPoints );

    /*!
       This method fills a buffer for the VTK collection (a *.pvd file).

       \param time the current time
       \param[out] vtkCollectionStringStream the stringstream object (a file buffer)
     */
    void composeVTKCollection ( const std::string& variableName,
                                std::stringstream& vtkCollectionStringStream);

    /*!
       This method fills a buffer for the header part of a *.vtu file.

       \param numPoints we need to specify the dimension of the data set written in the VTU file
       \param[out] vtuHeaderStringStream the stringstream object (a file buffer)
     */
    void composeVTUHeaderStream ( UInt numPoints,
                                  std::stringstream& vtuHeaderStringStream );

    /*!
       This method fills a buffer for the footer part of a *.vtu file.

       \param[out] vtuHeaderStringStream the stringstream object (a file buffer)
     */
    void composeVTUFooterStream ( std::stringstream& vtuFooterStringStream );

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
    void composeVTUGeoStream ( const feSpacePtr_Type& _feSpacePtr,
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

    void composeTypeDataHeaderStream (where_Type where,
                                      std::stringstream& typeDataHeaderStringStream);

    /*!
       \sa composeTypeDataHeaderStream

       \param where Node or Cell
       \param typeDataHeaderStringStream  the stringstream object (a file buffer)
     */
    void composeTypeDataFooterStream (where_Type where,
                                      std::stringstream& typeDataFooterStringStream);
    /*!
       This method writes in a buffer (according to the VTK XML format) the content
       of the Data section of the file.

       \param dvar the ExporterData object
       \param localToGlobalMap a map to query the dvar object for global IDs
       \param dataArraysStringStream the stringstream object (a file buffer)
     */
    void composeDataArrayStream (const exporterData_Type& dvar,
                                 const std::map<UInt, UInt>& localToGlobalMap,
                                 std::stringstream& dataArraysStringStream);

    // to be checked - do not use for now
    void composeDataArrayStream (where_Type where,
                                 std::stringstream& dataArraysStringStream);

    //! The scalar reader (specialization of the parent class method)
    /*!
      @param dvar the ExporterData object
    */
    virtual void readScalar ( ExporterData<mesh_Type>& /*dvar*/ )
    {
        ERROR_MSG ( "ExporterVTK::readScalar has not yet been implemented.");
    }
    //! The vector reader (specialization of the parent class method)
    /*!
      @param dvar the ExporterData object
    */
    virtual void readVector ( ExporterData<mesh_Type>& /*dvar*/ )
    {
        ERROR_MSG ( "ExporterVTK::readVector has not yet been implemented.");
    }
    //! The reader for VTU files
    /*!
      @param dvar the ExporterData object
    */
    void readVTUFiles ( exporterData_Type& dvar );
    //! A routine for loading values stored in binary format in a VTU file
    /*!
      @param line a line read from file
      @param values the list of values extracted from the line
      @param numBits the size of data to be read
    */
    void readBinaryData ( const std::string& line, std::vector<Real>& values, const UInt& numBits );
    //! A routine for loading values stored in ASCII format in a VTU file
    /*!
      @param line a line read from file
      @param values the list of values extracted from the line
    */
    void readASCIIData ( const std::string& line, std::vector<Real>& values );
    //@}

    //! @name Private members
    //@{
    EXPORT_MODE M_exportMode;

    FLOAT_PRECISION M_floatPrecision;

    std::map< std::string, std::list<std::string> > M_pvtuFiles;
    //@}

};


// ==============
// Implementation
// ==============

// ==============
// Constructors
// ==============

template<typename MeshType>
ExporterVTK<MeshType>::ExporterVTK() :
    super(),
    M_exportMode (ASCII_EXPORT),
    M_floatPrecision ( DOUBLE_PRECISION )
{
}


template<typename MeshType> ExporterVTK<MeshType >::ExporterVTK (
    const GetPot& data_file,
    const std::string prefix)
    :
    super (data_file, prefix)
{
    this->setDataFromGetPot (data_file);
}


template<typename MeshType>
void ExporterVTK<MeshType >::setDataFromGetPot (
    const GetPot& data_file,
    const std::string& section )
{
    super::setDataFromGetPot ( data_file, section );

    switch ( data_file ( (section + "/exportMode").c_str(), 1) )
    {
        case 1:
            M_exportMode = ASCII_EXPORT;
            break;
        case 2:
            M_exportMode = BINARY_EXPORT;
            break;
        default:
            ERROR_MSG ( "Unsupported export mode!" );
            break;
    }

    switch ( data_file ( (section + "/floatPrecision").c_str(), 2) )
    {
        case 1:
            M_floatPrecision = SINGLE_PRECISION;
            break;
        case 2:
            M_floatPrecision = DOUBLE_PRECISION;
            break;
        default:
            ERROR_MSG ( "Unsupported float precision requirement!" );
            break;
    }

}
// ==============
// Destructor
// ==============
template<typename MeshType>
ExporterVTK<MeshType>::~ExporterVTK()
{
    // Right before dying, the exporter will write a VTK collection in a pvd file
    // This allows ParaView to load a single file for the entire time series,
    // and each frame will be associated to the correct time value.
    // This will happen only if at least one time step was actually post-processed
    if ( this->M_timeSteps.size() )
    {
        std::stringstream buffer ("");

        // a unique time collection is produced by the leader process
        if (this->M_procId == 0)
        {
            for (typename super::dataVectorIterator_Type iData = this->M_dataVector.begin();
                    iData != this->M_dataVector.end(); ++iData)
            {
                composeVTKCollection ( iData->variableName(), buffer );

                std::string filename ( this->M_postDir + this->M_prefix + "_" + iData->variableName() + ".pvd" );
                std::ofstream vtkCollectionFile;
                vtkCollectionFile.open ( filename.c_str() );
                ASSERT (vtkCollectionFile.is_open(), "There is an error while opening " + filename );
                ASSERT (vtkCollectionFile.good(), "There is an error while writing to " + filename );
                vtkCollectionFile << buffer.str();
                vtkCollectionFile.close();

                buffer.str ("");
            }
        }
    }

}

// =====================
// Public methods
// =====================

template<typename MeshType>
void ExporterVTK<MeshType>::postProcess (const Real& time)
{
    // typedef std::list< ExporterData >::iterator Iterator;

    this->computePostfix();

    std::size_t found ( this->M_postfix.find ( "*" ) );
    if ( found == string::npos )
    {
        if (this->M_procId == 0)
        {
            std::cout << "  X-  ExporterVTK post-processing" << std::flush;
        }

        LifeChrono chrono;
        chrono.start();

        this->M_timeSteps.push_back (time);

        for (typename super::dataVectorIterator_Type iData = this->M_dataVector.begin();
                iData != this->M_dataVector.end(); ++iData)
        {
            std::ofstream vtkFile;
            std::stringstream buffer ("");

            // a unique PVTU file + a time collection is produced by the leader process
            if (this->M_procId == 0)
            {
                composePVTUStream (*iData, buffer);

                std::string vtkPFileName ( this->M_prefix + "_" + iData->variableName() +
                                           this->M_postfix + ".pvtu" );
                std::string vtkPFileNameWithDir (this->M_postDir + vtkPFileName);

                std::ofstream vtkPFile;
                vtkPFile.open ( vtkPFileNameWithDir.c_str() );
                ASSERT (vtkPFile.is_open(), "There is an error while opening " + vtkPFileName );
                ASSERT (vtkPFile.good(), "There is an error while writing to " + vtkPFileName );
                vtkPFile << buffer.str();
                vtkPFile.close();

                buffer.str ("");

                this->M_pvtuFiles[iData->variableName()].push_back (vtkPFileName);

            }


            // redundant. should be done just the first time
            std::map<UInt, UInt> ltgPointsMap, gtlPointsMap;
            std::vector<Vector> pointsCoords;
            createPointsMaps ( iData->feSpacePtr(), gtlPointsMap, ltgPointsMap, pointsCoords );

            composeVTUHeaderStream ( gtlPointsMap.size(), buffer );
            composeVTUGeoStream ( iData->feSpacePtr(), gtlPointsMap, pointsCoords, buffer );

            composeTypeDataHeaderStream (iData->where(), buffer);
            composeDataArrayStream (*iData, ltgPointsMap, buffer);
            composeTypeDataFooterStream (iData->where(), buffer);

            composeVTUFooterStream ( buffer );

            // each process writes its own file
            std::string filename ( this->M_postDir + this->M_prefix + "_" + iData->variableName() +
                                   this->M_postfix + "." + this->M_procId + ".vtu" );
            vtkFile.open ( filename.c_str() );
            ASSERT (vtkFile.is_open(), "There is an error while opening " + filename );
            ASSERT (vtkFile.good(), "There is an error while writing to " + filename );
            vtkFile << buffer.str();
            vtkFile.close();

            buffer.str ("");

        }
        chrono.stop();
        if (!this->M_procId)
        {
            std::cout << "...done in " << chrono.diff() << " s." << std::endl;
        }
    }
}


template<typename MeshType>
void ExporterVTK<MeshType>::import (const Real& /*time*/)
{
    this->computePostfix();

    assert ( this->M_postfix != "*****" );

    if (this->M_procId == 0)
    {
        std::cout << "  X-  ExporterVTK importing ..." << std::endl;
    }

    LifeChrono chrono;
    chrono.start();
    for (typename super::dataVectorIterator_Type iData = this->M_dataVector.begin();
            iData != this->M_dataVector.end(); ++iData)
    {
        this->readVTUFiles (*iData);
    }
    chrono.stop();
    if (this->M_procId == 0)
    {
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

template <typename MeshType>
UInt ExporterVTK<MeshType>::whichCellType ( const feSpacePtr_Type& _feSpacePtr )
{
    ASSERT ( _feSpacePtr.get(), "\nA pointer to a valid FE object is required!");

    UInt vtkCellType (0);

    switch ( _feSpacePtr->fe().refFE().type() )
    {
        case FE_P1_2D:
            // case FE_P1bubble_2D:
            vtkCellType = VTK_TRIANGLE;
            break;
        case FE_P2_2D:
            vtkCellType = VTK_QUADRATIC_TRIANGLE;
            break;
        case FE_P0_3D:
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
            ERROR_MSG ( "WARNING: the element is not yet implemented in ExporterVTK\n" )
            break;
    }

    return vtkCellType;
}


template <typename MeshType>
void
ExporterVTK<MeshType>::composeDataArrayStream (const exporterData_Type& dvar,
                                               const std::map<UInt, UInt>& localToGlobalMap,
                                               std::stringstream& dataArraysStringStream)
{
    const UInt start        ( dvar.start() );
    const UInt numGlobalDOF ( dvar.numDOF() );
    const UInt numMyDOF     ( localToGlobalMap.size() );

    std::stringstream dataToBeEncoded;
    dataToBeEncoded.str ("");
    std::string encodedDataString;

    std::string formatString;
    std::stringstream nComponents;
    nComponents << dvar.fieldDim();

    int32_type sizeOfFloat = 4;
    std::string floatTypeString;
    switch ( M_floatPrecision )
    {
        case SINGLE_PRECISION:
            ASSERT ( sizeof (float) == 4, "\nThis piece of code assumes sizeof(float) == 4" )
            sizeOfFloat = sizeof (float);
            floatTypeString = "Float32";
            break;
        case DOUBLE_PRECISION:
            ASSERT ( sizeof (Real) == 8, "\nThis piece of code assumes sizeof(float) == 4" )
            sizeOfFloat = sizeof (Real);
            floatTypeString = "Float64";
            break;
        default:
            ERROR_MSG ( "WARNING: this float precision cannot be handled in ExporterVTK\n" )
            break;
    }

    int32_type lengthOfRawData ( dvar.fieldDim() *numMyDOF * sizeOfFloat );

    switch ( M_exportMode )
    {
        case ASCII_EXPORT:
            formatString = "ascii";
            break;
        case BINARY_EXPORT:
            formatString = "binary";
            dataToBeEncoded.write ( reinterpret_cast<char*> ( &lengthOfRawData ),
                                    sizeof (int32_type) );
            lengthOfRawData += sizeof (int32_type);
            break;
        default:
            ERROR_MSG ( "WARNING: this export mode cannot be handled in ExporterVTK\n" )
            break;
    }

    dataArraysStringStream << "\t\t\t\t<DataArray type=\"" << floatTypeString << "\" Name=\""
                           << dvar.variableName() << "\" NumberOfComponents=\""
                           << nComponents.str() << "\" format=\"" << formatString << "\">\n";

    switch ( M_exportMode )
    {
        case ASCII_EXPORT:
            dataArraysStringStream.setf (std::ios_base::fixed);
            dataArraysStringStream.precision (5);
            dataArraysStringStream.width (12);

            for (UInt iDOF = 0; iDOF < numMyDOF; ++iDOF)
            {
                const Int id = localToGlobalMap.find (iDOF)->second;
                for (UInt iCoor = 0; iCoor < dvar.fieldDim(); ++iCoor)
                {
                    dataArraysStringStream << dvar ( start + id +
                                                     iCoor * numGlobalDOF ) << " ";
                }
            }
            break;
        case BINARY_EXPORT:
            for (UInt iDOF = 0; iDOF < numMyDOF; ++iDOF)
            {
                const Int id = localToGlobalMap.find (iDOF)->second;
                for (UInt iCoor = 0; iCoor < dvar.fieldDim(); ++iCoor)
                {
                    if ( M_floatPrecision == SINGLE_PRECISION )
                    {
                        const float value ( dvar ( start + id + iCoor * numGlobalDOF ) );
                        dataToBeEncoded.write ( reinterpret_cast<const char*> (&value), sizeof (float) );
                    }
                    else
                    {
                        const Real value ( dvar ( start + id + iCoor * numGlobalDOF ) );
                        dataToBeEncoded.write ( reinterpret_cast<const char*> (&value), sizeof (Real) );
                    }
                }
            }

            encodedDataString = base64_encode (reinterpret_cast<const unsigned char*> ( dataToBeEncoded.str().c_str() ),
                                               lengthOfRawData );
            dataArraysStringStream << encodedDataString;

            break;
        default:
            ERROR_MSG ( "WARNING: this export mode cannot be handled in ExporterVTK\n" )
            break;
    }

    dataArraysStringStream << "\n\t\t\t\t</DataArray>\n";

    dataArraysStringStream << "\t\t\t\t<DataArray type=\"Int32\" NumberOfComponents=\"1\" "
                           << "Name=\"GlobalId\" format=\"ascii\">\n";
    for ( UInt iDOF = 0; iDOF < numMyDOF; ++iDOF )
    {
        dataArraysStringStream << localToGlobalMap.find (  iDOF )->second << " ";
    }
    dataArraysStringStream << "\n\t\t\t\t</DataArray>\n";
}


template <typename MeshType>
void
ExporterVTK<MeshType>::readVTUFiles ( exporterData_Type& dvar )
{
    ASSERT ( this->M_numImportProc, "The number of pieces to be loaded was not specified." );

    UInt numPoints, numCells;
    std::vector<Real> inputValues;
    std::vector<Real> globalIDs;

    const UInt start        ( dvar.start() );
    const UInt numGlobalDOF ( dvar.numDOF() );

    // Each processor will read all the files, and fill just its own component of the vectors
    for ( UInt iProc = 0; iProc < this->M_numImportProc; ++iProc )
    {
        std::string filename ( this->M_postDir + this->M_prefix + "_" + dvar.variableName() +
                               this->M_postfix + "." + iProc + ".vtu" );
        std::ifstream inputFile ( filename.c_str() );

        if (this->M_procId == 0)
        {
            std::cout << "\tfile " << filename << std::endl;
        }

        ASSERT (inputFile.is_open(), "There is an error while opening " + filename );

        // file parsing: line by line
        std::string line;
        size_t found;
        std::stringstream parseLine;

        while ( inputFile.good() && getline ( inputFile, line ) )
        {
            // this is essentially a consistency check: the number of DOF is explicitly
            // written in the VTK files. We will check that the number of values read
            // from file matches this number
            found = line.find ( "NumberOfPoints" );
            if ( found != std::string::npos )
            {
                // place the get pointer at the first " after NumberOfPoints
                found = line.find ( "\"", found, 1 );
                // after the " we'll find the number: parse the substring
                parseLine.str ( line.substr (found + 1) );
                parseLine >> numPoints;

                // do the same for NumberOfCells
                found = line.find ( "NumberOfCells" );
                found = line.find ( "\"", found, 1 );
                parseLine.str ( line.substr (found + 1) );
                parseLine >> numCells;
            }

            // load all PointData arrays
            if ( line.find ( "<PointData" ) != std::string::npos )
            {
                while ( inputFile.good() && getline ( inputFile, line ) )
                {
                    if ( line.find ( dvar.variableName() ) != std::string::npos )
                    {
                        inputValues.resize ( dvar.fieldDim() *numPoints );

                        if ( line.find ( "binary" ) != std::string::npos )
                        {
                            UInt numBitsFloat;
                            found = line.find ( "Float" );
                            parseLine.str ( line.substr (found + 5) );
                            parseLine >> numBitsFloat;
                            ASSERT (inputFile.good(), "There is an error while opening " + filename );
                            getline ( inputFile, line );
                            readBinaryData ( line, inputValues, numBitsFloat );
                        }
                        else
                        {
                            ASSERT (inputFile.good(), "There is an error while opening " + filename );
                            getline ( inputFile, line );
                            readASCIIData ( line, inputValues );
                        }
                    }
                    if ( line.find ( "GlobalId" ) != std::string::npos )
                    {
                        globalIDs.resize ( numPoints );
                        ASSERT (inputFile.good(), "There is an error while opening " + filename );
                        getline ( inputFile, line );
                        readASCIIData ( line, globalIDs );
                    }
                }

                for (UInt iPoint = 0; iPoint < numPoints; ++iPoint)
                {
                    const Int id = globalIDs[iPoint];
                    if ( dvar.feSpacePtr()->map().map (Repeated)->MyGID ( id ) )
                    {
                        for (UInt iCoor = 0; iCoor < dvar.fieldDim(); ++iCoor)
                        {
                            dvar ( start + id + iCoor * numGlobalDOF ) =
                                inputValues[ iPoint * dvar.fieldDim() + iCoor ];
                        }
                    }
                }
            }
        }
        inputFile.close();
    }
}


template <typename MeshType>
void
ExporterVTK<MeshType>::readBinaryData ( const std::string& line, std::vector<Real>& values, const UInt& numBits )
{
    std::stringstream decodedData, dataToBeDecoded;
    dataToBeDecoded.str ("");
    std::string decodedDataString;
    UInt sizeOfFloat ( 0 ), sizeOfVector ( values.size() ), lengthOfRawData;

    switch ( numBits )
    {
        case 32:
            ASSERT ( sizeof (float) == 4, "\nThis piece of code assumes sizeof(float) == 4" );
            sizeOfFloat = sizeof (float);
            break;
        case 64:
            ASSERT ( sizeof (Real) == 8, "\nThis piece of code assumes sizeof(Real) == 8" );
            sizeOfFloat = sizeof (Real);
            break;
        default:
            ERROR_MSG ( "unmanaged float type" );
            break;
    }

    lengthOfRawData = sizeOfVector * sizeOfFloat + sizeof (int32_type);

    // assign the block of char to a stringstream (to convert it into a string)
    dataToBeDecoded.write ( line.c_str(), line.size() );

    // perform the string decoding and store the string in a stringstream
    // (to access it "bitwise")
    decodedDataString = base64_decode ( dataToBeDecoded.str() );
    decodedData.str ( decodedDataString );

    ASSERT ( lengthOfRawData == decodedDataString.size(), "unexpected line length" );


    // the first value in the string is the size of the subsequent bunch of data
    int32_type* inputInt = new int32_type;
    decodedData.read ( reinterpret_cast< char*> ( inputInt ),
                       sizeof (int32_type) );

    ASSERT ( *inputInt - sizeOfVector, "Inconsistent size of data!" );

    switch ( numBits )
    {
        case 32:
        {
            float* inputFloat = new float[sizeOfVector];
            decodedData.read ( reinterpret_cast< char*> ( inputFloat ),
                               sizeOfVector * sizeOfFloat );
            std::vector<float> inputValuesTemp (sizeOfVector);
            inputValuesTemp.assign ( inputFloat, inputFloat + sizeOfVector );
            values.assign ( inputValuesTemp.begin(), inputValuesTemp.end() );
            //for (UInt i = 0; i < inputValues.size(); i++)
            //  inputValues[i] = inputValuesTemp[i];
            break;
        }
        case 64:
        {
            Real* inputReal = new Real[sizeOfVector];
            decodedData.read ( reinterpret_cast< char*> ( inputReal ),
                               sizeOfVector * sizeOfFloat );
            values.assign ( inputReal, inputReal + sizeOfVector );
            break;
        }
    }
}


template <typename MeshType>
void
ExporterVTK<MeshType>::readASCIIData ( const std::string& line, std::vector<Real>& values )
{
    std::stringstream readData ( line );

    // simply parse the line to fill the vector of values
    for (UInt i = 0; i < values.size(); i++)
    {
        readData >> values[i];
    }

}


/*
    preliminary attempt at managing simultaneously all data associated to Nodes
    as opposed to data associated to Cells.

    Untested for now.
 */
template <typename MeshType>
void
ExporterVTK<MeshType>::composeDataArrayStream (where_Type where,
                                               std::stringstream& dataArraysStringStream)
{


    typename super::iterator_Type it;
    std::pair<typename super::iterator_Type, typename super::iterator_Type> rangeFound;

    rangeFound = this->M_whereToDataMap.equal_range (where);

    if ( rangeFound.first != rangeFound.second )
    {
        debugStream (8000) << "\n[ExporterVTK::composeDataArrayStream] found data\n";

        dataArraysStringStream.setf (std::ios_base::fixed);
        dataArraysStringStream.precision (5);
        dataArraysStringStream.width (12);

        for (it = rangeFound.first; it != rangeFound.second; ++it)
        {
            switch ( it->second.fieldType() )
            {
                case exporterData_Type::ScalarField:

                    dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
                                           << it->second.variableName() << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";

                    for (UInt iValue = 1; iValue <= it->second.size(); ++iValue)
                    {
                        dataArraysStringStream << it->second ( iValue ) << " ";
                    }
                    break;
                case exporterData_Type::VectorField:

                    dataArraysStringStream << "\t\t\t\t<DataArray type=\"Float32\" Name=\""
                                           << it->second.variableName() << "\" NumberOfComponents=\""
                                           << nDimensions << "\" format=\"ascii\">\n";

                    for (UInt iValue = 1; iValue <= it->second.size(); ++iValue)
                    {
                        for (UInt iCoor = 0; iCoor < nDimensions; ++iCoor)
                        {
                            dataArraysStringStream << it->second ( iValue + iCoor * it->second.size() ) << " ";
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
                default:
                    ERROR_MSG ( "Unknown field type" );
                    break;
            }
            dataArraysStringStream << "\n\t\t\t\t</DataArray>\n";
        }
        // return dataArraysStringStream;
    }
}


template <typename MeshType>
void ExporterVTK<MeshType>::composePVTUStream (const exporterData_Type& dvar,
                                               std::stringstream& pVTUStringStream)
{

    std::string floatTypeString;
    switch ( M_floatPrecision )
    {
        case SINGLE_PRECISION:
            floatTypeString = "Float32";
            break;
        case DOUBLE_PRECISION:
            floatTypeString = "Float64";
            break;
        default:
            ERROR_MSG ( "unmanaged float type" );
            break;
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
        default:
            ERROR_MSG ( "Cannot manage this data location" );
            break;
    }
    pVTUStringStream << "\t\t<" << whereString << " ";

    pVTUStringStream << ">\n";

    std::string formatString;
    std::stringstream nComponents;
    nComponents << dvar.fieldDim();
    switch ( M_exportMode )
    {
        case ASCII_EXPORT:
            formatString = "ascii";
            break;
        case BINARY_EXPORT:
            formatString = "binary";
            break;
        default:
            ERROR_MSG ( "Cannot manage this export mode" );
            break;
    }
    pVTUStringStream << "\t\t\t<PDataArray type=\"" << floatTypeString << "\" Name=\""
                     << dvar.variableName() << "\" NumberOfComponents=\""
                     << nComponents.str() << "\" format=\"" << formatString << "\">\n";

    pVTUStringStream << "\t\t\t</PDataArray>\n";

    pVTUStringStream << "\t\t\t<PDataArray type=\"Int32\" Name=\"GlobalId\" NumberOfComponents=\"1\" "
                     << "format=\"ascii\">\n";
    pVTUStringStream << "\t\t\t</PDataArray>\n";

    pVTUStringStream << "\t\t</" << whereString << ">\n";

    for ( Int iProc = 0; iProc < dvar.feSpacePtr()->map().comm().NumProc(); ++iProc )
    {
        std::stringstream fileName ( ( this->M_prefix + "_" + dvar.variableName() +
                                       this->M_postfix + "." + iProc + ".vtu").c_str() );

        //footer part of the file
        pVTUStringStream << "\t\t<Piece Source=\"" << fileName.str() << "\"/>\n";
    }

    pVTUStringStream << "\t</PUnstructuredGrid>\n";
    pVTUStringStream << "</VTKFile>\n";

}


template <typename MeshType>
void ExporterVTK<MeshType>::createPointsMaps ( const feSpacePtr_Type& _feSpacePtr,
                                               std::map<UInt, UInt>& globalToLocalPointsMap,
                                               std::map<UInt, UInt>& localToGlobalPointsMap,
                                               std::vector<Vector>& coordinatesOfPoints )
{
    ASSERT ( this->M_mesh.get(), "\nA pointer to a valid mesh object is required!");
    ASSERT ( _feSpacePtr.get(), "\nA pointer to a valid FESpace object is required!");

    const ID numVertices = this->M_mesh->numVertices();
    const ID numElements = this->M_mesh->numElements();

    // careful: the vertex map in the mesh is repeated. to know how many non vertex dofs I have
    // in the partitioned mesh I need to look at repeated maps
    const UInt numPoints = _feSpacePtr->map().map (Repeated)->NumMyElements() / _feSpacePtr->fieldDim();

    coordinatesOfPoints.resize ( nDimensions, ZeroVector (numPoints) );

    Real x, y, z;

    // The global ID of the considered DOF
    UInt globalPointId (0);
    // The local ID of the considered DOF
    UInt positionInPartitionedMesh (0);
    // Helper iterator
    std::pair< std::map<UInt, UInt>::iterator, bool > returnType;

    // Vertex based Dof: the coordinates are available from the Point List
    for (ID iVertex = 0; iVertex < numVertices; ++iVertex)
    {
        globalPointId = this->M_mesh->point (iVertex).id();
        localToGlobalPointsMap.insert ( std::pair<UInt, UInt> ( positionInPartitionedMesh, globalPointId ) );
        globalToLocalPointsMap.insert ( std::pair<UInt, UInt> ( globalPointId, positionInPartitionedMesh ) );

        for (ID jCoor = 0; jCoor < nDimensions; ++jCoor)
        {
            coordinatesOfPoints[jCoor][positionInPartitionedMesh] =
                this->M_mesh->point (iVertex).coordinate (jCoor);
        }
        ++positionInPartitionedMesh;
    }
    ASSERT ( positionInPartitionedMesh == numVertices, "didn't store all vertices in the maps");
    ASSERT ( localToGlobalPointsMap.size() == globalToLocalPointsMap.size(),
             "problem in storing the local to global and global to local maps" );

    // Now I store the coordinates of the supplementary nodes in a temporary vector
    for ( UInt iElement = 0; iElement < numElements; ++iElement )
    {
        _feSpacePtr->fe().update ( this->M_mesh->element ( iElement ), UPDATE_ONLY_CELL_NODES );
        for ( UInt iPoint = _feSpacePtr->dof().numLocalVertices();
                iPoint < _feSpacePtr->dof().numLocalDof(); ++iPoint )
        {
            _feSpacePtr->fe().coorMap ( x, y, z,
                                        _feSpacePtr->fe().refFE().xi ( iPoint ),
                                        _feSpacePtr->fe().refFE().eta ( iPoint ),
                                        _feSpacePtr->fe().refFE().zeta ( iPoint ) );

            globalPointId = _feSpacePtr->dof().localToGlobalMap ( iElement, iPoint );

            returnType = globalToLocalPointsMap.insert ( std::pair<UInt, UInt> ( globalPointId,
                                                                                 positionInPartitionedMesh ) );

            // by looping over mesh elements I get repetition in the list of points
            // update the maps only when finding a new point
            if ( returnType.second )
            {
                localToGlobalPointsMap.insert ( std::pair<UInt, UInt> ( positionInPartitionedMesh, globalPointId ) );
                coordinatesOfPoints[0][positionInPartitionedMesh] = x;
                coordinatesOfPoints[1][positionInPartitionedMesh] = y;
                coordinatesOfPoints[2][positionInPartitionedMesh] = z;
                ++positionInPartitionedMesh;
            }

        }
    }
    ASSERT ( positionInPartitionedMesh == numPoints, "didn't store all points in the maps" );
    ASSERT ( localToGlobalPointsMap.size() == globalToLocalPointsMap.size(),
             "problem in storing the local to global and global to local maps" );
    ASSERT ( localToGlobalPointsMap.size() == numPoints,
             "problem in storing the local to global and global to local maps" );

}


template <typename MeshType>
void ExporterVTK<MeshType>::composeVTKCollection ( const std::string& variableName,
                                                   std::stringstream& vtkCollectionStringStream )
{
    //    ASSERT( this->M_timeSteps.size(), "No time values to be saved in the VTK collection!");
    //    ASSERT( this->M_pvtuFiles[variableName].size(), "No file names to be saved in the VTK collection!");
    //    ASSERT( this->M_pvtuFiles[variableName].size() == this->M_timeSteps.size(),
    //            "The number of post-processed files does not match the number of time steps in the list!" );

    //header part of the file
    vtkCollectionStringStream << "<?xml version=\"1.0\"?>\n";
    vtkCollectionStringStream << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    vtkCollectionStringStream << "\t<Collection>\n";

    typedef std::list<Real>::const_iterator realIterator;
    realIterator iTime = this->M_timeSteps.begin();
    typedef std::list<std::string>::const_iterator stringIterator;
    stringIterator iFileName = this->M_pvtuFiles[variableName].begin();

    // Someone might have added a variable but never post processed it (so that no pvtu files are generated)
    if (iFileName != this->M_pvtuFiles[variableName].end() )
        for ( ; iTime != this->M_timeSteps.end(); ++iTime, ++iFileName)
            vtkCollectionStringStream << "\t\t<DataSet timestep=\"" << *iTime << "\" group=\"\" part=\"0\" "
                                      << "file=\"" << *iFileName << "\" />\n";

    vtkCollectionStringStream << "\t</Collection>\n";
    vtkCollectionStringStream << "</VTKFile>\n";

}


template <typename MeshType>
void ExporterVTK<MeshType>::composeVTUGeoStream ( const feSpacePtr_Type& _feSpacePtr,
                                                  const std::map<UInt, UInt>& globalToLocalPointsMap,
                                                  //const std::map<UInt, UInt>& localToGlobalPointsMap,
                                                  const std::vector<Vector>& coordinatesOfPoints,
                                                  std::stringstream& vtuGeoStringStream )
{
    ASSERT ( this->M_mesh.get(), "\nA pointer to a valid mesh object is required!");
    ASSERT ( _feSpacePtr.get(), "\nA pointer to a valid FESpace object is required!");

    debugStream (8000) << "\n[ExporterVTK::composeVTUHeaderStream]\n";

    const UInt numPoints = globalToLocalPointsMap.size();
    const UInt numElements = this->M_mesh->numElements();
    const UInt numLocalDof = _feSpacePtr->dof().numLocalDof();

    vtuGeoStringStream << "\t\t\t<Points>\n";
    vtuGeoStringStream << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"" << nDimensions
                       << "\" format=\"ascii\">\n";

    for ( UInt iPoint = 0; iPoint < numPoints; ++iPoint )
    {
        for ( UInt iCoor = 0; iCoor < nDimensions; ++iCoor )
        {
            vtuGeoStringStream << coordinatesOfPoints[iCoor][ iPoint ] << " ";
        }
    }
    vtuGeoStringStream << "\n\t\t\t\t</DataArray>\n";

    vtuGeoStringStream << "\t\t\t</Points>\n";

    // connectivity
    // cells_size = nldof*(nV+1);
    vtuGeoStringStream << "\t\t\t<Cells>\n";
    vtuGeoStringStream << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    for (UInt iElement = 0; iElement < numElements; ++iElement)
    {
        for ( UInt jPoint = 0; jPoint < numLocalDof; ++jPoint)
        {
            // UInt globalElementId( this->M_mesh->element(iElement).id() );
            const UInt globalPointId ( _feSpacePtr->dof().localToGlobalMap ( iElement, jPoint ) );
            ASSERT ( globalToLocalPointsMap.find ( globalPointId ) != globalToLocalPointsMap.end(),
                     "didn't find a local ID for global point" );
            const UInt localId ( globalToLocalPointsMap.find ( globalPointId )->second );
            vtuGeoStringStream << localId << " ";
        }
    }
    vtuGeoStringStream << "\n\t\t\t\t</DataArray>\n";

    vtuGeoStringStream << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";

    for (UInt iElement = 1; iElement <= numElements; ++iElement)
    {
        vtuGeoStringStream << iElement* numLocalDof << " ";
    }

    vtuGeoStringStream << "\n\t\t\t\t</DataArray>\n";

    vtuGeoStringStream << "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";

    // definition of cells types
    UInt vtkCellType = whichCellType (_feSpacePtr);
    for (UInt iElement = 1; iElement <= numElements; ++iElement)
    {
        vtuGeoStringStream << vtkCellType << " ";
    }

    vtuGeoStringStream << "\n\t\t\t\t</DataArray>\n";
    vtuGeoStringStream << "\t\t\t</Cells>\n";

}


template <typename MeshType>
void ExporterVTK<MeshType>::composeVTUHeaderStream ( UInt numPoints,
                                                     std::stringstream& vtuHeaderStringStream )
{

    //header part of the file
    vtuHeaderStringStream << "<?xml version=\"1.0\"?>\n";
    vtuHeaderStringStream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    vtuHeaderStringStream << "\t<UnstructuredGrid>\n";
    vtuHeaderStringStream << "\t\t<Piece NumberOfPoints=\"" << numPoints << "\""
                          << " NumberOfCells=\"" << this->M_mesh->numElements() << "\">\n";


}


template <typename MeshType>
void ExporterVTK<MeshType>::composeVTUFooterStream ( std::stringstream& vtuFooterStringStream )
{
    debugStream (8000) << "\n[ExporterVTK::composeFooter]\n";

    //footer part of the file
    vtuFooterStringStream << "\t\t</Piece>\n";
    vtuFooterStringStream << "\t</UnstructuredGrid>\n";
    vtuFooterStringStream << "</VTKFile>\n";


}


template <typename MeshType>
void
ExporterVTK<MeshType>::composeTypeDataHeaderStream (where_Type where,
                                                    std::stringstream& dataHeaderStringStream)
{
    debugStream (8000) << "\n[ExporterVTK::composeTypeDataHeaderStream] where = " << where << "\n";

    std::string whereString;
    switch ( where )
    {
        case exporterData_Type::Node:
            whereString = "PointData";
            break;
        case exporterData_Type::Cell:
            whereString = "CellData";
            break;
        default:
            ERROR_MSG ( "Cannot manage this data location")
            break;
    }
    dataHeaderStringStream << "\t\t\t<" << whereString << " ";

    dataHeaderStringStream << ">\n";
}


template <typename MeshType>
void
ExporterVTK<MeshType>::composeTypeDataFooterStream (where_Type where,
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
        default:
            ERROR_MSG ( "Cannot manage this data location")
            break;
    }
    dataFooterStringStream << "\t\t\t</" << whereString << ">\n";
}

}
#endif // define EXPORTERVTK_H
