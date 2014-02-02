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
  @brief This file provides the class  ExporterHDF5 for post-processing with hdf5

  @date 11-11-2008
  @author Simone Deparis <simone.deparis@epfl.ch>

  @contributor Radu Popescu <radu.popescu@epfl.ch>
  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef EXPORTER_HDF5_H
#define EXPORTER_HDF5_H 1

#include <sstream>


#include <Epetra_ConfigDefs.h>
#include <EpetraExt_DistArray.h>
#include <EpetraExt_HDF5.h>
#include <Epetra_Comm.h>
#include <Epetra_IntVector.h>
#include <Epetra_MultiVector.h>

#include <boost/algorithm/string.hpp>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>

#ifndef HAVE_HDF5

#warning warning you should reconfigure Trilinos with  -D TPL_ENABLE_HDF5:BOOL=ON

#else

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/util/StringUtility.hpp>
#include <lifev/core/fem/ReferenceFE.hpp>
#include <lifev/core/fem/ReferenceFEScalar.hpp>
#include <lifev/core/filter/Exporter.hpp>

namespace LifeV
{

//! Hdf5 data exporter, implementation of Exporter
/*!
  @author Simone Deparis <simone.deparis@epfl.ch>
  @author Radu Popescu <radu.popescu@epfl.ch>

  Usage: two steps
  <ol>
  <li> first: add the variables using addVariable
  <li> second: call postProcess( time );
  </ol>
*/
template<typename MeshType>
class ExporterHDF5 : public Exporter<MeshType>
{

public:
    //! @name Public typedefs
    //@{
    typedef MeshType mesh_Type;
    typedef Exporter<MeshType> super;
    typedef typename super::meshPtr_Type meshPtr_Type;
    typedef typename super::vector_Type    vector_Type;
    typedef typename super::vectorPtr_Type vectorPtr_Type;
    typedef typename super::exporterData_Type exporterData_Type;

    typedef EpetraExt::HDF5 hdf5_Type;
    typedef boost::shared_ptr<hdf5_Type> hdf5Ptr_Type;
    typedef std::vector<std::vector<Int> > graph_Type;
    typedef boost::shared_ptr<graph_Type> graphPtr_Type;
    typedef boost::shared_ptr<std::vector<meshPtr_Type> > serial_meshPtr_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor for ExporterHDF5
    ExporterHDF5();

    //! Constructor for ExporterHDF5
    /*!
      @param dfile the GetPot data file where you must provide an [exporter] section with:
      "start"     (start index for sections in the hdf5 data structure 0 for 000, 1 for 001 etc.),
      "save"      (how many time steps per postprocessing)
      "multimesh" ( = true if the mesh has to be saved at each post-processing step)
      @param mesh the mesh
      @param the prefix for the case file (ex. "test" for test.case)
      @param the procId determines de CPU id. if negative, it ussemes there is only one processor
    */
    ExporterHDF5 (const GetPot& dfile, meshPtr_Type mesh, const std::string& prefix, const Int& procId);

    //! Constructor for ExporterHDF5 without prefix and procID
    /*!
      @param dfile the GetPot data file where you must provide an [exporter] section with:
      "start"     (start index for sections in the hdf5 data structure 0 for 000, 1 for 001 etc.),
      "save"      (how many time steps per postprocessing)
      "multimesh" ( = true if the mesh has to be saved at each post-processing step)
      @param mesh the mesh
    */
    ExporterHDF5 (const GetPot& dfile, const std::string& prefix);

    //! Destructor for ExporterHDF5
    virtual ~ExporterHDF5() {}

    //@}

    //! @name Public Methods
    //@{
    virtual void postProcess (const Real& time);

    //! Import data from previous simulations at a certain time
    /*!
      @param Time the time of the data to be imported
      @return number of iteration corresponding at the time step
    */
    UInt importFromTime ( const Real& Time );

    //! Import data from previous simulations at a certain time
    /*!
      @param Time the time of the data to be imported
      @return the simulation time corresponding to the iteration
    */
    Real importFromIter ( const UInt& );

    //! Import data from previous simulations
    /*!
      @param time the solver time
    */
    void import (const Real& startTime, const Real& dt); // dt is used to rebuild the history up to now

    //! Import data from previous simulations
    /*!
      @param time the solver time
    */
    void import (const Real& time);

    //! Close the Hdf5 file
    /*!
      Close the HDF5 file.
    */
    void closeFile()
    {
        M_HDF5->Close();
    }

    //! Read variable
    void readVariable ( exporterData_Type& dvar);

    //@}

    //! @name Set Methods
    //@{

    //! Set data from file.
    /*!
     * @param dataFile data file.
     * @param section section in the data file.
     */
    void setDataFromGetPot ( const GetPot& dataFile, const std::string& section = "exporter" );

    //@}

    //! @name Get Methods
    //@{

    //! returns the type of the map to use for the VectorEpetra
    MapEpetraType mapType() const;

    //@}

protected:

    //! @name Protected Methods
    //@{
    //! Define the shape of the elements
    void defineShape();
    //! write empty xdmf file
    void writeInitXdmf();
    //! append to xdmf file
    void writeXdmf (const Real& time);
    //! save position and write closing lines
    void writeCloseLinesXdmf();
    //! remove closing lines
    void removeCloseLinesXdmf();

    void writeTopology  ( std::ofstream& xdmf );
    void writeGeometry  ( std::ofstream& xdmf );
    void writeAttributes ( std::ofstream& xdmf );
    void writeScalarDatastructure  ( std::ofstream& xdmf, const exporterData_Type& dvar );
    void writeVectorDatastructure  ( std::ofstream& xdmf, const exporterData_Type& dvar );

    void writeVariable (const exporterData_Type& dvar);
    void writeScalar (const exporterData_Type& dvar);
    void writeVector (const exporterData_Type& dvar);

    void writeGeometry();

    void readScalar ( exporterData_Type& dvar);
    void readVector ( exporterData_Type& dvar);
    //@}

    //! @name Protected data members
    //@{
    hdf5Ptr_Type      M_HDF5;
    std::ofstream     M_xdmf;

    const std::string M_closingLines;
    std::streampos    M_closingLinesPosition;
    std::string       M_outputFileName;

    //! do we want to write on file the connectivity?
    bool                        M_printConnectivity;
    //@}

};



// ===================================================
// Constructors
// ===================================================
template<typename MeshType>
ExporterHDF5<MeshType>::ExporterHDF5() :
    super               (),
    M_HDF5              (),
    M_closingLines      ( "\n    </Grid>\n\n  </Domain>\n</Xdmf>\n"),
    M_outputFileName    ( "noninitialisedFileName" ),
    M_printConnectivity ( true )
{
}

template<typename MeshType>
ExporterHDF5<MeshType>::ExporterHDF5 (const GetPot& dfile, meshPtr_Type mesh, const std::string& prefix,
                                      const Int& procId) :
    super               ( dfile, prefix ),
    M_HDF5              (),
    M_closingLines      ( "\n    </Grid>\n\n  </Domain>\n</Xdmf>\n"),
    M_outputFileName    ( "noninitialisedFileName" )
{
    M_printConnectivity = dfile ( ( prefix + "/printConnectivity" ).data(), 1);
    this->setMeshProcId ( mesh, procId );
}

template<typename MeshType>
ExporterHDF5<MeshType>::ExporterHDF5 (const GetPot& dfile, const std::string& prefix) :
    super               ( dfile, prefix ),
    M_HDF5              (),
    M_closingLines      ( "\n    </Grid>\n\n  </Domain>\n</Xdmf>\n"),
    M_outputFileName    ( "noninitialisedFileName" )
{
    M_printConnectivity = dfile ( ( prefix + "/printConnectivity" ).data(), 1);
}

// ===================================================
// Methods
// ===================================================

template<typename MeshType>
void ExporterHDF5<MeshType>::postProcess (const Real& time)
{
    if ( M_HDF5.get() == 0)
    {
        M_HDF5.reset (new hdf5_Type (this->M_dataVector.begin()->storedArrayPtr()->comm() ) );
        M_outputFileName = this->M_prefix + ".h5";
        M_HDF5->Create (this->M_postDir + M_outputFileName);

        // write empty xdmf file
        writeInitXdmf();

        if (!this->M_multimesh)
        {
            writeGeometry(); // see also writeGeometry
            M_HDF5->Flush();
        }
    }

    // typedef std::list< exporterData_Type >::const_iterator Iterator;

    this->computePostfix();

    std::size_t found ( this->M_postfix.find ( "*" ) );
    if ( found == std::string::npos )
    {
        if (!this->M_procId)
        {
            std::cout << "  X-  HDF5 post-processing ...                 " << std::flush;
        }
        LifeChrono chrono;
        chrono.start();
        for (typename super::dataVectorIterator_Type i = this->M_dataVector.begin(); i != this->M_dataVector.end(); ++i)
        {
            writeVariable (*i);
        }
        // pushing time
        this->M_timeSteps.push_back (time);

        writeXdmf (time);

        if (this->M_multimesh)
        {
            writeGeometry(); // see also writeGeometry
        }

        chrono.stop();

        // Write to file without closing the file
        M_HDF5->Flush();

        if (!this->M_procId)
        {
            std::cout << "done in " << chrono.diff() << " s." << std::endl;
        }
    }
}

template<typename MeshType>
UInt ExporterHDF5<MeshType>::importFromTime ( const Real& Time )
{
    // Container for the time and the postfix
    std::pair< Real, Int > SelectedTimeAndPostfix;
    if ( !this->M_procId )
    {
        // Open the xmf file
        std::ifstream xmfFile;
        xmfFile.open ( ( this->M_postDir + this->M_prefix + ".xmf" ).c_str(), std::ios::in );

        // Vector of TimeStep
        std::vector< std::pair< Real, Int > > TimeAndPostfix;
        if ( xmfFile.is_open() )
        {
            // Define some variables
            std::string line;
            std::vector<std::string> stringsVector;

            // Read one-by-one all the lines of the file
            while ( !xmfFile.eof() )
            {
                std::getline ( xmfFile, line, '\n' );

                // If the line begin with "<!-- Time " it is the beginning of a new block
                if ( !line.compare ( 0, 10, "<!-- Time " ) )
                {
                    boost::split ( stringsVector, line, boost::is_any_of ( " " ) );
                    TimeAndPostfix.push_back ( std::make_pair ( string2number ( stringsVector[2] ),
                                                                string2number ( stringsVector[4] ) ) );
                }
            }
        }
        xmfFile.close();

        // Find the closest time step
        SelectedTimeAndPostfix = TimeAndPostfix.front();
        for ( std::vector< std::pair< Real, Int > >::const_iterator i = TimeAndPostfix.begin();
                i < TimeAndPostfix.end() ; ++i )
            if ( std::fabs ( SelectedTimeAndPostfix.first - Time ) >= std::fabs ( (*i).first - Time ) )
            {
                SelectedTimeAndPostfix = *i;
            }
    }
    this->M_dataVector.begin()->storedArrayPtr()->comm().Broadcast ( &SelectedTimeAndPostfix.second, 1, 0 );
    this->M_timeIndex = SelectedTimeAndPostfix.second;
    this->computePostfix();

    // Importing
    if ( !this->M_procId )
    {
        std::cout << "  X-  HDF5 importing ...                       " << std::flush;
    }

    LifeChrono chrono;
    chrono.start();
    for (typename super::dataVectorIterator_Type i = this->M_dataVector.begin(); i != this->M_dataVector.end(); ++i)
    {
        this->readVariable (*i);
    }

    chrono.stop();
    if ( !this->M_procId )
        std::cout << "done in " << chrono.diff() << " s. (Time " << SelectedTimeAndPostfix.first
                  << ", Iteration " << SelectedTimeAndPostfix.second << " )" << std::endl;

    return static_cast <UInt> ( SelectedTimeAndPostfix.second );
}

template<typename MeshType>
Real ExporterHDF5<MeshType>::importFromIter ( const UInt& iter )
{
    // Container for the time and the postfix
    std::pair< Real, Int > SelectedTimeAndPostfix;
    if ( !this->M_procId )
    {
        // Open the xmf file
        std::ifstream xmfFile;
        xmfFile.open ( ( this->M_postDir + this->M_prefix + ".xmf" ).c_str(), std::ios::in );

        // Vector of TimeStep
        std::vector< std::pair< Real, UInt > > TimeAndPostfix;

        if ( xmfFile.is_open() )
        {
            // Define some variables
            std::string line;
            std::vector<std::string> stringsVector;

            // Read one-by-one all the lines of the file
            while ( !xmfFile.eof() )
            {
                std::getline ( xmfFile, line, '\n' );

                // If the line begin with "<!-- Time " it is the beginning of a new block
                if ( !line.compare ( 0, 10, "<!-- Time " ) )
                {
                    boost::split ( stringsVector, line, boost::is_any_of ( " " ) );
                    TimeAndPostfix.push_back ( std::make_pair ( string2number ( stringsVector[2] ),
                                                                string2number ( stringsVector[4] ) ) );
                }
            }
        }

        xmfFile.close();

        // Find the closest time step
        SelectedTimeAndPostfix = TimeAndPostfix.front();
        bool found             = false;

        for ( std::vector< std::pair< Real, UInt > >::const_iterator i = TimeAndPostfix.begin();
                i < TimeAndPostfix.end()  ; ++i )
        {
            if ( i->second == iter )
            {
                SelectedTimeAndPostfix = *i;
                found = true;
                break;
            }
        }

        ASSERT (found, "Selected iteration not found");

    }

    this->M_dataVector.begin()->storedArrayPtr()->Comm().Broadcast ( &SelectedTimeAndPostfix.second, 1, 0 );
    this->M_timeIndex = SelectedTimeAndPostfix.second;

    std::ostringstream index;
    index.fill ('0');

    index << std::setw (5) << this->M_timeIndex;
    this->M_postfix = "." + index.str();

    // Importing
    if ( !this->M_procId )
        std::cout << "  X-  HDF5 importing iteration "
                  << index.str()
                  << " at time " << SelectedTimeAndPostfix.first
                  << " ... " << std::flush;

    LifeChrono chrono;
    chrono.start();
    for (typename super::dataVectorIterator_Type i = this->M_dataVector.begin(); i != this->M_dataVector.end(); ++i)
    {
        this->readVariable (*i);
    }
    chrono.stop();

    if ( !this->M_procId )
        std::cout << "done in " << chrono.diff() << " s. (Time " << SelectedTimeAndPostfix.first
                  << ", Iteration " << SelectedTimeAndPostfix.second << " )" << std::endl;

    return SelectedTimeAndPostfix.first;
}


template<typename MeshType>
void ExporterHDF5<MeshType>::import (const Real& Tstart, const Real& dt)
{
    // dt is used to rebuild the history up to now
    Real time (Tstart - this->M_timeIndex * dt);

    for ( UInt count (0); count < this->M_timeIndex; ++count )
    {
        this->M_timeSteps.push_back (time);
        time += dt;
    }

    time += dt;

    import (time);
}

template<typename MeshType>
void ExporterHDF5<MeshType>::import (const Real& time)
{
    if ( M_HDF5.get() == 0)
    {
        M_HDF5.reset (new hdf5_Type (this->M_dataVector.begin()->storedArrayPtr()->comm() ) );
        M_HDF5->Open (this->M_postDir + this->M_prefix + ".h5"); //!! Simone
    }

    this->M_timeSteps.push_back (time);

    this->computePostfix();

    assert ( this->M_postfix != "*****" );

    if (!this->M_procId)
    {
        std::cout << "  X-  HDF5 importing ...      " << std::endl;
    }

    LifeChrono chrono;
    chrono.start();
    for (typename super::dataVectorIterator_Type i = this->M_dataVector.begin(); i != this->M_dataVector.end(); ++i)
    {
        this->readVariable (*i); ///!!! Simone
    }
    chrono.stop();
    if (!this->M_procId)
    {
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }
}

template <typename MeshType>
void ExporterHDF5<MeshType>::readVariable (exporterData_Type& dvar)
{
    if ( M_HDF5.get() == 0)
    {
        M_HDF5.reset (new hdf5_Type (dvar.storedArrayPtr()->blockMap().Comm() ) );
        M_HDF5->Open (this->M_postDir + this->M_prefix + ".h5"); //!! Simone
    }
    super::readVariable (dvar);
}

// ===================================================
// Set Methods
// ===================================================
template<typename MeshType>
void ExporterHDF5<MeshType>::setDataFromGetPot ( const GetPot& dataFile, const std::string& section )
{
    super::setDataFromGetPot ( dataFile, section );
    M_printConnectivity = dataFile ( ( section + "/printConnectivity" ).data(), 1);
}

// ===================================================
// Get Methods
// ===================================================
template<typename MeshType>
MapEpetraType ExporterHDF5<MeshType>::mapType() const
{
    return Unique;
}

// ===================================================
// Protected Methods
// ===================================================
template<typename MeshType>
void ExporterHDF5<MeshType>::defineShape()
{
}

// write empty xdmf file
template <typename MeshType>
void ExporterHDF5<MeshType>::writeInitXdmf()
{
    if (this->M_procId == 0)
    {
        M_xdmf.open ( (this->M_postDir + this->M_prefix + ".xmf").c_str(), std::ios_base::out );

        M_xdmf << "<?xml version=\"1.0\" ?>\n"
               << "<!DOCTYPE Xdmf SYSTEM \""
               << this->M_prefix
               << ".xdmf\" [\n"
               << "<!ENTITY DataFile \""
               << this->M_prefix
               << ".h5\">\n"
               << "]>\n"
               << "<!-- "
               << this->M_prefix
               << ".h5 is generated by LifeV -->\n"
               << "<Xdmf>\n"
               << "  <Domain Name=\""
               << this->M_prefix
               << "\">\n"
               << "    <Grid Name=\""
               << this->M_prefix
               << "Grid\" GridType=\"Collection\" CollectionType=\"Temporal\">\n"
               << "\n";

        writeCloseLinesXdmf();
    }
}

template <typename MeshType>
void ExporterHDF5<MeshType>::writeXdmf (const Real& time)
{
    /*
      strategy: write the topology,
      <Topology
      Type="Tetrahedron"
      NumberOfElements="183"
      BaseOffset="1">
      <DataStructure Format="HDF"
      Dimensions="602  4" ???
      DataType="Int"
      Precision="8">
      &DataFile;:/Connections/Values
      </DataStructure>
      </Topology>
      and then the geometry
      <Geometry Type="X_Y_Z">
      <DataStructure Format="HDF"
      Dimensions="183"
      DataType="Float"
      Precision="8">
      &DataFile;:/PointsX/Values
      </DataStructure>
      <DataStructure Format="HDF"
      Dimensions="183"
      DataType="Float"
      Precision="8">
      &DataFile;:/PointsY/Values
      </DataStructure>
      <DataStructure Format="HDF"
      Dimensions="183"
      DataType="Float"
      Precision="8">
      &DataFile;:/PointsY/Values
      </DataStructure>
      </Geometry>


      In this aim we create
      two Epetra_IntVector, one with a repeated map and the second with a unique map
      two Epetra_MultiVector, one with a repeated map and the second with a unique map

      The Int vectors are for the connections, the Real ones for the Vertices
    */

    if (this->M_procId == 0)
    {
        removeCloseLinesXdmf();

        // write grid with time, topology, geometry and attributes
        // NOTE: The first line (<!-- Time t Iteration i -->) is used in function importFromTime.
        //       Check compatibility after any change on it!
        M_xdmf <<
               "<!-- Time " << time << " Iteration " << this->M_postfix.substr (1, 5) << " -->\n" <<
               "    <Grid Name=\"Mesh " << time << "\">\n" <<
               "      <Time TimeType=\"Single\" Value=\"" << time << "\" />\n";
        writeTopology (M_xdmf);
        writeGeometry (M_xdmf);
        writeAttributes (M_xdmf);

        M_xdmf << "\n"
               "    </Grid>\n\n";



        // write closing lines
        writeCloseLinesXdmf();
    }
}

// save position and write closing lines
template <typename MeshType>
void ExporterHDF5<MeshType>::writeCloseLinesXdmf()
{
    // save position
    M_closingLinesPosition = M_xdmf.tellp();

    // write closing lines
    M_xdmf << M_closingLines;
    M_xdmf.flush();

}

// remove closing lines
template <typename MeshType>
void ExporterHDF5<MeshType>::removeCloseLinesXdmf()
{
    M_xdmf.seekp (M_closingLinesPosition);
}

template <typename MeshType>
void ExporterHDF5<MeshType>::writeTopology  ( std::ofstream& xdmf )
{
    std::string FEstring;

    switch ( MeshType::elementShape_Type::S_shape )
    {
        case TETRA:
            FEstring = "Tetrahedron";
            break;
        case HEXA:
            FEstring = "Hexahedron";
            break;
        case TRIANGLE:
            FEstring = "Triangle";
            break;
        case QUAD:
            FEstring = "Quadrilateral";
            break;
        case LINE:
            FEstring = "Polyline";
            break;
        default:
            ERROR_MSG ( "FE not allowed in HDF5 Exporter" );
    }

    xdmf << "      <Topology\n"
         << "         Type=\""
         << FEstring
         << "\"\n"
         << "         NumberOfElements=\""
         << this->M_mesh->numGlobalElements()
         << "\"\n"
         << "         BaseOffset=\""
         << 0
         << "\">\n"
         << "         <DataStructure Format=\"HDF\"\n"
         << "                        Dimensions=\""
         << this->M_mesh->numGlobalElements()
         << " "
         << this->M_mesh->numLocalVertices()
         << "\"\n" << "                        DataType=\"Int\"\n"
         << "                        Precision=\"8\">\n" <<  "             "
         << M_outputFileName << ":/Connections/Values\n"
         << "         </DataStructure>\n" << "      </Topology>\n";
}

template <typename MeshType>
void ExporterHDF5<MeshType>::writeGeometry  ( std::ofstream& xdmf )
{

    std::string postfix_string;

    // see also in postProcess
    if (this->M_multimesh)
    {
        postfix_string = this->M_postfix;
    }
    else
    {
        postfix_string = "";
    }


    xdmf <<
         "      <Geometry Type=\"X_Y_Z\">\n" <<
         "         <DataStructure Format=\"HDF\"\n" <<
         "                        Dimensions=\"" << this->M_mesh->numGlobalVertices() << "\"\n" <<
         "                        DataType=\"Float\"\n" <<
         "                        Precision=\"8\">\n" <<
         "             " << M_outputFileName << ":/" << "PointsX" << postfix_string << "/Values\n" <<
         "         </DataStructure>\n" <<
         "         <DataStructure Format=\"HDF\"\n" <<
         "                        Dimensions=\"" << this->M_mesh->numGlobalVertices() << "\"\n" <<
         "                        DataType=\"Float\"\n" <<
         "                        Precision=\"8\">\n" <<
         "             " << M_outputFileName << ":/" << "PointsY" << postfix_string << "/Values\n" <<
         "         </DataStructure>\n" <<
         "         <DataStructure Format=\"HDF\"\n" <<
         "                        Dimensions=\"" << this->M_mesh->numGlobalVertices() << "\"\n" <<
         "                        DataType=\"Float\"\n" <<
         "                        Precision=\"8\">\n" <<
         "             " << M_outputFileName << ":/" << "PointsZ" << postfix_string << "/Values\n" <<
         "         </DataStructure>\n" <<
         "      </Geometry>\n" <<
         "\n";
}

template <typename MeshType>
void ExporterHDF5<MeshType>::writeAttributes  ( std::ofstream& xdmf )
{

    // Loop on the variables to output
    for (typename super::dataVectorIterator_Type i = this->M_dataVector.begin(); i != this->M_dataVector.end(); ++i)
    {
        xdmf <<
             "\n      <Attribute\n" <<
             "         Type=\"" << i->typeName() << "\"\n" <<
             "         Center=\"" << i->whereName() << "\"\n" <<
             "         Name=\"" << i->variableName() << "\">\n";

        switch ( i->fieldType() )
        {
            case exporterData_Type::ScalarField:
                writeScalarDatastructure (xdmf, *i);
                break;
            case exporterData_Type::VectorField:
                writeScalarDatastructure (xdmf, *i);
                //writeVectorFieldstructure(xdmf, *i);
                break;
        }

        xdmf <<
             "      </Attribute>\n";
    }
}

template <typename MeshType>
void ExporterHDF5<MeshType>::writeScalarDatastructure  ( std::ofstream& xdmf, const exporterData_Type& dvar )
{

    Int globalUnknowns (0);
    switch ( dvar.where() )
    {
        case exporterData_Type::Node:
            globalUnknowns = this->M_mesh->numGlobalVertices();
            break;
        case exporterData_Type::Cell:
            globalUnknowns = this->M_mesh->numGlobalElements();
            break;
    }

    //exported vectors are threedimensional also in 2D (check)
    UInt dim = (dvar.fieldDim() == 1) ? 1 : nDimensions;

    // First: hyperslab definition, then description of the data
    xdmf <<

         "         <DataStructure ItemType=\"HyperSlab\"\n" <<
         "                        Dimensions=\"" << globalUnknowns << " " << dim << "\"\n" <<
         "                        Type=\"HyperSlab\">\n" <<
         "           <DataStructure  Dimensions=\"3 2\"\n" <<
         "                           Format=\"XML\">\n" <<
         "               0    0\n" <<
         "               1    1\n" <<
         "               " << globalUnknowns << " " << dim << "\n" <<
         "           </DataStructure>\n" <<

         "           <DataStructure  Format=\"HDF\"\n" <<
         "                           Dimensions=\"" << dvar.numDOF() << " " << dim << "\"\n" <<
         "                           DataType=\"Float\"\n" <<
         "                           Precision=\"8\">\n" <<
         "               " << M_outputFileName << ":/" << dvar.variableName()
         << this->M_postfix  << "/Values\n" << // see also in writeVector/scalar
         "           </DataStructure>\n" <<
         "         </DataStructure>\n";

}

template <typename MeshType>
void ExporterHDF5<MeshType>::writeVectorDatastructure  ( std::ofstream& xdmf, const exporterData_Type& dvar )
{


    std::string coord[3] = {"X", "Y", "Z"}; // see also wr_vector

    xdmf << "         <DataStructure ItemType=\"Function\"\n"
         << "                        Dimensions=\""
         << this->M_mesh->numGlobalVertices()
         << " "
         << nDimensions
         << "\"\n"
         << "                        Function=\"JOIN($0 , $1, $2)\">\n";

    for ( UInt i (0); i < nDimensions; ++i )
    {
        xdmf << "           <DataStructure  Format=\"HDF\"\n"
             << "                           Dimensions=\""
             << this->M_mesh->numGlobalVertices()
             << " 1\"\n"
             << "                           DataType=\"Float\"\n"
             << "                           Precision=\"8\">\n"
             << "               "
             << M_outputFileName
             << ":/"
             << dvar.variableName()
             << coord[i]
             << this->M_postfix
             << "/Values\n"
             << "           </DataStructure>\n";
    }

    xdmf << "         </DataStructure>\n";
}

template <typename MeshType>
void ExporterHDF5<MeshType>::writeVariable (const exporterData_Type& dvar)
{

    switch ( dvar.fieldType() )
    {
        case exporterData_Type::ScalarField:
            writeScalar (dvar);
            break;
        case exporterData_Type::VectorField:
            writeVector (dvar);
            break;
    }
}

template <typename MeshType>
void ExporterHDF5<MeshType>::writeScalar (const exporterData_Type& dvar)
{
    /* Examples:
       M_HDF5->Write("map-" + toString(Comm.NumProc()), Map);
       M_HDF5->Write("matrix", Matrix);
       M_HDF5->Write("LHS", LHS);
       M_HDF5->Write("RHS", RHS);
    */

    UInt size  = dvar.numDOF();
    UInt start = dvar.start();

    MapEpetra subMap (dvar.storedArrayPtr()->blockMap(), start, size);
    vector_Type subVar (subMap);
    subVar.subset (*dvar.storedArrayPtr(), start);

    std::string varname (dvar.variableName() + this->M_postfix); // see also in writeAttributes
    bool writeTranspose (true);
    M_HDF5->Write (varname, subVar.epetraVector(), writeTranspose );
}

template <typename MeshType>
void ExporterHDF5<MeshType>::writeVector (const exporterData_Type& dvar)
{

    UInt size  = dvar.numDOF();
    UInt start = dvar.start();

    // solution array has to be reordered and stored in a Multivector.
    // Using auxiliary arrays:
    Real**                                  ArrayOfPointers (new Real*[nDimensions]);
    boost::shared_array< boost::shared_ptr<vector_Type> >
    		ArrayOfVectors (new boost::shared_ptr<vector_Type>[nDimensions]);

    Int MyLDA;


    // Building subsets (new vectors) of the original data and than taking a view of them to
    // build a multivector.
    // Note: the contents of ArrayOfPointers[0,1,2] must not be deleted explicitly, since their
    // content belongs to ArrayOfVectors[0,1,2].
    // ArrayOfVectors[0,1,2] are deleted when ArrayOfVectors is destroyed

    for (UInt d ( 0 ); d < dvar.fieldDim(); ++d)
    {
        MapEpetra subMap (dvar.storedArrayPtr()->blockMap(), start + d * size, size);
        ArrayOfVectors[d].reset (new  vector_Type (subMap) );
        ArrayOfVectors[d]->subset (*dvar.storedArrayPtr(), start + d * size);

        ArrayOfVectors[d]->epetraVector().ExtractView (&ArrayOfPointers[d], &MyLDA);
    }

    for (UInt d ( dvar.fieldDim() ); d < nDimensions; ++d)
    {
        MapEpetra subMap (dvar.storedArrayPtr()->blockMap(), start, size);
        ArrayOfVectors[d].reset (new  vector_Type (subMap) );

        ArrayOfVectors[d]->epetraVector().ExtractView (&ArrayOfPointers[d], &MyLDA);
    }

    MapEpetra subMap (dvar.storedArrayPtr()->blockMap(), start, size);
    Epetra_MultiVector multiVector (View, *subMap.map (Unique), ArrayOfPointers, nDimensions);


    bool writeTranspose (true);
    std::string varname (dvar.variableName() + this->M_postfix); // see also in writeAttributes
    M_HDF5->Write (varname, multiVector, writeTranspose);

    delete[] ArrayOfPointers;
}

template <typename MeshType>
void ExporterHDF5<MeshType>::writeGeometry()
{
    /*
      2 variables:
      &DataFile;:/Connections/Values
      &DataFile;:/PointsX/Values
      &DataFile;:/PointsY/Values
      &DataFile;:/PointsZ/Values
      note:
      if (this->M_multimesh)
      &DataFile;:/Points+ this->M_postfix + X/Values
      see also writeGeometry
    */

    /* Steps:
       generate local  connections and local coordinates (not repeated)
       write out
    */

    // We need a map ,but it's not always possible to use that from the variables
    // (if we write out a P0 variable)
    // We build a map for the connections based on the element numbers and for the points we fake a P1 map

    ASSERT (this->M_dataVector.size() > 0 , "hdf5exporter: ListData is empty");

    // Connections
    // Need to use elements not dofs for this map. Recover local element lists
    UInt numberOfPoints = MeshType::elementShape_Type::S_numPoints;

    std::vector<Int> elementList;
    UInt ownedElements = this->M_mesh->elementList().countElementsWithFlag ( EntityFlags::GHOST, &Flag::testOneNotSet );
    elementList.reserve ( ownedElements * numberOfPoints );
    UInt elementCount = 0;
    for ( ID i = 0; i < this->M_mesh->numElements(); ++i )
    {
        typename MeshType::element_Type const& element (this->M_mesh->element (i) );
        if ( element.isOwned() )
        {
            UInt lid = elementCount * numberOfPoints;
            for (ID j = 0; j < numberOfPoints; ++j, ++lid)
            {
                elementList[lid] = element.id() * numberOfPoints + j;
            }
            elementCount++;
        }
    }

    Epetra_Map connectionsMap ( this->M_mesh->numGlobalElements() *numberOfPoints,
                                ownedElements * numberOfPoints,
                                &elementList[0],
                                0, this->M_dataVector.begin()->storedArrayPtr()->comm() );

    Epetra_IntVector connections (connectionsMap);
    elementCount = 0;
    for (ID i = 0; i < this->M_mesh->numElements(); ++i)
    {
        typename MeshType::element_Type const& element (this->M_mesh->element (i) );
        if ( element.isOwned() )
        {
            UInt lid = elementCount * numberOfPoints;
            for (ID j = 0; j < numberOfPoints; ++j, ++lid)
            {
                connections[lid] = element.point (j).id();
            }
            elementCount++;
        }
    }

    this->M_dataVector.begin()->storedArrayPtr()->comm().Barrier();

    // Points

    // Build a map for linear elements, even though the original FE might be P0
    // This gives the right map for the coordinate arrays

    MapEpetra subMap;
    switch ( MeshType::elementShape_Type::S_shape )
    {
        case TETRA:
        {
            const ReferenceFE& refFEP1 = feTetraP1;
            DOF tmpDof ( *this->M_mesh, refFEP1 );
            std::vector<Int> myGlobalElements ( tmpDof.globalElements ( *this->M_mesh ) );
            // Create the map
            MapEpetra tmpMapP1 ( -1, myGlobalElements.size(), &myGlobalElements[0],
                                 this->M_dataVector.begin()->storedArrayPtr()->mapPtr()->commPtr() );
            subMap = tmpMapP1;
            break;
        }
        case TRIANGLE:
        {
            const ReferenceFE& refFEP1 = feTriaP1;
            DOF tmpDof ( *this->M_mesh, refFEP1 );
            std::vector<Int> myGlobalElements ( tmpDof.globalElements ( *this->M_mesh ) );
            // Create the map
            MapEpetra tmpMapP1 ( -1, myGlobalElements.size(), &myGlobalElements[0],
                                 this->M_dataVector.begin()->storedArrayPtr()->mapPtr()->commPtr() );
            subMap = tmpMapP1;
            break;
        }
        case HEXA:
        {
            const ReferenceFE& refFEQ1 = feHexaQ1;
            DOF tmpDof ( *this->M_mesh, refFEQ1 );
            std::vector<Int> myGlobalElements ( tmpDof.globalElements ( *this->M_mesh ) );
            // Create the map
            MapEpetra tmpMapQ1 ( -1, myGlobalElements.size(), &myGlobalElements[0],
                                 this->M_dataVector.begin()->storedArrayPtr()->mapPtr()->commPtr() );
            subMap = tmpMapQ1;
            break;
        }
        case LINE:
        {
            const ReferenceFE& refFEP11D = feSegP1;
            DOF tmpDof ( *this->M_mesh, refFEP11D );
            std::vector<Int> myGlobalElements ( tmpDof.globalElements ( *this->M_mesh ) );
            // Create the map
            MapEpetra tmpMapQ11D ( -1, myGlobalElements.size(), &myGlobalElements[0],
                                   this->M_dataVector.begin()->storedArrayPtr()->mapPtr()->commPtr() );
            subMap = tmpMapQ11D;
            break;
        }
        default:
            ERROR_MSG ( "FE not allowed in HDF5 Exporter" );

    }

    VectorEpetra pointsX (subMap);
    VectorEpetra pointsY (subMap);
    VectorEpetra pointsZ (subMap);

    Int gid;
    for (ID i = 0; i < this->M_mesh->numVertices(); ++i)
    {
        typename MeshType::point_Type point;
        if ( this->M_multimesh )
        {
            point = this->M_mesh->point (i);
        }
        else
        {
            point = this->M_mesh->meshTransformer().pointInitial (i);
        }

        if ( point.isOwned() )
        {

            gid = point.id();

            bool insertedX (true);
            bool insertedY (true);
            bool insertedZ (true);

            insertedX = insertedX && pointsX.setCoefficient (gid, point.x() );
            insertedY = insertedY && pointsY.setCoefficient (gid, point.y() );
            insertedZ = insertedZ && pointsZ.setCoefficient (gid, point.z() );
        }
    }

    // Now we are ready to export the vectors to the hdf5 file

    std::string pointsXVarname ("PointsX");
    std::string pointsYVarname ("PointsY");
    std::string pointsZVarname ("PointsZ");
    std::string connectionsVarname ("Connections");

    if (this->M_multimesh)
    {
        pointsXVarname      += this->M_postfix; // see also in writeGeometry
        pointsYVarname      += this->M_postfix; // see also in writeGeometry
        pointsZVarname      += this->M_postfix; // see also in writeGeometry
    }

    if ( this->M_printConnectivity )
    {
        M_HDF5->Write (connectionsVarname, connections);
        this->M_printConnectivity = false;
    }

    // bool writeTranspose (true);
    M_HDF5->Write (pointsXVarname, pointsX.epetraVector(), true);
    M_HDF5->Write (pointsYVarname, pointsY.epetraVector(), true);
    M_HDF5->Write (pointsZVarname, pointsZ.epetraVector(), true);
}

template <typename MeshType>
void ExporterHDF5<MeshType>::readScalar (exporterData_Type& dvar)
{

    UInt size  = dvar.numDOF();
    UInt start = dvar.start();

    MapEpetra subMap (dvar.storedArrayPtr()->blockMap(), start, size);
    Epetra_MultiVector* subVar (0);

    std::string varname (dvar.variableName() ); // see also in writeAttributes
    if (this->M_postfix != "")
    {
        varname += this->M_postfix;
    }
    bool readTranspose (true);
    M_HDF5->Read (varname, *subMap.map (this->mapType() ), subVar, readTranspose);

    dvar.storedArrayPtr()->subset (*subVar, subMap, 0, start );

    delete subVar;

}

template <typename MeshType>
void ExporterHDF5<MeshType>::readVector ( exporterData_Type& dvar)
{
    UInt size  = dvar.numDOF();
    UInt start = dvar.start();

    using namespace boost;

    // solution array has first to be read has Multivector.

    // first read the multivector:
    MapEpetra subMap (dvar.storedArrayPtr()->blockMap(), start, size);
    Epetra_MultiVector* subVar (0);

    bool readTranspose (true);
    std::string varname (dvar.variableName() ); // see also in writeAttributes

    if (this->M_postfix != "")
    {
        varname += this->M_postfix;
    }

    M_HDF5->Read (varname, *subMap.map (this->mapType() ), subVar, readTranspose);


    // then put back value in our VectorEpetra

    for (UInt d ( 0 ); d < dvar.fieldDim(); ++d)
    {
        dvar.storedArrayPtr()->subset (*subVar, subMap,  0, start + d * size, d );
    }

    delete subVar;
}

} // Namespace LifeV

#endif // HAVE_HDF5

#endif // EXPORTER_HDF5_H
