//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
  @file
  @brief This file provides the class  Hdf5exporter for post-processing with hdf5

  @author Simone Deparis <simone.deparis@epfl.ch>
  @author Radu Popescu <radu.popescu@epfl.ch>
  @date 11-11-2008
 */

#ifndef HDF5EXPORTER_H
#define HDF5EXPORTER_H 1

#include <vector>

#include <lifeconfig.h>

#ifndef HAVE_HDF5
#warning warning you should reconfigure with --with-hdf5=... flag

#else
#include <life/lifecore/util_string.hpp>
#include <life/lifefilters/exporter.hpp>
#include <EpetraExt_HDF5.h>
#include <Epetra_MultiVector.h>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include <EpetraExt_DistArray.h>
#include <Epetra_IntVector.h>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/refFEScalar.hpp>
#include <Epetra_Comm.h>
#include <sstream>
#include <string>
#include <map>

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
template<typename Mesh>
class Hdf5exporter : public Exporter<Mesh>
{

public:

    typedef Exporter<Mesh> super;
    typedef typename super::mesh_ptrtype mesh_ptrtype;
    typedef typename super::vector_rawtype vector_type;
    typedef typename super::vector_ptrtype vector_ptrtype;

    typedef EpetraExt::HDF5 hdf5_type;
    typedef boost::shared_ptr<hdf5_type> hdf5_ptrtype;
    typedef std::vector<std::vector<int> > graph_type;
    typedef boost::shared_ptr<graph_type> graph_ptrtype;
    typedef boost::shared_ptr<std::vector<mesh_ptrtype> > serial_mesh_ptrtype;


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor for Hdf5exporter
    Hdf5exporter();

    //! Constructor for Hdf5exporter
    /*!
       @param dfile the GetPot data file where you must provide an [exporter] section with:
          "start"     (start index for sections in the hdf5 data structure 0 for 000, 1 for 001 etc.),
          "save"      (how many time steps per postprocessing)
          "multimesh" ( = true if the mesh has to be saved at each post-processing step)
       @param mesh the mesh
       @param the prefix for the case file (ex. "test" for test.case)
       @param the procId determines de CPU id. if negative, it ussemes there is only one processor
    */
    Hdf5exporter(const GetPot& dfile, mesh_ptrtype mesh, const std::string& prefix, const int& procId);

    //! Constructor for Hdf5exporter without prefix and procID
    /*!
       @param dfile the GetPot data file where you must provide an [exporter] section with:
          "start"     (start index for sections in the hdf5 data structure 0 for 000, 1 for 001 etc.),
          "save"      (how many time steps per postprocessing)
          "multimesh" ( = true if the mesh has to be saved at each post-processing step)
       @param mesh the mesh
    */
    Hdf5exporter(const GetPot& dfile, const std::string& prefix);

    //! Destructor for Hdf5exporter
    ~Hdf5exporter() {}

    //@}


    //! @name Methods
    //@{
    virtual void postProcess(const Real& time);

    //! Import data from previous simulations at a certain time
    /*!
       @param Time the time of the data to be imported
       @return number of iteration corresponding at the time step
     */
    UInt importFromTime( const Real& Time );

    //! Import data from previous simulations at a certain time
    /*!
       @param Time the time of the data to be imported
       @return the simulation time corresponding to the iteration
     */
    double importFromIter( const UInt& );

    //! Import data from previous simulations
    /*!
       @param time the solver time
    */
    void import(const Real& /*Tstart*/, const Real& /*dt*/); // dt is used to rebuild the history up to now

    //! Import data from previous simulations
    /*!
       @param time the solver time
    */
    void import(const Real& /*Tstart*/);

    //! Close the Hdf5 file
    /*!
         Close the HDF5 file.
     */
    void CloseFile() {M_HDF5->Close();}

    //@}


    //! @name Get Methods
    //@{

    //! returns the type of the map to use for the EpetraVector
    EpetraMapType mapType() const;

    void rd_var( ExporterData& dvar);

    //@}

protected:

    //! @name Private Methods
    //@{

    //! Define the shape of the elements
    void defineShape();
    //! write empty xdmf file
    void M_wr_initXdmf();
    //! append to xdmf file
    void M_wr_Xdmf(const Real& time);
    //! save position and write closing lines
    void M_wr_closeLinesXdmf();
    //! remove closing lines
    void M_wr_removeCloseLinesXdmf();

    void M_wr_topology  ( std::ofstream& xdmf );
    void M_wr_geometry  ( std::ofstream& xdmf );
    void M_wr_attributes( std::ofstream& xdmf );
    void M_wr_scalar_datastructure  ( std::ofstream& xdmf, const ExporterData& dvar );
    void M_wr_vector_datastructure  ( std::ofstream& xdmf, const ExporterData& dvar );

    void M_wr_var(const ExporterData& dvar);
    void M_wr_scalar(const ExporterData& dvar);
    void M_wr_vector(const ExporterData& dvar);

    void M_wr_geo();

    void readScalar( ExporterData& dvar);
    void readVector( ExporterData& dvar);

    //@}

    hdf5_ptrtype      M_HDF5;
    std::ofstream     M_xdmf;

    const std::string M_closingLines;
    std::streampos    M_closingLinesPosition;
    std::string       M_outputFileName;

};



// ===================================================
// Constructors
// ===================================================
template<typename Mesh>
Hdf5exporter<Mesh>::Hdf5exporter():
        super               (),
        M_HDF5              (),
        M_closingLines      ( "\n    </Grid>\n\n  </Domain>\n</Xdmf>\n"),
        M_outputFileName    ( "noninitialisedFileName" )
{
}

template<typename Mesh>
Hdf5exporter<Mesh>::Hdf5exporter(const GetPot& dfile, mesh_ptrtype mesh, const std::string& prefix,
                                 const int& procId) :
        super               ( dfile, prefix ),
        M_HDF5              (),
        M_closingLines      ( "\n    </Grid>\n\n  </Domain>\n</Xdmf>\n"),
        M_outputFileName    ( "noninitialisedFileName" )
{
    setMeshProcId( mesh, procId );
}

template<typename Mesh>
Hdf5exporter<Mesh>::Hdf5exporter(const GetPot& dfile, const std::string& prefix):
        super               ( dfile, prefix ),
        M_HDF5              (),
        M_closingLines      ( "\n    </Grid>\n\n  </Domain>\n</Xdmf>\n"),
        M_outputFileName    ( "noninitialisedFileName" )
{
}

// ===================================================
// Methods
// ===================================================

template<typename Mesh>
void Hdf5exporter<Mesh>::postProcess(const Real& time)
{
    if ( M_HDF5.get() == 0)
    {
        M_HDF5.reset(new hdf5_type(this->M_listData.begin()->storedArray()->Comm()));
        M_outputFileName=this->M_prefix+".h5";
        M_HDF5->Create(this->M_postDir+M_outputFileName);

        // write empty xdmf file
        M_wr_initXdmf();

        if (!this->M_multimesh)
        {
            M_wr_geo(); // see also M_wr_geometry
            M_HDF5->Flush();
        }
    }

    typedef std::list< ExporterData >::const_iterator Iterator;

    this->computePostfix();

    if ( this->M_postfix != "*****" )
    {
        if (!this->M_procId) std::cout << "  x-  HDF5 post-processing ...        " << std::flush;
        Chrono chrono;
        chrono.start();
        for (Iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i)
        {
            M_wr_var(*i);
        }
        // pushing time
        this->M_timeSteps.push_back(time);

        M_wr_Xdmf(time);

        if (this->M_multimesh)
        {
            M_wr_geo(); // see also M_wr_geometry
        }

        chrono.stop();

        // Write to file without closing the file
        M_HDF5->Flush();

        if (!this->M_procId) std::cout << "         done in " << chrono.diff() << " s." << std::endl;
    }
}

template<typename Mesh>
UInt Hdf5exporter<Mesh>::importFromTime( const Real& Time )
{
    // Container for the time and the postfix
    std::pair< Real, int > SelectedTimeAndPostfix;
    if ( !this->M_procId )
    {
        // Open the xmf file
        std::ifstream xmfFile;
        xmfFile.open( ( this->M_postDir + this->M_prefix + ".xmf" ).c_str(), std::ios::in );

        // Vector of TimeStep
        std::vector< std::pair< Real, int > > TimeAndPostfix;
        if ( xmfFile.is_open() )
        {
            // Define some variables
            std::string line;
            std::vector<std::string> stringsVector;

            // Read one-by-one all the lines of the file
            while ( !xmfFile.eof() )
            {
                std::getline( xmfFile, line, '\n' );

                // If the line begin with "<!-- Time " it is the beginning of a new block
                if ( !line.compare( 0, 10, "<!-- Time " ) )
                {
                    boost::split( stringsVector, line, boost::is_any_of( " " ) );
                    TimeAndPostfix.push_back( std::make_pair( string2number( stringsVector[2] ), string2number( stringsVector[4] ) ) );
                }
            }
        }
        xmfFile.close();

        // Find the closest time step
        SelectedTimeAndPostfix = TimeAndPostfix.front();
        for ( std::vector< std::pair< Real, int > >::const_iterator i = TimeAndPostfix.begin(); i < TimeAndPostfix.end() ; ++i )
            if ( std::abs( SelectedTimeAndPostfix.first - Time ) >= std::abs( (*i).first - Time ) )
                SelectedTimeAndPostfix = *i;

        //std::cout << "  x-  HDF5 import from time " << SelectedTimeAndPostfix.first << " iteration " << SelectedTimeAndPostfix.second << std::endl;
    }

    this->M_listData.begin()->storedArray()->Comm().Broadcast( &SelectedTimeAndPostfix.second, 1, 0 );
    this->M_count = SelectedTimeAndPostfix.second;
    this->computePostfix();

    // Importing
    if ( !this->M_procId )
        std::cout << "  x-  HDF5 importing ...                       "<< std::flush;

    Chrono chrono;
    chrono.start();
    for ( std::list< ExporterData >::iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i )
        this->rd_var(*i);

    chrono.stop();
    if ( !this->M_procId )
        std::cout << "done in " << chrono.diff() << " s. (Time " << SelectedTimeAndPostfix.first
                  << ", Iteration " << SelectedTimeAndPostfix.second << " )" << std::endl;

    return static_cast <UInt> ( SelectedTimeAndPostfix.second );
}



template<typename Mesh>
double Hdf5exporter<Mesh>::importFromIter( const UInt& iter )
{
    // Container for the time and the postfix
    std::pair< Real, int > SelectedTimeAndPostfix;
    if ( !this->M_procId )
    {
        // Open the xmf file
        std::ifstream xmfFile;
        xmfFile.open( ( this->M_postDir + this->M_prefix + ".xmf" ).c_str(), std::ios::in );

        // Vector of TimeStep
        std::vector< std::pair< Real, int > > TimeAndPostfix;

        if ( xmfFile.is_open() )
        {
            // Define some variables
            std::string line;
            std::vector<std::string> stringsVector;

            // Read one-by-one all the lines of the file
            while ( !xmfFile.eof() )
            {
                std::getline( xmfFile, line, '\n' );

                // If the line begin with "<!-- Time " it is the beginning of a new block
                if ( !line.compare( 0, 10, "<!-- Time " ) )
                {
                    boost::split( stringsVector, line, boost::is_any_of( " " ) );
                    TimeAndPostfix.push_back( std::make_pair( string2number( stringsVector[2] ), string2number( stringsVector[4] ) ) );
                }
            }
        }

        xmfFile.close();

        // Find the closest time step
        SelectedTimeAndPostfix = TimeAndPostfix.front();
        bool found             = false;

        for ( std::vector< std::pair< Real, int > >::const_iterator i = TimeAndPostfix.begin(); i < TimeAndPostfix.end()  ; ++i )
        {
            if ( i->second == iter )
            {
                SelectedTimeAndPostfix = *i;
                found = true;
                break;
            }
        }

        ASSERT(found, "Selected iteration not found");

    }



    //std::cout << "  x-  HDF5 import from time " << SelectedTimeAndPostfix.first << " iteration " << SelectedTimeAndPostfix.second << std::endl;
    this->M_listData.begin()->storedArray()->Comm().Broadcast( &SelectedTimeAndPostfix.second, 1, 0 );
    this->M_count = SelectedTimeAndPostfix.second;

    std::ostringstream index;
    index.fill('0');

    index << std::setw(5) << this->M_count;
    this->M_postfix = "." + index.str();

    //    this->computePostfix();

    // Importing
    if ( !this->M_procId )
        std::cout << "  x-  HDF5 importing iteration "
                  << index.str()
                  << " at time " << SelectedTimeAndPostfix.first
                  << " ... " << std::flush;

    Chrono chrono;
    chrono.start();
    for ( std::list< ExporterData >::iterator i = this->M_listData.begin(); i != this->M_listData.end(); ++i )
    {
        this->rd_var(*i);
    }
    chrono.stop();

    if ( !this->M_procId )
        std::cout << "done in " << chrono.diff() << " s. (Time " << SelectedTimeAndPostfix.first
                  << ", Iteration " << SelectedTimeAndPostfix.second << " )" << std::endl;

    return SelectedTimeAndPostfix.first;
}


template<typename Mesh>
void Hdf5exporter<Mesh>::import(const Real& Tstart, const Real& dt)
{
    // dt is used to rebuild the history up to now
    Real time(Tstart - this->M_count*dt);

    for ( UInt count(0); count < this->M_count; ++count )
    {
        this->M_timeSteps.push_back(time);
        time += dt;
    }

    time += dt;

    import(time);
}

template<typename Mesh>
void Hdf5exporter<Mesh>::import(const Real& time)
{
    if ( M_HDF5.get() == 0)
    {
        M_HDF5.reset(new hdf5_type(this->M_listData.begin()->storedArray()->Comm()));
        M_HDF5->Open(this->M_postDir+this->M_prefix+".h5"); //!! Simone
    }

    this->M_timeSteps.push_back(time);

    this->computePostfix();

    assert( this->M_postfix != "*****" );

    if (!this->M_procId) std::cout << "  x-  HDF5 importing ..."<< std::endl;

    Chrono chrono;
    chrono.start();
    for (std::list< ExporterData >::iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i)
    {
        this->rd_var(*i); ///!!! Simone
    }
    chrono.stop();
    if (!this->M_procId) std::cout << "      done in " << chrono.diff() << " s." << std::endl;
}

// ===================================================
// Get Methods
// ===================================================
template<typename Mesh>
EpetraMapType Hdf5exporter<Mesh>::mapType() const
{
    return Unique;
}

// ===================================================
// Private Methods
// ===================================================
template<typename Mesh>
void Hdf5exporter<Mesh>::defineShape()
{
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_var(const ExporterData& dvar)
{

    switch ( dvar.type() )
    {
    case ExporterData::Scalar:
        M_wr_scalar(dvar);
        break;
    case ExporterData::Vector:
        M_wr_vector(dvar);
        break;
    }
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_scalar(const ExporterData& dvar)
{
    /* Examples:
    M_HDF5->Write("map-" + toString(Comm.NumProc()), Map);
    M_HDF5->Write("matrix", Matrix);
    M_HDF5->Write("LHS", LHS);
    M_HDF5->Write("RHS", RHS);
    */

    UInt size  = dvar.size();
    UInt start = dvar.start();

    EpetraMap subMap(dvar.storedArray()->BlockMap(), start, size);
    vector_type subVar(subMap);
    subVar.subset(*dvar.storedArray(),start);

    std::string varname (dvar.variableName()+ this->M_postfix); // see also in M_wr_attributes
    bool writeTranspose (true);
    M_HDF5->Write(varname, subVar.getEpetraVector(), writeTranspose );
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_vector(const ExporterData& dvar)
{

    UInt size  = dvar.size();
    UInt start = dvar.start();

    using namespace boost;

    // solution array has to be reordered and stored in a Multivector.
    // Using auxiliary arrays:
    //shared_array<double*>                   ArrayOfPointers(new double*[nDimensions]);
    double **                                 ArrayOfPointers(new double*[nDimensions]);
    shared_array< shared_ptr<vector_type> > ArrayOfVectors (new shared_ptr<vector_type>[nDimensions]);

    int MyLDA;


    // Building subsets (new vectors) of the original data and than taking a view of them to
    // build a multivector.
    // Note: the contents of ArrayOfPointers[0,1,2] must not be deleted explicitly, since their
    // content belongs to ArrayOfVectors[0,1,2].
    // ArrayOfVectors[0,1,2] are deleted when ArrayOfVectors is destroyed

    for (UInt d ( 0 ); d < nDimensions; ++d)
    {
        EpetraMap subMap(dvar.storedArray()->BlockMap(), start+d*size, size);
        ArrayOfVectors[d].reset(new  vector_type(subMap));
        ArrayOfVectors[d]->subset(*dvar.storedArray(),start+d*size);

        ArrayOfVectors[d]->getEpetraVector().ExtractView(&ArrayOfPointers[d], &MyLDA);
    }

    EpetraMap subMap(dvar.storedArray()->BlockMap(), start, size);
    Epetra_MultiVector multiVector(View, *subMap.getMap(Unique), ArrayOfPointers, nDimensions);


    bool writeTranspose (true);
    std::string varname (dvar.variableName() + this->M_postfix); // see also in M_wr_attributes
    M_HDF5->Write(varname, multiVector, writeTranspose);

    delete[] ArrayOfPointers;
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_geo()
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
      see also M_wr_geometry
    */

    /* Steps:
       generate local  connections and local coordinates (not repeated)
       write out
    */

    // We need a map ,but it's not always possible to use that from the variables
    // (if we write out a P0 variable)
    // We build a map for the connections based on the element numbers and for the points we fake a P1 map

    ASSERT (this->M_listData.size() > 0 , "hdf5exporter: ListData is empty");

    // Connections
    // Need to use elements not dofs for this map. Recover local element lists

    std::vector<int> elementList;
    elementList.reserve(this->M_mesh->numElements()*Mesh::ElementShape::numPoints);
    for (ID i=1; i <= this->M_mesh->numElements(); ++i)
    {
        typename Mesh::ElementType const& element (this->M_mesh->element(i));
        UInt lid=(i-1)*Mesh::ElementShape::numPoints;
        for (ID j=1; j<= Mesh::ElementShape::numPoints; ++j, ++lid)
        {
            elementList[lid] = (element.id()-1)*Mesh::ElementShape::numPoints+(j-1);
        }
    }

    Epetra_Map connectionsMap(this->M_mesh->numGlobalElements()*Mesh::ElementShape::numPoints,
                              this->M_mesh->numElements()*Mesh::ElementShape::numPoints,
                              &elementList[0],
                              0, this->M_listData.begin()->storedArray()->Comm());

    Epetra_IntVector connections(connectionsMap);
    for (ID i=1; i <= this->M_mesh->numElements(); ++i)
    {
        typename Mesh::ElementType const& element (this->M_mesh->element(i));
        UInt lid=(i-1)*Mesh::ElementShape::numPoints;
        for (ID j=1; j<= Mesh::ElementShape::numPoints; ++j, ++lid)
        {
            connections[lid] = element.point(j).id();
        }
    }

    this->M_listData.begin()->storedArray()->Comm().Barrier();

    // this offset is needed by hdf5 since it starts numbering from 0
    //int const hdf5Offset(this->M_listData.begin()->storedArray()->BlockMap().IndexBase());
    int const hdf5Offset(0);

    // Points

    // Build a map for linear elements, even though the origianl FE might be P0
    // This gives the right map for the coordinate arrays

    EpetraMap subMap;
    switch ( Mesh::ElementShape::Shape )
    {
    case TETRA:
    {
        const RefFE & refFEP1 = feTetraP1;
        EpetraMap tmpMapP1(refFEP1, *this->M_mesh,
                           this->M_listData.begin()->storedArray()->getMap_ptr()->CommPtr());
        subMap = tmpMapP1;
        break;
    }
    case HEXA:
    {
        const RefFE & refFEQ1 = feHexaQ1;
        EpetraMap tmpMapQ1(refFEQ1, *this->M_mesh,
                           this->M_listData.begin()->storedArray()->getMap_ptr()->CommPtr());
        subMap = tmpMapQ1;
        break;
    }
    case LINE:
    {
        const RefFE & refFEP11D = feSegP1;
        EpetraMap tmpMapQ11D(refFEP11D, *this->M_mesh,
                             this->M_listData.begin()->storedArray()->getMap_ptr()->CommPtr());
        subMap = tmpMapQ11D;
        break;
    }
    default:
        ERROR_MSG( "FE not allowed in HDF5 Exporter" );

    }
    //    EpetraMap subMap(refFE, *this->M_mesh,
    //		     const_cast<Epetra_Comm&>(this->M_listData.begin()->storedArray()->Comm()));

    EpetraVector pointsX(subMap);
    EpetraVector pointsY(subMap);
    EpetraVector pointsZ(subMap);

    int gid;
    for (ID i=1; i <= this->M_mesh->numVertices(); ++i)
    {
        typename Mesh::PointType const& point (this->M_mesh->pointList(i));
        gid = point.id() - hdf5Offset;

        bool insertedX(true);
        bool insertedY(true);
        bool insertedZ(true);

        insertedX = insertedX && pointsX.checkAndSet(gid, point.x());
        insertedY = insertedY && pointsY.checkAndSet(gid, point.y());
        insertedZ = insertedZ && pointsZ.checkAndSet(gid, point.z());
    }

    // Now we are ready to export the vectors to the hdf5 file

    std::string pointsXVarname("PointsX");
    std::string pointsYVarname("PointsY");
    std::string pointsZVarname("PointsZ");
    std::string connectionsVarname("Connections");

    if (this->M_multimesh)
    {
        connectionsVarname  += this->M_postfix; // see also in M_wr_topology
        pointsXVarname      += this->M_postfix; // see also in M_wr_geometry
        pointsYVarname      += this->M_postfix; // see also in M_wr_geometry
        pointsZVarname      += this->M_postfix; // see also in M_wr_geometry
    }

    M_HDF5->Write(connectionsVarname, connections);
    // bool writeTranspose (true);
    M_HDF5->Write(pointsXVarname, pointsX.getEpetraVector(), true);
    M_HDF5->Write(pointsYVarname, pointsY.getEpetraVector(), true);
    M_HDF5->Write(pointsZVarname, pointsZ.getEpetraVector(), true);
}


// write empty xdmf file
template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_initXdmf()
{
    if (this->M_procId == 0)
    {
        M_xdmf.open( (this->M_postDir+this->M_prefix+".xmf").c_str(), std::ios_base::out );

        M_xdmf <<
        "<?xml version=\"1.0\" ?>\n" <<
        "<!DOCTYPE Xdmf SYSTEM \"" << this->M_prefix << ".xdmf\" [\n" <<
        "<!ENTITY DataFile \"" << this->M_prefix << ".h5\">\n" <<
        "]>\n" <<
        "<!-- " << this->M_prefix << ".h5 is generated by LifeV -->\n" <<
        "<Xdmf>\n" <<
        "  <Domain Name=\"" << this->M_prefix << "\">\n" <<
        "    <Grid Name=\"" << this->M_prefix << "Grid\" GridType=\"Collection\" CollectionType=\"Temporal\">\n" <<
        "\n";

        M_wr_closeLinesXdmf();
    }
}

// save position and write closing lines
template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_closeLinesXdmf()
{
    // save position
    M_closingLinesPosition = M_xdmf.tellp();

    // write closing lines
    M_xdmf << M_closingLines;
    M_xdmf.flush();

}

// remove closing lines
template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_removeCloseLinesXdmf()
{
    M_xdmf.seekp(M_closingLinesPosition);
}


template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_Xdmf(const Real& time)
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

      The Int vectors are for the connections, the double ones for the Vertices
    */

    if (this->M_procId == 0)
    {
        M_wr_removeCloseLinesXdmf();

        // write grid with time, topology, geometry and attributes
        // NOTE: The first line (<!-- Time t Iteration i -->) is used in function importFromTime.
        //       Check compatibility after any change on it!
        M_xdmf <<
        "<!-- Time " << time << " Iteration " << this->M_postfix.substr(1,5) << " -->\n" <<
        "    <Grid Name=\"Mesh " << time << "\">\n" <<
        "      <Time TimeType=\"Single\" Value=\"" << time << "\" />\n";
        M_wr_topology(M_xdmf);
        M_wr_geometry(M_xdmf);
        M_wr_attributes(M_xdmf);

        M_xdmf << "\n"
        "    </Grid>\n\n";



        // write closing lines
        M_wr_closeLinesXdmf();
    }
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_topology  ( std::ofstream& xdmf )
{
    std::string FEstring;

    switch ( Mesh::ElementShape::Shape )
    {
    case TETRA:
        FEstring = "Tetrahedron";
        break;
    case HEXA:
        FEstring = "Hexahedron";
        break;
    case LINE:
        FEstring = "Polyline";
        break;
    default:
        ERROR_MSG( "FE not allowed in HDF5 Exporter" );
    }

    xdmf <<
    "      <Topology\n" <<
    "         Type=\"" << FEstring <<"\"\n" <<
    "         NumberOfElements=\"" << this->M_mesh->numGlobalElements() << "\"\n" <<
    "         BaseOffset=\"1\">\n" <<
    "         <DataStructure Format=\"HDF\"\n" <<
    "                        Dimensions=\""<< this->M_mesh->numGlobalElements() << " " << this->M_mesh->numLocalVertices() << "\"\n" <<
    "                        DataType=\"Int\"\n" <<
    "                        Precision=\"8\">\n" <<
    "             " << M_outputFileName << ":/Connections/Values\n" <<
    "         </DataStructure>\n" <<
    "      </Topology>\n";
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_geometry  ( std::ofstream& xdmf )
{

    std::string postfix_string;

    // see also in postProcess
    if (this->M_multimesh)
        postfix_string = this->M_postfix;
    else
        postfix_string = "";


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

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_attributes  ( std::ofstream& xdmf )
{

    // Loop on the variables to output
    for (std::list< ExporterData >::const_iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i)
    {
        xdmf <<
        "\n      <Attribute\n" <<
        "         Type=\"" << i->typeName() << "\"\n" <<
        "         Center=\"" << i->whereName() << "\"\n" <<
        "         Name=\"" << i->variableName()<<"\">\n";

        switch ( i->type() )
        {
        case ExporterData::Scalar:
            M_wr_scalar_datastructure(xdmf, *i);
            break;
        case ExporterData::Vector:
            M_wr_scalar_datastructure(xdmf, *i);
            //M_wr_vector_datastructure(xdmf, *i);
            break;
        }

        xdmf <<
        "      </Attribute>\n";
    }
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_scalar_datastructure  ( std::ofstream& xdmf, const ExporterData& dvar )
{

    Int globalUnknowns (0);
    switch ( dvar.where() )
    {
    case ExporterData::Node:
        globalUnknowns = this->M_mesh->numGlobalVertices();
        break;
    case ExporterData::Cell:
        globalUnknowns = this->M_mesh->numGlobalElements();
        break;
    }

    // First: hyperslab definition, then description of the data
    xdmf <<

    "         <DataStructure ItemType=\"HyperSlab\"\n" <<
    "                        Dimensions=\"" << globalUnknowns << " " << dvar.typeDim() << "\"\n" <<
    "                        Type=\"HyperSlab\">\n" <<
    "           <DataStructure  Dimensions=\"3 2\"\n" <<
    "                           Format=\"XML\">\n" <<
    "               0    0\n" <<
    "               1    1\n" <<
    "               " << globalUnknowns << " " << dvar.typeDim() << "\n" <<
    "           </DataStructure>\n" <<

    "           <DataStructure  Format=\"HDF\"\n" <<
    "                           Dimensions=\"" << dvar.size() << " " << dvar.typeDim() << "\"\n" <<
    "                           DataType=\"Float\"\n" <<
    "                           Precision=\"8\">\n" <<
    "               " << M_outputFileName << ":/" << dvar.variableName() << this->M_postfix  <<"/Values\n" << // see also in M_wr_vector/scalar
    "           </DataStructure>\n" <<
    "         </DataStructure>\n";

}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_vector_datastructure  ( std::ofstream& xdmf, const ExporterData& dvar )
{


    string coord[3]={"X","Y","Z"}; // see also wr_vector

    xdmf <<
    "         <DataStructure ItemType=\"Function\"\n" <<
    "                        Dimensions=\"" << this->M_mesh->numGlobalVertices() << " " << dvar.typeDim() << "\"\n" <<
    "                        Function=\"JOIN($0 , $1, $2)\">\n";

    for (int i(0); i < dvar.typeDim(); ++i)
    {
        xdmf <<
        "           <DataStructure  Format=\"HDF\"\n" <<
        "                           Dimensions=\"" << this->M_mesh->numGlobalVertices() << " 1\"\n" <<
        "                           DataType=\"Float\"\n" <<
        "                           Precision=\"8\">\n" <<
        "               " << M_outputFileName << ":/" << dvar.variableName()<< coord[i] << this->M_postfix  <<"/Values\n" << // see also in M_wr_vector/scalar
        "           </DataStructure>\n";
    }

    xdmf <<
    "         </DataStructure>\n";


}

template <typename Mesh>
void Hdf5exporter<Mesh>::rd_var(ExporterData& dvar)
{
    if ( M_HDF5.get() == 0)
    {
        M_HDF5.reset(new hdf5_type(dvar.storedArray()->BlockMap().Comm()));
        M_HDF5->Open(this->M_postDir+this->M_prefix+".h5"); //!! Simone
    }
    super::rd_var(dvar);
}

template <typename Mesh>
void Hdf5exporter<Mesh>::readScalar(ExporterData& dvar)
{

    UInt size  = dvar.size();
    UInt start = dvar.start();

    EpetraMap subMap(dvar.storedArray()->BlockMap(), start, size);
    Epetra_MultiVector* subVar(0);

    std::string varname (dvar.variableName()); // see also in M_wr_attributes
    if (this->M_postfix!="")
    {
        varname += this->M_postfix;
    }
    bool readTranspose (true);
    M_HDF5->Read(varname, *subMap.getMap(this->mapType()), subVar, readTranspose);

    dvar.storedArray()->subset(*subVar, subMap, 0, start );

    delete subVar;

}

template <typename Mesh>
void Hdf5exporter<Mesh>::readVector( ExporterData& dvar)
{
    UInt size  = dvar.size();
    UInt start = dvar.start();

    using namespace boost;

    // solution array has first to be read has Multivector.

    // first read the multivector:
    EpetraMap subMap(dvar.storedArray()->BlockMap(), start, size);
    Epetra_MultiVector* subVar(0);

    bool readTranspose (true);
    std::string varname (dvar.variableName()); // see also in M_wr_attributes

    if (this->M_postfix!="")
    {
        varname += this->M_postfix;
    }

    M_HDF5->Read(varname, *subMap.getMap(this->mapType()), subVar, readTranspose);


    // then put back value in our EpetraVector

    for (UInt d ( 0 ); d < nDimensions; ++d)
    {
        dvar.storedArray()->subset(*subVar, subMap,  0, start+d*size, d );
    }

    delete subVar;
}

}
#endif

#endif
