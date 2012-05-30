//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2009-2010 EPFL, Politecnico di Milano

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
 *  @file
 *  @brief This file provides an interface for post-processing with ensight
 *
 *  @author M.A. Fernandez
 *  @author C. Prud'homme
 *  @author S. Deparis
 *  @date 1-10-2005
 *
 *  @maintainer Radu Popescu <radu.popescu@epfl.ch>
 */

#ifndef EXPORTER_ENSIGHT_H
#define EXPORTER_ENSIGHT_H

#include <lifev/core/filter/Exporter.hpp>

namespace LifeV
{

const int ensightOffset = 1; //the offset of the IDs in ensight files

/**
 * @class ExporterEnsight
 * @brief ExporterEnsight data exporter
 */
template<typename MeshType>
class ExporterEnsight : public Exporter<MeshType>
{

public:
    //! @name Public typedefs
    //@{
    typedef MeshType mesh_Type;
    typedef Exporter<MeshType> super;
    typedef typename super::meshPtr_Type  meshPtr_Type;
    typedef typename super::vectorPtr_Type vectorPtr_Type;
    typedef typename super::exporterData_Type exporterData_Type;
    //@}

    //! @name Constructors and destructor
    //@{
    //! Default constructor
    ExporterEnsight();


    //! Constructor for ExporterEnsight
    /*!
      @param dfile the GetPot data file where you must provide and [ensight] section with:
      "start" (start index for filenames 0 for 000, 1 for 001 etc.),
      "save" (how many time steps per posptrocessing)
      "multimesh" (=true if the mesh has to be saved at each post-processing step)

      @param mesh the mesh

      @param the prefix for the case file (ex. "test" for test.case)

      @param the procId determines de CPU id. if negative, it ussemes there is only one processor
    */
    ExporterEnsight(const GetPot& dfile, meshPtr_Type mesh, const std::string& prefix, const Int& procId );

    //! Constructor for ExporterEnsight
    ExporterEnsight(const GetPot& dfile, const std::string& prefix);
    //@}

    //! @name Public methods
    //@{

    //! Post-process the variables added to the list
    /*!
      @param time the solver time
    */
    void postProcess(const Real& time);

    //! Import data from previous simulations at a certain time
    /*!
      @param Time the time of the data to be imported

      Not yet implemented for ExporterEnsight
    */
    UInt importFromTime( const Real& /*time*/ ) { assert(false); return 0; }

    //! Import data from previous simulations
    /*!
      @param time the solver time

      dt is used to rebuild the history up to now
    */
    void import(const Real& startTime, const Real& dt);

    //! Read  only last timestep
    void import(const Real& startTime);

    //! Read variable
    void readVariable(exporterData_Type& dvar) {super::readVariable(dvar);}

    //! Set the mesh and the processor id
    void setMeshProcId( const meshPtr_Type mesh, const Int& procId );

    //! temporary: the method should work form the Exporter class
    void exportPID ( MeshPartitioner<MeshType> const & /*meshPart*/ )
    {
        std::cerr << "  X-  exportPID is not working with VTK" << std::endl;
    }
    //@}

    //! @name Get methods
    //@{

    //! returns the type of the map to use for the VectorEpetra
    MapEpetraType mapType() const;

    //@}

private:
    //! @name Private methods
    //@{
    void defineShape();

    void writeCase(const Real& time);
    void writeAsciiGeometry( const std::string geoFile );

    void writeAscii(const exporterData_Type& dvar);
    void writeAsciiValues(const exporterData_Type& dvar, const std::string& suffix);
    void caseMeshSection(std::ofstream& casef);
    void caseVariableSection(std::ofstream& casef);
    void caseTimeSection(std::ofstream& casef, const Real& time);
    void writeGlobalIDs(const std::string& filename);

    void readAscii( exporterData_Type& dvar );
    void readAsciiValues( exporterData_Type& dvar, const std::string& suffix );
    void readGlobalIDs( const std::string& filename,
                        std::vector<Real>& globalDOF );

    void readScalar( exporterData_Type& dvar ) { readAsciiValues(dvar, ".scl"); }
    void readVector( exporterData_Type& dvar ) { readAsciiValues(dvar, ".vct"); }

    void initProcId();
    void setNodesMap( std::vector<Int> ltGNodesMap );
    void initNodesMap();
    //@}

    //! @name Private members
    //@{
    UInt                        M_steps;
    std::vector<Int>            M_ltGNodesMap;
    std::string                 M_me;

    std::string                 M_FEstr;
    std::string                 M_bdFEstr;
    UInt                        M_nbLocalDof;
    UInt                        M_nbLocalBdDof;

    bool                        M_firstTimeStep;
    //@}
};


// ==============
// Implementation
// ==============

// ==============
// Constructors
// ==============

template<typename MeshType>
ExporterEnsight<MeshType>::ExporterEnsight():
        super(),
        M_steps(0),
        M_ltGNodesMap(),
        M_me(),
        M_firstTimeStep(true)
{
}

template<typename MeshType>
ExporterEnsight<MeshType>::ExporterEnsight(const GetPot& dfile, meshPtr_Type mesh, const std::string& prefix,
                                           const Int& procId)
        :
        super(dfile, prefix),
        M_steps(0),
        M_ltGNodesMap(),
        M_me(),
        M_firstTimeStep(true)
{
    this->setDataFromGetPot(dfile,"exporter");
    this->setMeshProcId(mesh,procId);
}

template<typename MeshType>
ExporterEnsight<MeshType>::ExporterEnsight(const GetPot& dfile, const std::string& prefix):
        super(dfile,prefix),
        M_steps(0),
        M_ltGNodesMap(),
        M_me(),
        M_firstTimeStep(true)
{
    this->setDataFromGetPot(dfile,"exporter");
}

// =====================
// Public methods
// =====================

template<typename MeshType>
void ExporterEnsight<MeshType>::postProcess(const Real& time)
{
    // typedef std::list< exporterData_Type >::iterator Iterator;
    if( M_firstTimeStep )
    {
        if (!this->M_multimesh)
            writeAsciiGeometry( this->M_postDir + this->M_prefix + this->M_me + ".geo" );

        writeGlobalIDs( this->M_postDir + super::M_prefix + "_globalIDs" +
                        this->M_me + ".scl" );

        M_firstTimeStep = false;
    }
    this->computePostfix();

    std::size_t found( this->M_postfix.find( "*" ) );
    if ( found == string::npos )
    {
        if (!this->M_procId) std::cout << "  X-  ExporterEnsight post-processing ...        " << std::flush;
        LifeChrono chrono;
        chrono.start();
        for (typename super::dataVectorIterator_Type i=this->M_dataVector.begin(); i != this->M_dataVector.end(); ++i)
        {
            if ( i->regime() != exporterData_Type::NullRegime )
                writeAscii(*i);
            if (i->regime() == exporterData_Type::SteadyRegime)
                i->setRegime( exporterData_Type::NullRegime );
        }
        writeCase(time);

        if (this->M_multimesh)
            writeAsciiGeometry( this->M_postDir + this->M_prefix + this->M_postfix + this->M_me+".geo" );
        chrono.stop();
        if (!this->M_procId) std::cout << "      done in " << chrono.diff() << " s." << std::endl;
    }

}

template<typename MeshType>
void ExporterEnsight<MeshType>::import(const Real& startTime, const Real& dt)
{
    // dt is used to rebuild the history up to now
    Real time(startTime - this->M_timeIndex * dt);

    for ( UInt count(0); count < this->M_timeIndex; ++count)
    {
        this->M_timeSteps.push_back(time);
        ++this->M_steps;
        time += dt;
    }

    time += dt;

    import(time);

}

template<typename MeshType>
void ExporterEnsight<MeshType>::import(const Real& time)
{
    this->M_timeSteps.push_back(time);
    ++this->M_steps;

    // typedef std::list< exporterData_Type >::iterator Iterator;

    this->computePostfix();

    assert( this->M_postfix.find( "*" ) == string::npos );

    if (!this->M_procId) std::cout << "  X-  ExporterEnsight importing ..."<< std::endl;

    LifeChrono chrono;
    chrono.start();
    for (typename super::dataVectorIterator_Type i=this->M_dataVector.begin(); i != this->M_dataVector.end(); ++i)
    {
        this->readVariable(*i);
    }
    chrono.stop();
    if (!this->M_procId) std::cout << "      done in " << chrono.diff() << " s." << std::endl;

}

template<typename MeshType>
void ExporterEnsight<MeshType>::setMeshProcId( const meshPtr_Type mesh, const Int& procId )
{
    super::setMeshProcId( mesh, procId );

    initNodesMap();
    initProcId();

    typedef typename MeshType::elementShape_Type elementShape_Type;

    switch ( elementShape_Type::S_shape )
    {
    case TETRA:
        M_FEstr = "tetra4";
        M_bdFEstr = "tria3";
        M_nbLocalBdDof = 3;
        M_nbLocalDof = 4;
        break;
    case HEXA:
        M_FEstr = "hexa8";
        M_bdFEstr = "quad4";
        M_nbLocalBdDof = 4;
        M_nbLocalDof = 8;
        break;
    case TRIANGLE:
        M_FEstr = "tria3";
        M_bdFEstr = "bar2";
        M_nbLocalBdDof = 2;
        M_nbLocalDof = 3;
        break;
    case QUAD:
        M_FEstr = "quad4";
        M_bdFEstr = "bar2";
        M_nbLocalBdDof = 4;
        M_nbLocalDof = 3;
        break;
    default:
        ERROR_MSG( "FE not allowed in ExporterEnsight writer" );
    }

}

// ===================
// Get methods
// ===================

template<typename MeshType>
MapEpetraType ExporterEnsight<MeshType>::mapType() const
{
    return Repeated;
}

// ===================
// Private methods
// ===================

template <typename MeshType>
void ExporterEnsight<MeshType>::writeCase(const Real& time)
{
    std::ofstream casef( (this->M_postDir + this->M_prefix + this->M_me + ".case").c_str() );
    casef << "FORMAT\n";
    casef << "type: ensight\n";
    caseMeshSection(casef);
    caseVariableSection(casef);
    caseTimeSection(casef,time);
}

template <typename MeshType>
void ExporterEnsight<MeshType>::writeAsciiGeometry(const std::string gFile)
{
    using std::setw;

    std::ofstream geoFile(gFile.c_str() );
    ID vertexNumber = this->M_mesh->numVertices();
    ID elementNumber = this->M_mesh->numElements();
    UInt part=0;
    geoFile << "Geometry file\n";
    geoFile << "Geometry file\n";
    geoFile << "node id given\n";
    geoFile << "element id given\n";
    geoFile << "coordinates\n";
    geoFile.setf(std::ios::right | std::ios_base::scientific);
    geoFile.precision(5);
    geoFile << setw(8) <<  vertexNumber << "\n";
    for (ID i=0; i < vertexNumber; ++i)
    {
        geoFile << setw(8) << i + ensightOffset;
        for (UInt icoor=0; icoor<nDimensions; icoor++)
            geoFile << setw(12) << float(this->M_mesh->pointList(i).coordinatesArray()[icoor]);
        geoFile << "\n";
    }

    geoFile<< "part";

    ++part;
    geoFile << setw(8) << part << "\n";
    geoFile << "full geometry\n";
    // elements
    geoFile << M_FEstr << "\n";
    geoFile << setw(8) << elementNumber << "\n";
    for (ID i=0; i < elementNumber; ++i)
    {
        geoFile << setw(8) << i + ensightOffset;
        for (ID j=0; j< M_nbLocalDof; ++j)
        {
            geoFile << setw(8) << this->M_mesh->element(i).point(j).localId() + ensightOffset;
        }
        geoFile << "\n";

    }
}

template <typename MeshType>
void ExporterEnsight<MeshType>::writeAscii(const exporterData_Type& dvar)
{

    switch ( dvar.fieldType() )
    {
    case exporterData_Type::ScalarField:
        writeAsciiValues(dvar,".scl");
        break;
    case exporterData_Type::VectorField:
        writeAsciiValues(dvar,".vct");
        break;
    }

}

template <typename MeshType>
void ExporterEnsight<MeshType>::writeAsciiValues(const exporterData_Type& dvar, const std::string& suffix)
{
    using std::setw;

    std::ofstream exportFile;

    if ( dvar.regime() == exporterData_Type::SteadyRegime )
        exportFile.open( (this->M_postDir + super::M_prefix + "_" + dvar.variableName() +
                          this->M_me + suffix).c_str() );
    else
        exportFile.open( (this->M_postDir + super::M_prefix + "_" + dvar.variableName() +
                          this->M_postfix + this->M_me + suffix).c_str() );

    UInt count=0;

    UInt size  = dvar.numDOF();
    UInt start = dvar.start();
    UInt vertexNumber = static_cast<UInt> (this->M_ltGNodesMap.size());

    if( suffix.compare(".vct") == 0 )
        exportFile<<"Vector per node\n";
    else if( suffix.compare(".scl") == 0 )
        exportFile<<"Scalar per node\n";

    exportFile.setf(std::ios::right | std::ios_base::scientific);
    exportFile.precision(5);

    for (UInt i=0; i<vertexNumber; ++i)
        for (UInt j=0; j< dvar.fieldDim(); ++j)
        {
            Int id = this->M_ltGNodesMap[i];
            exportFile << setw(12) << float(dvar(start + j * size + id)) ;
            ++count;
            if ( count == 6 )
            {
                exportFile << "\n";
                count=0;
            }
        }
    exportFile << std::endl;
}

template <typename MeshType>
void ExporterEnsight<MeshType>::writeGlobalIDs(const std::string& filename)
{
    using std::setw;

    std::ofstream globalIDsFile;

    globalIDsFile.open( filename.c_str() );

    UInt count=0;

    UInt vertexNumber = static_cast<UInt> (this->M_ltGNodesMap.size());
    globalIDsFile<<"Node global ID "<<vertexNumber<<"\n";

    for (UInt i=0; i<vertexNumber; ++i)
    {
        Int id = this->M_ltGNodesMap[i];
        globalIDsFile << setw(12) << id ;
        ++count;
        if ( count == 6 )
        {
            globalIDsFile << "\n";
            count=0;
        }
    }
    globalIDsFile << std::endl;
    globalIDsFile.close();
}

template <typename MeshType>
void ExporterEnsight<MeshType>::caseMeshSection(std::ofstream& casef)
{
    casef << "GEOMETRY\n";
    if ( this->M_multimesh )
    {
        std::string stars(".");
        for (UInt cc(0); cc<this->M_timeIndexWidth; ++cc) stars+="*";

        casef << "model: 1 " + this->M_prefix + stars << this->M_me << ".geo change_coords_only\n";
    }
    else
        casef << "model: 1 " + this->M_prefix + this->M_me + ".geo\n";
}

template <typename MeshType>
void ExporterEnsight<MeshType>::caseVariableSection(std::ofstream& casef)
{
    // typedef typename std::list< exporterData_Type >::const_iterator Iterator;
    casef << "VARIABLE\n";
    std::string aux, str;
    for (typename super::dataVectorIterator_Type i=this->M_dataVector.begin(); i != this->M_dataVector.end(); ++i)
    {
        if ( i->regime() == exporterData_Type::SteadyRegime )
            str = "";
        else
        {
            std::string stars(".");
            for (UInt cc(0); cc<this->M_timeIndexWidth; ++cc) stars+="*";

            str = stars;
        }
        aux = i->variableName() + " " + super::M_prefix + "_" + i->variableName();
        switch ( i->fieldType() )
        {
        case exporterData_Type::ScalarField:
            casef << "scalar per node: 1 " +  aux + str << this->M_me << ".scl\n";
            break;
        case exporterData_Type::VectorField:
            casef << "vector per node: 1 " +  aux +  str << this->M_me << ".vct\n";
            break;
        }
    }
}

template <typename MeshType>
void ExporterEnsight<MeshType>::caseTimeSection(std::ofstream& casef, const Real& time)
{
    this->M_timeSteps.push_back(time);
    ++this->M_steps;
    casef << "TIME\n";
    casef << "time set: 1\n";
    casef << "number of steps: " <<  this->M_steps << "\n";
    casef << "filename start number: " << this->M_timeIndexStart << "\n";
    casef << "filename increment: 1\n";
    casef << "time values:\n";

    UInt count=0;

    typedef std::list<Real>::const_iterator Iterator;
    for (Iterator i=this->M_timeSteps.begin(); i != this->M_timeSteps.end(); ++i)
    {
        casef << *i << " " ;
        ++count;
        if ( count == 6)
        {
            casef <<"\n";
            count = 0;
        }
    }
}

template <typename MeshType>
void ExporterEnsight<MeshType>::readAscii(exporterData_Type& dvar)
{

    switch ( dvar.fieldType() )
    {
    case exporterData_Type::ScalarField:
        writeAsciiValues(dvar,".scl");
        break;
    case exporterData_Type::VectorField:
        writeAsciiValues(dvar,".vct");
        break;
    }

}

template <typename MeshType> void ExporterEnsight<MeshType>::readAsciiValues(exporterData_Type& dvar, const std::string& suffix)
{
    ASSERT( this->M_numImportProc, "The number of pieces to be loaded was not specified." );

    // this vector lists the global IDs of DOFs listed in each piece
    std::vector<Real> globalDOF;

    // Each processor will read all the files, and fill just its own component of the vectors
    for( UInt iProc = 0; iProc < this->M_numImportProc; ++iProc )
    {
        // build the postfix for the file corresponding to part iProc
        std::ostringstream index;
        index.fill( '0' );
        index << std::setw(1) << "." ;
        index << std::setw(3) << iProc;

        // fill globalDOF, the list of DOF IDs on which we are going to operate
        readGlobalIDs( this->M_postDir + super::M_prefix + "_globalIDs" + index.str() + ".scl", globalDOF );

        // open the file with the vector field to be imported
        std::string filename( this->M_postDir + super::M_prefix + "_" + dvar.variableName() +
                              this->M_postfix + index.str() + suffix );
        std::ifstream importFile( filename.c_str() );

        // debugging step
        if (!this->M_procId) std::cout << "\tfile "<< filename << std::endl;

        ASSERT(importFile.good(), std::stringstream("There is an error while reading " +
                                                    filename).str().c_str() );
        // parameters to access ExporterData structures
        UInt size  = dvar.numDOF();
        UInt start = dvar.start();

        // discard the header of the file
        std::string trashcan;
        importFile >> trashcan >> trashcan >> trashcan;

        // loop over the DOFs listed in the current piece
        for (UInt i=0; i<globalDOF.size(); ++i)
        {
            // extract the global ID of each DOF
            Int id = globalDOF[i];
            // we are working with vectorial fields
            for (UInt j=0; j<dvar.fieldDim(); ++j)
            {
                // this is the value of the field, to be imported
                Real value(0);
                ASSERT(importFile.good(), std::stringstream("There is an error while reading " +
                                                            filename).str().c_str() );
                importFile >> value;

                // do the actual import only if the global ID belongs to the current process
                if( dvar.feSpacePtr()->map().map(Repeated)->MyGID( id ) )
                {
                    dvar(start + j*size + id) = value;
                }
            }
        }

        ASSERT(!importFile.fail(), std::stringstream("There is an error while reading " +
                                                     filename).str().c_str() );
        importFile.close();
    }
}

template<typename MeshType>
void ExporterEnsight<MeshType>::initProcId()
{
    std::ostringstream index;
    index.fill( '0' );
    if (this->M_procId >=0)
    {
        index << std::setw(1) << "." ;
        index << std::setw(3) << this->M_procId;
    }
    M_me = index.str();
}

template<typename MeshType>
void ExporterEnsight<MeshType>::setNodesMap( std::vector<Int> ltGNodesMap )
{
    M_ltGNodesMap = ltGNodesMap;
}

template<typename MeshType>
void ExporterEnsight<MeshType>::initNodesMap()
{
    UInt vertexNumber = this->M_mesh->numVertices();
    M_ltGNodesMap.resize(vertexNumber);
    for (UInt i=0; i<vertexNumber; ++i)
    {
        M_ltGNodesMap[i] = this->M_mesh->pointList( i ).id();
    }
}

template<typename MeshType>
void ExporterEnsight<MeshType>::readGlobalIDs( const std::string& filename,
                                               std::vector<Real>& globalDOF )
{
    std::ifstream globalIDsFile( filename.c_str() );
    globalDOF.resize(0);

    if (!this->M_procId) std::cout << "\tfile "<< filename << std::endl;

    ASSERT(globalIDsFile.good(), std::stringstream("There is an error while reading " +
                                                   filename).str().c_str() );

    UInt vertexNumber; // = static_cast<UInt> (this->M_ltGNodesMap.size()); // will become this->M_numDOF, TP 5/2011


    // file parsing: line by line
    std::string line;
    std::stringstream parseLine;

    getline( globalIDsFile, line ); // "Node global ID N"
    parseLine.str(line);
    std::string trashcan;
    parseLine >> trashcan >> trashcan >> trashcan >> vertexNumber;

    for( UInt iNode = 0; iNode < vertexNumber; ++iNode )
    {
        Real gID(0);
        globalIDsFile >> gID;
        globalDOF.push_back( gID );
        // std::cout << globalDOF[iNode] << " " << std::flush;
    }
    // std::cout << std::endl;

    ASSERT(!globalIDsFile.fail(), std::stringstream("There is an error while reading " +
                                                    filename).str().c_str() );
    globalIDsFile.close();
}

} // Namespace LifeV

#endif // EXPORTER_ENSIGHT_H
