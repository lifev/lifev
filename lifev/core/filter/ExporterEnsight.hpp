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
 *  @contributor Tiziano Passerini <tiziano@mathcs.emory.edu>
 *  @maintainer Radu Popescu <radu.popescu@epfl.ch>
 */

#ifndef EXPORTER_ENSIGHT_H
#define EXPORTER_ENSIGHT_H

#include <lifev/core/filter/Exporter.hpp>

namespace LifeV
{

const UInt ensightOffset = 1; //the offset of the IDs in ensight files

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
    typedef MeshType                          mesh_Type;
    typedef Exporter<mesh_Type>               super;
    typedef typename super::meshPtr_Type      meshPtr_Type;
    typedef typename super::commPtr_Type      commPtr_Type;
    typedef typename super::vectorPtr_Type    vectorPtr_Type;
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
    ExporterEnsight (const GetPot& dfile, meshPtr_Type mesh, const std::string& prefix, const Int& procId );

    //! Constructor for ExporterEnsight
    ExporterEnsight (const GetPot& dfile, const std::string& prefix);
    //@}

    //! @name Public methods
    //@{

    //! Post-process the variables added to the list
    /*!
      @param time the solver time
    */
    void postProcess (const Real& time);

    //! Import data from previous simulations at a certain time
    /*!
      @param Time the time of the data to be imported

      Not yet implemented for ExporterEnsight
    */
    UInt importFromTime ( const Real& /*time*/ )
    {
        ERROR_MSG ( "ExporterEnsight::importFromTime has not yet been implemented.");
        return -1;
    }

    //! Import data from previous simulations and rebuild the internal time counters
    /*!
      @param importTime the time of the snapshot to be imported
      @param dt the time step, is used to rebuild the history up to now
    */
    void import (const Real& importTime, const Real& dt);

    //! Import data from previous simulations
    /*!
      @param importTime the time of the snapshot to be imported
    */
    void import (const Real& importTime);

    //! Read variable
    /*!
      @param dvar the ExporterData object
    */
    void readVariable (exporterData_Type& dvar)
    {
        super::readVariable (dvar);
    }

    //! Set the mesh and the processor id
    /*!
      @param mesh a pointer to the mesh
      @param procId the ID of the current process
    */
    void setMeshProcId ( const meshPtr_Type mesh, const Int& procId );

    //! temporary: the method should work form the Exporter class
    void exportPID (  meshPtr_Type /*meshPart*/, commPtr_Type comm, const bool /*binaryFormat*/ = false )
    {
        if ( !comm->MyPID() )
        {
            std::cerr << "  X-  exportPID is not working with Ensight" << std::endl;
        }
    }
    //@}

    //! @name Get methods
    //@{

    //! returns the type of the map to use for the VectorEpetra
    MapEpetraType mapType() const
    {
        return Repeated;
    }

    //@}

private:
    //! @name Private methods
    //@{
    //! Compose the .case file
    /*!
      @param Time the time of the snapshot
    */
    void writeCase (const Real& time);
    //! Compose the .geo file
    /*!
      @param geoFile the name of the file to be produced
    */
    void writeAsciiGeometry ( const std::string geoFile );
    //! The ASCII writer
    /*!
      @param dvar the ExporterData object
    */
    void writeAscii ( const exporterData_Type& dvar);
    //! The ASCII writer
    /*!
      @param dvar the ExporterData object
      @param suffix the file suffix (.scl or .vct)
    */
    void writeAsciiValues (const exporterData_Type& dvar, const std::string& suffix);
    //! Compose the "mesh" section of the .case file
    /*!
      @param casef the file object

      Specify whether or not the mesh is changing for each snapshot
    */
    void caseMeshSection (std::ofstream& casef);
    //! Compose the "variable" section of the .case file
    /*!
      @param casef the file object

      The file name will be different based on the regime and the kind of
      field (scalar / vector)
    */
    void caseVariableSection (std::ofstream& casef);
    //! Compose the "time" section of the .case file
    /*!
      @param casef the file object
      @param time the current time

      The file will contain the updated number of time steps and the complete
      list of time values
    */
    void caseTimeSection (std::ofstream& casef, const Real& time);
    //! Dump on file the IDs of the global DOFs associated to this process
    /*!
      @param filename the file name
    */
    void writeGlobalIDs (const std::string& filename);
    //! The ASCII reader
    /*!
      @param dvar the ExporterData object
    */
    void readAscii ( exporterData_Type& dvar );
    //! The ASCII reader
    /*!
      @param dvar the ExporterData object
      @param suffix the file suffix (.scl or .vct)
    */
    void readAsciiValues ( exporterData_Type& dvar, const std::string& suffix );
    //! Read from file and store in a vector a list of global IDs
    /*!
      @param filename the file name
      @param globalDOF the stored list of IDs
    */
    void readGlobalIDs ( const std::string& filename,
                         std::vector<Real>& globalDOF );
    //! The generic reader (specialization of the parent class method)
    /*!
      @param dvar the ExporterData object
    */
    void readScalar ( exporterData_Type& dvar )
    {
        readAsciiValues (dvar, ".scl");
    }
    //! The generic reader (specialization of the parent class method)
    /*!
      @param dvar the ExporterData object
    */
    void readVector ( exporterData_Type& dvar )
    {
        readAsciiValues (dvar, ".vct");
    }
    //! initialize the internal data structures storing the ID of the current process
    void initProcId();
    //! Set the local-to-global map of DOFs
    /*!
      @param ltGNodesMap the local-to-global map

      ltGNodesMap[i] is the global ID of the i-th DOF owned by the current process
    */
    void setNodesMap ( std::vector<Int> ltGNodesMap )
    {
        M_ltGNodesMap = ltGNodesMap;
    }
    //! Build the local-to-global map of DOFs
    /*!
      @param ltGNodesMap the local-to-global map

      ltGNodesMap[i] is the global ID of the i-th DOF owned by the current process
    */
    void initNodesMap();
    //@}

    //! @name Private members
    //@{
    //! the counter of the number of steps processed
    UInt                        M_steps;
    //! the local-to-global map: ltGNodesMap[i] is the global ID of the i-th DOF
    //! owned by the current process
    std::vector<Int>            M_ltGNodesMap;
    //! the ID of the current process
    std::string                 M_me;
    //! a string label for the FE space
    std::string                 M_FEstr;
    //! a string label for the boundary FE space
    std::string                 M_bdFEstr;
    //! the number of local DOFs (per element)
    UInt                        M_nbLocalDof;
    //! the number of local boundary DOFs (per element)
    UInt                        M_nbLocalBdDof;
    //! are we performing the first post-processing operation?
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
ExporterEnsight<MeshType>::ExporterEnsight() :
    super(),
    M_steps (0),
    M_ltGNodesMap(),
    M_me(),
    M_firstTimeStep (true)
{
}

template<typename MeshType>
ExporterEnsight<MeshType>::ExporterEnsight (const GetPot& dfile, meshPtr_Type mesh, const std::string& prefix,
                                            const Int& procId)
    :
    super (dfile, prefix),
    M_steps (0),
    M_ltGNodesMap(),
    M_me(),
    M_firstTimeStep (true)
{
    this->setDataFromGetPot (dfile, "exporter");
    this->setMeshProcId (mesh, procId);
}

template<typename MeshType>
ExporterEnsight<MeshType>::ExporterEnsight (const GetPot& dfile, const std::string& prefix) :
    super (dfile, prefix),
    M_steps (0),
    M_ltGNodesMap(),
    M_me(),
    M_firstTimeStep (true)
{
    this->setDataFromGetPot (dfile, "exporter");
}

// =====================
// Public methods
// =====================

template<typename MeshType>
void ExporterEnsight<MeshType>::postProcess (const Real& time)
{
    // writing the geo file and the list of global IDs, but only upon the first instance
    if ( M_firstTimeStep )
    {
        if (!this->M_multimesh)
        {
            writeAsciiGeometry ( this->M_postDir + this->M_prefix + this->M_me + ".geo" );
        }

        writeGlobalIDs ( this->M_postDir + super::M_prefix + "_globalIDs" +
                         this->M_me + ".scl" );

        M_firstTimeStep = false;
    }
    // prepare the file postfix
    this->computePostfix();

    // the postfix will be full of stars, if this time step is not going to generate a snapshot
    std::size_t found ( this->M_postfix.find ( "*" ) );
    if ( found == std::string::npos )
    {
        if (!this->M_procId)
        {
            std::cout << "  X-  ExporterEnsight post-processing ...        " << std::flush;
        }
        LifeChrono chrono;
        chrono.start();
        for (typename super::dataVectorIterator_Type i = this->M_dataVector.begin();
                i != this->M_dataVector.end(); ++i)
        {
            // the "regime" attribute needs to be valid
            if ( i->regime() != exporterData_Type::NullRegime )
            {
                writeAscii (*i);
            }
            // if the solution is steady, we do not need to export it at each time step
            if (i->regime() == exporterData_Type::SteadyRegime)
            {
                i->setRegime ( exporterData_Type::NullRegime );
            }
        }
        // write an updated case file
        writeCase (time);

        // write an updated geo file, if needed
        if (this->M_multimesh)
        {
            writeAsciiGeometry ( this->M_postDir + this->M_prefix + this->M_postfix + this->M_me + ".geo" );
        }
        chrono.stop();
        if (!this->M_procId)
        {
            std::cout << "      done in " << chrono.diff() << " s." << std::endl;
        }
    }

}

template<typename MeshType>
void ExporterEnsight<MeshType>::import (const Real& startTime, const Real& dt)
{
    // dt is used to rebuild the history up to now
    Real time (startTime - this->M_timeIndex * dt);

    for ( UInt count (0); count < this->M_timeIndex; ++count)
    {
        this->M_timeSteps.push_back (time);
        ++this->M_steps;
        time += dt;
    }

    time += dt;

    import (time);

}

template<typename MeshType>
void ExporterEnsight<MeshType>::import (const Real& time)
{
    this->M_timeSteps.push_back (time);
    ++this->M_steps;

    // typedef std::list< exporterData_Type >::iterator Iterator;

    this->computePostfix();

    assert ( this->M_postfix.find ( "*" ) == std::string::npos );

    if (!this->M_procId)
    {
        std::cout << "  X-  ExporterEnsight importing ..." << std::endl;
    }

    LifeChrono chrono;
    chrono.start();
    for (typename super::dataVectorIterator_Type i = this->M_dataVector.begin(); i != this->M_dataVector.end(); ++i)
    {
        this->readVariable (*i);
    }
    chrono.stop();
    if (!this->M_procId)
    {
        std::cout << "      done in " << chrono.diff() << " s." << std::endl;
    }

}

template<typename MeshType>
void ExporterEnsight<MeshType>::setMeshProcId ( const meshPtr_Type mesh, const Int& procId )
{
    super::setMeshProcId ( mesh, procId );

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
            ERROR_MSG ( "FE not allowed in ExporterEnsight writer" );
            break;
    }

}

// ===================
// Get methods
// ===================

// ===================
// Private methods
// ===================

template <typename MeshType>
void ExporterEnsight<MeshType>::writeCase (const Real& time)
{
    std::string filename ( this->M_postDir + this->M_prefix + this->M_me + ".case" );
    std::ofstream casef ( filename.c_str() );
    ASSERT (casef.is_open(), "There is an error while opening " + filename );
    ASSERT (casef.good(), "There is an error while writing to " + filename );
    casef << "FORMAT\ntype: ensight\n";
    caseMeshSection (casef);
    caseVariableSection (casef);
    caseTimeSection (casef, time);
    casef.close();
}

template <typename MeshType>
void ExporterEnsight<MeshType>::writeAsciiGeometry (const std::string gFile)
{
    using std::setw;

    std::ofstream geoFile (gFile.c_str() );
    ASSERT (geoFile.is_open(), "There is an error while opening " + gFile );
    ID vertexNumber = this->M_mesh->numVertices();
    ID elementNumber = this->M_mesh->numElements();
    UInt part = 0;
    ASSERT (geoFile.good(), "There is an error while writing to " + gFile );
    geoFile << "Geometry file\nGenerated by LifeV\nnode id given\nelement id given\ncoordinates\n";
    geoFile.setf (std::ios::right | std::ios_base::scientific);
    geoFile.precision (5);
    ASSERT (geoFile.good(), "There is an error while writing to " + gFile );
    geoFile << setw (8) <<  vertexNumber << "\n";
    for (ID i = 0; i < vertexNumber; ++i)
    {
        ASSERT (geoFile.good(), "There is an error while writing to " + gFile );
        geoFile << setw (8) << i + ensightOffset;
        for (UInt icoor = 0; icoor < nDimensions; icoor++)
        {
            ASSERT (geoFile.good(), "There is an error while writing to " + gFile );
            geoFile << setw (12) << static_cast<float> (this->M_mesh->pointList (i).coordinatesArray() [icoor]);
        }
        ASSERT (geoFile.good(), "There is an error while opening " + gFile );
        geoFile << "\n";
    }

    ASSERT (geoFile.good(), "There is an error while writing to " + gFile );
    ++part;
    geoFile << "part" << setw (8) << part << "\nfull geometry\n"
            // elements
            << M_FEstr << "\n" << setw (8) << elementNumber << "\n";
    for (ID i = 0; i < elementNumber; ++i)
    {
        ASSERT (geoFile.good(), "There is an error while writing to " + gFile );
        geoFile << setw (8) << i + ensightOffset;
        for (ID j = 0; j < M_nbLocalDof; ++j)
        {
            ASSERT (geoFile.good(), "There is an error while writing to " + gFile );
            geoFile << setw (8) << this->M_mesh->element (i).point (j).localId() + ensightOffset;
        }
        ASSERT (geoFile.good(), "There is an error while writing to " + gFile );
        geoFile << "\n";

    }
    geoFile.close();
}

template <typename MeshType>
void ExporterEnsight<MeshType>::writeAscii (const exporterData_Type& dvar)
{

    switch ( dvar.fieldType() )
    {
        case exporterData_Type::ScalarField:
            writeAsciiValues (dvar, ".scl");
            break;
        case exporterData_Type::VectorField:
            writeAsciiValues (dvar, ".vct");
            break;
        default:
            ERROR_MSG ( "Unknown field type" )
            break;
    }

}

template <typename MeshType>
void ExporterEnsight<MeshType>::writeAsciiValues (const exporterData_Type& dvar, const std::string& suffix)
{
    using std::setw;

    std::ofstream exportFile;
    std::string filename;

    if ( dvar.regime() == exporterData_Type::SteadyRegime )
        filename = this->M_postDir + super::M_prefix + "_" + dvar.variableName() +
                   this->M_me + suffix;
    else
        filename = this->M_postDir + super::M_prefix + "_" + dvar.variableName() +
                   this->M_postfix + this->M_me + suffix;

    exportFile.open ( filename.c_str() );

    ASSERT (exportFile.is_open(), "There is an error while opening " + filename );

    UInt count = 0;

    const UInt size  = dvar.numDOF();
    const UInt start = dvar.start();
    const UInt vertexNumber = static_cast<UInt> (this->M_ltGNodesMap.size() );

    ASSERT (exportFile.good(), "There is an error while writing to " + filename );
    if ( suffix.compare (".vct") == 0 )
    {
        exportFile << "Vector per node\n";
    }
    else if ( suffix.compare (".scl") == 0 )
    {
        exportFile << "Scalar per node\n";
    }

    exportFile.setf (std::ios::right | std::ios_base::scientific);
    exportFile.precision (5);

    for (UInt i = 0; i < vertexNumber; ++i)
        for (UInt j = 0; j < dvar.fieldDim(); ++j)
        {
            const Int id = this->M_ltGNodesMap[i];
            ASSERT (exportFile.good(), "There is an error while writing to " + filename );
            exportFile << setw (12) << static_cast<float> (dvar (start + j * size + id) ) ;
            ++count;
            if ( count == 6 )
            {
                ASSERT (exportFile.good(), "There is an error while writing to " + filename );
                exportFile << "\n";
                count = 0;
            }
        }
    ASSERT (exportFile.good(), "There is an error while writing to " + filename );
    exportFile << std::endl;

    exportFile.close();
}

template <typename MeshType>
void ExporterEnsight<MeshType>::writeGlobalIDs (const std::string& filename)
{
    using std::setw;

    std::ofstream globalIDsFile;

    globalIDsFile.open ( filename.c_str() );
    ASSERT (globalIDsFile.is_open(), "There is an error while opening " + filename );

    UInt count = 0;

    const UInt vertexNumber = static_cast<UInt> (this->M_ltGNodesMap.size() );
    ASSERT (globalIDsFile.good(), "There is an error while writing to " + filename );
    globalIDsFile << "Node global ID " << vertexNumber << "\n";

    for (UInt i = 0; i < vertexNumber; ++i)
    {
        const Int id = this->M_ltGNodesMap[i];
        ASSERT (globalIDsFile.good(), "There is an error while writing to " + filename );
        globalIDsFile << setw (12) << id ;
        ++count;
        if ( count == 6 )
        {
            ASSERT (globalIDsFile.good(), "There is an error while writing to " + filename );
            globalIDsFile << "\n";
            count = 0;
        }
    }
    ASSERT (globalIDsFile.good(), "There is an error while writing to " + filename );
    globalIDsFile << std::endl;
    globalIDsFile.close();
}

template <typename MeshType>
void ExporterEnsight<MeshType>::caseMeshSection (std::ofstream& casef)
{
    ASSERT (casef.good(), "There is an error while writing to file" );
    casef << "GEOMETRY\n";
    if ( this->M_multimesh )
    {
        std::string stars (".");
        for (UInt cc (0); cc < this->M_timeIndexWidth; ++cc)
        {
            stars += "*";
        }

        ASSERT (casef.good(), "There is an error while writing to file" );
        casef << "model: 1 " + this->M_prefix + stars << this->M_me << ".geo change_coords_only\n";
    }
    else
    {
        ASSERT (casef.good(), "There is an error while writing to file" );
        casef << "model: 1 " + this->M_prefix + this->M_me + ".geo\n";
    }
}

template <typename MeshType>
void ExporterEnsight<MeshType>::caseVariableSection (std::ofstream& casef)
{
    ASSERT (casef.good(), "There is an error while writing to file" );
    casef << "VARIABLE\n";
    std::string aux, str;
    for (typename super::dataVectorIterator_Type i = this->M_dataVector.begin(); i != this->M_dataVector.end(); ++i)
    {
        if ( i->regime() == exporterData_Type::SteadyRegime )
        {
            str = "";
        }
        else
        {
            std::string stars (".");
            for (UInt cc (0); cc < this->M_timeIndexWidth; ++cc)
            {
                stars += "*";
            }

            str = stars;
        }
        aux = i->variableName() + " " + super::M_prefix + "_" + i->variableName();
        switch ( i->fieldType() )
        {
            case exporterData_Type::ScalarField:
                ASSERT (casef.good(), "There is an error while writing to file" );
                casef << "scalar per node: 1 " +  aux + str << this->M_me << ".scl\n";
                break;
            case exporterData_Type::VectorField:
                ASSERT (casef.good(), "There is an error while writing to file" );
                casef << "vector per node: 1 " +  aux +  str << this->M_me << ".vct\n";
                break;
        }
    }
}

template <typename MeshType>
void ExporterEnsight<MeshType>::caseTimeSection (std::ofstream& casef, const Real& time)
{
    ASSERT (casef.good(), "There is an error while writing to file" );
    this->M_timeSteps.push_back (time);
    ++this->M_steps;
    casef << "TIME\ntime set: 1\nnumber of steps: " <<  this->M_steps << "\n"
          << "filename start number: " << this->M_timeIndexStart << "\n"
          << "filename increment: 1\ntime values:\n";

    UInt count = 0;

    typedef std::list<Real>::const_iterator Iterator;
    for (Iterator i = this->M_timeSteps.begin(); i != this->M_timeSteps.end(); ++i)
    {
        ASSERT (casef.good(), "There is an error while writing to file" );
        casef << *i << " " ;
        ++count;
        if ( count == 6)
        {
            ASSERT (casef.good(), "There is an error while writing to file" );
            casef << "\n";
            count = 0;
        }
    }
}

template <typename MeshType>
void ExporterEnsight<MeshType>::readAscii (exporterData_Type& dvar)
{

    switch ( dvar.fieldType() )
    {
        case exporterData_Type::ScalarField:
            writeAsciiValues (dvar, ".scl");
            break;
        case exporterData_Type::VectorField:
            writeAsciiValues (dvar, ".vct");
            break;
        default:
            ERROR_MSG ("Unsupported field type!");
            break;
    }

}

template <typename MeshType> void ExporterEnsight<MeshType>::readAsciiValues (exporterData_Type& dvar,
                                                                              const std::string& suffix)
{
    ASSERT ( this->M_numImportProc, "The number of pieces to be loaded was not specified." );

    // this vector lists the global IDs of DOFs listed in each piece
    std::vector<Real> globalDOF;

    // Each processor will read all the files, and fill just its own component of the vectors
    for ( UInt iProc = 0; iProc < this->M_numImportProc; ++iProc )
    {
        // build the postfix for the file corresponding to part iProc
        std::ostringstream index;
        index.fill ( '0' );
        index << std::setw (1) << "." ;
        index << std::setw (3) << iProc;

        // fill globalDOF, the list of DOF IDs on which we are going to operate
        readGlobalIDs ( this->M_postDir + super::M_prefix + "_globalIDs" + index.str() + ".scl", globalDOF );

        // open the file with the vector field to be imported
        const std::string filename ( this->M_postDir + super::M_prefix + "_" + dvar.variableName() +
                                     this->M_postfix + index.str() + suffix );
        std::ifstream importFile ( filename.c_str() );

        // debugging step
        if (!this->M_procId)
        {
            std::cout << "\tfile " << filename << std::endl;
        }

        ASSERT (importFile.is_open(), "There is an error while reading " + filename );
        // discard the header of the file
        std::string line;
        ASSERT (importFile.good(), "There is an error while reading from " + filename );
        getline ( importFile, line ); // the line looks like this: "Vector/Scalar per node"

        // parameters to access ExporterData structures
        const UInt size  = dvar.numDOF();
        const UInt start = dvar.start();

        // loop over the DOFs listed in the current piece
        for (UInt i = 0; i < globalDOF.size(); ++i)
        {
            // extract the global ID of each DOF
            const Int id = globalDOF[i];
            // we are working with vectorial fields
            for (UInt j = 0; j < dvar.fieldDim(); ++j)
            {
                // this is the value of the field, to be imported
                Real value (0);

                ASSERT (importFile.good(), "There is an error while reading from " + filename );
                importFile >> value;

                // do the actual import only if the global ID belongs to the current process
                if ( dvar.feSpacePtr()->map().map (Repeated)->MyGID ( id ) )
                {
                    dvar (start + j * size + id) = value;
                }
            }
        }

        ASSERT (!importFile.fail(), "There is an error while reading " + filename );
        importFile.close();
    }
}

template<typename MeshType>
void ExporterEnsight<MeshType>::initProcId()
{
    std::ostringstream index;
    index.fill ( '0' );
    if (this->M_procId >= 0)
    {
        index << std::setw (1) << "." ;
        index << std::setw (3) << this->M_procId;
    }
    M_me = index.str();
}

template<typename MeshType>
void ExporterEnsight<MeshType>::initNodesMap()
{
    UInt vertexNumber = this->M_mesh->numVertices();
    M_ltGNodesMap.resize (vertexNumber);
    for (UInt i = 0; i < vertexNumber; ++i)
    {
        M_ltGNodesMap[i] = this->M_mesh->pointList ( i ).id();
    }
}

template<typename MeshType>
void ExporterEnsight<MeshType>::readGlobalIDs ( const std::string& filename,
                                                std::vector<Real>& globalDOF )
{
    std::ifstream globalIDsFile ( filename.c_str() );

    if (!this->M_procId)
    {
        std::cout << "\tfile " << filename << std::endl;
    }

    ASSERT (globalIDsFile.is_open(), "There is an error while opening " + filename );

    UInt vertexNumber;
    globalDOF.resize (0);

    // file parsing: line by line
    std::string line;
    std::stringstream parseLine;

    ASSERT (globalIDsFile.good(), "There is an error while reading " + filename );
    getline ( globalIDsFile, line ); // the line looks like this: "Node global ID N"
    parseLine.str (line);
    std::string trashcan;
    parseLine >> trashcan >> trashcan >> trashcan >> vertexNumber;

    for ( UInt iNode = 0; iNode < vertexNumber; ++iNode )
    {
        Real gID (0);
        ASSERT (globalIDsFile.good(), "There is an error while reading " + filename );
        globalIDsFile >> gID;
        globalDOF.push_back ( gID );
    }

    ASSERT (!globalIDsFile.fail(), "There is an error while reading " + filename );
    globalIDsFile.close();
}

} // Namespace LifeV

#endif // EXPORTER_ENSIGHT_H
