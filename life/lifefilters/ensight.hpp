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

#ifndef ENSIGHT_H
#define ENSIGHT_H

#include <life/lifefilters/exporter.hpp>

namespace LifeV
{

/**
 * @class Ensight
 * @brief Ensight data exporter
 */
template<typename MeshType>
class Ensight : public Exporter<MeshType>
{

public:
    //! @name Public typedefs
    //@{
    typedef MeshType mesh_Type;
    typedef Exporter<MeshType> super;
    typedef typename super::meshPtr_Type  meshPtr_Type;
    typedef typename super::vectorPtr_Type vectorPtr_Type;
    //@}

    //! @name Constructors and destructor
    //@{
    //! Default constructor
    Ensight();


    //! Constructor for Ensight
    /*!
      @param dfile the GetPot data file where you must provide and [ensight] section with:
      "start" (start index for filenames 0 for 000, 1 for 001 etc.),
      "save" (how many time steps per posptrocessing)
      "multimesh" (=true if the mesh has to be saved at each post-processing step)

      @param mesh the mesh

      @param the prefix for the case file (ex. "test" for test.case)

      @param the procId determines de CPU id. if negative, it ussemes there is only one processor
    */
    Ensight(const GetPot& dfile, meshPtr_Type mesh, const std::string& prefix, const Int& procId );

    //! Constructor for Ensight
    Ensight(const GetPot& dfile, const std::string& prefix);
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

      Not yet implemented for Ensight
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
    void readVariable(ExporterData& dvar) {super::readVariable(dvar);}

    //! Set the mesh and the processor id
    void setMeshProcId( const meshPtr_Type mesh, const Int& procId );
    //@}

    //! @name Get methods
    //@{

    //! returns the type of the map to use for the EpetraVector
    EpetraMapType mapType() const;

    // DEPRECATED
    void __attribute__((__deprecated__)) rd_var(ExporterData& dvar) {super::rd_var(dvar);}

    //@}

private:
    //! @name Private methods
    //@{
    void defineShape();

    void writeCase(const Real& time);
    void writeAsciiGeometry( const std::string geoFile );

    void writeAscii(const ExporterData& dvar);
    void writeAsciiScalar(const ExporterData& dvar);
    void writeAsciiVector(const ExporterData& dvar);
    void caseMeshSection(std::ofstream& casef);
    void caseVariableSection(std::ofstream& casef);
    void caseTimeSection(std::ofstream& casef, const Real& time);

    void readScalar( ExporterData& dvar );
    void readVector( ExporterData& dvar );

    void initProcId();
    void setNodesMap( std::vector<Int> ltGNodesMap );
    void initNodesMap();
    //@}

    //! @name Private members
    //@{
    std::string M_importDir;
    UInt M_steps;
    std::vector<Int> M_ltGNodesMap;
    std::string M_me;

    std::string                 M_FEstr;
    std::string                 M_bdFEstr;
    UInt                        M_nbLocalDof;
    UInt                        M_nbLocalBdDof;
    //@}
};


// ==============
// Implementation
// ==============

// ==============
// Constructors
// ==============

template<typename MeshType>
Ensight<MeshType>::Ensight():
    super(),
    M_importDir("./"),
    M_steps(0),
    M_ltGNodesMap(),
    M_me()
{
}

template<typename MeshType>
Ensight<MeshType>::Ensight(const GetPot& dfile, meshPtr_Type mesh, const std::string& prefix,
                       const Int& procId)
    :
    super(dfile, prefix),
    M_importDir(dfile("exporter/import_dir", "./")),
    M_steps(0),
    M_ltGNodesMap(),
    M_me()
{
    this->setMeshProcId(mesh,procId);
}

template<typename MeshType>
Ensight<MeshType>::Ensight(const GetPot& dfile, const std::string& prefix):
    super(dfile,prefix),
    M_importDir(dfile("exporter/import_dir", "./")),
    M_steps(0),
    M_ltGNodesMap(),
    M_me()
{
}

// =====================
// Public methods
// =====================

template<typename MeshType>
void Ensight<MeshType>::postProcess(const Real& time)
{
    typedef std::list< ExporterData >::iterator Iterator;

    this->computePostfix();

    if ( this->M_postfix != "*****" )
    {
        if (!this->M_procId) std::cout << "  x-  Ensight post-processing ...        " << std::flush;
        Chrono chrono;
        chrono.start();
        for (Iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i)
        {
            if (i->steady() < 2 )
                writeAscii(*i);
            if (i->steady() == 1) i->setSteady(2);
        }
        writeCase(time);

        if (this->M_multimesh)
            writeAsciiGeometry( this->M_postDir + this->M_prefix + this->M_postfix + this->M_me+".geo" );
        chrono.stop();
        if (!this->M_procId) std::cout << "      done in " << chrono.diff() << " s." << std::endl;
    }
}

template<typename MeshType>
void Ensight<MeshType>::import(const Real& startTime, const Real& dt)
{
    // dt is used to rebuild the history up to now
    Real time(startTime - this->M_startIndex * dt);

    for ( UInt count(0); count < this->M_startIndex; ++count)
    {
        this->M_timeSteps.push_back(time);
        ++this->M_steps;
        time += dt;
    }

    time += dt;

    import(time);

}

template<typename MeshType>
void Ensight<MeshType>::import(const Real& time)
{
    this->M_timeSteps.push_back(time);
    ++this->M_steps;

    typedef std::list< ExporterData >::iterator Iterator;

    this->computePostfix();

    assert( this->M_postfix != "*****" );

    if (!this->M_procId) std::cout << "  x-  Ensight importing ..."<< std::endl;

    Chrono chrono;
    chrono.start();
    for (Iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i)
    {
        this->readVariable(*i);
    }
    chrono.stop();
    if (!this->M_procId) std::cout << "      done in " << chrono.diff() << " s." << std::endl;

}

template<typename MeshType>
void Ensight<MeshType>::setMeshProcId( const meshPtr_Type mesh, const Int& procId )
{
    super::setMeshProcId( mesh, procId );

    initNodesMap();
    initProcId();

    typedef typename MeshType::ElementShape ElementShape;

    switch ( ElementShape::S_shape )
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
        ERROR_MSG( "FE not allowed in Ensight writer" );
    }

    if (!this->M_multimesh)
        writeAsciiGeometry( this->M_postDir + this->M_prefix + this->M_me + ".geo" );
}

// ===================
// Get methods
// ===================

template<typename MeshType>
EpetraMapType Ensight<MeshType>::mapType() const
{
    return Repeated;
}

// ===================
// Private methods
// ===================

template <typename MeshType>
void Ensight<MeshType>::writeCase(const Real& time)
{
    std::ofstream casef( (this->M_postDir + this->M_prefix + this->M_me + ".case").c_str() );
    casef << "FORMAT\n";
    casef << "type: ensight\n";
    caseMeshSection(casef);
    caseVariableSection(casef);
    caseTimeSection(casef,time);
}

template <typename MeshType>
void Ensight<MeshType>::writeAsciiGeometry(const std::string gFile)
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
    for (ID i=1; i <= vertexNumber; ++i)
    {
        geoFile << setw(8) << i ;
        for (UInt icoor=0; icoor<nDimensions; icoor++)
            geoFile << setw(12) << float(this->M_mesh->pointList(i).coor()[icoor]);
        geoFile << "\n";
    }

    geoFile<< "part";

    ++part;
    geoFile << setw(8) << part << "\n";
    geoFile << "full geometry\n";
    // elements
    geoFile << M_FEstr << "\n";
    geoFile << setw(8) << elementNumber << "\n";
    for (ID i=1; i <= elementNumber; ++i)
    {
        geoFile << setw(8) << i ;
        for (ID j=1; j<= M_nbLocalDof; ++j)
        {
            geoFile << setw(8) << this->M_mesh->element(i).point(j).localId();
        }
        geoFile << "\n";

    }
}

template <typename MeshType>
void Ensight<MeshType>::writeAscii(const ExporterData& dvar)
{

    switch ( dvar.type() )
    {
    case ExporterData::Scalar:
        writeAsciiScalar(dvar);
        break;
    case ExporterData::Vector:
        writeAsciiVector(dvar);
        break;
    }

}

template <typename MeshType>
void Ensight<MeshType>::writeAsciiScalar(const ExporterData& dvar)
{
    using std::setw;

    std::ofstream scalarFile;

    if (dvar.steady() )
        scalarFile.open( (this->M_postDir + super::M_prefix + "_" + dvar.variableName() +
                          this->M_me + ".scl").c_str() );
    else
        scalarFile.open( (this->M_postDir + super::M_prefix + "_" + dvar.variableName() +
                          this->M_postfix + this->M_me + ".scl").c_str() );

    UInt count=0;

    UInt start = dvar.start();
    UInt vertexNumber = static_cast<UInt> (this->M_ltGNodesMap.size());
    scalarFile<<"Scalar per node\n";

    scalarFile.setf(std::ios::right | std::ios_base::scientific);
    scalarFile.precision(5);

    for (UInt i=0; i<vertexNumber; ++i)
    {
        Int id = this->M_ltGNodesMap[i];
        scalarFile << setw(12) << float(dvar(start + id)) ;
        ++count;
        if ( count == 6 )
        {
            scalarFile << "\n";
            count=0;
        }
    }
    scalarFile << std::endl;
}

template <typename MeshType> void Ensight<MeshType>::writeAsciiVector(const ExporterData& dvar)
{
    using std::setw;

    std::ofstream vectorFile;

    if (dvar.steady() )
        vectorFile.open( (this->M_postDir + super::M_prefix + "_" + dvar.variableName() +
                          this->M_me + ".vct").c_str() );
    else
        vectorFile.open( (this->M_postDir + super::M_prefix + "_" + dvar.variableName() +
                          this->M_postfix + this->M_me + ".vct").c_str() );

    UInt count=0;

    UInt size   = dvar.size();
    UInt start = dvar.start();
    UInt vertexNumber = static_cast<UInt> (this->M_ltGNodesMap.size());

    vectorFile<<"Vector per node\n";

    vectorFile.setf(std::ios::right | std::ios_base::scientific);
    vectorFile.precision(5);

    for (UInt i=0; i<vertexNumber; ++i)
        for (UInt j=0; j< nDimensions; ++j)
        {
            Int id = this->M_ltGNodesMap[i];
            vectorFile << setw(12) << float(dvar(start + j * size + id)) ;
            ++count;
            if ( count == 6 )
            {
                vectorFile << "\n";
                count=0;
            }
        }
    vectorFile << std::endl;
}

template <typename MeshType>
void Ensight<MeshType>::caseMeshSection(std::ofstream& casef)
{
    casef << "GEOMETRY\n";
    if ( this->M_multimesh )
        casef << "model: 1 " + this->M_prefix + ".*****"<< this->M_me << ".geo change_coords_only\n";
    else
        casef << "model: 1 " + this->M_prefix + this->M_me + ".geo\n";
}

template <typename MeshType>
void Ensight<MeshType>::caseVariableSection(std::ofstream& casef)
{
    typedef std::list< ExporterData >::const_iterator Iterator;
    casef << "VARIABLE\n";
    std::string aux, str;
    for (Iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i)
    {
        if (i-> steady() )
            str = "";
        else
            str = ".*****";
        aux = i->variableName() + " " + super::M_prefix + "_" + i->variableName();
        switch ( i->type() )
        {
        case ExporterData::Scalar:
            casef << "scalar per node: 1 " +  aux + str << this->M_me << ".scl\n";
            break;
        case ExporterData::Vector:
            casef << "vector per node: 1 " +  aux +  str << this->M_me << ".vct\n";
            break;
        }
    }
}

template <typename MeshType>
void Ensight<MeshType>::caseTimeSection(std::ofstream& casef, const Real& time)
{
    this->M_timeSteps.push_back(time);
    ++this->M_steps;
    casef << "TIME\n";
    casef << "time set: 1\n";
    casef << "number of steps: " <<  this->M_steps << "\n";
    casef << "filename start number: 0\n";
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
void Ensight<MeshType>::readScalar( ExporterData& dvar )
{

    std::string filename( M_importDir + super::M_prefix + "_" + dvar.variableName() +
                          this->M_postfix + this->M_me + ".scl" );
    std::ifstream scalarFile( filename.c_str() );

    if (!this->M_procId) std::cout << "\tfile "<< filename << std::endl;

    ASSERT(scalarFile.good(), std::stringstream("There is an error while reading " +
                                                filename).str().c_str() );

    UInt start = dvar.start();

    UInt vertexNumber = static_cast<UInt> (this->M_ltGNodesMap.size());

    std::string trashcan;

    scalarFile >> trashcan >> trashcan >> trashcan;

    scalarFile.setf(std::ios::right | std::ios_base::scientific);
    scalarFile.precision(5);

    std::map<Int,Int>::iterator iter;

    for (UInt i=0; i<vertexNumber; ++i)
    {
        ASSERT(scalarFile.good(), std::stringstream("There is an error while reading " +
                                                    filename).str().c_str() );

        Int id = this->M_ltGNodesMap[i];
        scalarFile.width(12);
        scalarFile >> dvar(start + id) ;
    }

    ASSERT(!scalarFile.fail(), std::stringstream("There is an error while reading " +
                                                 filename).str().c_str() );
}

template <typename MeshType> void Ensight<MeshType>::readVector(ExporterData& dvar)
{

    std::string filename( M_importDir + super::M_prefix + "_" + dvar.variableName() +
                          this->M_postfix + this->M_me + ".vct" );
    std::ifstream vectorFile( filename.c_str() );

    if (!this->M_procId) std::cout << "\tfile "<< filename << std::endl;

    ASSERT(vectorFile.good(), std::stringstream("There is an error while reading " +
                                                filename).str().c_str() );

    UInt size   = dvar.size();
    UInt start = dvar.start();
    UInt vertexNumber = static_cast<UInt> (this->M_ltGNodesMap.size());

    std::string trashcan;

    vectorFile >> trashcan >> trashcan >> trashcan;

    vectorFile.setf(std::ios::right | std::ios_base::scientific);
    vectorFile.precision(5);

    for (UInt i=0; i<vertexNumber; ++i)
        for (UInt j=0; j< nDimensions; ++j)
        {
            ASSERT(vectorFile.good(), std::stringstream("There is an error while reading " +
                                                        filename).str().c_str() );

            Int id = this->M_ltGNodesMap[i];
            vectorFile.width(12);
            vectorFile >> dvar(start + j*size + id) ;
        }

    ASSERT(!vectorFile.fail(), std::stringstream("There is an error while reading " +
                                                 filename).str().c_str() );
}

template<typename MeshType>
void Ensight<MeshType>::initProcId()
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
void Ensight<MeshType>::setNodesMap( std::vector<Int> ltGNodesMap )
{
    M_ltGNodesMap = ltGNodesMap;
}

template<typename MeshType>
void Ensight<MeshType>::initNodesMap()
{
    UInt vertexNumber = this->M_mesh->numVertices();
    M_ltGNodesMap.resize(vertexNumber);
    for (UInt i=0; i<vertexNumber; ++i)
    {
        Int id = this->M_mesh->pointList( i + 1 ).id();
        M_ltGNodesMap[i] = id;
    }
}

} // Namespace LifeV

#endif // ENSIGHT_H
