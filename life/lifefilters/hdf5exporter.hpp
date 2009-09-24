/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

/*!
  \file hdf5exporter.hpp
  \author S. Deparis
  \date 11/2007
  \version 1.0

  \brief This file provides an interface for post-processing with hdf5

  Usage: two steps
  - first: add the variables using addVariable
  - second: call postProcess( time );
*/

#ifndef _HDF5EXPORTER_H_
#define _HDF5EXPORTER_H_


#include <lifeconfig.h>
#ifndef HAVE_HDF5
#warning warning you should reconfigure with --with-hdf5=... flag
#else

#include <life/lifefilters/exporter.hpp>
#include <hdf5.h>
#include <EpetraExt_HDF5.h>
#include <Epetra_MultiVector.h>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <EpetraExt_DistArray.h>
#include <Epetra_IntVector.h>



namespace LifeV
{

/**
 * \class Hdf5exporter
 * \brief Hdf5exporter data exporter
 */
template<typename Mesh>
class Hdf5exporter : public Exporter<Mesh> {

    typedef Exporter<Mesh> super;
    typedef typename super::mesh_ptrtype  mesh_ptrtype;
    typedef typename super::vector_rawtype vector_type;
    typedef typename super::vector_ptrtype vector_ptrtype;

    typedef EpetraExt::HDF5         hdf5_type;
    typedef boost::shared_ptr<hdf5_type> hdf5_ptrtype;

public:

    /**
       Constructor for Hdf5exporter

       \param dfile the GetPot data file where you must provide and [hdf5exporter] section with:
       "start" (start index for filenames 0 for 000, 1 for 001 etc.),
       "save" (how many time steps per posptrocessing)
       "multimesh" (=true if the mesh has to be saved at each post-processing step)

       \param mesh the mesh

       \param the prefix for the case file (ex. "test" for test.case)

       \param the procId determines de CPU id. if negative, it ussemes there is only one processor
    */
    Hdf5exporter(const GetPot& dfile, mesh_ptrtype mesh, const std::string prefix, const int procId);

    Hdf5exporter(const GetPot& dfile, const std::string prefix);


    /**
       setters
    */

    void setMeshProcId( mesh_ptrtype mesh, int const procId );

    /**
       Post-porcess the variables added to the list

       \param time the solver time
    */
    void postProcess(const Real& time);

    /**
       Import data from previous simulations

       \param time the solver time
    */
    virtual void import(const Real& Tstart, const Real& dt) // dt is used to rebuild the history up to now
    {
        ASSERT(false,"Hdf5exporter::import not yet coded");
    }

    //! Read  only last timestep
    void import(const Real& Tstart)
    {
        ASSERT(false,"Hdf5exporter::import not yet coded");
    }

private:

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

    hdf5_ptrtype M_HDF5;
    std::ofstream M_xdmf;

    std::string const M_closingLines;
    std::streampos M_closingLinesPosition;
};


//
// Implementation
//
template<typename Mesh>
Hdf5exporter<Mesh>::Hdf5exporter(const GetPot& dfile, mesh_ptrtype mesh, const std::string prefix,
                                 const int procId)
    :
    super(dfile, mesh, prefix,procId),
    M_HDF5(),
    M_closingLines ( "\n    </Grid>\n\n  </Domain>\n</Xdmf>\n")
{
    setMeshProcId(mesh,procId);
}

template<typename Mesh>
Hdf5exporter<Mesh>::Hdf5exporter(const GetPot& dfile, const std::string prefix):
    super(dfile,prefix),
    M_HDF5(),
    M_closingLines ( "\n    </Grid>\n\n  </Domain>\n</Xdmf>\n")
{
}


template<typename Mesh>
void Hdf5exporter<Mesh>::setMeshProcId( mesh_ptrtype mesh , int const procId )
{

    initMeshProcId( mesh, procId );

    typedef typename Mesh::VolumeShape ElementShape;

    switch ( ElementShape::numPoints )
        {
        case 4:
            this->M_FEstr = "tetra4";
            this->M_bdFEstr = "tria3";
            this->M_nbLocalBdDof = 3;
            this->M_nbLocalDof = 4;
            break;
        case 8:
            this->M_FEstr = "hexa8";
            this->M_bdFEstr = "quad4";
            this->M_nbLocalBdDof = 4;
            this->M_nbLocalDof = 8;
            break;
        default:
            ERROR_MSG( "FE not allowed in Hdf5exporter writer" );
        }
}

template<typename Mesh>
void Hdf5exporter<Mesh>::postProcess(const Real& time)
{
    if ( M_HDF5.get() == 0)
        {
            M_HDF5.reset(new hdf5_type(this->M_listData.begin()->getPtr()->Comm()));
            M_HDF5->Create(this->M_prefix+".h5");

            // write empty xdmf file
            M_wr_initXdmf();


            if (!this->M_multimesh)
                M_wr_geo(); // see also M_wr_geometry
        }


    typedef std::list< ExporterData >::const_iterator Iterator;


    this->getPostfix();

    if ( this->M_postfix != "***" )
        {
            if (!this->M_procId) std::cout << "  x-  Hdf5exporter post-processing..."<< std::endl;
            Chrono chrono;
            chrono.start();
            for (Iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i)
                {
                    M_wr_var(*i);
                }
            // pushing time
            this->M_timeSteps.push_back(time);
            ++this->M_steps;

            M_wr_Xdmf(time);

            if (this->M_multimesh)
                M_wr_geo(); // see also M_wr_geometry
            chrono.stop();

            if (!this->M_procId) std::cout << "      done in " << chrono.diff() << " s." << std::endl;
        }
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_var(const ExporterData& dvar)
{

    switch( dvar.type() )
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

    UInt dim   = dvar.dim();
    UInt start = dvar.start();

    EpetraMap subMap(dvar.getPtr()->BlockMap(), start, dim);
    vector_type subVar(subMap);
    subVar.subset(*dvar.getPtr(),start);

    std::string varname (this->M_prefix+"_"+dvar.variableName()+ this->M_postfix); // see also in M_wr_attributes
    bool writeTranspose (true);
    M_HDF5->Write(varname, subVar.getEpetraVector(), writeTranspose );
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_vector(const ExporterData& dvar)
{


    /*
    UInt dim   = dvar.dim();
    UInt start = dvar.start();

    EpetraMap subMapDist(dvar.getPtr()->BlockMap(), start, nDimensions*dim);

    // Unfortunately, we need row-oriented write-out in hdf5/xmf (at least for now in paraview)
    // EpetraExt::DistArray< double > distVar(*subMapDist.getMap(Unique), nDimensions);
    EpetraVector varT(subMapDist);

    double* rowMajor;
    int     lda (nDimensions);
    int*    rowGID = varT.BlockMap().MyGlobalElements ();

    for (int i(0); i<nDimensions; ++i)
    {
        rowMajor = varT.getEpetraVector().Values() + i;

        EpetraMap subMap(dvar.getPtr()->BlockMap(), start+i*dim, dim);

        vector_type subVar(subMap);
        subVar.subset(*dvar.getPtr(), start+i*dim);

        double* colMajor (subVar.getEpetraVector().Values());

        for (int j(0); j < varT.getEpetraVector().MyLength () / 3  ; ++j, ++colMajor, rowMajor += lda)
            *rowMajor = *colMajor;

    }

    std::string varname (this->M_prefix+"_"+dvar.variableName() + this->M_postfix); // see also in M_wr_attributes
    M_HDF5->Write(varname, varT.getEpetraVector());


    UInt dim   = dvar.dim();
    UInt start = dvar.start();

    string coord[3]={"X","Y","Z"}; // see also wr_vector_datastructure
    for (int i(0); i<nDimensions; ++i)
    {
        EpetraMap subMap(dvar.getPtr()->BlockMap(), start+i*dim, dim);
        vector_type subVar(subMap);
        subVar.subset(*dvar.getPtr(), start+i*dim);

        std::string varname (this->M_prefix+"_"+dvar.variableName() + coord[i] + this->M_postfix); // see also in M_wr_attributes
        M_HDF5->Write(varname, subVar.getEpetraVector());
    }


    */

    UInt dim   = dvar.dim();
    UInt start = dvar.start();

    using namespace boost;

    shared_array<double*>                   ArrayOfPointers(new double*[nDimensions]);
    shared_array< shared_ptr<vector_type> > ArrayOfVectors (new shared_ptr<vector_type>[nDimensions]);
    int MyLDA;


    // Building subsets (new vectors) of the original data and than taking a view of them to
    // build a multivector.
    // Note: the contents of ArrayOfPointers[0,1,2] must not be deleted explicitly, since their
    // content belongs to ArrayOfVectors[0,1,2].
    // ArrayOfVectors[0,1,2] are deleted when ArrayOfVectors is destroyed

    for (UInt d ( 0 ); d < nDimensions; ++d)
    {
        EpetraMap subMap(dvar.getPtr()->BlockMap(), start+d*dim, dim);
        ArrayOfVectors[d].reset(new  vector_type(subMap));
        ArrayOfVectors[d]->subset(*dvar.getPtr(),start+d*dim);

        ArrayOfVectors[d]->getEpetraVector().ExtractView(&ArrayOfPointers[d], &MyLDA);
    }

    EpetraMap subMap(dvar.getPtr()->BlockMap(), start, dim);
    Epetra_MultiVector multiVector(View, *subMap.getMap(Unique), &ArrayOfPointers[0], nDimensions);


    bool writeTranspose (true);
    std::string varname (this->M_prefix+"_"+dvar.variableName() + this->M_postfix); // see also in M_wr_attributes
    M_HDF5->Write(varname, multiVector, writeTranspose);


}

template <typename Mesh>
void Hdf5exporter<Mesh>::Hdf5exporter<Mesh>::M_wr_geo()
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

    // I need a map, I suppose var[0] has one:
    //     dvar.getPtr()->BlockMap()
    // but this map probably includes more dimensions or P2 elements ...
    // ok: I assume that anyway it starts with the P1 elements of the first dimensions,
    //  I will discart the following entries.

    ASSERT (this->M_listData.size() > 0 , "hdf5exporter: ListData is empty");

    // Connections
    // linear map for volumes. we Will have to insert via local ids.
    Epetra_Map connectionsMap( this->M_nbLocalDof * this->M_mesh->numGlobalVolumes(),
                this->M_nbLocalDof * this->M_mesh->numVolumes(),
                0, this->M_listData.begin()->getPtr()->Comm());

    Epetra_IntVector connections(connectionsMap);

    int* connPtr;
    connections.ExtractView(&connPtr);

    for (ID i=1; i <= this->M_mesh->numVolumes(); ++i)
        {
            typename Mesh::VolumeType const& vol (this->M_mesh->volume(i));
            for (ID j=1; j<= this->M_nbLocalDof; ++j, ++connPtr)
                {
                    *connPtr = vol.point(j).id();
                }
        }


    this->M_listData.begin()->getPtr()->Comm().Barrier();

    // this offset is needed by hdf5 since it starts numbering from 0
    //int const hdf5Offset(this->M_listData.begin()->getPtr()->BlockMap().IndexBase());
    int const hdf5Offset(0);

    // Points
    EpetraMap subMap(this->M_listData.begin()->getPtr()->BlockMap(), hdf5Offset, this->M_mesh->numGlobalVertices(),
                     this->M_listData.begin()->getPtr()->BlockMap().IndexBase() - hdf5Offset );

    EpetraVector pointsX(subMap);
    EpetraVector pointsY(subMap);
    EpetraVector pointsZ(subMap);


    int gid;
    for(ID i=1; i <= this->M_mesh->numVertices(); ++i)
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
            connectionsVarname += this->M_postfix; // see also in M_wr_topology
            pointsXVarname      += this->M_postfix; // see also in M_wr_geometry
            pointsYVarname      += this->M_postfix; // see also in M_wr_geometry
            pointsZVarname      += this->M_postfix; // see also in M_wr_geometry
        }


    M_HDF5->Write(connectionsVarname, connections);
    bool writeTranspose (true);
    M_HDF5->Write(pointsXVarname, pointsX.getEpetraVector(), true);
    M_HDF5->Write(pointsYVarname, pointsY.getEpetraVector(), true);
    M_HDF5->Write(pointsZVarname, pointsZ.getEpetraVector(), true);



}


// write empty xdmf file
template
<typename Mesh> void Hdf5exporter<Mesh>::M_wr_initXdmf()
{
    if (this->M_procId == 0)
    {
        M_xdmf.open( (this->M_prefix+".xmf").c_str(), std::ios_base::out );

        M_xdmf <<
            "<?xml version=\"1.0\" ?>\n" <<
            "<!DOCTYPE Xdmf SYSTEM \"prova.xdmf\" [\n" <<
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
template
<typename Mesh> void Hdf5exporter<Mesh>::M_wr_closeLinesXdmf()
{
    // save position
    M_closingLinesPosition = M_xdmf.tellp();

    // write closing lines
    M_xdmf << M_closingLines;
    M_xdmf.flush();

}


// remove closing lines
template
<typename Mesh> void Hdf5exporter<Mesh>::M_wr_removeCloseLinesXdmf()
{
    M_xdmf.seekp(M_closingLinesPosition);
}


template
<typename Mesh> void Hdf5exporter<Mesh>::M_wr_Xdmf(const Real& time)
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
        M_xdmf <<
            "    <!-- Time " << time << " -->\n"
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

template
<typename Mesh> void Hdf5exporter<Mesh>::M_wr_topology  ( std::ofstream& xdmf )
{

    xdmf <<
        "      <Topology\n" <<
        "         Type=\"Tetrahedron\"\n" <<
        "         NumberOfElements=\"" << this->M_mesh->numGlobalVolumes() << "\"\n" <<
        "         BaseOffset=\"1\">\n" <<
        "         <DataStructure Format=\"HDF\"\n" <<
        "                        Dimensions=\""<< this->M_mesh->numGlobalVolumes() << " 4\"\n" <<
        "                        DataType=\"Int\"\n" <<
        "                        Precision=\"8\">\n" <<
        "             &DataFile;:/Connections/Values\n" <<
        "         </DataStructure>\n" <<
        "      </Topology>\n";
}

template
<typename Mesh> void Hdf5exporter<Mesh>::M_wr_geometry  ( std::ofstream& xdmf )
{

    std::string geoVarName;

    // see also in postProcess
    if (this->M_multimesh)
        geoVarName = "Points" + this->M_postfix;
    else
        geoVarName = "Points";


    xdmf <<
        "      <Geometry Type=\"X_Y_Z\">\n" <<
        "         <DataStructure Format=\"HDF\"\n" <<
        "                        Dimensions=\"" << this->M_mesh->numGlobalVertices() << "\"\n" <<
        "                        DataType=\"Float\"\n" <<
        "                        Precision=\"8\">\n" <<
        "             &DataFile;:/" << geoVarName << "X/Values\n" <<
        "         </DataStructure>\n" <<
        "         <DataStructure Format=\"HDF\"\n" <<
        "                        Dimensions=\"" << this->M_mesh->numGlobalVertices() << "\"\n" <<
        "                        DataType=\"Float\"\n" <<
        "                        Precision=\"8\">\n" <<
        "             &DataFile;:/" << geoVarName << "Y/Values\n" <<
        "         </DataStructure>\n" <<
        "         <DataStructure Format=\"HDF\"\n" <<
        "                        Dimensions=\"" << this->M_mesh->numGlobalVertices() << "\"\n" <<
        "                        DataType=\"Float\"\n" <<
        "                        Precision=\"8\">\n" <<
        "             &DataFile;:/" << geoVarName << "Z/Values\n" <<
        "         </DataStructure>\n" <<
        "      </Geometry>\n" <<
        "\n";
}


template
<typename Mesh> void Hdf5exporter<Mesh>::M_wr_attributes  ( std::ofstream& xdmf )
{

    // Loop on the variables to output
    for (std::list< ExporterData >::const_iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i)
        {
            xdmf <<
                "\n      <Attribute\n" <<
                "         Type=\"" << i->typeName() << "\"\n" <<
                "         Center=\"Node\"\n" <<
                "         Name=\"" << this->M_prefix+"_"+i->variableName()<<"\">\n";

            switch( i->type() )
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

template
<typename Mesh> void Hdf5exporter<Mesh>::M_wr_scalar_datastructure  ( std::ofstream& xdmf, const ExporterData& dvar )
{
    xdmf <<
                "           <DataStructure  Format=\"HDF\"\n" <<
                "                           Dimensions=\"" << this->M_mesh->numGlobalVertices()  << " " << dvar.typeDim() << "\"\n" <<
                "                           DataType=\"Float\"\n" <<
                "                           Precision=\"8\">\n" <<
                "               &DataFile;:/" << this->M_prefix+"_"+dvar.variableName() << this->M_postfix  <<"/Values\n" << // see also in M_wr_vector/scalar
                "           </DataStructure>\n";

}

template
<typename Mesh> void Hdf5exporter<Mesh>::M_wr_vector_datastructure  ( std::ofstream& xdmf, const ExporterData& dvar )
{


    string coord[3]={"X","Y","Z"}; // see also wr_vector

    xdmf <<
                "         <DataStructure ItemType=\"Function\"\n" <<
                "                        Dimensions=\"" << this->M_mesh->numGlobalVertices() << " " << dvar.typeDim() << "\"\n" <<
                "                        Function=\"JOIN($0 , $1, $2)\">\n";

            for(int i(0); i < dvar.typeDim(); ++i)
                {
                    xdmf <<
                "           <DataStructure  Format=\"HDF\"\n" <<
                "                           Dimensions=\"" << this->M_mesh->numGlobalVertices() << " 1\"\n" <<
                "                           DataType=\"Float\"\n" <<
                "                           Precision=\"8\">\n" <<
                "               &DataFile;:/" << this->M_prefix+"_"+dvar.variableName()<< coord[i] << this->M_postfix  <<"/Values\n" << // see also in M_wr_vector/scalar
                "           </DataStructure>\n";
                }

    xdmf <<
                "         </DataStructure>\n";


}


}
#endif
#endif
