/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

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
#error you should reconfigure with --with-hdf5=... flag
#endif

#include <life/lifefilters/exporter.hpp>
#include <hdf5.h>
#include <EpetraExt_HDF5.h>

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



private:

    void M_wr_Xdmf(const Real& time);
    void M_wr_geo( const std::string geo_file );
    //void M_getPostfix();
    void M_wr_var(const ExporterData& dvar);
    void M_wr_scalar(const ExporterData& dvar);
    void M_wr_vector(const ExporterData& dvar);
    void M_Xdmf_mesh_section(std::ofstream& casef);
    void M_Xdmf_variable_section(std::ofstream& casef);
    void M_Xdmf_time_section(std::ofstream& casef, const Real& time);

    hdf5_ptrtype HDF5;

};


//
// Implementation
//
template<typename Mesh>
Hdf5exporter<Mesh>::Hdf5exporter(const GetPot& dfile, mesh_ptrtype mesh, const std::string prefix,
                                 const int procId)
    :
    super(dfile, mesh, prefix,procId),
    HDF5()
{
    setMeshProcId(mesh,procId);
}

template<typename Mesh>
Hdf5exporter<Mesh>::Hdf5exporter(const GetPot& dfile, const std::string prefix):
    super(dfile,prefix)
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
    if ( HDF5.get() == 0)
        {
            HDF5.reset(new hdf5_type(this->M_listData.begin()->getPtr()->Comm()));
            HDF5->Create(this->M_prefix+".h5");

            if (!this->M_multimesh)
                M_wr_geo( "geo" );
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
            M_wr_Xdmf(time);

            if (this->M_multimesh)
                M_wr_geo( "geo." + this->M_postfix);
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
    HDF5->Write("map-" + toString(Comm.NumProc()), Map);
    HDF5->Write("matrix", Matrix);
    HDF5->Write("LHS", LHS);
    HDF5->Write("RHS", RHS);
    */

    std::string scalars ("Scalars");
    std::string varname (dvar.prefix()+"."+this->M_postfix);


    HDF5->Write(varname, dvar.getPtr()->getEpetraVector());
}

template <typename Mesh> void Hdf5exporter<Mesh>::M_wr_vector(const ExporterData& dvar)
{

    UInt dim   = dvar.dim();
    UInt start = dvar.start();
    UInt nVert = this->M_mesh->numVertices();

    for (UInt d ( 0 ); d < nDimensions; ++d)
        {
            EpetraMap subMap(dvar.getPtr()->Map(), start, dim);
            vector_type subVar(subMap);
            subVar.subset(*dvar.getPtr(),start);

            std::string vectors ("Vectors");
            std::ostringstream varname;
            varname << dvar.prefix() << "." << d << "." << this->M_postfix;
            HDF5->Write(varname.str(), subVar.getEpetraVector());

        }

}

template <typename Mesh> void Hdf5exporter<Mesh>::M_wr_geo(const std::string geo_name)
{
#if 0
    std::ofstream geof( geo_file.c_str() );
    ID nV = this->M_mesh->numVertices();
    ID nE = this->M_mesh->numVolumes();
    UInt part=0;
    geof << "Geometry file\n";
    geof << "Geometry file\n";
    geof << "node id given\n";
    geof << "element id given\n";
    geof << "coordinates\n";
    geof.setf(std::ios::right | std::ios_base::scientific);
    geof.precision(5);
    geof.width(8);
    geof << nV << "\n";

    for(ID i=1; i <= nV; ++i)
        {
            geof.width(8);
            geof << i ;
            geof.width(12);
            geof << this->M_mesh->pointList(i).x();
            geof.width(12);
            geof << this->M_mesh->pointList(i).y();
            geof.width(12);
            geof << this->M_mesh->pointList(i).z();
            geof << "\n";
        }

    geof<< "part";

    geof.width(8);
    ++part;
    geof << part << "\n";
    geof << "full geometry\n";
    // elements
    geof << this->M_FEstr << "\n";
    geof.width(8);
    geof << nE << "\n";
    for (ID i=1; i <= nE; ++i)
        {
            geof.width(8);
            geof << i ;
            for (ID j=1; j<= this->M_nbLocalDof; ++j)
                {
                    geof.width(8);
                    //geof << M_mesh->volume(i).point(j).id();
                    geof << this->M_mesh->volume(i).point(j).localId();
                }
            geof << "\n";

        }

    geof.close();
#endif
}

template
<typename Mesh> void Hdf5exporter<Mesh>::M_wr_Xdmf(const Real& time)
{
  std::ofstream casef( (this->M_prefix+this->M_me+".case").c_str() );
    casef << "FORMAT\n";
    casef << "type: hdf5exporter\n";
    M_Xdmf_mesh_section(casef);
    M_Xdmf_variable_section(casef);
    M_Xdmf_time_section(casef,time);
    casef.close();
}

template <typename Mesh> void Hdf5exporter<Mesh>::M_Xdmf_mesh_section(std::ofstream& casef)
{
    casef << "GEOMETRY\n";
    if ( this->M_multimesh )
      casef << "model: 1 "+this->M_prefix+".***"<< this->M_me << ".geo change_coords_only\n";
    else
        casef << "model: 1 "+this->M_prefix+this->M_me+".geo\n";
}




template <typename Mesh> void Hdf5exporter<Mesh>::M_Xdmf_time_section(std::ofstream& casef, const Real& time)
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

template <typename Mesh> void Hdf5exporter<Mesh>::M_Xdmf_variable_section(std::ofstream& casef)
{
    typedef std::list< ExporterData >::const_iterator Iterator;
    casef << "VARIABLE\n";
    std::string aux;
    for (Iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i)
        {
            aux = i->prefix()+" "+i->prefix();
            switch( i->type() )
                {
                case ExporterData::Scalar:
		  casef << "scalar per node: 1 "+aux+".***" << this->M_me << ".scl\n";
                    break;
                case ExporterData::Vector:
                    casef << "vector per node: 1 "+aux+".***" << this->M_me << ".vct\n";
                    break;
                }
        }
}


}


#endif
