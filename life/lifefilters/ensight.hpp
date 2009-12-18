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
  \file ensight.hpp
  \author M.A. Fernandez, C. Prud'homme, S. Deparis
  \date 10/2005 08/2007
  \version 1.0

  \brief This file provides an interface for post-processing with ensight

  Usage: two steps
  - first: add the variables using addVariable
  - second: call postProcess( time );
*/
#ifndef _ENSIGHT_H_
#define _ENSIGHT_H_

#include <life/lifefilters/exporter.hpp>


namespace LifeV
{

/**
 * \class Ensight
 * \brief Ensight data exporter
 */
template<typename Mesh>
class Ensight : public Exporter<Mesh> {

public:

    typedef Exporter<Mesh> super;
    typedef typename super::mesh_ptrtype  mesh_ptrtype;
    typedef typename super::vector_ptrtype vector_ptrtype;

    /**
       Constructor for Ensight

       \param dfile the GetPot data file where you must provide and [ensight] section with:
       "start" (start index for filenames 0 for 000, 1 for 001 etc.),
       "save" (how many time steps per posptrocessing)
       "multimesh" (=true if the mesh has to be saved at each post-processing step)

       \param mesh the mesh

       \param the prefix for the case file (ex. "test" for test.case)

       \param the procId determines de CPU id. if negative, it ussemes there is only one processor
    */
    Ensight(const GetPot& dfile, mesh_ptrtype mesh, const std::string prefix, const int procId );

    Ensight(const GetPot& dfile, const std::string prefix);


    /**
       setters
    */

    void setMeshProcId( mesh_ptrtype mesh, int const procId );

    //! returns the type of the map to use for the EpetraVector
    EpetraMapType mapType() const;

    /**
       Post-porcess the variables added to the list

       \param time the solver time
    */
    void postProcess(const Real& time);

    /**
       Import data from previous simulations

       \param time the solver time
    */
    void import(const Real& Tstart, const Real& dt); // dt is used to rebuild the history up to now

    //! Read  only last timestep
    void import(const Real& Tstart);

    void rd_var(ExporterData& dvar)
    {super::rd_var(dvar);}

private:

    void M_wr_case(const Real& time);
    void M_wr_ascii_geo( const std::string geo_file );

    void M_wr_ascii(const ExporterData& dvar);
    void M_wr_ascii_scalar(const ExporterData& dvar);
    void M_wr_ascii_vector(const ExporterData& dvar);
    void M_case_mesh_section(std::ofstream& casef);
    void M_case_variable_section(std::ofstream& casef);
    void M_case_time_section(std::ofstream& casef, const Real& time);

    void M_rd_scalar( ExporterData& dvar );
    void M_rd_vector( ExporterData& dvar );



};


//
// Implementation
//
template<typename Mesh>
Ensight<Mesh>::Ensight(const GetPot& dfile, mesh_ptrtype mesh, const std::string prefix,
                       int const procId)
    :
    super(dfile, mesh, prefix,procId)
{
    setMeshProcId(mesh,procId);
}

template<typename Mesh>
Ensight<Mesh>::Ensight(const GetPot& dfile, const std::string prefix):
    super(dfile,prefix)
{
}

template<typename Mesh>
void Ensight<Mesh>::setMeshProcId( mesh_ptrtype mesh , int const procId )
{

    initMeshProcId( mesh, procId );

    typedef typename Mesh::ElementShape ElementShape;

    switch ( ElementShape::Shape )
        {
        case TETRA:
            this->M_FEstr = "tetra4";
            this->M_bdFEstr = "tria3";
            this->M_nbLocalBdDof = 3;
            this->M_nbLocalDof = 4;
            break;
        case HEXA:
            this->M_FEstr = "hexa8";
            this->M_bdFEstr = "quad4";
            this->M_nbLocalBdDof = 4;
            this->M_nbLocalDof = 8;
            break;
        case TRIANGLE:
            this->M_FEstr = "tria3";
            this->M_bdFEstr = "bar2";
            this->M_nbLocalBdDof = 2;
            this->M_nbLocalDof = 3;
            break;
        case QUAD:
            this->M_FEstr = "quad4";
            this->M_bdFEstr = "bar2";
            this->M_nbLocalBdDof = 4;
            this->M_nbLocalDof = 3;
            break;
        default:
            ERROR_MSG( "FE not allowed in Ensight writer" );
        }

    if (!this->M_multimesh)
      M_wr_ascii_geo( this->M_post_dir+this->M_prefix+this->M_me+".geo" );


}
template<typename Mesh>
EpetraMapType Ensight<Mesh>::mapType() const
{
    return Repeated;
}

template<typename Mesh>
void Ensight<Mesh>::postProcess(const Real& time)
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
            	    		M_wr_ascii(*i);
            	if (i->steady() == 1) i->set_steady(2);
                }
            M_wr_case(time);

            if (this->M_multimesh)
                M_wr_ascii_geo( this->M_post_dir+this->M_prefix+this->M_postfix+this->M_me+".geo" );
            chrono.stop();
            if (!this->M_procId) std::cout << "      done in " << chrono.diff() << " s." << std::endl;
        }
}

template <typename Mesh>
void Ensight<Mesh>::M_wr_ascii(const ExporterData& dvar)
{

    switch( dvar.type() )
        {
        case ExporterData::Scalar:
            M_wr_ascii_scalar(dvar);
            break;
        case ExporterData::Vector:
            M_wr_ascii_vector(dvar);
            break;
        }

}



template <typename Mesh>
void Ensight<Mesh>::M_wr_ascii_scalar(const ExporterData& dvar)
{
    using std::setw;

	std::ofstream sclf;

	if (dvar.steady() )
		sclf.open( (this->M_post_dir+super::M_prefix+"_"+dvar.variableName()+this->M_me+".scl").c_str() );
	else
		sclf.open( (this->M_post_dir+super::M_prefix+"_"+dvar.variableName()+this->M_postfix+this->M_me+".scl").c_str() );

    UInt count=0;

    UInt start = dvar.start();
    //    UInt nVert = this->M_mesh->numVertices();
    UInt nVert = static_cast<UInt> (this->M_LtGNodesMap.size());
    sclf<<"Scalar per node\n";
    //for (UInt i=start;i<dim;++i)

    sclf.setf(std::ios::right | std::ios_base::scientific);
    sclf.precision(5);

    for (UInt i=0; i<nVert; ++i)
        {
    	    // int id = this->M_mesh->pointList( i ).id();
            int id = this->M_LtGNodesMap[i];
            sclf << setw(12) << float(dvar(start+id)) ;
            ++count;
            if ( count == 6 )
                {
                    sclf << "\n";
                    count=0;
                }
        }
    sclf << std::endl;
}

template <typename Mesh> void Ensight<Mesh>::M_wr_ascii_vector(const ExporterData& dvar)
{
    using std::setw;

	std::ofstream vctf;

	if (dvar.steady() )
		vctf.open( (this->M_post_dir+super::M_prefix+"_"+dvar.variableName()+this->M_me+".vct").c_str() );
	else
		vctf.open( (this->M_post_dir+super::M_prefix+"_"+dvar.variableName()+this->M_postfix+this->M_me+".vct").c_str() );

    UInt count=0;

    UInt dim   = dvar.size();
    UInt start = dvar.start();
    //    UInt nVert = this->M_mesh->numVertices();
    UInt nVert = static_cast<UInt> (this->M_LtGNodesMap.size());

    vctf<<"Vector per node\n";
    //for (UInt i=start;i<dim;++i)

    vctf.setf(std::ios::right | std::ios_base::scientific);
    vctf.precision(5);

    for (UInt i=0; i<nVert; ++i)
        for (UInt j=0; j< nDimensions;++j)
            {
                // int id = this->M_mesh->pointList( i ).id();
                int id = this->M_LtGNodesMap[i];
                vctf << setw(12) << float(dvar(start+j*dim+id)) ;
                ++count;
                if ( count == 6 )
                    {
                        vctf << "\n";
                        count=0;
                    }
            }
    vctf << std::endl;
}


template <typename Mesh> void Ensight<Mesh>::M_wr_ascii_geo(const std::string geo_file)
{
    using std::setw;

    std::ofstream geof(geo_file.c_str() );
    ID nV = this->M_mesh->numVertices();
    ID nE = this->M_mesh->numElements();
    UInt part=0;
    geof << "Geometry file\n";
    geof << "Geometry file\n";
    geof << "node id given\n";
    geof << "element id given\n";
    geof << "coordinates\n";
    geof.setf(std::ios::right | std::ios_base::scientific);
    geof.precision(5);
    geof << setw(8) <<  nV << "\n";
    for(ID i=1; i <= nV; ++i)
        {
            geof << setw(8) << i ;
            for(UInt icoor=0; icoor<nDimensions; icoor++)
            	geof << setw(12) << float(this->M_mesh->pointList(i).coor()[icoor]);
            geof << "\n";
        }

    geof<< "part";

    ++part;
    geof << setw(8) << part << "\n";
    geof << "full geometry\n";
    // elements
    geof << this->M_FEstr << "\n";
    geof << setw(8) << nE << "\n";
    for (ID i=1; i <= nE; ++i)
        {
            geof << setw(8) << i ;
            for (ID j=1; j<= this->M_nbLocalDof; ++j)
                {
                    //geof << M_mesh->volume(i).point(j).id();
                    geof << setw(8) << this->M_mesh->element(i).point(j).localId();
                }
            geof << "\n";

        }

    // boundary parts
#if 0
    std::set<EntityFlag> flags;
    std::map< EntityFlag, __gnu_cxx::slist<ID> > faces;
    EntityFlag marker;
    typedef std::set<EntityFlag>::const_iterator Iterator_flag;
    typedef __gnu_cxx::slist<ID>::const_iterator Iterator_face;

    Iterator_flag result;

    ID nBF = M_mesh->numBFaces();
    for (ID i=1; i <= nBF; ++i)
        {
            marker = M_mesh->boundaryFace(i).marker();
            flags.insert(marker);
            faces[marker].push_front(i);
        }

    for (Iterator_flag i=flags.begin(); i!= flags.end(); ++i)
        {
            marker = *i;
            geof<< "part";
            ++part;
            geof<< setw(8) << part << "\n";
            geof<<"boundary ref "<< marker << "\n";
            __gnu_cxx::slist<ID>& faceList= faces[marker];
            geof<< M_bdFEstr << "\n";
            geof<< setw(8) << faceList.size() << "\n";
            for (Iterator_face j=faceList.begin(); j!= faceList.end(); ++j)
                {
                    for( ID k=1; k <= M_nbLocalBdDof; ++k)
                        {
                            geof << setw(8) << M_mesh->boundaryFace(*j).point(k).id();
                        }
                    geof << "\n";
                }
        }
#endif
}

template
<typename Mesh> void Ensight<Mesh>::M_wr_case(const Real& time)
{
  std::ofstream casef( (this->M_post_dir+this->M_prefix+this->M_me+".case").c_str() );
    casef << "FORMAT\n";
    casef << "type: ensight\n";
    M_case_mesh_section(casef);
    M_case_variable_section(casef);
    M_case_time_section(casef,time);
}

template <typename Mesh> void Ensight<Mesh>::M_case_mesh_section(std::ofstream& casef)
{
    casef << "GEOMETRY\n";
    if ( this->M_multimesh )
      casef << "model: 1 "+this->M_prefix+".*****"<< this->M_me << ".geo change_coords_only\n";
    else
        casef << "model: 1 "+this->M_prefix+this->M_me+".geo\n";
}




template <typename Mesh> void Ensight<Mesh>::M_case_time_section(std::ofstream& casef, const Real& time)
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

template <typename Mesh> void Ensight<Mesh>::M_case_variable_section(std::ofstream& casef)
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
            aux = i->variableName()+" "+super::M_prefix+"_"+i->variableName();
            switch( i->type() )
                {
                case ExporterData::Scalar:
		  casef << "scalar per node: 1 "+ aux + str << this->M_me << ".scl\n";
                    break;
                case ExporterData::Vector:
                    casef << "vector per node: 1 "+ aux+ str << this->M_me << ".vct\n";
                    break;
                }
        }
}

template<typename Mesh>
void Ensight<Mesh>::import(const Real& Tstart, const Real& dt)
{
    // dt is used to rebuild the history up to now
    Real time(Tstart - this->M_count*dt);

    for ( UInt count(0); count < this->M_count; ++count)
    {
        this->M_timeSteps.push_back(time);
        ++this->M_steps;
        time += dt;
    }

    time += dt;

    import(time);

}

template<typename Mesh>
void Ensight<Mesh>::import(const Real& time)
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
            this->rd_var(*i);
        }
    chrono.stop();
    if (!this->M_procId) std::cout << "      done in " << chrono.diff() << " s." << std::endl;

}

template <typename Mesh>
void Ensight<Mesh>::M_rd_scalar( ExporterData& dvar )
{

    std::string filename( this->M_import_dir+super::M_prefix+"_"+dvar.variableName()+this->M_postfix+this->M_me+".scl" );
    std::ifstream sclf( filename.c_str() );

    if (!this->M_procId) std::cout << "\tfile "<< filename << std::endl;

    ASSERT(sclf.good(), std::stringstream("There is an error while reading " + filename).str().c_str() );

    // UInt count=0;

//    UInt dim = dvar.size();
    UInt start = dvar.start();
    //    UInt nVert = this->M_mesh->numVertices();
    UInt nVert = static_cast<UInt> (this->M_LtGNodesMap.size());

    std::string trashcan;

    sclf >> trashcan >> trashcan >> trashcan; // <<"Scalar per node\n";

    sclf.setf(std::ios::right | std::ios_base::scientific);
    sclf.precision(5);

	std::map<int,int>::iterator iter;

    for (UInt i=0; i<nVert; ++i)
        {
    	    ASSERT(sclf.good(), std::stringstream("There is an error while reading " + filename).str().c_str() );

            // int id = this->M_mesh->pointList( i ).id();
            int id = this->M_LtGNodesMap[i];
            sclf.width(12);
            sclf >> dvar(start+id) ;
            /*
            ++count;
            if ( count == 6 )
                {
                    sclf << "\n";
                    count=0;
                }
            */
        }
    ASSERT(!sclf.fail(), std::stringstream("There is an error while reading " + filename).str().c_str() );
}

template <typename Mesh> void Ensight<Mesh>::M_rd_vector(ExporterData& dvar)
{

    std::string filename( this->M_import_dir+super::M_prefix+"_"+dvar.variableName()+this->M_postfix+this->M_me+".vct" );
    std::ifstream vctf( filename.c_str() );

    if (!this->M_procId) std::cout << "\tfile "<< filename << std::endl;

    ASSERT(vctf.good(), std::stringstream("There is an error while reading " + filename).str().c_str() );

//    UInt count=0;

    UInt dim   = dvar.size();
    UInt start = dvar.start();
    //    UInt nVert = this->M_mesh->numVertices();
    UInt nVert = static_cast<UInt> (this->M_LtGNodesMap.size());

    std::string trashcan;

    vctf >> trashcan >> trashcan >> trashcan;  // <<"Vector per node\n";

    vctf.setf(std::ios::right | std::ios_base::scientific);
    vctf.precision(5);

    for (UInt i=0; i<nVert; ++i)
        for (UInt j=0; j< nDimensions;++j)
            {
                ASSERT(vctf.good(), std::stringstream("There is an error while reading " + filename).str().c_str() );

                // int id = this->M_mesh->pointList( i ).id();
                int id = this->M_LtGNodesMap[i];
                vctf.width(12);
                vctf >> dvar(start+j*dim+id) ;
                /*
                ++count;
                if ( count == 6 )
                    {
                        vctf << "\n";
                        count=0;
                    }
                */
            }
    ASSERT(!vctf.fail(), std::stringstream("There is an error while reading " + filename).str().c_str() );
}



}


#endif
