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
  \file ensight.h
  \author M.A. Fernandez
  \date 10/2005
  \version 1.0

  \brief This file provides an interface for post-processing with ensight

  Usage: two steps
  - first: add the variables using addVariable
  - second: call postProcess( time );
*/
#ifndef _ENSIGHT_H_
#define _ENSIGHT_H_

#include <life/lifecore/life.hpp>
#include <life/lifemesh/markers.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifearray/tab.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <list>
#include <ext/slist>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/chrono.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

namespace LifeV
{
  
  typedef boost::numeric::ublas::vector_range< Vector >  VectorRange;

  class EnsightData {
    
  public:   
    enum Type{Scalar,Vector};
  
    EnsightData(const Type type, const std::string prefix, const VectorRange& vr, const UInt dim);
    
    std::string prefix() const;
    Real operator()(const UInt i) const;
    UInt dim() const;
    Type type() const;

  private:
    std::string M_prefix;
    const VectorRange M_vr;
    UInt M_dim;
    Type M_type;

  };



  template<typename Mesh> 
  class Ensight {

  public:

    Ensight(const GetPot& dfile, Mesh& mesh, const std::string prefix);
   
    /**
       Constructor for Ensight

       \param dfile the GetPot data file where you must provide and [ensight] section with: 
       "start" (start index for filenames 0 for 000, 1 for 001 etc.), 
       "verbose" (how many time steps per posptrocessing) 
       "multimesh" (=1 if the mesh has to be saved at each post-processing step)

       \param mesh the mesh

       \param the prefix for the case file (ex. "test" for test.case)
    */


    void addVariable(const EnsightData::Type type, const std::string prefix, const VectorRange& vr, const ID dim);
    /**
       Adds a new variable to be post-processed 
    
       \param type the type fo the variable Ensight::Scalar or Ensight::Vector

       \param prefix the prefix of the files storing the variable (ex: "velocity" for velocity.***)

       \param vr an ublas::vector_range type given a view of the varialbe (ex: subrange(fluid.u(),0,3*dimU) )

       \param dim the number of Dof for that variable
    */


    void postProcess(const Real& time);
    /**
       Post-porcess the variables added to the list
       
       \param time the solver time
    */

    
    
  private:

    void M_wr_case(const Real& time);
    void M_wr_ascii_geo( const std::string geo_file );
    void M_getPostfix();
    void M_wr_ascii(const EnsightData& dvar);  
    void M_wr_ascii_scalar(const EnsightData& dvar); 
    void M_wr_ascii_vector(const EnsightData& dvar);  
    void M_case_mesh_section(std::ofstream& casef);    
    void M_case_variable_section(std::ofstream& casef);
    void M_case_time_section(std::ofstream& casef, const Real& time);

    Mesh& M_mesh;
    std::string M_prefix;
    UInt M_count;
    UInt M_verbose;
    bool M_multimesh;
    UInt M_steps;
    std::list<Real> M_timeSteps;
    std::string M_FEstr;
    std::string M_bdFEstr;
    UInt M_nbLocalDof;
    UInt M_nbLocalBdDof;
    std::string M_postfix;
    std::list<EnsightData> M_listData;
  };

  
  //
  // Implementation
  //
  template<typename Mesh> Ensight<Mesh>::Ensight(const GetPot& dfile, Mesh& mesh, const std::string prefix):
    M_mesh(mesh),
    M_prefix(prefix),
    M_count(dfile("ensight/start",0)),
    M_verbose(dfile("ensight/verbose",1)),
    M_multimesh(dfile("ensight/multimesh",0)),
    M_steps(0)
  {
    typedef typename Mesh::VolumeShape ElementShape;
    
    switch ( ElementShape::numPoints )
      {
      case 4:
        M_FEstr = "tetra4";
        M_bdFEstr = "tria3";
        M_nbLocalBdDof = 3; 
        M_nbLocalDof = 4;
        break;
      case 8:
        M_FEstr = "hexa8";
        M_bdFEstr = "quad4";    
        M_nbLocalBdDof = 4;
        M_nbLocalDof = 8;
        break;
      default:
        ERROR_MSG( "FE not allowed in Ensight writer" );
      }	
    if (!M_multimesh)
      M_wr_ascii_geo( M_prefix+".geo" );
  
  }
  
  template <typename Mesh> void Ensight<Mesh>::M_getPostfix()
  {
    
    std::ostringstream index;
  
    
    if ( fmod( float( M_count ), float( M_verbose ) ) == 0.0 )
      {
        index << ( M_count / M_verbose );
        
        switch ( index.str().size() )
          {
          case 1:
            M_postfix = "00" + index.str();
            break;
          case 2:
            M_postfix = "0" + index.str();
            break;
          case 3:
            M_postfix = index.str();
            break;
          }

      }
    else
      M_postfix = "***";
    
    ++M_count;
  }

  
  
  template <typename Mesh> void Ensight<Mesh>::M_wr_ascii_vector(const EnsightData& dvar)
  {
    std::ofstream vctf( (dvar.prefix()+"."+M_postfix+".vct").c_str() );
   
    UInt count=0;

    UInt dim = dvar.dim();
    
    vctf<<"Vector per node"<<std::endl;
    for (UInt i=0;i<dim;++i)
      for (UInt j=0; j< nDimensions;++j)
        {
          vctf.setf(std::ios_base::scientific); 
          vctf.precision(5);    
          vctf.width(12);
          vctf << dvar(j*dim+i);
          ++count;
          if ( count == 6 ) 
            {
              vctf << std::endl;
              count=0;
            }
        } 
    vctf << std::endl;
    vctf.close();
  }

  template<typename Mesh> 
  void Ensight<Mesh>::addVariable(const EnsightData::Type type, const std::string prefix, const VectorRange& vr, const ID dim) 
  {
    M_listData.push_back( EnsightData(type,prefix,vr,dim) );
  }
  
  template<typename Mesh> 
  void Ensight<Mesh>::postProcess(const Real& time)
  {
    typedef std::list< EnsightData >::const_iterator Iterator;

    M_getPostfix();

    if ( M_postfix != "***" )
      {
	std::cout << "  x-  Ensight post-processing..."<< std::endl;
	Chrono chrono;   
	chrono.start();
	for (Iterator i=M_listData.begin(); i != M_listData.end(); ++i)
	  {
	    M_wr_ascii(*i);
	  }
	M_wr_case(time);
	
	if (M_multimesh)
	  M_wr_ascii_geo( M_prefix+"."+M_postfix+".geo" );
	chrono.stop();
	std::cout << "      done in " << chrono.diff() << " s." << std::endl;
      }
  }
  
  template <typename Mesh> 
  void Ensight<Mesh>::M_wr_ascii(const EnsightData& dvar) 
  {

    switch( dvar.type() )
      {
      case EnsightData::Scalar:
	M_wr_ascii_scalar(dvar);
	break;
      case EnsightData::Vector:
	M_wr_ascii_vector(dvar);
	break;
      }
    
  }


 
  template <typename Mesh> 
  void Ensight<Mesh>::M_wr_ascii_scalar(const EnsightData& dvar) 
  {

    std::ofstream sclf( (dvar.prefix()+"."+M_postfix+".scl").c_str() );
   
    UInt count=0;
       
    UInt dim = dvar.dim();
    
    sclf<<"Scalar per node"<<std::endl;
    for (UInt i=0;i<dim;++i) 
      {
        sclf.setf(std::ios_base::scientific);   
        sclf.precision(5);      
        sclf.width(12);
        sclf << dvar(i);
        ++count;
        if ( count == 6 ) 
          {
            sclf << std::endl;
            count=0;
          }
      }
    sclf << std::endl;
    sclf.close();
  }


  template <typename Mesh> void Ensight<Mesh>::M_wr_ascii_geo(const std::string geo_file) 
  {
    std::ofstream geof( geo_file.c_str() );
    ID nV = M_mesh.numVertices();
    ID nE = M_mesh.numVolumes();
    ID nBF = M_mesh.numBFaces();
    UInt part=0;

    geof<<"Geometry file"<<std::endl;
    geof<<"Geometry file"<<std::endl;
    geof<<"node id assign"<<std::endl;
    geof<<"element id assign"<<std::endl;
    geof<< "coordinates" << std::endl;
    geof.width(8);
    geof<< nV << std::endl;
   
    for(ID i=1; i <= nV; ++i) 
      {
        for (ID j=1; j <= nDimensions; ++j) 
          {
            geof.setf(std::ios_base::scientific);
            geof.precision(5);  
            geof.width(12);
            geof <<  M_mesh.point(i).coordinate(j);
          }
        
        geof<< std::endl;
      } 

    // volume parts
    EntityFlag marker;
    std::set<EntityFlag> flags;
    std::map< EntityFlag, __gnu_cxx::slist<ID> > geo;

    typedef std::set<EntityFlag>::const_iterator Iterator_flag;
    typedef __gnu_cxx::slist<ID>::const_iterator Iterator_geo;
    Iterator_flag result;

    for (ID i=1; i <= nE; ++i)
      {
        marker = M_mesh.volume(i).marker();
        flags.insert(marker);
        geo[marker].push_front(i);
      }
    for (Iterator_flag i=flags.begin(); i!= flags.end(); ++i)
      {
        marker = *i;
        geof<< "part ";
        geof.width(8);
        ++part;
        geof<< part << std::endl;
        geof<<"subdomain ref "<< marker << std::endl;
        __gnu_cxx::slist<ID>& volumeList= geo[marker];
        geof<< M_FEstr << std::endl;
        geof.width(8);
        geof<< volumeList.size() << std::endl;
        for (Iterator_geo j=volumeList.begin(); j!= volumeList.end(); ++j) 
          {
            for( ID k=1; k <= M_nbLocalDof; ++k) 
              {
                geof.width(8);
                geof << M_mesh.volume(*j).point(k).id();
              }
            geof<<"\n";
          }
      }
    
    geo.clear();
    flags.clear();
    
    
    // boundary parts
       
   
    for (ID i=1; i <= nBF; ++i)
      {
        marker = M_mesh.boundaryFace(i).marker();
        flags.insert(marker);
        geo[marker].push_front(i);
      }
      
    for (Iterator_flag i=flags.begin(); i!= flags.end(); ++i)
      {
        marker = *i;
        geof<< "part";
        geof.width(8);
        ++part;
        geof<< part << std::endl;
        geof<<"boundary ref "<< marker << std::endl;
        __gnu_cxx::slist<ID>& faceList= geo[marker];
        geof<< M_bdFEstr << std::endl;
        geof.width(8);
        geof<< faceList.size() << std::endl;
        for (Iterator_geo j=faceList.begin(); j!= faceList.end(); ++j) 
          {
            for( ID k=1; k <= M_nbLocalBdDof; ++k) 
              {
                geof.width(8);
                geof << M_mesh.boundaryFace(*j).point(k).id();
              }
            geof<<"\n";
          }
      }
     
    geo.clear();
    flags.clear();
    geof.close();
  }
  
  template 
  <typename Mesh> void Ensight<Mesh>::M_wr_case(const Real& time) 
  {
    std::ofstream casef( (M_prefix+".case").c_str() );
    casef << "FORMAT\n";
    casef << "type: ensight\n";
    M_case_mesh_section(casef);
    M_case_variable_section(casef);
    M_case_time_section(casef,time);
    casef.close();
  }

  template <typename Mesh> void Ensight<Mesh>::M_case_mesh_section(std::ofstream& casef)
  {
    casef << "GEOMETRY\n";
    if ( M_multimesh )
      casef << "model: 1 "+M_prefix+".***.geo change_coords_only\n";
    else
      casef << "model: 1 "+M_prefix+".geo\n";    
  }




  template <typename Mesh> void Ensight<Mesh>::M_case_time_section(std::ofstream& casef, const Real& time)
  {
    M_timeSteps.push_back(time);
    ++M_steps;
    casef << "TIME\n";
    casef << "time set: 1\n"; 
    casef << "number of steps: " <<  M_steps << "\n";
    casef << "filename start number: 0\n";
    casef << "filename increment: 1\n";
    casef << "time values:\n";
   
    UInt count=0;
   
    typedef std::list<Real>::const_iterator Iterator;
    for (Iterator i=M_timeSteps.begin(); i != M_timeSteps.end(); ++i)
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
    typedef std::list< EnsightData >::const_iterator Iterator;
    casef << "VARIABLE\n";
    std::string aux;
    for (Iterator i=M_listData.begin(); i != M_listData.end(); ++i)
      {
	aux = i->prefix()+" "+i->prefix();
	switch( i->type() )
	  {
	  case EnsightData::Scalar:
	    casef << "scalar per node: 1 "+aux+".***.scl\n";
	    break;
	  case EnsightData::Vector:
	    casef << "vector per node: 1 "+aux+".***.vct\n";
	    break;
	  }
      }
  }


}


#endif
