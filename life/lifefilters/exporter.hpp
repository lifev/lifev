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
  \file exporter.hpp
  \author S. Deparis
  \date 11/2008
  \version 1.0

  \brief This file provides an interface for post-processing

  Usage: two steps
  - first: add the variables using addVariable
  - second: call postProcess( time );
*/
#ifndef _EXPORTER_H_
#define _EXPORTER_H_

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
#include <life/lifearray/EpetraVector.hpp>

#include <boost/shared_ptr.hpp>


namespace LifeV
{

//typedef boost::numeric::ublas::vector_range< Vector >  VectorRange;
//typedef EpetraVector VectorRange;
//typedef double* VectorRange;

class ExporterData {

    typedef EpetraVector vector_rawtype;
    typedef boost::shared_ptr<vector_rawtype> vector_ptrtype;

public:
    enum Type{Scalar,Vector};

    ExporterData(const Type type, const std::string prefix, vector_ptrtype const& vec, UInt start, UInt size, UInt steady);

    std::string prefix() const;
    Real operator()(const UInt i) const;
    Real& operator()(const UInt i);
    UInt dim() const;
    UInt start() const { return M_start; }
    Type type() const;
    const vector_ptrtype getPtr() const { return M_vr; }
    UInt steady() const {return M_steady; }
    void set_steady(UInt i) {M_steady = i;}

private:
    std::string M_prefix;
    const vector_ptrtype M_vr;
    UInt M_dim;
    UInt M_start;
    Type M_type;
    UInt M_steady;

};

/**
 * \class Exporter
 * \brief generic data exporter
 */
template<typename Mesh>
class Exporter {

public:
    typedef boost::shared_ptr<Mesh> mesh_ptrtype;
    typedef ExporterData::vector_rawtype vector_rawtype;
    typedef ExporterData::vector_ptrtype vector_ptrtype;

    /**
       Constructor for Exporter

       \param dfile the GetPot data file where you must provide and [ensight] section with:
       "start" (start index for filenames 0 for 000, 1 for 001 etc.),
       "save" (how many time steps per posptrocessing)
       "multimesh" (=true if the mesh has to be saved at each post-processing step)

       \param mesh the mesh

       \param the prefix for the case file (ex. "test" for test.case)

       \param the procId determines de CPU id. if negative, it ussemes there is only one processor
    */
    Exporter(const GetPot& dfile, mesh_ptrtype mesh, const std::string prefix, const int procId );

    Exporter(const GetPot& dfile, const std::string prefix);

    virtual ~Exporter() {};

    /**
       setters
    */

    virtual void setMeshProcId( mesh_ptrtype mesh, int const procId );
    void initMeshProcId       ( mesh_ptrtype mesh, int const procId );


    /**
       Adds a new variable to be post-processed

       \param type the type fo the variable Ensight::Scalar or Ensight::Vector

       \param prefix the prefix of the files storing the variable (ex: "velocity" for velocity.***)

       \param vr an ublas::vector_range type given a view of the varialbe (ex: subrange(fluid.u(),0,3*dimU) )

       \param dim the number of Dof for that variable
    */
    void addVariable(const ExporterData::Type type, const std::string prefix, vector_ptrtype const& map, UInt start, UInt size, UInt steady =0 );


    /**
       Post-porcess the variables added to the list

       \param time the solver time
    */
    virtual void postProcess(const Real& time) = 0;

protected:

    void getPostfix();

    mesh_ptrtype M_mesh;
    std::string M_prefix;
    std::string M_post_dir;
    UInt M_count;
    UInt M_save;
    bool M_multimesh;
    UInt M_steps;
    std::list<Real> M_timeSteps;
    std::string M_FEstr;
    std::string M_bdFEstr;
    UInt M_nbLocalDof;
    UInt M_nbLocalBdDof;
    std::string M_postfix;
    std::list<ExporterData> M_listData;

    int M_procId;
    std::string M_me;
};


//
// Implementation
//
template<typename Mesh>
Exporter<Mesh>::Exporter(const GetPot& dfile, mesh_ptrtype mesh, const std::string prefix,
						 int const procId)
    :
    M_prefix(prefix),
    M_post_dir(dfile("exporter/post_dir", "./")),
    M_count(dfile("exporter/start",0)),
    M_save(dfile("exporter/save",1)),
    M_multimesh(dfile("exporter/multimesh",true)),
    M_steps(0)
{
    setMeshProcId(mesh,procId);
}

template<typename Mesh>
Exporter<Mesh>::Exporter(const GetPot& dfile, const std::string prefix):
    M_prefix(prefix),
    M_post_dir(dfile("exporter/post_dir", "./")),
    M_count(dfile("exporter/start",0)),
    M_save(dfile("exporter/save",1)),
    M_multimesh(dfile("exporter/multimesh",true)),
    M_steps(0)
{
}

template<typename Mesh>
void Exporter<Mesh>::setMeshProcId( mesh_ptrtype mesh , int const procId )
{
    initMeshProcId(mesh,procId);
}

template<typename Mesh>
void Exporter<Mesh>::initMeshProcId( mesh_ptrtype mesh , int const procId )
{
    M_mesh = mesh;

    M_procId = procId;

    std::ostringstream index;
    index.fill( '0' );
    if (M_procId >=0)
        {
            index << std::setw(1) << "." ;
            index << std::setw(3) << M_procId;
        }
    M_me = index.str();

}

template <typename Mesh> void Exporter<Mesh>::getPostfix()
{

    std::ostringstream index;
    index.fill( '0' );

    if (M_count % M_save == 0)
        {
            index << std::setw(3) << ( M_count / M_save );

            M_postfix = "." + index.str();

        }
    else
        M_postfix = "***";

    ++M_count;
}


template<typename Mesh>
void Exporter<Mesh>::addVariable(const ExporterData::Type type, const std::string prefix, vector_ptrtype const& vr, UInt start, UInt dim, UInt steady)
{
    M_listData.push_back( ExporterData(type,prefix,vr,start, dim, steady) );
}


}


#endif
