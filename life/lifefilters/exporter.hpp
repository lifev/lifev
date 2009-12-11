/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2008 EPFL

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

//@HEADER

/*!
  @file exporter.hpp
  @brief This file provides an interface for post-processing

  @author Simone Deparis <simone.deparis@epfl.ch>
  @date 11-11-2008

  Usage: two steps
  <ol>
    <li> first: add the variables using addVariable
    <li> second: call postProcess( time );
  </ol>
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

//! ExporterData - Holds the datastructure of the an array to import/export
/*!
  @author Simone Deparis <simone.deparis@epfl.ch>
  @date 11-11-2008

  This class holds the datastructure of one datum
  to help the importer/exporter
  like:   variable name, shared pointer to data, size, and few others.
 */
class ExporterData
{
public:

    //! @name Public Types
    //@{

    typedef EpetraVector                      vector_rawtype;
    typedef boost::shared_ptr<vector_rawtype> vector_ptrtype;

    //! Type of data stored.
    enum Type
    {
        Scalar, /*!< Scalar stands for scalarfield */
        Vector  /*!< Vector stands for vectorfield */
    };

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructor with all the data to be stored
    /*!
        Constructor with all the data to be stored
        @param type         - scalar or vector field
        @param variableName - name assigned to this variable in output file
        @param vec          - shared pointer to array
        @param start        - address of first datum in the array to be stored.
                              Useful in case you want to define a subrange of the vector *vec
        @param size         - size of the array to store (not yet multiplied by the dimension of the vector)
        @param steady       - if  file name for postprocessing has to include time dependency
    */
    ExporterData(const Type type, const std::string variableName, vector_ptrtype const& vec, UInt start, UInt size, UInt steady);

    //@}

    //! @name Operators
    //@{

    //! Accessor to component (const)
    /*!
        @param i which component to access
    */
    Real operator()(const UInt i) const;

    //! Accessor to component (reference)
    /*!
        @param i which component to access
    */
    Real& operator()(const UInt i);

    //@}

    //! @name Set Methods
    //@{

    //! file name for postprocessing has to include time dependency
    void set_steady(UInt i) {M_steady = i;}

    //@}

    //! @name Get Methods
    //@{

    //! name assigned to this variable in output file
    const std::string&  variableName() const;

    //! size of the stored array. was: dim()
    const UInt& size() const;

    //! address of first datum in the array
    const UInt& start() const { return M_start; }

    //! scalar or vector field
    const Type& type() const;

    //! shared pointer to array
    const vector_ptrtype storedArray() const { return M_vr; }

    //! returns 0 if file name for postprocessing has to include time dependency
    UInt steady() const {return M_steady; }

    //! returns Scalar or Vector strings
    std::string typeName() const;

    //! returns 1 (if Scalar) or 3 (if Vector)
    UInt typeDim() const;

    //@}


private:

    //! @name Private Methods
    //@{

    //! dim has been replaced by size() and will be removed
    const UInt dim() const { assert(false); }

    //! shared pointer to array, replaced by storedArray() and will be removed
    const vector_ptrtype getPtr() const { assert(false); return M_vr; }

    //@}

    //! name assigned to this variable in output file
    std::string M_variableName;

    //! pointer to storedArray
    const vector_ptrtype M_vr;

    //! address of first datum in the array
    UInt M_size;

    //! address of first datum in the array
    UInt M_start;

    //! scalar or vector field
    Type M_type;

    //! equal to 0 if file name for postprocessing has to include time dependency
    UInt M_steady;

};


//! Exporter - Pure virtual class that describes a generic exporter
/*!
  @author Simone Deparis <simone.deparis@epfl.ch>
  @date 11-11-2008

  This class is pure virtual and describes a generic exporter that can
  also do import
 */
template<typename Mesh>
class Exporter {

public:

    typedef boost::shared_ptr<Mesh>      mesh_ptrtype;
    typedef ExporterData::vector_rawtype vector_rawtype;
    typedef ExporterData::vector_ptrtype vector_ptrtype;

    //! @name Constructor & Destructor
    //@{

    //! Constructor for Exporter
    /*!
        \param dfile the GetPot data file where you must provide an [exporter] section with:
          "start"     (start index for filenames 0 for 000, 1 for 001 etc.),
          "save"      (how many time steps per postprocessing)
          "multimesh" ( = true if the mesh has to be saved at each post-processing step)
       \param mesh the mesh
       \param the prefix for the case file (ex. "test" for test.case)
       \param the procId determines de CPU id. if negative, it ussemes there is only one processor
    */
    Exporter(const GetPot& dfile, mesh_ptrtype mesh, const std::string prefix, const int procId );

    //! Constructor for Exporter without prefix and procID
    /*!
        In this case prefix and procID should be set separately
        \param dfile the GetPot data file where you must provide an [exporter] section with:
          "start"     (start index for filenames 0 for 000, 1 for 001 etc.),
          "save"      (how many time steps per postprocessing)
          "multimesh" ( = true if the mesh has to be saved at each post-processing step)
       \param the prefix for the case file (ex. "test" for test.case)
    */
    Exporter(const GetPot& dfile, const std::string prefix);

    //! Destructor
    virtual ~Exporter() {};
    //@}

    //! @name Methods
    //@{

    //! Adds a new variable to be post-processed
    /*!
        \param type the type fo the variable Ensight::Scalar or Ensight::Vector
        \param prefix the prefix of the files storing the variable (ex: "velocity" for velocity.***)
        \param vr an ublas::vector_range type given a view of the varialbe (ex: subrange(fluid.u(),0,3*dimU) )
        \param size the number of Dof for that variable
    */
    void addVariable(const ExporterData::Type type, const std::string variableName, vector_ptrtype const& map, UInt start, UInt size, UInt steady =0 );

    //! Post-process the variables added to the list
    /*!
        \param time the solver time
    */
    virtual void postProcess(const Real& time) = 0;

    //! Import data from previous simulations
    /*!
       \param time the solver time
    */
    virtual void import(const Real& Tstart, const Real& dt) = 0; // dt is used to rebuild the history up to now

    //! Read  only last timestep
    virtual void import(const Real& Tstart) = 0;

    //! Set the output folder for postprocessing
    /*!
     * @param outputDirectory output folder
     */
    void setOutputDirectory( const std::string& outputDirectory );

    void initMeshProcId       ( mesh_ptrtype mesh, int const procId );


protected:

    void setNodesMap       ( std::vector<int> LtGNodesMap );

    //! compute postfix
    void computePostfix();

    //@}

public:

    //! @name Set Methods
    //@{

    virtual void setMeshProcId( mesh_ptrtype mesh, int const procId );

    //@}

    //! @name Get Methods
    //@{

    //! returns the type of the map to use for the EpetraVector
    virtual EpetraMapType mapType() const = 0;

    //@}

private:
    void initProcId       ( int const procId );

    void initNodesMap       ();

protected:

    mesh_ptrtype M_mesh;
    std::string M_prefix;
    std::string M_post_dir;
    std::string M_import_dir;
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

    std::vector<int> M_LtGNodesMap;
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
    M_import_dir(dfile("exporter/import_dir", "./")),
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
    M_import_dir(dfile("exporter/import_dir", "./")),
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

    initProcId(procId);

    initNodesMap();
}

template<typename Mesh>
void Exporter<Mesh>::initProcId( int const procId )
{
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

template<typename Mesh>
void Exporter<Mesh>::setNodesMap( std::vector<int> LtGNodesMap )
{
    M_LtGNodesMap = LtGNodesMap;
}

template<typename Mesh>
void Exporter<Mesh>::setOutputDirectory( const std::string& outputDirectory )
{
    M_post_dir = outputDirectory;
}

template<typename Mesh>
void Exporter<Mesh>::initNodesMap()
{
    UInt nVert = M_mesh->numVertices();
    M_LtGNodesMap.resize(nVert);
    for (UInt i=0; i<nVert; ++i)
        {
            int id = this->M_mesh->pointList( i+1 ).id();
            M_LtGNodesMap[i] = id;
        }

}

template <typename Mesh> void Exporter<Mesh>::computePostfix()
{

    std::ostringstream index;
    index.fill( '0' );

    if (M_count % M_save == 0)
        {
            index << std::setw(5) << ( M_count / M_save );

            M_postfix = "." + index.str();

        }
    else
        M_postfix = "*****";

    ++M_count;
}


template<typename Mesh>
void Exporter<Mesh>::addVariable(const ExporterData::Type type, const std::string variableName, vector_ptrtype const& vr, UInt start, UInt size, UInt steady)
{
    M_listData.push_back( ExporterData(type,variableName,vr,start, size, steady) );
}


}


#endif
