//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

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

    //! Where is data centered /
    enum Where
    {
        Node, /*!< Node centered */
        Cell  /*!< Cell centered */
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
        @param where        - where is variable located (Node or Cell)
    */
  ExporterData(const Type type, const std::string variableName, vector_ptrtype& vec, UInt start, UInt size, UInt steady, const Where where=ExporterData::Node);

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

    //! size of the stored array
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

    //! Node or Cell centered ?
    const Where& where() const;

    //! returns Node or Cell centered string
    std::string whereName() const;

    //@}


private:

    //! @name Private Methods
    //@{

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

    //! Node or Cell centered
    Where M_where;
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

    //! Empty constructor for Exporter
    Exporter();

    //! Constructor for Exporter without prefix and procID
    /*!
        In this case prefix and procID should be set separately
        @param dfile the GetPot data file where you must provide an [exporter] section with:
          "start"     (start index for filenames 0 for 000, 1 for 001 etc.),
          "save"      (how many time steps per postprocessing)
          "multimesh" ( = true if the mesh has to be saved at each post-processing step)
       @param the prefix for the case file (ex. "test" for test.case)
    */
    Exporter(const GetPot& dfile, const std::string& prefix);

    //! Destructor
    virtual ~Exporter() {};
    //@}


    //! @name Methods
    //@{

    //! Adds a new variable to be post-processed
    /*!
        @param type the type fo the variable Ensight::Scalar or Ensight::Vector
        @param prefix the prefix of the files storing the variable (ex: "velocity" for velocity.***)
        @param vr an ublas::vector_range type given a view of the varialbe (ex: subrange(fluid.u(),0,3*dimU) )
        @param size size of the stored array
    */
  void addVariable(const ExporterData::Type type, const std::string variableName, vector_ptrtype& vector, UInt start, UInt size,
                   UInt steady = 0, ExporterData::Where where = ExporterData::Node );

    //! Post-process the variables added to the list
    /*!
        @param time the solver time
    */
    virtual void postProcess(const Real& time) = 0;

    //! Import data from previous simulations at a certain time
    /*!
       @param Time the time of the data to be imported
     */
    virtual UInt importFromTime( const Real& Time ) = 0;

    //! Import data from previous simulations
    /*!
       @param time the solver time
       @param dt time step used to rebuild the history up to now
    */
    virtual void import(const Real& Tstart, const Real& dt) = 0;

    //! Read  only last timestep
    virtual void import(const Real& Tstart) = 0;

    virtual void rd_var(ExporterData& dvar);

protected:

    //! compute postfix
    void computePostfix();

    virtual void M_rd_scalar( ExporterData& dvar ) = 0;
    virtual void M_rd_vector( ExporterData& dvar ) = 0;

    //@}


    //! @name Set Methods
    //@{

public:

    //! Set data from file.
    /*!
     * @param dataFile data file.
     * @param section section in the data file.
     */
    void setDataFromGetPot( const GetPot& dataFile, const std::string& section = "exporter" );

    //! Set prefix.
    /*!
     * @param prefix prefix.
     */
    void setPrefix( const std::string& prefix );

    //! Set the folder for pre/postprocessing
    /*!
     * @param Directory output folder
     */
    void setDirectory( const std::string& Directory );

    //! Set the folder for pre/postprocessing
    /*!
     * @param Directory output folder
     */
    void setStartIndex( const UInt& StartIndex );

    //! Set how many time step between two saves.
    /*!
     * @param save steps
     */
    void setSave( const UInt& save );

    //! Set if to save the mesh at each time step.
    /*!
     * @param multimesh multimesh
     */
    void setMultimesh( const bool& multimesh );

    virtual void setMeshProcId( const mesh_ptrtype mesh, const int& procId );

    //@}


    //! @name Get Methods
    //@{

    const UInt& getStartIndex( void );

    //! returns the type of the map to use for the EpetraVector
    virtual EpetraMapType mapType() const = 0;

    //@}

protected:

    std::string                 M_prefix;
    std::string                 M_post_dir;
    UInt                        M_count;
    UInt                        M_save;
    bool                        M_multimesh;
    mesh_ptrtype                M_mesh;
    int                         M_procId;
    std::string                 M_postfix;

    std::list<ExporterData>     M_listData;
    std::list<Real>             M_timeSteps;
};



// ===================================================
// Constructors
// ===================================================
template<typename Mesh>
Exporter<Mesh>::Exporter():
    M_prefix        ( "output"),
    M_post_dir      ( "./" ),
    M_count         ( 0 ),
    M_save          ( 1 ),
    M_multimesh     ( true )
{}

template<typename Mesh>
Exporter<Mesh>::Exporter( const GetPot& dfile, const std::string& prefix ):
    M_prefix        ( prefix ),
    M_post_dir      ( dfile("exporter/post_dir", "./") ),
    M_count         ( dfile("exporter/start",0) ),
    M_save          ( dfile("exporter/save",1) ),
    M_multimesh     ( dfile("exporter/multimesh",true) )
{}

// ===================================================
// Methods
// ===================================================
template<typename Mesh>
void Exporter<Mesh>::addVariable(const ExporterData::Type type,
                                 const std::string variableName,
                                 vector_ptrtype& vr,
                                 UInt start,
                                 UInt size,
                                 UInt steady,
                                 ExporterData::Where where)
{
  M_listData.push_back( ExporterData(type,variableName,vr,start, size, steady, where) );
}

template <typename Mesh>
void Exporter<Mesh>::rd_var(ExporterData& dvar)
{
    switch( dvar.type() )
    {
    case ExporterData::Scalar:
        M_rd_scalar(dvar);
        break;
    case ExporterData::Vector:
        M_rd_vector(dvar);
        break;
    }
}

template <typename Mesh>
void Exporter<Mesh>::computePostfix()
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

// ===================================================
// Set Methods
// ===================================================
template<typename Mesh>
void Exporter<Mesh>::setDataFromGetPot( const GetPot& dataFile, const std::string& section )
{
    M_post_dir      = dataFile( ( section + "/post_dir"  ).data(), "./" );
    M_count         = dataFile( ( section + "/start"     ).data(), 0 );
    M_save          = dataFile( ( section + "/save"      ).data(), 1 );
    M_multimesh     = dataFile( ( section + "/multimesh" ).data(), true );
}

template<typename Mesh>
void Exporter<Mesh>::setPrefix( const std::string& prefix )
{
    M_prefix = prefix;
}

template<typename Mesh>
void Exporter<Mesh>::setDirectory( const std::string& Directory )
{
    M_post_dir = Directory;
}

template<typename Mesh>
void Exporter<Mesh>::setStartIndex( const UInt& StartIndex )
{
    M_count = StartIndex;
}

template<typename Mesh>
void Exporter<Mesh>::setSave( const UInt& save )
{
    M_save = save;
}

template<typename Mesh>
void Exporter<Mesh>::setMultimesh( const bool& multimesh )
{
    M_multimesh = multimesh;
}

template<typename Mesh>
void Exporter<Mesh>::setMeshProcId( const mesh_ptrtype mesh , const int& procId )
{
    M_mesh   = mesh;
    M_procId = procId;
}

// ===================================================
// Get Methods
// ===================================================
template<typename Mesh>
const UInt& Exporter<Mesh>::getStartIndex( void )
{
    return M_count;
}

}

#endif
