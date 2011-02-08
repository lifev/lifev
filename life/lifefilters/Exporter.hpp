//@HEADER
/*
*******************************************************************************

Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

This file is part of LifeV.

LifeV is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LifeV is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER
/*!
  @file
  @brief exporter and ExporterData classes provide interfaces for post-processing

  @date 11-11-2008
  @author Simone Deparis <simone.deparis.epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
  @contributor Tiziano Passerini <tiziano@mathcs.emory.edu>

    Usage: two steps
    <ol>
        <li> first: add the variables using addVariable
        <li> second: call postProcess( time );
    </ol>
*/

#ifndef EXPORTER_H
#define EXPORTER_H 1

#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifearray/VectorEpetra.hpp>
#include <life/lifefilters/GetPot.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <life/lifecore/LifeV.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifemesh/MarkerDefinitions.hpp>

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
template< typename MeshType >
class ExporterData
{
public:

    //! @name Public Types
    //@{

    typedef MeshType                          mesh_Type;
    typedef VectorEpetra                      vector_Type;
    typedef boost::shared_ptr<vector_Type>    vectorPtr_Type;
    typedef FESpace< mesh_Type, MapEpetra >   feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>   feSpacePtr_Type;

    //! FieldTypeEnum of data stored.
    enum FieldTypeEnum
    {
        ScalarField, /*!< ScalarField stands for scalar field */
        VectorField, /*!< VectorField stands for vector field */
    };

    //! Where is data centered? /
    enum WhereEnum
    {
        Node, /*!< Node centered */
        Cell  /*!< Cell centered */
    };

    //! Time regime of the field /
    enum FieldRegimeEnum
    {
        UnsteadyRegime = 0, /*!< The field is in unsteady regime */
        SteadyRegime = 1, /*!< The field is in steady regime */
        NullRegime = 2 /*!< DEPRECATED */
    };

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Constructor with all the data to be stored
    /*!
        Constructor with all the data to be stored
        @param fieldType    - scalar or vector field
        @param variableName - name assigned to this variable in output file
        @param feSpacePtr   - shared pointer to variable FESpace
        @param vectorPtr    - shared pointer to variable array
        @param start        - address of first datum in the array to be stored.
                              Useful in case you want to define a subrange of the vector *vec
        @param regime       - if UnsteadyRegime, then the file name for postprocessing has to include time dependency
        @param where        - where is variable located (Node or Cell)
    */
    ExporterData(const FieldTypeEnum&   fieldType,
                 const std::string&     variableName,
                 const feSpacePtr_Type& feSpacePtr,
                 const vectorPtr_Type&  vectorPtr,
                 const UInt&            start,
                 const FieldRegimeEnum& regime,
                 const WhereEnum&       where = Node);

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
    void setRegime(FieldRegimeEnum regime) {M_regime = regime;}

    //@}


    //! @name Get Methods
    //@{

    //! name assigned to this variable in output file
    const std::string& variableName() const
    {
        return M_variableName;
    }

    //! number of (scalar) DOFs of the variable
    const UInt& __attribute__ (( deprecated )) size() const
    {
        return numDOF();
    }
    const UInt& numDOF() const
    {
        return M_numDOF;
    }

    //! address of first datum in the array
    const UInt& start() const { return M_start; }

    //! scalar or vector field
    const FieldTypeEnum& __attribute__ (( deprecated )) type() const
    {
        return fieldType();
    }
    const FieldTypeEnum& fieldType() const
    {
        return M_fieldType;
    }

    //! shared pointer to array
    const vectorPtr_Type storedArrayPtr() const
    {
        return M_storedArrayPtr;
    }

    //! returns UnsteadyRegime (=0) if file name for postprocessing has to include time dependency
    FieldRegimeEnum regime() const {return M_regime; }

    //! returns Scalar or Vector strings
    std::string typeName() const;

    //! returns 1 (if Scalar) or 3 (if Vector)
    UInt typeDim() const;

    //! Node or Cell centered ?
    const WhereEnum& where() const
    {
        return M_where;
    }

    //! returns Node or Cell centered string
    std::string whereName() const;

    const feSpacePtr_Type& feSpacePtr() const {return M_feSpacePtr;};
    //@}

private:
    //! @name Private data members
    //@{
    //! name assigned to this variable in output file
    std::string M_variableName;

    //! pointer to the FESpace of the variable
    feSpacePtr_Type M_feSpacePtr;

    //! pointer to storedArray
    vectorPtr_Type M_storedArrayPtr;

    //! number of (scalar) DOFs of the variable
    UInt M_numDOF;

    //! address of first datum in the array
    UInt M_start;

    //! scalar or vector field
    FieldTypeEnum M_fieldType;

    //! equal to UnsteadyRegime if file name for postprocessing has to include time dependency
    FieldRegimeEnum M_regime;

    //! Node or Cell centered
    WhereEnum M_where;
    //@}
};


//! Exporter - Pure virtual class that describes a generic exporter
/*!
    @author Simone Deparis <simone.deparis@epfl.ch>
    @date 11-11-2008

    This class is pure virtual and describes a generic exporter that can
    also do import
 */
template<typename MeshType>
class Exporter
{

public:
    //! @name Public typedefs
    //@{
    typedef MeshType                                    mesh_Type;
    typedef boost::shared_ptr<MeshType>                 meshPtr_Type;
    typedef ExporterData<mesh_Type>                     exporterData_Type;
    typedef typename exporterData_Type::vector_Type     vector_Type;
    typedef typename exporterData_Type::vectorPtr_Type  vectorPtr_Type;
    typedef typename exporterData_Type::feSpacePtr_Type feSpacePtr_Type;
    typedef typename exporterData_Type::WhereEnum       WhereEnum;
    typedef typename exporterData_Type::FieldTypeEnum   FieldTypeEnum;
    typedef typename exporterData_Type::FieldRegimeEnum FieldRegimeEnum;
    typedef typename std::vector<exporterData_Type >    dataVector_Type;
    typedef typename dataVector_Type::iterator          dataVectorIterator_Type;
    typedef typename std::multimap<WhereEnum, UInt >    whereToDataIdMap_Type;
    typedef typename std::multimap<FE_TYPE, UInt >      feTypeToDataIdMap_Type;
    //@}

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
        @param type the type of the variable exporterData_Type::FieldTypeEnum
        @param variableName the name of the variable (to be used in file prefix)
        @param feSpacePtr a pointer to the FESpace of the variable
        @param vectorPtr a pointer to the vector containing the values of the variable
        @param start location in the vector where the storing of the variable starts
        @param regime if UnsteadyRegime the filename should change at each time step
        @param where choose whether the variable is defined on Nodes of Elements
    */
    void addVariable(const FieldTypeEnum& type,
                     const std::string& variableName,
                     const feSpacePtr_Type& feSpacePtr,
                     const vectorPtr_Type& vectorPtr,
                     const UInt& start,
                     const FieldRegimeEnum& regime = exporterData_Type::UnsteadyRegime,
                     const WhereEnum& where = exporterData_Type::Node );

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
    virtual void import(const Real& startTime, const Real& dt) = 0;

    //! Read  only last timestep
    virtual void import(const Real& startTime) = 0;

    virtual void readVariable(exporterData_Type& dvar);
    //@}

    //! @name Set Methods
    //@{

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
    void setPrefix( const std::string& prefix )
    {
        M_prefix = prefix;
    }

    //! Set the folder for pre/postprocessing
    /*!
     * @param Directory output folder
     */
    void setPostDir( const std::string& Directory )
    {
        M_postDir = Directory;
    }

    //! Set the folder for pre/postprocessing
    /*!
     * @param Directory output folder
     */
    void setTimeIndexStart( const UInt& StartIndex )
    {
        M_timeIndexStart = StartIndex;
    }

    //! Set how many time step between two saves.
    /*!
     * @param save steps
     */
    void setSave( const UInt& save )
    {
        M_save = save;
    }

    //! Set if to save the mesh at each time step.
    /*!
     * @param multimesh multimesh
     */
    void setMultimesh( const bool& multimesh )
    {
        M_multimesh = multimesh;
    }

    virtual void setMeshProcId( const meshPtr_Type mesh, const int& procId );

    //! Close the output file
    /*!
         This method is only used by  some of the exporter which derived from this class.
     */
    virtual void closeFile() {}
    //@}

    //! @name Get Methods
    //@{

    const UInt& timeIndexStart() const { return M_timeIndexStart; }
    const UInt& timeIndex() const { return M_timeIndex; }

    //! returns the type of the map to use for the VectorEpetra
    virtual MapEpetraType mapType() const = 0;
    //@}

protected:

    //! @name Protected methods
    //@{
    //! compute postfix
    void computePostfix();

    virtual void readScalar( exporterData_Type& dvar ) = 0;
    virtual void readVector( exporterData_Type& dvar ) = 0;

    //@}

    //! @name Protected data members
    //@{
    std::string                 M_prefix;
    std::string                 M_postDir;
    UInt                        M_timeIndexStart;
    UInt                        M_timeIndex;
    UInt                        M_save;
    bool                        M_multimesh;
    UInt                        M_timeIndexWidth;
    meshPtr_Type                M_mesh;
    int                         M_procId;
    std::string                 M_postfix;

    whereToDataIdMap_Type       M_whereToDataIdMap;
    feTypeToDataIdMap_Type      M_feTypeToDataIdMap;
    dataVector_Type             M_dataVector;
    std::list<Real>             M_timeSteps;
    //@}
};

// ==================================================
// EXPORTERDATA: IMPLEMENTATION
// ==================================================

// =================
// Constructor
// =================

template< typename MeshType >
ExporterData<MeshType>::ExporterData( const FieldTypeEnum&   type,
                                      const std::string&     variableName,
                                      const feSpacePtr_Type& feSpacePtr,
                                      const vectorPtr_Type&  vectorPtr,
                                      const UInt&            start,
                                      const FieldRegimeEnum& regime,
                                      const WhereEnum&       where ):
        M_variableName      ( variableName ),
        M_feSpacePtr        ( feSpacePtr ),
        M_storedArrayPtr    ( vectorPtr ),
        M_numDOF            ( feSpacePtr->dim() ),
        M_start             ( start ),
        M_fieldType         ( type ),
        M_regime            ( regime ),
        M_where             ( where )
{}

// ==============
// Operators
// ==============

template< typename MeshType >
Real ExporterData<MeshType>::operator()( const UInt i ) const
{
    return (*M_storedArrayPtr)[i];
}

template< typename MeshType >
Real& ExporterData<MeshType>::operator()( const UInt i )
{
    return (*M_storedArrayPtr)[i];
}

template< typename MeshType >
std::string ExporterData<MeshType>::typeName() const
{
    switch (M_fieldType)
    {
    case ScalarField:
        return "Scalar";
    case VectorField:
        return "Vector";
    }

    return "ERROR string";
}

template< typename MeshType >
UInt ExporterData<MeshType>::typeDim() const
{
    /*switch ( M_fieldType )
    {
    case ScalarField:
        return 1;
    case VectorField:
        return 3;
    }*/

    return M_feSpacePtr->fieldDim();
}

template< typename MeshType >
std::string ExporterData<MeshType>::whereName() const
{
    switch (M_where)
    {
    case Node:
        return "Node";
    case Cell:
        return "Cell";
    }

    return "ERROR string";
}

// ==================================================
// EXPORTER: IMPLEMENTATION
// ==================================================

// ===================================================
// Constructors
// ===================================================
template<typename MeshType>
Exporter<MeshType>::Exporter():
        M_prefix        ( "output"),
        M_postDir       ( "./" ),
        M_timeIndexStart( 0 ),
        M_timeIndex     ( M_timeIndexStart ),
        M_save          ( 1 ),
        M_multimesh     ( true ),
        M_timeIndexWidth( 5 )
{}

template<typename MeshType>
Exporter<MeshType>::Exporter( const GetPot& dfile, const std::string& prefix ):
        M_prefix        ( prefix ),
        M_postDir       ( dfile("exporter/post_dir", "./") ),
        M_timeIndexStart( dfile("exporter/start",0) ),
        M_timeIndex     ( M_timeIndexStart ),
        M_save          ( dfile("exporter/save",1) ),
        M_multimesh     ( dfile("exporter/multimesh",true) ),
        M_timeIndexWidth( dfile("exporter/time_id_width",5) )
{}

// ===================================================
// Methods
// ===================================================
template<typename MeshType>
void Exporter<MeshType>::addVariable(const FieldTypeEnum& type,
                                     const std::string& variableName,
                                     const feSpacePtr_Type& feSpacePtr,
                                     const vectorPtr_Type& vectorPtr,
                                     const UInt& start,
                                     const FieldRegimeEnum& regime,
                                     const WhereEnum& where )
{
    M_dataVector.push_back( exporterData_Type(type, variableName, feSpacePtr, vectorPtr, start, regime, where) );
    M_whereToDataIdMap.insert( std::pair<WhereEnum,UInt >(where, M_dataVector.size()-1 ) );
    M_feTypeToDataIdMap.insert( std::pair<FE_TYPE,UInt >(feSpacePtr->fe().refFE().type(), M_dataVector.size()-1 ) );
}

template <typename MeshType>
void Exporter<MeshType>::readVariable(exporterData_Type& dvar)
{
    switch ( dvar.type() )
    {
    case exporterData_Type::ScalarField:
        readScalar(dvar);
        break;
    case exporterData_Type::VectorField:
        readVector(dvar);
        break;
    }
}

template <typename MeshType>
void Exporter<MeshType>::computePostfix()
{
    std::ostringstream index;
    index.fill( '0' );

    if (M_timeIndex % M_save == 0)
    {
        index << std::setw(M_timeIndexWidth) << ( M_timeIndex / M_save );

        M_postfix = "." + index.str();
    }
    else
    {
        // M_postfix = "*****";
        std::string stars("");
        for (UInt cc(0); cc<M_timeIndexWidth; ++cc) stars+="*";
        M_postfix = stars;
    }

    ++M_timeIndex;
}

// ===================================================
// Set Methods
// ===================================================
template<typename MeshType>
void Exporter<MeshType>::setDataFromGetPot( const GetPot& dataFile, const std::string& section )
{
    M_postDir        = dataFile( ( section + "/post_dir"  ).data(), "./" );
    M_timeIndexStart = dataFile( ( section + "/start"     ).data(), 0 );
    M_timeIndex      = M_timeIndexStart;
    M_save           = dataFile( ( section + "/save"      ).data(), 1 );
    M_multimesh      = dataFile( ( section + "/multimesh" ).data(), true );
    M_timeIndexWidth = dataFile( ( section + "/time_id_width" ).data(), 5);
}

template<typename MeshType>
void Exporter<MeshType>::setMeshProcId( const meshPtr_Type mesh , const int& procId )
{
    M_mesh   = mesh;
    M_procId = procId;
}

} // Namespace LifeV

#endif // EXPORTER_H
