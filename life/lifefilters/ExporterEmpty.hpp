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
    @brief This file provides a dummy interface for post-processing

    @date 27-12-2008
    @author Simone Deparis <simone.deparis@epfl.ch>

    @maintainer Radu Popescu <radu.popescu@epfl.ch>
 */

#ifndef EXPORTER_EMPTY_H
#define EXPORTER_EMPTY_H 1

#include <life/lifefilters/Exporter.hpp>

namespace LifeV
{

//! @name Class that provides a dummy interface for post-processing
template<typename MeshType>
class ExporterEmpty : public Exporter<MeshType>
{

public:
    //! @name Public typedefs
    //@{
    typedef MeshType mesh_Type;
    typedef Exporter<MeshType> super;
    typedef typename super::meshPtr_Type  meshPtr_Type;
    typedef typename super::vectorPtr_Type vectorPtr_Type;
    //@}

    //! @name Constructors
    //@{
    ExporterEmpty();
    ExporterEmpty(const GetPot& dfile, meshPtr_Type mesh, const std::string& prefix, const int& procId);
    ExporterEmpty(const GetPot& dfile, const std::string& prefix);
    //@}

    //! @name Public methods
    //@{
    void postProcess(const Real& /*time*/) {}
    UInt importFromTime( const Real& /*time*/ ) { assert(false); return 0; }
    void import(const Real& /*startTime*/, const Real& /*dt*/) {}
    void import(const Real& /*startTime*/) {}
    //@}

    //! @name Get methods
    //@{
    EpetraMapType mapType() const;
    //@}

private:
    //! @name Private methods
    //@{
    virtual void readScalar( ExporterData& /*dvar*/) {}
    virtual void readVector( ExporterData& /*dvar*/) {}
    //@}

};

// ======================
// IMPLEMENTATION
// ======================

// ======================
// Constructors
// ======================

template<typename MeshType>
ExporterEmpty<MeshType>::ExporterEmpty():
        super()
{
}

template<typename MeshType>
ExporterEmpty<MeshType>::ExporterEmpty(const GetPot& dfile, meshPtr_Type mesh, const std::string& prefix,
                         const int& procId):
        super(dfile, prefix)
{
    setMeshProcId(mesh, procId);
}

template<typename MeshType>
ExporterEmpty<MeshType>::ExporterEmpty(const GetPot& dfile, const std::string& prefix):
        super(dfile, prefix)
{
}

// ====================
// Get methods
// ====================

template<typename MeshType>
EpetraMapType ExporterEmpty<MeshType>::mapType() const
{
    return Unique;
}

} // Namespace LifeV

#endif // EXPORTER_EMPTY_H
