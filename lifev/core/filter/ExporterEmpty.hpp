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

#include <lifev/core/filter/Exporter.hpp>

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
    typedef typename super::commPtr_Type commPtr_Type;


    //@}

    //! @name Static members

    //! returns the type of the map to use for the VectorEpetra
    static MapEpetraType const MapType;
    //@}

    //! @name Constructors
    //@{
    ExporterEmpty();
    ExporterEmpty (const GetPot& dfile, meshPtr_Type mesh, const std::string& prefix, const int& procId);
    ExporterEmpty (const GetPot& dfile, const std::string& prefix);
    //@}

    //! @name Public methods
    //@{
    void postProcess (const Real& /*time*/) {}
    void exportPID ( meshPtr_Type /*meshPart*/, commPtr_Type /*comm*/, const bool /*binaryFormat*/ = false ) {}
    void exportRegionMarkerID ( std::shared_ptr<MeshType> /*mesh*/, std::shared_ptr<Epetra_Comm> /*comm*/  ) {}
    UInt importFromTime ( const Real& /*time*/ )
    {
        assert (false);
        return 0;
    }
    void import (const Real& /*startTime*/, const Real& /*dt*/) {}
    void import (const Real& /*startTime*/) {}
    //@}

    //! @name Get methods
    //@{
    MapEpetraType mapType() const;
    //@}

private:
    //! @name Private methods
    //@{
    virtual void readScalar ( ExporterData<MeshType>& /*dvar*/) {}
    virtual void readVector ( ExporterData<MeshType>& /*dvar*/) {}
    //@}

};

// ======================
// IMPLEMENTATION
// ======================

template<typename MeshType>
MapEpetraType const ExporterEmpty<MeshType>::MapType (Unique);


// ======================
// Constructors
// ======================

template<typename MeshType>
ExporterEmpty<MeshType>::ExporterEmpty() :
    super()
{
}

template<typename MeshType>
ExporterEmpty<MeshType>::ExporterEmpty (const GetPot& dfile, meshPtr_Type mesh, const std::string& prefix,
                                        const int& procId) :
    super (dfile, prefix)
{
    this->setMeshProcId (mesh, procId);
}

template<typename MeshType>
ExporterEmpty<MeshType>::ExporterEmpty (const GetPot& dfile, const std::string& prefix) :
    super (dfile, prefix)
{
}

// ====================
// Get methods
// ====================

template<typename MeshType>
MapEpetraType ExporterEmpty<MeshType>::mapType() const
{
    return MapType;
}

} // Namespace LifeV

#endif // EXPORTER_EMPTY_H
