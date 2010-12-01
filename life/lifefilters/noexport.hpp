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
    @date 2008-12-27

    This file provides an interface as a fake post-processing.
 */

#ifndef _NOEXPORT_H_
#define _NOEXPORT_H_

#include <life/lifefilters/exporter.hpp>

namespace LifeV
{

template<typename Mesh>
class NoExport : public Exporter<Mesh>
{

public:

    typedef Exporter<Mesh> super;
    typedef typename super::mesh_ptrtype  mesh_ptrtype;
    typedef typename super::vector_ptrtype vector_ptrtype;

    NoExport();
    NoExport(const GetPot& dfile, mesh_ptrtype mesh, const std::string& prefix, const int& procId );
    NoExport(const GetPot& dfile, const std::string& prefix);

    EpetraMapType mapType() const;
    void postProcess(const Real& /*time*/) {}
    UInt importFromTime( const Real& /*Time*/ ) { assert(false); return 0; }
    void import(const Real& /*Tstart*/, const Real& /*dt*/) {} // dt is used to rebuild the history up to now
    void import(const Real& /*Tstart*/) {}

private:
    virtual void M_rd_scalar( ExporterData& /*dvar*/ ) {}
    virtual void M_rd_vector( ExporterData& /*dvar*/ ) {}

};

template<typename Mesh>
NoExport<Mesh>::NoExport():
        super()
{
}

template<typename Mesh>
NoExport<Mesh>::NoExport(const GetPot& dfile, mesh_ptrtype mesh, const std::string& prefix,
                         const int& procId)
        :
        super(dfile,prefix)
{
    setMeshProcId(mesh,procId);
}

template<typename Mesh>
NoExport<Mesh>::NoExport(const GetPot& dfile, const std::string& prefix):
        super(dfile,prefix)
{
}

template<typename Mesh>
EpetraMapType NoExport<Mesh>::mapType() const
{
    return Unique;
}

}

#endif
