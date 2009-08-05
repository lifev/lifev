/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): simone <simone@mac-simao.local>
       Date: 2008-12-27

  Copyright (C) 2008 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file noexport.hpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2008-12-27
  \brief This file provides an interface as a fake post-processing.

 */


#ifndef _NOEXPORT_H_
#define _NOEXPORT_H_

#include <life/lifefilters/exporter.hpp>


namespace LifeV
{

/**
 * \class Ensight
 * \brief Ensight data exporter
 */
template<typename Mesh>
class NoExport : public Exporter<Mesh> {

public:

    typedef Exporter<Mesh> super;
    typedef typename super::mesh_ptrtype  mesh_ptrtype;
    typedef typename super::vector_ptrtype vector_ptrtype;

    /**
       Constructor for NoExport

       \param dfile the GetPot data file where you must provide and [ensight] section with:
       "start" (start index for filenames 0 for 000, 1 for 001 etc.),
       "save" (how many time steps per posptrocessing)
       "multimesh" (=true if the mesh has to be saved at each post-processing step)

       \param mesh the mesh

       \param the prefix for the case file (ex. "test" for test.case)

       \param the procId determines de CPU id. if negative, it ussemes there is only one processor
    */
    NoExport(const GetPot& dfile, mesh_ptrtype mesh, const std::string prefix, const int procId );

    NoExport(const GetPot& dfile, const std::string prefix);


    /**
       Post-porcess the variables added to the list

       \param time the solver time
    */
    void postProcess(const Real& time) {}

    /**
       Import data from previous simulations

       \param time the solver time
    */
    void import(const Real& Tstart, const Real& dt) {} // dt is used to rebuild the history up to now

    //! Read  only last timestep
    void import(const Real& Tstart) {}
};

//
// Implementation
//
template<typename Mesh>
NoExport<Mesh>::NoExport(const GetPot& dfile, mesh_ptrtype mesh, const std::string prefix,
                       int const procId)
    :
    super(dfile, mesh, prefix,procId)
{
    setMeshProcId(mesh,procId);
}

template<typename Mesh>
NoExport<Mesh>::NoExport(const GetPot& dfile, const std::string prefix):
    super(dfile,prefix)
{
}


}

#endif
