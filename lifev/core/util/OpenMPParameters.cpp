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
  @brief Class that handles setting OpenMP parameters

  @date 2013-03
  @author Radu Popescu <radu.popescu@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#include <lifev/core/util/OpenMPParameters.hpp>

namespace LifeV
{

OpenMPParameters::OpenMPParameters()
    : numThreads (1), chunkSize (0)
{
#ifdef _OPENMP
    scheduler = omp_sched_static;
#endif
}

void OpenMPParameters::apply()
{
#ifdef _OPENMP
    omp_set_num_threads (numThreads);
    omp_set_schedule (scheduler, chunkSize);
#endif
}

} // namespace LifeV
