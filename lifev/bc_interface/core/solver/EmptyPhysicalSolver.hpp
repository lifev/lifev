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
 *  @file
 *  @brief File containing  a default physical solver in orde to use the bc interface
 *
 *  @date 04 - 2014
 *  @author Simone Rossi <simone.rossi@epfl.ch>
 *
 *  @maintainer Simone Rossi <simone.rossi@epfl.ch>
 */
#ifndef EMPTYPHYSICALSOLVER_HPP_
#define EMPTYPHYSICALSOLVER_HPP_

#include <boost/shared_ptr.hpp>
//#include <lifev/core/array/VectorEpetra.hpp>

template < class vector_Type >
class EmptyPhysicalSolver
{
public:
	typedef boost::shared_ptr<vector_Type> solutionPtr_Type;
};


#endif /* DEFAULTPHYSICALSOLVER_HPP_ */
