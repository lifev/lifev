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
  @brief Cardiac Electrophysiology Test
  @author Lucia Mirabella <lucia.mirabella@mail.polimi.it> and Mauro Perego <mauro.perego@polimi.it>
  @date 11-2007
  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
  @last update 11-2010
 */

#ifndef __HEART_H
#define __HEART_H

#define MONODOMAIN

#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifemesh/partitionMesh.hpp>

#ifdef MONODOMAIN
#include <life/lifesolver/HeartMonodomainSolver.hpp>
#else
#include <life/lifesolver/HeartBidomainSolver.hpp>
#endif
#include <life/lifesolver/HeartIonicSolver.hpp>
#include <life/lifefilters/ensight.hpp>
#include <life/lifefilters/hdf5exporter.hpp>
#include <life/lifefilters/noexport.hpp>



namespace LifeV
{
/*!
  \class Heart

  3D Action potential propagation class


 */
class Heart
{
public:

    //! @name Typedefs
    //@{

#ifdef MONODOMAIN
    typedef HeartMonodomainSolver< RegionMesh3D<LinearTetra> >::vector_Type  	vector_Type;
    typedef HeartMonodomainSolver<RegionMesh3D<LinearTetra> >::matrix_Type      	matrix_Type;
#else
    typedef HeartBidomainSolver< RegionMesh3D<LinearTetra> >::vector_Type  	vector_Type;
    typedef HeartBidomainSolver<RegionMesh3D<LinearTetra> >::matrix_Type    	matrix_Type;
#endif
    typedef boost::shared_ptr<vector_Type> 					vectorPtr_Type;
    typedef boost::shared_ptr<matrix_Type>     				matrixPtr_Type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Heart( Int argc,
           char** argv );

    virtual ~Heart() {}

    //@}

    /** @name  Methods
     */
    //@{

    //! To build the system and iterate
    void run();

    //! To compute the righthand side of the system
#ifdef MONODOMAIN
    void computeRhs( vector_Type& rhs,
                     HeartMonodomainSolver< RegionMesh3D<LinearTetra> >& electricModel,
                     boost::shared_ptr< HeartIonicSolver< RegionMesh3D<LinearTetra> > > ionicModel,
                     HeartMonodomainData& dataMonodomain );
#else
    void computeRhs( vector_Type& rhs,
                     HeartBidomainSolver< RegionMesh3D<LinearTetra> >& electricModel,
                     boost::shared_ptr< HeartIonicSolver< RegionMesh3D<LinearTetra> > > ionicModel,
                     HeartBidomainData& dataBidomain );
#endif
    //@}


private:
    UInt ion_model;
    UInt nbeq;
    boost::shared_ptr<HeartFunctors> M_heart_fct;
};
}
#endif /* __HEART_H */
