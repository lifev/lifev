//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

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
 * @file heart.hpp
 * @brief Cardiac Electrophysiology Test
 * @author Lucia Mirabella <lucia.mirabella@mail.polimi.it> and Mauro Perego <mauro.perego@polimi.it>
 * @date 11-2007
 * @contributors Ricardo Ruiz-Baier
 * @last update 11-2010
 */

#ifndef __HEART_H
#define __HEART_H

#define BIDOMAIN

#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>
#ifdef MONODOMAIN
	#include <life/lifesolver/monodomainSolver.hpp>
#else
	#include <life/lifesolver/bidomainSolver.hpp>
#endif
#include <life/lifesolver/ionicSolver.hpp>
#include <life/lifefilters/ensight.hpp>
#include <life/lifefilters/hdf5exporter.hpp>
#include <life/lifefilters/noexport.hpp>

//#include <life/lifesolver/heartFunctors.hpp>




namespace LifeV
{
/*!
 * \class Heart
 * \brief 3D Action potential propagation class
 *
 *  @author Lucia Mirabella
 *  @see
 */
class Heart
//     :
//     public LifeV::Application
{
public:


    /** @name Typedefs
     */
    //@{

#ifdef MONODOMAIN
	typedef MonodomainSolver< RegionMesh3D<LinearTetra> >::vector_type  	vector_type;
	typedef MonodomainSolver<RegionMesh3D<LinearTetra> >::matrix_type      	matrix_type;
#else
    typedef BidomainSolver< RegionMesh3D<LinearTetra> >::vector_type  	vector_type;
    typedef BidomainSolver<RegionMesh3D<LinearTetra> >::matrix_type      	matrix_type;
#endif
    typedef boost::shared_ptr<vector_type> 					vector_ptrtype;
    typedef boost::shared_ptr<matrix_type>        				matrix_ptrtype;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Heart( int argc,
              char** argv );

    ~Heart()
        {}

    //@}

    /** @name  Methods
     */
    //@{

    //! To build the system and iterate
    void run();

    //! To compute the righthand side of the system
#ifdef MONODOMAIN
    void computeRhs( vector_type& rhs,
    		MonodomainSolver< RegionMesh3D<LinearTetra> >& electricModel,
    		boost::shared_ptr< IonicSolver< RegionMesh3D<LinearTetra> > > ionicModel,
    		DataMonodomain& dataMonodomain );
#else
    void computeRhs( vector_type& rhs,
    		BidomainSolver< RegionMesh3D<LinearTetra> >& electricModel,
    		boost::shared_ptr< IonicSolver< RegionMesh3D<LinearTetra> > > ionicModel,
    		DataBidomain& dataBidomain );
#endif
    //@}


private:
	UInt ion_model;
	UInt nbeq;
    boost::shared_ptr<HeartFunctors> M_heart_fct;

};
}
#endif /* __HEART_H */
