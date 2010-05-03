/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Miguel A. Fernandez <miguel.fernandez@inria.fr>
            Christoph Winkelmann <christoph.winkelmann@epfl.ch>
      Date: 2003-06-09

 Copyright (C) 2004 EPFL

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
  \file ionicSolver.hpp
  \author L. Mirabella M. Perego
  \date 11/2007

 \file ionicSolver.hpp
  \modif J.Castelneau
  \date 01/2009

  \brief This file contains a Oseen equation solver class
*/
#ifndef _IONICSOLVER_H_
#define _IONICSOLVER_H_

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>
#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifefem/sobolevNorms.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifesolver/dataIonic.hpp>
#include <boost/shared_ptr.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdf_template.hpp>
#include <fstream>

namespace LifeV
{
/*!
  \class IonicSolver

  This class implements a ionic model solver.
*/


template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class IonicSolver
{

public:

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;

    typedef DataIonic<Mesh> data_type;

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&,
                                 const Real&, const ID& );

    typedef boost::function<Real ( Real const&, Real const&, Real const&,
                                   Real const&, ID const& )> source_type;

    typedef boost::function<Real (const EntityFlag&, const Real&, const Real&, const Real&, const ID&)> fct_TauClose;

    typedef Mesh mesh_type;

    typedef typename SolverType::matrix_type      matrix_type;
    typedef boost::shared_ptr<matrix_type>        matrix_ptrtype;
    typedef typename SolverType::vector_type      vector_type;

    typedef typename SolverType::prec_raw_type    prec_raw_type;
    typedef typename SolverType::prec_type        prec_type;

    Real zero_scalar( const Real& /* t */,
    const Real& /* x */,
    const Real& /* y */,
    const Real& /* z */,
    const ID& /* i */ )
    {
        return 0.0;
    };


    //! Constructor
    /*!
      \param dataType
      \param recovery variable FE space
	  \param Epetra communicator
    */
    IonicSolver( const data_type&          dataType,
           FESpace<Mesh, EpetraMap>& uFEspace,
           Epetra_Comm&              comm );


    //! virtual destructor
    virtual ~IonicSolver();


    //! Returns the recovery variable FE space
    FESpace<Mesh, EpetraMap>& recoveryFESpace()   {return M_uFESpace;}

    //! Return maps
    Epetra_Map const& getRepeatedEpetraMap() const { return *M_localMap.getMap(Repeated); }

    EpetraMap const& getMap() const { return M_localMap; }

    virtual void updateRepeated( )=0;

    //! Update the ionic model elvecs
    virtual void updateElvec( UInt eleIDw )=0;

    //! Solves the ionic model
    virtual void ionModelSolve( const vector_type& u, const Real timestep )=0;

    //! Computes the term Iion
    virtual void computeIion(ElemVec& elvec, ElemVec& elvec_u, FESpace<Mesh, EpetraMap>& uFESpace )=0;

    //! Initialize
    virtual void initialize( )=0;

    const data_type&               M_data;

    // FE space
    FESpace<Mesh, EpetraMap>&      M_uFESpace;

    //! MPI communicator
    Epetra_Comm*                   M_comm;

    int                            M_me;

    //! Map
    EpetraMap                      M_localMap;

    //! Boolean that indicates if output is sent to cout
    bool                           M_verbose;


protected:

    UInt dim_u() const           { return M_uFESpace.dim(); }

}; // class IonicSolver


//
// IMPLEMENTATION
//

//! Constructor
template<typename Mesh, typename SolverType>
IonicSolver<Mesh, SolverType>::
IonicSolver( const data_type&          dataType,
       FESpace<Mesh, EpetraMap>& uFEspace,
       Epetra_Comm&              comm ):
    M_data                   ( dataType ),
    M_uFESpace               ( uFEspace ),
    M_comm                   ( &comm ),
    M_me                     ( M_comm->MyPID() ),
    M_localMap               ( M_uFESpace.map() ),
    M_verbose                ( M_me == 0)
{
}

template<typename Mesh, typename SolverType>
IonicSolver<Mesh, SolverType>::
~IonicSolver()
{
}


///////////////////////

template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class Rogers_McCulloch : public virtual IonicSolver<Mesh, SolverType>
{
public:
	typedef typename IonicSolver<Mesh, SolverType>::data_type data_type;
	typedef typename IonicSolver<Mesh, SolverType>::vector_type vector_type;
	typedef typename IonicSolver<Mesh, SolverType>::Function Function;

    Rogers_McCulloch( const data_type&          dataType,
           FESpace<Mesh, EpetraMap>& uFEspace,
           Epetra_Comm&              comm );
    virtual ~Rogers_McCulloch();

    //! Update the ionic model elvecs
    void updateRepeated( );

    //! Update the ionic model elvecs
    void updateElvec( UInt eleID );

    //! Solves the ionic model
    void ionModelSolve( const vector_type& u, const Real timestep );

    //! Computes the term  (G1 (u-u0)(u-u0-aA)(u-u0-A)+G2(u-u0)w) for the PDE righthand side
    void computeIion(ElemVec& elvec, ElemVec& elvec_u, FESpace<Mesh, EpetraMap>& uFESpace );

    //! Returns the local solution vector
    const vector_type& solution_w() const {return M_sol_w;}

    //! Initialize
    void initialize( );
//    void initialize( const vector_type& );

protected:

    //! Global solution _w
    vector_type                    	M_sol_w;
    vector_type						M_wVecRep;
    ElemVec 						M_elvec;


private:
};


//
// IMPLEMENTATION
//

//! Constructor
template<typename Mesh, typename SolverType>
Rogers_McCulloch<Mesh, SolverType>::
Rogers_McCulloch( const data_type&          dataType,
       FESpace<Mesh, EpetraMap>& uFEspace,
       Epetra_Comm&              comm ):
    	   IonicSolver<Mesh, SolverType>( dataType, uFEspace, comm),
    	   M_sol_w                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
    	   M_wVecRep( M_sol_w, Repeated ),
    	   M_elvec ( IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbNode, 1 )
{
}




template<typename Mesh, typename SolverType>
Rogers_McCulloch<Mesh, SolverType>::
~Rogers_McCulloch()
{
}


template<typename Mesh, typename SolverType>
void Rogers_McCulloch<Mesh, SolverType>::updateRepeated( )
{

	M_wVecRep=M_sol_w;

}

template<typename Mesh, typename SolverType>
void Rogers_McCulloch<Mesh, SolverType>::updateElvec( UInt eleID )
{
	M_elvec.zero();
	UInt ig;
		//! Filling local elvec_w with recovery variable values in the nodes
		for ( UInt iNode = 0 ; iNode < ( UInt ) IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbNode ; iNode++ )
		{
			ig = IonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobal( eleID, iNode + 1 );
			M_elvec.vec()[ iNode ] = M_wVecRep[ig];
		}

}

template<typename Mesh, typename SolverType>
void Rogers_McCulloch<Mesh, SolverType>::ionModelSolve( const vector_type& u, const Real timestep )
{
	//! Solving dw/dt= b/(A*T) (u - u0 - A d w)
	Real G = this->M_data.b/this->M_data.A/this->M_data.T;

	Real alpha = 1/timestep + G*this->M_data.A*this->M_data.d;

	IonicSolver<Mesh, SolverType>::M_comm->Barrier(); //    MPI_Barrier(MPI_COMM_WORLD);

	M_sol_w*=1.0/timestep;
	vector_type temp(u);
	temp.getEpetraVector().PutScalar (G*this->M_data.u0);
	M_sol_w+= G*u;
	M_sol_w-= temp;
	M_sol_w*=1/alpha;
	M_sol_w.GlobalAssemble();

	IonicSolver<Mesh, SolverType>::M_comm->Barrier(); //    MPI_Barrier(MPI_COMM_WORLD);

}

template<typename Mesh, typename SolverType>
void Rogers_McCulloch<Mesh, SolverType>::computeIion( ElemVec& elvec, ElemVec& elvec_u,
                                                      FESpace<Mesh, EpetraMap>& uFESpace )
{
    //Computing  G1 (u-u0)(u-u0-aA)(u-u0-A)+G2(u-u0)w

 	Real u_ig, w_ig, s;
    Real G1 = this->M_data.c1/this->M_data.T/pow(this->M_data.A,2.0);
	Real G2 = this->M_data.c2/this->M_data.T;

    for ( int i = 0;i < uFESpace.fe().nbNode;i++ )
    {
        s = 0.0;
        for ( int ig = 0; ig < uFESpace.fe().nbQuadPt();ig++ )
        {
            u_ig = 0.0; w_ig = 0.0;
            for ( int i = 0;i < uFESpace.fe().nbNode;i++ ){
                u_ig += elvec_u( i ) * uFESpace.fe().phi( i, ig );
                w_ig += M_elvec( i ) * uFESpace.fe().phi( i, ig );
            }
            s += (G1*(u_ig - this->M_data.u0) *
                  (u_ig - this->M_data.u0 - this->M_data.a*this->M_data.A)*
                  (u_ig - this->M_data.u0 - this->M_data.A)
                  + G2 * (u_ig - this->M_data.u0) * w_ig) * uFESpace.fe().phi( i, ig ) *
                uFESpace.fe().weightDet( ig );
        }
        elvec( i ) +=s;
    }

// 	Real u_ig, w_ig;

// 	Real G1 = this->M_data.c1/this->M_data.T/pow(this->M_data.A,2.0);
// 	Real G2 = this->M_data.c2/this->M_data.T;

//     for ( int ig = 0; ig < uFESpace.fe().nbQuadPt;ig++ )
//     {
//         u_ig = w_ig = 0.;
//         for ( int i = 0;i < uFESpace.fe().nbNode;i++ )
//             u_ig += elvec_u( i ) * uFESpace.fe().phi( i, ig );
//         for ( int i = 0;i < uFESpace.fe().nbNode;i++ )
//             w_ig += M_elvec( i ) * uFESpace.fe().phi( i, ig );

//         for ( int i = 0;i < uFESpace.fe().nbNode;i++ )
//         {
//             elvec( i ) += (G1*(u_ig- this->M_data.u0)*(u_ig - this->M_data.u0 - this->M_data.a*this->M_data.A)*(u_ig - this->M_data.u0 - this->M_data.A) + G2 * (u_ig - this->M_data.u0) * w_ig) * uFESpace.fe().phi( i, ig ) * uFESpace.fe().weightDet( ig );
//         }
//     }


}

template<typename Mesh, typename SolverType>
void Rogers_McCulloch<Mesh, SolverType>::initialize( )
{
	M_sol_w.getEpetraVector().PutScalar (0.0);
}


} // namespace LifeV


#endif //_IONICSOLVER_H_
