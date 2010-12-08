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
 License along with LifeV. If not, see <http://www.gnu.org/licences/>.

************************************************************************
*/
//@HEADER
/*!
 * @file
 * @brief Class for choosing and solving ionic models in electrophysiology. Mitchel-Schaeffer, Rogers-McCulloch, Luo-Rudy I
 * @author Lucia Mirabella <lucia.mirabella@mail.polimi.it> and Mauro Perego <mauro.perego@polimi.it>
 * @date 11-2007
 * @contributors J.Castelneau (INRIA), Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 * @mantainer Lucia Mirabella <lucia.mirabella@mail.polimi.it>
 * @last update 11-2010
 */

#ifndef _IONICSOLVER_H_
#define _IONICSOLVER_H_

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
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

namespace LifeV
{
//! IonicSolver - This class implements a ionic model solver.


template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class IonicSolver
{

public:
    //! @name Type definitions
    //@{

    typedef DataIonic data_Type;

    typedef Real ( *Function ) ( const Real&,
                                 const Real&,
                                 const Real&,
                                 const Real&,
                                 const ID& );

    typedef boost::function<Real (const EntityFlag& ref,
                                  const Real& x,
                                  const Real& y,
                                  const Real& z,
                                  const ID& i)> fct_TauClose;

    typedef Mesh mesh_type;

    typedef typename SolverType::matrix_type      matrix_Type;
    typedef typename SolverType::vector_type      vector_Type;

    typedef typename SolverType::prec_raw_type    precRaw_Type;
    typedef typename SolverType::prec_type        prec_Type;

 //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     * @param dataType
     * @param mesh
     * @param recovery FE space
     * @param Epetra communicator
     */
    IonicSolver( const data_Type& dataType,
                 const Mesh& mesh,
                 FESpace<Mesh, EpetraMap>& uFEspace,
                 Epetra_Comm& comm );

    //! Destructor
    virtual ~IonicSolver();

 //@}

    //! @name Methods
    //@{


    //! Returns the recovery variable FE space
    FESpace<Mesh, EpetraMap>& recoveryFESpace()   {return M_uFESpace;}

    //! Return maps
    Epetra_Map const& getRepeatedEpetraMap() const { return *M_localMap.getMap(Repeated); }

    EpetraMap const& getMap() const { return M_localMap; }

    virtual void updateRepeated( )=0;

    //! Update the ionic model elvecs
    virtual void updateElvec( UInt eleIDw )=0;

    //! Solves the ionic model
    virtual void ionModelSolve( const vector_Type& u,
                                const Real timeStep )=0;

    //! Computes the term -1/ \Cm u^n (G (1-u^n/vp) (1-u^n/v_th) + eta_1 v^{n+1})
    //! for the PDE righthand side
    virtual void computeIion( Real Capacitance,
                              ElemVec& elvec,
                              ElemVec& elvec_u,
                              FESpace<Mesh, EpetraMap>& uFESpace )=0;

    //! Initialize
    virtual void initialize( )=0;

    const data_Type&               M_data;

    const Mesh&               M_mesh;

    // FE space
    FESpace<Mesh, EpetraMap>&      M_uFESpace;

    //! MPI communicator
    Epetra_Comm*                   M_comm;

    Int                            M_me;

    //! Map
    EpetraMap                      M_localMap;

    //! Boolean that indicates if output is sent to cout
    bool                           M_verbose;

    //@}

protected:

    UInt dim_u() const           { return M_uFESpace.dim(); }

}; // class IonicSolver


//
// IMPLEMENTATION
//

//! Constructor
template<typename Mesh, typename SolverType>
IonicSolver<Mesh, SolverType>::
IonicSolver( const data_Type& dataType,
             const Mesh& mesh,
             FESpace<Mesh, EpetraMap>& uFEspace,
             Epetra_Comm& comm ):
    M_data                   ( dataType ),
    M_mesh                   ( mesh ),
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

template< typename Mesh,
	  typename SolverType = LifeV::SolverTrilinos >
class MitchellSchaeffer : public virtual IonicSolver<Mesh, SolverType>
{
public:
	typedef typename IonicSolver<Mesh, SolverType>::data_Type	data_Type;
	typedef typename IonicSolver<Mesh, SolverType>::vector_type	vector_Type;
	typedef typename IonicSolver<Mesh, SolverType>::Function 	Function;
	typedef typename IonicSolver<Mesh, SolverType>::fct_TauClose	fct_TauClose;

    MitchellSchaeffer( const data_Type& dataType,
                        const Mesh& mesh,
                        FESpace<Mesh, EpetraMap>& uFEspace,
                        Epetra_Comm& comm );

    virtual ~MitchellSchaeffer();

    void updateRepeated( );

    void updateElvec( UInt eleID );

    void ionModelSolve( const vector_Type& u, const Real timeStep );

    void computeIion( Real Capacitance,
                      ElemVec& elvec,
                      ElemVec& elvec_u,
                      FESpace<Mesh, EpetraMap>& uFESpace );

    const vector_Type& solution_w() const {return M_sol_w;}

    void initialize( );

    void setHeteroTauClose(fct_TauClose);

    Real fct_Tau_Close(const EntityFlag& ref,
                       const Real& x,
                       const Real& y,
                       const Real& z,
                       const ID& i) const;

protected:

    //! Global solution _w
    vector_Type                    	M_sol_w;
    vector_Type				M_wVecRep;
    ElemVec 				M_elvec;
    UInt				order_bdf;
    BdfT<vector_Type> 			bdf_w;
    fct_TauClose			M_TauClose;

private:
};


//
// IMPLEMENTATION
//

//! Constructor
template<typename Mesh, typename SolverType>
MitchellSchaeffer<Mesh, SolverType>::
MitchellSchaeffer( const data_Type& dataType,
                    const Mesh& mesh,
                    FESpace<Mesh, EpetraMap>& uFEspace,
                    Epetra_Comm& comm ):
    IonicSolver<Mesh, SolverType>( dataType, mesh, uFEspace, comm),
    M_sol_w ( IonicSolver<Mesh, SolverType>::M_localMap ),
    M_wVecRep( M_sol_w, Repeated ),
    M_elvec ( IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbNode, 1 ),
    order_bdf ( IonicSolver<Mesh, SolverType>::M_data.M_BDForder ),
    bdf_w( order_bdf )
{
	setHeteroTauClose(this->M_data.M_ShdPtr->get_heterotauclose());
}

template<typename Mesh, typename SolverType>
MitchellSchaeffer<Mesh, SolverType>::
~MitchellSchaeffer()
{
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::updateRepeated( )
{
	M_wVecRep=M_sol_w;
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::updateElvec( UInt eleID )
{
	M_elvec.zero();
	UInt ig;
    //! Filling local elvec_w with recovery variable values in the nodes
    for ( UInt iNode = 0 ; iNode < IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbNode ; iNode++ )
    {
        ig = IonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobal( eleID, iNode + 1 );
        M_elvec.vec()[ iNode ] = M_wVecRep[ig];
    }
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::setHeteroTauClose(fct_TauClose fct)
{
    M_TauClose = fct;
}

template<typename Mesh, typename SolverType>
Real MitchellSchaeffer<Mesh, SolverType>::fct_Tau_Close(const EntityFlag& ref, const Real& x, const Real& y, const Real& z, const ID& i) const
{
    return M_TauClose(ref, x, y, z, i);
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::ionModelSolve( const vector_Type& u, const Real timeStep )
{
    //! Solving :
	//!           ((v_max-v_min)^{-2}  - w )/tau_open  if u < vcrit
	//! dw/dt ={
	//!            -w/tau_close   if u > vcrit

	Real aux1 = 1.0 / (bdf_w.coeff_der(0)/timeStep  + 1.0/this->M_data.M_tau_open );
	Real aux = 1.0/((this->M_data.M_potentialMaximum - this->M_data.M_potentialMinimum)*(this->M_data.M_potentialMaximum- this->M_data.M_potentialMinimum)* this->M_data.M_tau_open);
	Real aux2 = 1.0 / (bdf_w.coeff_der(0)/timeStep  + 1.0/this->M_data.M_tau_close);
	vector_Type M_time_der=bdf_w.time_der(timeStep);

	IonicSolver<Mesh, SolverType>::M_comm->Barrier();
	Real x, y, z;
	EntityFlag ref;
	UInt ID;
	for ( Int i = 0 ; i < u.getEpetraVector().MyLength() ; ++i )
	{
        Int ig=u.BlockMap().MyGlobalElements()[i];
		ID 	= ig;
		ref = this->M_mesh->point(ig).marker();
		x 	= this->M_mesh->point(ig).x();
        y 	= this->M_mesh->point(ig).y();
        z 	= this->M_mesh->point(ig).z();
   		if (u[ig] < this->M_data.M_criticalPotential)
            M_sol_w[ig] = aux1 * (aux + M_time_der[ig]);
        else if (this->M_data.M_hasHeterogeneousTauClose)
            M_sol_w[ig] = (1.0 / (bdf_w.coeff_der(0)/timeStep  + 1.0/fct_Tau_Close(ref,x,y,z,ID))) *  M_time_der[ig];//aux2 * M_time_der[ig];
        else
            M_sol_w[ig] = aux2 *  M_time_der[ig];
	}
	bdf_w.shift_right(M_sol_w);

	IonicSolver<Mesh, SolverType>::M_comm->Barrier();
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::computeIion(  Real,
                                                         ElemVec& elvec,
                                                         ElemVec& elvec_u,
                                                         FESpace<Mesh, EpetraMap>& uFESpace )
{
	for ( Int i = 0;i < uFESpace.fe().nbNode;i++ )
    {
        elvec( i ) =  this->M_reactionAmplitude*(((M_elvec( i ) / this->M_data.M_tau_in)
                                              * (elvec_u( i ) - this->M_data.M_potentialMinimum)
                                              * (elvec_u( i ) - this->M_data.M_potentialMinimum)
                                              * (this->M_data.M_potentialMaximum - elvec_u( i ))
                                              / (this->M_data.M_potentialMaximum - this->M_data.M_potentialMinimum) )
                                             - (( elvec_u( i ) - this->M_data.M_potentialMinimum )
                                                /(  this->M_data.M_tau_out * (this->M_data.M_potentialMaximum
                                                                              - this->M_data.M_potentialMinimum)))) ;
    }
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::
initialize( )
{
	M_sol_w.getEpetraVector().PutScalar (1.0/((this->M_data.M_potentialMaximum-this->M_data.M_potentialMinimum)*(this->M_data.M_potentialMaximum-this->M_data.M_potentialMinimum)));
	bdf_w.initialize_unk(M_sol_w);
	bdf_w.showMe();
}

template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class RogersMcCulloch : public virtual IonicSolver<Mesh, SolverType>
{
public:
	typedef typename IonicSolver<Mesh, SolverType>::data_Type data_Type;
	typedef typename IonicSolver<Mesh, SolverType>::vector_Type vector_Type;
	typedef typename IonicSolver<Mesh, SolverType>::Function Function;

    RogersMcCulloch( const data_Type& dataType,
                      const Mesh& mesh,
                      FESpace<Mesh, EpetraMap>& uFEspace,
                      Epetra_Comm& comm );
    virtual ~RogersMcCulloch();

    void updateRepeated( );

    void updateElvec( UInt eleID );

    void ionModelSolve( const vector_Type& u, const Real timeStep );

    void computeIion( Real Capacitance,
                      ElemVec& elvec,
                      ElemVec& elvec_u,
                      FESpace<Mesh, EpetraMap>& uFESpace );

    const vector_Type& solution_w() const {return M_sol_w;}

    void initialize( );

protected:

    //! Global solution _w
    vector_Type                    	M_sol_w;
    vector_Type						M_wVecRep;
    ElemVec 						M_elvec;
private:
};


//
// IMPLEMENTATION
//

//! Constructor
template<typename Mesh, typename SolverType>
RogersMcCulloch<Mesh, SolverType>::
RogersMcCulloch( const data_Type& dataType,
                  const Mesh& mesh,
                  FESpace<Mesh, EpetraMap>& uFEspace,
                  Epetra_Comm& comm ):
    IonicSolver<Mesh, SolverType>( dataType, mesh, uFEspace, comm),
    M_sol_w ( IonicSolver<Mesh, SolverType>::M_localMap ),
    M_wVecRep ( M_sol_w, Repeated ),
    M_elvec ( IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbNode, 1 )
{
}

template<typename Mesh, typename SolverType>
RogersMcCulloch<Mesh, SolverType>::
~RogersMcCulloch()
{
}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::updateRepeated( )
{
	M_wVecRep=M_sol_w;
}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::updateElvec( UInt eleID )
{
	M_elvec.zero();
	UInt ig;
		//! Filling local elvec_w with recovery variable values in the nodes
		for ( UInt iNode = 0 ; iNode < IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbNode ; iNode++ )
		{
			ig = IonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobal( eleID, iNode + 1 );
			M_elvec.vec()[ iNode ] = M_wVecRep[ig];
		}

}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::ionModelSolve( const vector_Type& u, const Real timeStep )
{
	//! Solving dw/dt= b/(A*T) (u - u0 - A d w)
	Real G = this->M_data.M_b/this->M_data.M_potentialAmplitude/this->M_data.M_timeUnit;

	Real alpha = 1/timeStep + G*this->M_data.M_potentialAmplitude*this->M_data.M_d;

	IonicSolver<Mesh, SolverType>::M_comm->Barrier();

	M_sol_w*=1/timeStep;
	vector_Type temp(u);
	temp.getEpetraVector().PutScalar (G*this->M_data.M_restPotential);
	M_sol_w+= G*u;
	M_sol_w-= temp;
	M_sol_w*=1/alpha;
	M_sol_w.GlobalAssemble();
	IonicSolver<Mesh, SolverType>::M_comm->Barrier();

}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::computeIion(  Real Capacitance,
                                                       ElemVec& elvec,
                                                       ElemVec& elvec_u,
                                                       FESpace<Mesh, EpetraMap>& uFESpace )
{
	Real u_ig, w_ig;

	Real G1 = this->M_data.M_c1/this->M_data.M_timeUnit/pow(this->M_data.M_potentialAmplitude,2.0);
	Real G2 = this->M_data.M_c2/this->M_data.M_timeUnit;

        for ( UInt ig = 0; ig < uFESpace.fe().nbQuadPt();ig++ )
    {
        u_ig = w_ig = 0.;
        for ( UInt i = 0;i < uFESpace.fe().nbNode;i++ )
            u_ig += elvec_u( i ) * uFESpace.fe().phi( i, ig );
        for ( UInt i = 0;i < uFESpace.fe().nbNode;i++ )
            w_ig += M_elvec( i ) * uFESpace.fe().phi( i, ig );

        for ( UInt i = 0;i < uFESpace.fe().nbNode;i++ )
        {
            elvec( i ) -= Capacitance*(G1*(u_ig- this->M_data.M_restPotential)*
(u_ig - this->M_data.M_restPotential - this->M_data.M_a*this->M_data.M_potentialAmplitude)*
(u_ig - this->M_data.M_restPotential - this->M_data.M_potentialAmplitude) + G2 * (u_ig - this->M_data.M_restPotential) * w_ig) *
 uFESpace.fe().phi( i, ig ) * uFESpace.fe().weightDet( ig );
        }
    }
}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::
initialize( )
{
	M_sol_w.getEpetraVector().PutScalar (0.);
}

template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class LuoRudy : public virtual IonicSolver<Mesh, SolverType>
{
public:
	typedef typename IonicSolver<Mesh, SolverType>::data_Type data_Type;
	typedef typename IonicSolver<Mesh, SolverType>::vector_Type vector_Type;
	typedef typename IonicSolver<Mesh, SolverType>::Function Function;

    LuoRudy( const data_Type&          dataType,
             const Mesh&          mesh,
             FESpace<Mesh, EpetraMap>& uFEspace,
             Epetra_Comm&              comm );
    virtual ~LuoRudy();

    void updateRepeated( );

    void updateElvec( UInt eleID);

    void compute_coeff( const Real& u_ig );

    void ionModelSolve( const vector_Type& u, const Real timeStep );

    void computeIion( Real Capacitance,
                      ElemVec& elvec,
                      ElemVec& elvec_u,
                      FESpace<Mesh, EpetraMap>& uFESpace );

    void initialize ( );

    //! Returns the local solution vector for each field
    const vector_Type& solution_h() const {return M_sol_h;}

    const vector_Type& solution_j() const {return M_sol_j;}

    const vector_Type& solution_m() const {return M_sol_m;}

    const vector_Type& solution_d() const {return M_sol_d;}

    const vector_Type& solution_f() const {return M_sol_f;}

    const vector_Type& solution_X() const {return M_sol_X;}

    const vector_Type& solution_Ca() const {return M_sol_Ca;}

    Real k_0, k_i, na_0, na_i, R, T, F, PR, c, e_na, g_k, e_k, g_k1, e_k1, e_kp, e_si;
    Real a_h, b_h, a_j, b_j, x_ii, a_m, b_m, a_d, b_d, a_f, b_f, a_X, b_X, ak1, bk1;
    Real k_p, k_1inf, h_inf, tau_h, j_inf, tau_j, m_inf, tau_m, d_inf, tau_d, f_inf,
    		tau_f, X_inf, tau_X;
	//fast sodium current
	Real M_Ina;
	//slow inward current
	Real M_Islow;
    //time dependent potassium current
    Real M_Ik;
    //time independent potassium current
    Real M_Ik1;
    //plateau potassium current
    Real M_Ikp;
    //background current
    Real M_Iback;
    //Total time independent potassium current
    Real M_Ik1t;

	vector_Type exp_vec_h;
	vector_Type exp_vec_j;
	vector_Type exp_vec_m;
	vector_Type exp_vec_d;
	vector_Type exp_vec_f;
	vector_Type exp_vec_X;
	vector_Type inf_vec_h;
	vector_Type inf_vec_j;
	vector_Type inf_vec_m;
	vector_Type inf_vec_d;
	vector_Type inf_vec_f;
	vector_Type inf_vec_X;
	vector_Type	dca_i_vec;

protected:
    //! Global solution h
    vector_Type                    	M_sol_h;
    //! Global solution j
    vector_Type                    	M_sol_j;
    //! Global solution m
    vector_Type                    	M_sol_m;
    //! Global solution d
    vector_Type                    	M_sol_d;
    //! Global solution f
    vector_Type                    	M_sol_f;
    //! Global solution X
    vector_Type                    	M_sol_X;
    //! Global solution Ca_i
    vector_Type                    	M_sol_Ca;

    vector_Type 					Iion;

    vector_Type						M_Iion_VecRep;

    ElemVec 						M_elvec_Iion;

private:
};


//
// IMPLEMENTATION
//


//! Constructor
template<typename Mesh, typename SolverType>
LuoRudy<Mesh, SolverType>::
LuoRudy( const data_Type& dataType,
		const Mesh& mesh,
		FESpace<Mesh, EpetraMap>& uFEspace,
		Epetra_Comm& comm ):
			IonicSolver<Mesh, SolverType>( dataType, mesh, uFEspace, comm),
			M_sol_h                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_sol_j                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_sol_m                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_sol_d                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_sol_f                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_sol_X                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_sol_Ca                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			Iion( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_Iion_VecRep( Iion, Repeated ),
			M_elvec_Iion ( IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbNode, 1 ),
			k_0(5.4),
			k_i(145.),
			na_0(140.),
			na_i(18.),
			R(8.314472), //% Joules/(Kelvin*mole)
			T(307.7532), //%kelvins
			F(96485.33838), //% coulumbs/mole
			PR(0.01833), //% Na/K permeability ratio
			c(1.), //% membrane capacitance set as 1 mui-F/cm^2
			e_na(1000.*(R*T/F)*log(na_0/na_i)),
			g_k(0.282*sqrt(k_0/5.4)),
			e_k(1000.*(R*T/F)*log((k_0+PR*na_0) / (k_i+PR*na_i))),
			g_k1(0.6047*sqrt(k_0/5.4)),
			e_k1(1000.*(R*T/F)*log(k_0/k_i)),
			e_kp(e_k1),
			exp_vec_h(IonicSolver<Mesh, SolverType>::M_localMap),
			exp_vec_j(IonicSolver<Mesh, SolverType>::M_localMap),
			exp_vec_m(IonicSolver<Mesh, SolverType>::M_localMap),
			exp_vec_d(IonicSolver<Mesh, SolverType>::M_localMap),
			exp_vec_f(IonicSolver<Mesh, SolverType>::M_localMap),
			exp_vec_X(IonicSolver<Mesh, SolverType>::M_localMap),
			inf_vec_h(IonicSolver<Mesh, SolverType>::M_localMap),
			inf_vec_j(IonicSolver<Mesh, SolverType>::M_localMap),
			inf_vec_m(IonicSolver<Mesh, SolverType>::M_localMap),
			inf_vec_d(IonicSolver<Mesh, SolverType>::M_localMap),
			inf_vec_f(IonicSolver<Mesh, SolverType>::M_localMap),
			inf_vec_X(IonicSolver<Mesh, SolverType>::M_localMap),
			dca_i_vec(IonicSolver<Mesh, SolverType>::M_localMap)
{
}

//! Destructor
template<typename Mesh, typename SolverType>
LuoRudy<Mesh, SolverType>::
~LuoRudy()
{
}

template<typename Mesh, typename SolverType>
void LuoRudy<Mesh, SolverType>::updateRepeated( )
{

	M_Iion_VecRep=Iion;
}

template<typename Mesh, typename SolverType>
void LuoRudy<Mesh, SolverType>::updateElvec( UInt eleID)
{
	M_elvec_Iion.zero();
	UInt ig;
		//! Filling local elvec with recovery variable values in the nodes
		for ( UInt iNode = 0 ; iNode < IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbNode ; iNode++ )
		{
			ig = IonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobal( eleID, iNode + 1 );
			M_elvec_Iion.vec()[ iNode ] = M_Iion_VecRep[ig];
		}

}


template<typename Mesh, typename SolverType>
void LuoRudy<Mesh, SolverType>::ionModelSolve( const vector_Type& u, const Real timeStep )
{
	//! Solving dw/dt=eta2 (u/vp -  eta3 w)
	Chrono chronoionmodelsolve;
	chronoionmodelsolve.start();
	IonicSolver<Mesh, SolverType>::M_comm->Barrier();

	for ( Int i = 0 ; i < u.getEpetraVector().MyLength() ; i++ )
	{
		Int ig=u.BlockMap().MyGlobalElements()[i];
		Real u_ig=u[ig];
		compute_coeff(u_ig);
        e_si = 7.7-13.0287*log(M_sol_Ca[ig]);
		//fast sodium current
		M_Ina = 23.*M_sol_m[ig]*M_sol_m[ig]*M_sol_m[ig]*M_sol_h[ig]*M_sol_j[ig]*(u_ig-e_na);
		//slow inward current
		M_Islow = 0.09*M_sol_d[ig]*M_sol_f[ig]*(u_ig-e_si);
        //change in ioniq concentration
        dca_i_vec.getEpetraVector().ReplaceGlobalValue(ig,0,-1e-4*M_Islow+0.07*(1e-4-M_sol_Ca[ig]));
        //time dependent potassium current
        M_Ik = g_k*M_sol_X[ig]*x_ii*(u_ig-e_k);
        //time independent potassium current
        M_Ik1 = g_k1*k_1inf*(u_ig-e_k1);
        //plateau potassium current
        M_Ikp = 0.0183*k_p*(u_ig-e_kp);
        //background current
        M_Iback = 0.03921*(u_ig+59.87);
        //Total time independent potassium current
        M_Ik1t = M_Ik1 + M_Ikp + M_Iback;
        // adding up the six ionic currents
        Iion.getEpetraVector().ReplaceGlobalValue(ig,0,M_Ina + M_Islow + M_Ik + M_Ik1t);
		exp_vec_h.getEpetraVector().ReplaceGlobalValue(ig,0,exp(-timeStep/tau_h));
		exp_vec_j.getEpetraVector().ReplaceGlobalValue(ig,0,exp(-timeStep/tau_j));
		exp_vec_m.getEpetraVector().ReplaceGlobalValue(ig,0,exp(-timeStep/tau_m));
		exp_vec_d.getEpetraVector().ReplaceGlobalValue(ig,0,exp(-timeStep/tau_d));
		exp_vec_f.getEpetraVector().ReplaceGlobalValue(ig,0,exp(-timeStep/tau_f));
		exp_vec_X.getEpetraVector().ReplaceGlobalValue(ig,0,exp(-timeStep/tau_X));
		inf_vec_h.getEpetraVector().ReplaceGlobalValue(ig,0,h_inf);
		inf_vec_j.getEpetraVector().ReplaceGlobalValue(ig,0,j_inf);
		inf_vec_m.getEpetraVector().ReplaceGlobalValue(ig,0,m_inf);
		inf_vec_d.getEpetraVector().ReplaceGlobalValue(ig,0,d_inf);
		inf_vec_f.getEpetraVector().ReplaceGlobalValue(ig,0,f_inf);
		inf_vec_X.getEpetraVector().ReplaceGlobalValue(ig,0,X_inf);
	}
	exp_vec_h.GlobalAssemble();
	exp_vec_j.GlobalAssemble();
	exp_vec_m.GlobalAssemble();
	exp_vec_d.GlobalAssemble();
	exp_vec_f.GlobalAssemble();
	exp_vec_X.GlobalAssemble();
	inf_vec_h.GlobalAssemble();
	inf_vec_j.GlobalAssemble();
	inf_vec_m.GlobalAssemble();
	inf_vec_d.GlobalAssemble();
	inf_vec_f.GlobalAssemble();
	inf_vec_X.GlobalAssemble();
	Iion.GlobalAssemble();
	dca_i_vec.GlobalAssemble();

	M_sol_h-=inf_vec_h;
	M_sol_h.getEpetraVector().Multiply(1., M_sol_h.getEpetraVector(), exp_vec_h.getEpetraVector(), 0.);
	M_sol_h+=inf_vec_h;

	M_sol_j-=inf_vec_j;
	M_sol_j.getEpetraVector().Multiply(1, M_sol_j.getEpetraVector(), exp_vec_j.getEpetraVector(), 0);
	M_sol_j+=inf_vec_j;

	M_sol_m-=inf_vec_m;
	M_sol_m.getEpetraVector().Multiply(1, M_sol_m.getEpetraVector(), exp_vec_m.getEpetraVector(), 0);
	M_sol_m+=inf_vec_m;

	M_sol_d-=inf_vec_d;
	M_sol_d.getEpetraVector().Multiply(1, M_sol_d.getEpetraVector(), exp_vec_d.getEpetraVector(), 0);
	M_sol_d+=inf_vec_d;

	M_sol_f-=inf_vec_f;
	M_sol_f.getEpetraVector().Multiply(1, M_sol_f.getEpetraVector(), exp_vec_f.getEpetraVector(), 0);
	M_sol_f+=inf_vec_f;

	M_sol_X-=inf_vec_X;
	M_sol_X.getEpetraVector().Multiply(1, M_sol_X.getEpetraVector(), exp_vec_X.getEpetraVector(), 0);
	M_sol_X+=inf_vec_X;

	M_sol_Ca+=timeStep*dca_i_vec;

	M_sol_h.GlobalAssemble();
	M_sol_j.GlobalAssemble();
	M_sol_m.GlobalAssemble();
	M_sol_d.GlobalAssemble();
	M_sol_f.GlobalAssemble();
	M_sol_X.GlobalAssemble();
	M_sol_Ca.GlobalAssemble();

	IonicSolver<Mesh, SolverType>::M_comm->Barrier();

	chronoionmodelsolve.stop();
    if (IonicSolver<Mesh, SolverType>::M_comm->MyPID()==0) std::cout << "Total ionmodelsolve time " << chronoionmodelsolve.diff() << " s." << std::endl;
}



template<typename Mesh, typename SolverType>
void LuoRudy<Mesh, SolverType>::compute_coeff( const Real& u_ig )
{
	   if (u_ig>=-40.)
	       {
	       a_h=0.;
	       b_h = 1./(0.13*(1.+exp( (u_ig+10.66)/(-11.1) )));
	       a_j=0.;
	       b_j = 0.3*exp(-2.535e-7*u_ig)/(1.+exp(-0.1*(u_ig+32.)));
	       }
	   else
	       {
	       a_h = 0.135*exp((80.+u_ig)/-6.8);
	       b_h = 3.56*exp(0.079*u_ig)+3.1e5*exp(0.35*u_ig);
	       a_j = (-1.2714e5*exp(0.2444*u_ig)-3.474e-5*exp(-0.04391*u_ig))*(u_ig+37.78)/(1+exp(0.311*(u_ig+79.23)));
	       b_j = 0.1212*exp(-0.01052*u_ig)/(1.+exp(-0.1378*(u_ig+40.14)));
	       }
	      a_m = 0.32*(u_ig+47.13)/(1.-exp(-0.1*(u_ig+47.13)));
	      b_m = 0.08*exp(-u_ig/11.);

	      //slow inward current
	      a_d = 0.095*exp(-0.01*(u_ig-5.))/(1.+exp(-0.072*(u_ig-5.)));
	      b_d = 0.07*exp(-0.017*(u_ig+44.))/(1.+exp(0.05*(u_ig+44.)));
	      a_f = 0.012*exp(-0.008*(u_ig+28.))/(1.+exp(0.15*(u_ig+28.)));
	      b_f = 0.0065*exp(-0.02*(u_ig+30.))/(1.+exp(-0.2*(u_ig+30.)));

	      //Time dependent potassium outward current
	      a_X = 0.0005*exp(0.083*(u_ig+50.))/(1.+exp(0.057*(u_ig+50.)));
	      b_X = 0.0013*exp(-0.06*(u_ig+20.))/(1.+exp(-0.04*(u_ig+20.)));

	      if(u_ig<=-100)
	    	  x_ii=1.;
	      else
	    	  x_ii=2.837*(exp(0.04*(u_ig+77.))-1.)/((u_ig+77.)*exp(0.04*(u_ig+35.)));


	      ak1 = 1.02/(1.+exp(0.2385*(u_ig-e_k1-59.215)));
	      bk1 = (0.49124*exp(0.08032*(u_ig-e_k1+5.476))+exp(0.06175*(u_ig-e_k1-594.31)))/(1.+exp(-0.5143*(u_ig-e_k1+4.753))); //%bk1

          //Plateau potassium outward current
	      k_p = 1./(1.+exp((7.488-u_ig)/5.98));

	      k_1inf = ak1/(ak1+bk1);
	      h_inf = a_h/(a_h + b_h);
	      tau_h = 1./(a_h + b_h);
	      j_inf = a_j/(a_j + b_j);
	      tau_j = 1./(a_j + b_j);
	      m_inf = a_m/(a_m + b_m);
	      tau_m = 1./(a_m + b_m);
	      d_inf = a_d/(a_d + b_d);
	      tau_d = 1./(a_d + b_d);
	      f_inf = a_f/(a_f + b_f);
	      tau_f = 1./(a_f + b_f);
	      X_inf = a_X/(a_X + b_X);
	      tau_X = 1./(a_X + b_X);
}

template<typename Mesh, typename SolverType>
void LuoRudy<Mesh, SolverType>::computeIion(  Real Capacitance,
                                               ElemVec& elvec,
                                               ElemVec& /*elvec_u*/,
                                               FESpace<Mesh, EpetraMap>& uFESpace )
{
   	Real Iion_ig;
    for ( UInt ig = 0; ig < uFESpace.fe().nbQuadPt();ig++ )
    {
        Iion_ig = 0.;
        for ( UInt i = 0;i < uFESpace.fe().nbNode;i++ )
            Iion_ig += M_elvec_Iion( i ) * uFESpace.fe().phi( i, ig );
        for ( UInt i = 0;i < uFESpace.fe().nbNode;i++ )
        {
            // divide by 1000 to convert microA in mA
            elvec( i ) -= Iion_ig/1000 * uFESpace.fe().phi( i, ig ) * uFESpace.fe().weightDet( ig );
        }
    }
}

template<typename Mesh, typename SolverType>
void LuoRudy<Mesh, SolverType>::
initialize( )
{
	 M_sol_h.getEpetraVector().PutScalar(1.);
     M_sol_j.getEpetraVector().PutScalar(1.);
     M_sol_m.getEpetraVector().PutScalar(0.);
     M_sol_d.getEpetraVector().PutScalar(0.);
     M_sol_f.getEpetraVector().PutScalar(1.);
     M_sol_X.getEpetraVector().PutScalar(0.);
     M_sol_Ca.getEpetraVector().PutScalar(0.0002);
	 M_sol_h.GlobalAssemble();
     M_sol_j.GlobalAssemble();
     M_sol_m.GlobalAssemble();
     M_sol_d.GlobalAssemble();
     M_sol_f.GlobalAssemble();
     M_sol_X.GlobalAssemble();
     M_sol_Ca.GlobalAssemble();
}

} // namespace LifeV


#endif //_IONICSOLVER_H_
