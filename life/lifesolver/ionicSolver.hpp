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
  @brief Class for choosing and solving ionic models in electrophysiology. Mitchel-Schaeffer, Rogers-McCulloch, Luo-Rudy I

  @date 11-2007
  @author Lucia Mirabella <lucia.mirabella@mail.polimi.it> and Mauro Perego <mauro.perego@polimi.it>

  @contributors J.Castelneau (INRIA), Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
  @mantainer Lucia Mirabella <lucia.mirabella@mail.polimi.it>
  @last update 11-2010
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

    typedef Real ( *function_Type ) ( const Real&,
                                 const Real&,
                                 const Real&,
                                 const Real&,
                                 const ID& );

    typedef boost::function<Real (const entityFlag_Type& ref,
                                  const Real& x,
                                  const Real& y,
                                  const Real& z,
                                  const ID& i)> functorTauClose_Type;

    typedef Mesh mesh_Type;

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
    Epetra_Map const& getRepeatedEpetraMap() const { return *M_localMap.map(Repeated); }

    EpetraMap const& getMap() const { return M_localMap; }

    virtual void updateRepeated( )=0;

    //! Update the ionic model elvecs
    virtual void updateElementSolution( UInt eleIDw )=0;

    //! Solves the ionic model
    virtual void solveIonicModel( const vector_Type& u,
                                const Real timeStep )=0;

    //! Computes the term -1/ \Cm u^n (G (1-u^n/vp) (1-u^n/v_th) + eta_1 v^{n+1})
    //! for the PDE righthand side
    virtual void computeIonicCurrent( Real Capacitance,
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

    UInt solutionUDimension() const           { return M_uFESpace.dim(); }

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
std::cout << "IonicSolver descrutor" << std::endl;
}

template< typename Mesh,
	  typename SolverType = LifeV::SolverTrilinos >
class MitchellSchaeffer : public virtual IonicSolver<Mesh, SolverType>
{
public:
	typedef typename IonicSolver<Mesh, SolverType>::data_Type	data_Type;
	typedef typename IonicSolver<Mesh, SolverType>::vector_Type	vector_Type;
	typedef typename IonicSolver<Mesh, SolverType>::function_Type 	function_Type;
	typedef typename IonicSolver<Mesh, SolverType>::functorTauClose_Type  functorTauClose_Type;

    MitchellSchaeffer( const data_Type& dataType,
                        const Mesh& mesh,
                        FESpace<Mesh, EpetraMap>& uFEspace,
                        Epetra_Comm& comm );

    virtual ~MitchellSchaeffer();

    void updateRepeated( );

    void updateElementSolution( UInt eleID );

    void solveIonicModel( const vector_Type& u, const Real timeStep );

    void computeIonicCurrent( Real Capacitance,
                      ElemVec& elvec,
                      ElemVec& elvec_u,
                      FESpace<Mesh, EpetraMap>& uFESpace );

    const vector_Type& solutionGatingW() const {return M_solutionGatingW;}

    void initialize( );

    void setHeteroTauClose(functorTauClose_Type);

    Real functorTauClose(const entityFlag_Type& ref,
                       const Real& x,
                       const Real& y,
                       const Real& z,
                       const ID& i) const;

protected:

    //! Global solution _w
    vector_Type                    	M_solutionGatingW;
    vector_Type				M_solutionGatingWRepeated;
    ElemVec 				M_elvec;
    UInt				BDForder;
    BdfT<vector_Type> 			bdf_w;
    functorTauClose_Type	M_TauClose;

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
    M_solutionGatingW ( IonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingWRepeated( M_solutionGatingW, Repeated ),
    M_elvec ( IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof(), 1 ),
    BDForder ( IonicSolver<Mesh, SolverType>::M_data.M_BDForder ),
    bdf_w( BDForder )
{

}

template<typename Mesh, typename SolverType>
MitchellSchaeffer<Mesh, SolverType>::
~MitchellSchaeffer()
{
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::updateRepeated( )
{
	M_solutionGatingWRepeated=M_solutionGatingW;
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::updateElementSolution( UInt eleID )
{
	M_elvec.zero();
	UInt ig;
    //! Filling local elvec_w with recovery variable values in the nodes
    for ( UInt iNode = 0 ; iNode < IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof() ; iNode++ )
    {
        ig = IonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobal( eleID, iNode + 1 );
        M_elvec.vec()[ iNode ] = M_solutionGatingWRepeated[ig];
    }
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::setHeteroTauClose(functorTauClose_Type fct)
{
    M_TauClose = fct;
}

template<typename Mesh, typename SolverType>
Real MitchellSchaeffer<Mesh, SolverType>::functorTauClose(const entityFlag_Type& ref,
                                                        const Real& x,
                                                        const Real& y,
                                                        const Real& z, const ID& i) const
{
    return M_TauClose(ref, x, y, z, i);
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::solveIonicModel( const vector_Type& u, const Real timeStep )
{
    //! Solving :
	//!           ((v_max-v_min)^{-2}  - w )/tau_open  if u < vcrit
	//! dw/dt ={
	//!            -w/tau_close   if u > vcrit

	Real aux1 = 1.0 / (bdf_w.coeff_der(0)/timeStep +
                       1.0/this->M_data.M_tau_open );
	Real aux = 1.0/((this->M_data.M_potentialMaximum -
                     this->M_data.M_potentialMinimum) *
                    (this->M_data.M_potentialMaximum -
                     this->M_data.M_potentialMinimum) *
                    this->M_data.M_tau_open);
	Real aux2 = 1.0 / (bdf_w.coeff_der(0)/timeStep +
                       1.0/this->M_data.M_tau_close);

	vector_Type M_time_der=bdf_w.time_der(timeStep);

	IonicSolver<Mesh, SolverType>::M_comm->Barrier();
	Real x, y, z;
	entityFlag_Type ref;
	UInt ID;
	for ( UInt i = 0 ; i < u.epetraVector().MyLength() ; ++i )
	{
        Int ig=u.blockMap().MyGlobalElements()[i];
		ID 	= ig;
		/*ref = this->M_mesh->point(ig).marker();
		x 	= this->M_mesh->point(ig).x();
        y 	= this->M_mesh->point(ig).y();
        z 	= this->M_mesh->point(ig).z();*/
   		if (u[ig] < this->M_data.M_criticalPotential)
            M_solutionGatingW[ig] = aux1 * (aux + M_time_der[ig]);
        else if (this->M_data.M_hasHeterogeneousTauClose)
            M_solutionGatingW[ig] = (1.0 / (bdf_w.coeff_der(0)/timeStep  +
                                  1.0/*fct_Tau_Close(ref,x,y,z,ID)*/)) *
                M_time_der[ig];//aux2 * M_time_der[ig];
        else
            M_solutionGatingW[ig] = aux2 *  M_time_der[ig];
	}
	bdf_w.shift_right(M_solutionGatingW);

	IonicSolver<Mesh, SolverType>::M_comm->Barrier();
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::computeIonicCurrent(  Real,
                                                         ElemVec& elvec,
                                                         ElemVec& elvec_u,
                                                         FESpace<Mesh, EpetraMap>& uFESpace )
{
	for ( UInt i = 0;i < uFESpace.fe().nbFEDof();i++ )
    {
        elvec( i ) =  this->M_data.M_reactionAmplitude*
            (((M_elvec( i ) / this->M_data.M_tau_in) *
              (elvec_u( i ) - this->M_data.M_potentialMinimum) *
              (elvec_u( i ) - this->M_data.M_potentialMinimum) *
              (this->M_data.M_potentialMaximum - elvec_u( i )) /
              (this->M_data.M_potentialMaximum - this->M_data.M_potentialMinimum) ) -
             (( elvec_u( i ) - this->M_data.M_potentialMinimum ) /
              (  this->M_data.M_tau_out * (this->M_data.M_potentialMaximum -
                                           this->M_data.M_potentialMinimum)))) ;
    }
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::
initialize( )
{
	M_solutionGatingW.epetraVector().PutScalar(1.0 /
                                        ((this->M_data.M_potentialMaximum-
                                          this->M_data.M_potentialMinimum)*
                                         (this->M_data.M_potentialMaximum-
                                          this->M_data.M_potentialMinimum)));
	bdf_w.initialize_unk(M_solutionGatingW);
	bdf_w.showMe();
}

template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class RogersMcCulloch : public virtual IonicSolver<Mesh, SolverType>
{
public:
	typedef typename IonicSolver<Mesh, SolverType>::data_Type data_Type;
	typedef typename IonicSolver<Mesh, SolverType>::vector_Type vector_Type;
	typedef typename IonicSolver<Mesh, SolverType>::function_Type function_Type;

    RogersMcCulloch( const data_Type& dataType,
                      const Mesh& mesh,
                      FESpace<Mesh, EpetraMap>& uFEspace,
                      Epetra_Comm& comm );
    virtual ~RogersMcCulloch();

    void updateRepeated( );

    void updateElementSolution( UInt eleID );

    void solveIonicModel( const vector_Type& u, const Real timeStep );

    void computeIonicCurrent( Real Capacitance,
                      ElemVec& elvec,
                      ElemVec& elvec_u,
                      FESpace<Mesh, EpetraMap>& uFESpace );

    const vector_Type& solutionGatingW() const {return M_solutionGatingW;}

    void initialize( );

protected:

    //! Global solution _w
    vector_Type                    	M_solutionGatingW;
    vector_Type						M_solutionGatingWRepeated;
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
    M_solutionGatingW ( IonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingWRepeated ( M_solutionGatingW, Repeated ),
    M_elvec ( IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof(), 1 )
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
	M_solutionGatingWRepeated=M_solutionGatingW;
}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::updateElementSolution( UInt eleID )
{
	M_elvec.zero();
	UInt ig;
    for ( UInt iNode = 0 ;
          iNode < IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof() ; iNode++ )
    {
        ig = IonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobal( eleID, iNode + 1 );
        M_elvec.vec()[ iNode ] = M_solutionGatingWRepeated[ig];
    }
}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::solveIonicModel( const vector_Type& u, const Real timeStep )
{
	//! Solving dw/dt= b/(A*T) (u - u0 - A d w)
	Real G = this->M_data.M_b / this->M_data.M_potentialAmplitude / this->M_data.M_timeUnit;

	Real alpha = 1 / timeStep + G*this->M_data.M_potentialAmplitude * this->M_data.M_d;

	IonicSolver<Mesh, SolverType>::M_comm->Barrier();

	M_solutionGatingW*=1/timeStep;
	vector_Type temp(u);
	temp.epetraVector().PutScalar (G*this->M_data.M_restPotential);
	M_solutionGatingW+= G*u;
	M_solutionGatingW-= temp;
	M_solutionGatingW*=1/alpha;
	M_solutionGatingW.globalAssemble();
	IonicSolver<Mesh, SolverType>::M_comm->Barrier();

}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::computeIonicCurrent(  Real Capacitance,
                                                              ElemVec& elvec,
                                                              ElemVec& elvec_u,
                                                              FESpace<Mesh, EpetraMap>& uFESpace )
{
	Real u_ig, w_ig;

	Real G1 = this->M_data.M_c1 / this->M_data.M_timeUnit / pow(this->M_data.M_potentialAmplitude,2.0);
	Real G2 = this->M_data.M_c2 / this->M_data.M_timeUnit;

        for ( UInt ig = 0; ig < uFESpace.fe().nbQuadPt();ig++ )
    {
        u_ig = w_ig = 0.;
        for ( UInt i = 0;i < uFESpace.fe().nbFEDof();i++ )
            u_ig += elvec_u( i ) * uFESpace.fe().phi( i, ig );
        for ( UInt i = 0;i < uFESpace.fe().nbFEDof();i++ )
            w_ig += M_elvec( i ) * uFESpace.fe().phi( i, ig );

        for ( UInt i = 0;i < uFESpace.fe().nbFEDof();i++ )
        {
            elvec( i ) -= Capacitance*(G1*(u_ig- this->M_data.M_restPotential)*
                                       (u_ig - this->M_data.M_restPotential -
                                        this->M_data.M_a*this->M_data.M_potentialAmplitude)*
                                       (u_ig - this->M_data.M_restPotential -
                                        this->M_data.M_potentialAmplitude) + G2 *
                                       (u_ig - this->M_data.M_restPotential) * w_ig) *
                uFESpace.fe().phi( i, ig ) * uFESpace.fe().weightDet( ig );
        }
    }
}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::
initialize( )
{
	M_solutionGatingW.epetraVector().PutScalar (0.);
}

template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class LuoRudy : public virtual IonicSolver<Mesh, SolverType>
{
public:
	typedef typename IonicSolver<Mesh, SolverType>::data_Type data_Type;
	typedef typename IonicSolver<Mesh, SolverType>::vector_Type vector_Type;
	typedef typename IonicSolver<Mesh, SolverType>::function_Type function_Type;

    LuoRudy( const data_Type&          dataType,
             const Mesh&          mesh,
             FESpace<Mesh, EpetraMap>& uFEspace,
             Epetra_Comm&              comm );
    virtual ~LuoRudy();

    void updateRepeated( );

    void updateElementSolution( UInt eleID);

    void computeODECoefficients( const Real& u_ig );

    void solveIonicModel( const vector_Type& u, const Real timeStep );

    void computeIonicCurrent( Real Capacitance,
                      ElemVec& elvec,
                      ElemVec& elvec_u,
                      FESpace<Mesh, EpetraMap>& uFESpace );

    void initialize ( );

    //! Returns the local solution vector for each field
    const vector_Type& solutionGatingH() const {return M_solutionGatingH;}

    const vector_Type& solutionGatingJ() const {return M_solutionGatingJ;}

    const vector_Type& solutionGatingM() const {return M_solutionGatingM;}

    const vector_Type& solutionGatingD() const {return M_solutionGatingD;}

    const vector_Type& solutionGatingF() const {return M_solutionGatingF;}

    const vector_Type& solutionGatingX() const {return M_solutionGatingX;}

    const vector_Type& solutionGatingCa() const {return M_solutionGatingCa;}

    Real M_K0, M_Ki, M_Na0, M_Nai, M_R, M_temperature, M_F, M_permeabilityRatio, M_c,
        M_Ena, M_Gk, M_Ek, M_Gk1, M_Ek1, M_Ekp, M_Esi, M_ah, M_bh, M_aj, M_bj, M_xii,
        M_am, M_bm, M_ad, M_bd, M_af, M_bf, M_aX, M_bX, M_ak1, M_bk1, M_Kp, M_K1inf,
        M_hinf, M_tauh, M_jinf, M_tauj, M_minf, M_taum, M_dinf, M_taud, M_finf, M_tauf,
        M_Xinf, M_tauX;

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

	vector_Type M_vectorExponentialh;
	vector_Type M_vectorExponentialj;
	vector_Type M_vectorExponentialm;
	vector_Type M_vectorExponentiald;
	vector_Type M_vectorExponentialf;
	vector_Type M_vectorExponentialX;
	vector_Type M_vectorInfimumh;
	vector_Type M_vectorInfimumj;
	vector_Type M_vectorInfimumm;
	vector_Type M_vectorInfimumd;
	vector_Type M_vectorInfimumf;
	vector_Type M_vectorInfimumX;
	vector_Type	M_vectorIonicChange;

protected:
    //! Global solution h
    vector_Type                    	M_solutionGatingH;
    //! Global solution j
    vector_Type                    	M_solutionGatingJ;
    //! Global solution m
    vector_Type                    	M_solutionGatingM;
    //! Global solution d
    vector_Type                    	M_solutionGatingD;
    //! Global solution f
    vector_Type                    	M_solutionGatingF;
    //! Global solution X
    vector_Type                    	M_solutionGatingX;
    //! Global solution Ca_i
    vector_Type                    	M_solutionGatingCa;

    vector_Type 					M_ionicCurrent;

    vector_Type						M_ionicCurrentRepeated;

    ElemVec 						M_elemVecIonicCurrent;

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
			M_solutionGatingH                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_solutionGatingJ                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_solutionGatingM                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_solutionGatingD                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_solutionGatingF                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_solutionGatingX                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_solutionGatingCa                  ( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_ionicCurrent( IonicSolver<Mesh, SolverType>::M_localMap ),
			M_ionicCurrentRepeated( M_ionicCurrent, Repeated ),
			M_elemVecIonicCurrent ( IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof(), 1 ),
			M_K0(5.4),
			M_Ki(145.),
			M_Na0(140.),
			M_Nai(18.),
			M_R(8.314472), //% Joules/(Kelvin*mole)
			M_temperature(307.7532), //%kelvins
			M_F(96485.33838), //% coulumbs/mole
			M_permeabilityRatio(0.01833), //% Na/K permeability ratio
			M_c(1.), //% membrane capacitance set as 1 mui-F/cm^2
			M_Ena(1000. * (M_R * M_temperature / M_F ) * log(M_Na0 / M_Nai) ),
			M_Gk(0.282 * sqrt(M_K0 / 5.4)),
			M_Ek(1000. * (M_R * M_temperature / M_F) * log((M_K0 + M_permeabilityRatio * M_Na0)
                                                           / (M_Ki + M_permeabilityRatio * M_Nai))),
			M_Gk1(0.6047 * sqrt(M_K0 / 5.4)),
			M_Ek1(1000.* (M_R * M_temperature / M_F)* log(M_K0 / M_Ki)),
			M_Ekp(M_Ek1),
			M_vectorExponentialh(IonicSolver<Mesh, SolverType>::M_localMap),
			M_vectorExponentialj(IonicSolver<Mesh, SolverType>::M_localMap),
			M_vectorExponentialm(IonicSolver<Mesh, SolverType>::M_localMap),
			M_vectorExponentiald(IonicSolver<Mesh, SolverType>::M_localMap),
			M_vectorExponentialf(IonicSolver<Mesh, SolverType>::M_localMap),
			M_vectorExponentialX(IonicSolver<Mesh, SolverType>::M_localMap),
			M_vectorInfimumh(IonicSolver<Mesh, SolverType>::M_localMap),
			M_vectorInfimumj(IonicSolver<Mesh, SolverType>::M_localMap),
			M_vectorInfimumm(IonicSolver<Mesh, SolverType>::M_localMap),
			M_vectorInfimumd(IonicSolver<Mesh, SolverType>::M_localMap),
			M_vectorInfimumf(IonicSolver<Mesh, SolverType>::M_localMap),
			M_vectorInfimumX(IonicSolver<Mesh, SolverType>::M_localMap),
			M_vectorIonicChange(IonicSolver<Mesh, SolverType>::M_localMap)
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

	M_ionicCurrentRepeated=M_ionicCurrent;
}

template<typename Mesh, typename SolverType>
void LuoRudy<Mesh, SolverType>::updateElementSolution( UInt eleID)
{
	M_elemVecIonicCurrent.zero();
	UInt ig;
		//! Filling local elvec with recovery variable values in the nodes
		for ( UInt iNode = 0 ; iNode < IonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof() ; iNode++ )
		{
			ig = IonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobal( eleID, iNode + 1 );
			M_elemVecIonicCurrent.vec()[ iNode ] = M_ionicCurrentRepeated[ig];
		}

}


template<typename Mesh, typename SolverType>
void LuoRudy<Mesh, SolverType>::solveIonicModel( const vector_Type& u, const Real timeStep )
{
	//! Solving dw/dt=eta2 (u/vp -  eta3 w)
	Chrono chronoionmodelsolve;
	chronoionmodelsolve.start();
	IonicSolver<Mesh, SolverType>::M_comm->Barrier();

	for ( Int i = 0 ; i < u.epetraVector().MyLength() ; i++ )
	{
		Int ig=u.blockMap().MyGlobalElements()[i];
		Real u_ig=u[ig];
		computeODECoefficients(u_ig);
        M_Esi = 7.7 - 13.0287 * log(M_solutionGatingCa[ig]);
		//fast sodium current
		M_Ina = 23.* M_solutionGatingM[ig] * M_solutionGatingM[ig] * M_solutionGatingM[ig] * M_solutionGatingH[ig] * M_solutionGatingJ[ig] * (u_ig - M_Ena);
		//slow inward current
		M_Islow = 0.09 * M_solutionGatingD[ig] * M_solutionGatingF[ig] * (u_ig - M_Esi);
        //change in ioniq concentration
        M_vectorIonicChange.epetraVector().ReplaceGlobalValue(ig,
                                                                 0,
                                                                 -1e-4 * M_Islow + 0.07*(1e-4 - M_solutionGatingCa[ig]));
        //time dependent potassium current
        M_Ik = M_Gk * M_solutionGatingX[ig] * M_xii*(u_ig - M_Ek);
        //time independent potassium current
        M_Ik1 = M_Gk1 * M_K1inf * (u_ig - M_Ek1);
        //plateau potassium current
        M_Ikp = 0.0183 * M_Kp * (u_ig - M_Ekp);
        //background current
        M_Iback = 0.03921 * (u_ig + 59.87);
        //Total time independent potassium current
        M_Ik1t = M_Ik1 + M_Ikp + M_Iback;
        // adding up the six ionic currents
        M_ionicCurrent.epetraVector().ReplaceGlobalValue(ig,
                                                  0,
                                                  M_Ina + M_Islow + M_Ik + M_Ik1t);
		M_vectorExponentialh.epetraVector().ReplaceGlobalValue(ig,
                                                                  0,
                                                                  exp(-timeStep / M_tauh));
		M_vectorExponentialj.epetraVector().ReplaceGlobalValue(ig,
                                                                  0,
                                                                  exp(-timeStep / M_tauj));
		M_vectorExponentialm.epetraVector().ReplaceGlobalValue(ig,
                                                                  0,
                                                                  exp(-timeStep / M_taum));
		M_vectorExponentiald.epetraVector().ReplaceGlobalValue(ig,
                                                                  0,
                                                                  exp(-timeStep / M_taud));
		M_vectorExponentialf.epetraVector().ReplaceGlobalValue(ig,
                                                                  0,
                                                                  exp(-timeStep / M_tauf));
		M_vectorExponentialX.epetraVector().ReplaceGlobalValue(ig,
                                                                  0,
                                                                  exp(-timeStep / M_tauX));
		M_vectorInfimumh.epetraVector().ReplaceGlobalValue(ig,
                                                              0,
                                                              M_hinf);
		M_vectorInfimumj.epetraVector().ReplaceGlobalValue(ig,
                                                              0,
                                                              M_jinf);
		M_vectorInfimumm.epetraVector().ReplaceGlobalValue(ig,
                                                              0,
                                                              M_minf);
		M_vectorInfimumd.epetraVector().ReplaceGlobalValue(ig,
                                                              0,
                                                              M_dinf);
		M_vectorInfimumf.epetraVector().ReplaceGlobalValue(ig,
                                                              0,
                                                              M_finf);
		M_vectorInfimumX.epetraVector().ReplaceGlobalValue(ig,
                                                              0,
                                                              M_Xinf);
	}
	M_vectorExponentialh.globalAssemble();
	M_vectorExponentialj.globalAssemble();
	M_vectorExponentialm.globalAssemble();
	M_vectorExponentiald.globalAssemble();
	M_vectorExponentialf.globalAssemble();
	M_vectorExponentialX.globalAssemble();
	M_vectorInfimumh.globalAssemble();
	M_vectorInfimumj.globalAssemble();
	M_vectorInfimumm.globalAssemble();
	M_vectorInfimumd.globalAssemble();
	M_vectorInfimumf.globalAssemble();
	M_vectorInfimumX.globalAssemble();
	M_ionicCurrent.globalAssemble();
	M_vectorIonicChange.globalAssemble();

	M_solutionGatingH-=M_vectorInfimumh;
	M_solutionGatingH.epetraVector().Multiply(1.,
                                       M_solutionGatingH.epetraVector(),
                                       M_vectorExponentialh.epetraVector(),
                                       0.);
	M_solutionGatingH+=M_vectorInfimumh;
	M_solutionGatingJ-=M_vectorInfimumj;

	M_solutionGatingJ.epetraVector().Multiply(1.,
                                       M_solutionGatingJ.epetraVector(),
                                       M_vectorExponentialj.epetraVector(),
                                       0.);
	M_solutionGatingJ+=M_vectorInfimumj;
	M_solutionGatingM-=M_vectorInfimumm;
	M_solutionGatingM.epetraVector().Multiply(1.,
                                       M_solutionGatingM.epetraVector(),
                                       M_vectorExponentialm.epetraVector(),
                                       0.);
	M_solutionGatingM+=M_vectorInfimumm;
	M_solutionGatingD-=M_vectorInfimumd;

	M_solutionGatingD.epetraVector().Multiply(1.,
                                       M_solutionGatingD.epetraVector(),
                                       M_vectorExponentiald.epetraVector(),
                                       0.);
	M_solutionGatingD+=M_vectorInfimumd;
	M_solutionGatingF-=M_vectorInfimumf;

	M_solutionGatingF.epetraVector().Multiply(1.,
                                       M_solutionGatingF.epetraVector(),
                                       M_vectorExponentialf.epetraVector(),
                                       0.);
	M_solutionGatingF+=M_vectorInfimumf;
	M_solutionGatingX-=M_vectorInfimumX;

	M_solutionGatingX.epetraVector().Multiply(1.,
                                       M_solutionGatingX.epetraVector(),
                                       M_vectorExponentialX.epetraVector(),
                                       0.);
	M_solutionGatingX+=M_vectorInfimumX;
	M_solutionGatingCa+=timeStep*M_vectorIonicChange;

	M_solutionGatingH.globalAssemble();
	M_solutionGatingJ.globalAssemble();
	M_solutionGatingM.globalAssemble();
	M_solutionGatingD.globalAssemble();
	M_solutionGatingF.globalAssemble();
	M_solutionGatingX.globalAssemble();
	M_solutionGatingCa.globalAssemble();

	IonicSolver<Mesh, SolverType>::M_comm->Barrier();

	chronoionmodelsolve.stop();
    if (IonicSolver<Mesh, SolverType>::M_comm->MyPID()==0)
        std::cout << "Total ionmodelsolve time " << chronoionmodelsolve.diff() << " s." << std::endl;
}



template<typename Mesh, typename SolverType>
void LuoRudy<Mesh, SolverType>::computeODECoefficients( const Real& u_ig )
{
    if (u_ig >= -40.)
    {
        M_ah=0.;
        M_bh = 1. / (0.13 * (1. + exp( (u_ig + 10.66) / (-11.1) )));
        M_aj=0.;
        M_bj = 0.3 * exp(-2.535e-7 * u_ig) / (1. + exp(-0.1 * (u_ig + 32.)));
    }
    else
    {
        M_ah = 0.135 * exp((80. + u_ig) / -6.8);
        M_bh = 3.56 * exp(0.079 *u_ig) + 3.1e5 * exp(0.35 * u_ig);
        M_aj = (-1.2714e5 * exp(0.2444 * u_ig)-3.474e-5 * exp(-0.04391 * u_ig))*
            (u_ig + 37.78) / (1 + exp(0.311 * (u_ig + 79.23)));
        M_bj = 0.1212 * exp(-0.01052 * u_ig) / (1. + exp(-0.1378 * (u_ig + 40.14)));
    }
    M_am = 0.32 * (u_ig + 47.13) / (1. - exp(-0.1 * (u_ig + 47.13)));
    M_bm = 0.08*exp(-u_ig/11.);

    //slow inward current
    M_ad = 0.095 * exp(-0.01 *(u_ig - 5.)) / (1. + exp(-0.072*(u_ig - 5.)));
    M_bd = 0.07  * exp(-0.017*(u_ig + 44.))/ (1. + exp( 0.05 *(u_ig + 44.)));
    M_af = 0.012 * exp(-0.008*(u_ig + 28.))/ (1. + exp( 0.15 *(u_ig + 28.)));
    M_bf = 0.0065* exp(-0.02 *(u_ig + 30.))/ (1. + exp( -0.2 *(u_ig + 30.)));

    //Time dependent potassium outward current
    M_aX = 0.0005 * exp(0.083 * (u_ig +50.)) / (1. + exp(0.057 * (u_ig + 50.)));
    M_bX = 0.0013 * exp(-0.06 * (u_ig +20.)) / (1. + exp(-0.04 * (u_ig + 20.)));

    if(u_ig<=-100)
        M_xii = 1.;
    else
        M_xii = 2.837 * (exp (0.04 * (u_ig + 77.)) -1.) / ((u_ig + 77.) * exp(0.04 * (u_ig + 35.)));
    M_ak1 = 1.02 / (1. + exp(0.2385 * (u_ig - M_Ek1 - 59.215)));
    M_bk1 = (0.49124 * exp(0.08032 * (u_ig - M_Ek1 + 5.476)) +
             exp(0.06175 * (u_ig - M_Ek1 - 594.31))) / (1. + exp(-0.5143 * (u_ig - M_Ek1 + 4.753)));
    //Plateau potassium outward current
    M_Kp = 1. / (1. + exp((7.488 - u_ig) / 5.98));

    M_K1inf = M_ak1 / (M_ak1 + M_bk1);
    M_hinf = M_ah   / (M_ah  + M_bh);
    M_tauh = 1.   / (M_ah + M_bh);
    M_jinf = M_aj / (M_aj + M_bj);
    M_tauj = 1.   / (M_aj + M_bj);
    M_minf = M_am / (M_am + M_bm);
    M_taum = 1.   / (M_am + M_bm);
    M_dinf = M_ad / (M_ad + M_bd);
    M_taud = 1.   / (M_ad + M_bd);
    M_finf = M_af / (M_af + M_bf);
    M_tauf = 1.   / (M_af + M_bf);
    M_Xinf = M_aX / (M_aX + M_bX);
    M_tauX = 1.   / (M_aX + M_bX);
}

template<typename Mesh, typename SolverType>
void LuoRudy<Mesh, SolverType>::computeIonicCurrent(  Real Capacitance,
                                              ElemVec& elvec,
                                              ElemVec& /*elvec_u*/,
                                              FESpace<Mesh, EpetraMap>& uFESpace )
{
   	Real Iion_ig;
    for ( UInt ig = 0; ig < uFESpace.fe().nbQuadPt();ig++ )
    {
        Iion_ig = 0.;
        for ( UInt i = 0;i < uFESpace.fe().nbFEDof();i++ )
            Iion_ig += M_elemVecIonicCurrent( i ) * uFESpace.fe().phi( i, ig );
        for ( UInt i = 0;i < uFESpace.fe().nbFEDof();i++ )
        {
            // divide by 1000 to convert microA in mA
            elvec( i ) -= Iion_ig * Capacitance *
                uFESpace.fe().phi( i, ig ) * uFESpace.fe().weightDet( ig );
        }
    }
}

template<typename Mesh, typename SolverType>
void LuoRudy<Mesh, SolverType>::
initialize( )
{
    M_solutionGatingH.epetraVector().PutScalar(1.);
    M_solutionGatingJ.epetraVector().PutScalar(1.);
    M_solutionGatingM.epetraVector().PutScalar(0.);
    M_solutionGatingD.epetraVector().PutScalar(0.);
    M_solutionGatingF.epetraVector().PutScalar(1.);
    M_solutionGatingX.epetraVector().PutScalar(0.);
    M_solutionGatingCa.epetraVector().PutScalar(0.0002);
    M_solutionGatingH.globalAssemble();
    M_solutionGatingJ.globalAssemble();
    M_solutionGatingM.globalAssemble();
    M_solutionGatingD.globalAssemble();
    M_solutionGatingF.globalAssemble();
    M_solutionGatingX.globalAssemble();
    M_solutionGatingCa.globalAssemble();
}

} // namespace LifeV


#endif //_IONICSOLVER_H_
