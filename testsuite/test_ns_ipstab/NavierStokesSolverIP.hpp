/*!
  \file NavierStokesSolverIP.h
  \author M.A. Fernandez
  \date 6/2003
  \version 1.0

  \brief This file contains a Navier-Stokes solver class which implements a
  stabilized implicit scheme with coupled solving.
*/
#ifndef _NAVIERSTOKESSOLVERIP_H_
#define _NAVIERSTOKESSOLVERIP_H_

#include "NavierStokesHandler.hpp"
#include "elemMat.hpp"
#include "elemVec.hpp"
#include "elemOper.hpp"
#include "values.hpp"
#include "pattern.hpp"
#include "assemb.hpp"
#include "bc_manage.hpp"
#include "SolverAztec.hpp"
//#include "lifeconfig.h"
//#if defined(HAVE_PETSC_H)
//#include "SolverPETSC.hpp"
//#endif /* HAVE_PETSC_H */
#include "bcCond.hpp"
#include "chrono.hpp"
#include "sobolevNorms.hpp"
#include "geoMap.hpp"

namespace LifeV
{
/*!
  \class NavierStokesSolverIP

  This class implements an NavierStokes solver via interior penalty
  stabilization. The resulting linear systems are solved by GMRES on the full
  matrix (u and p coupled).

*/
template<typename Mesh>
class NavierStokesSolverIP:
        public NavierStokesHandler<Mesh> {

public:
    typedef  typename  NavierStokesHandler<Mesh>::Function Function;

    //! Constructor
    /*!
      \param data_file GetPot data file
      \param refFE reference FE
      \param Qr volumic quadrature rule
      \param bdQr surface quadrature rule
      \param BCh boundary conditions for the velocity
    */
    NavierStokesSolverIP(const GetPot& data_file, const RefFE& refFE,
                         const QuadRule& Qr, const QuadRule& bdQr,
                         BC_Handler& BCh);

    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
    void timeAdvance(const Function source, const Real& time);

    //! Update convective term, bc treatment and solve the linearized ns system
    void iterate(const Real& time);

    void eval(Vector& fx0,Vector& gx0,Vector x0, int status);

private:

    //! Block pattern of M_u
    MSRPatt _pattM_u_block;

    //! Pattern for M
    MixedPattern<nDimensions,nDimensions,MSRPatt> _pattM_u;

    //! Block pattern of C
    MSRPatt _pattC;

    //! Matrix M_u: Vmass
    MixedMatr<nDimensions,nDimensions,MSRPatt,double> _M_u;

    //! Matrix CStokes
    MSRMatr<double> _CStokes;

    //! Matrix C: 1/dt*Vmass + mu*Vstiff operator + Convective_term
    MSRMatr<double> _C;

    //! Elementary matrices and vectors
    ElemMat _elmatC; //velocity stiffnes
    ElemMat _elmatM_u; //velocity mass
    ElemMat _elmatD;
    ElemMat _elmatDtr;
    ElemMat _elmatP;
    ElemVec _elvec; // Elementary right hand side

    //! Right  hand  side for the velocity
    PhysVectUnknown<Vector> _f_u;

    //! Right  hand  side global
    Vector _b;

    //! Global solution _u and _p
    Vector _x;

    //! Global solution _u and _p
    Vector _diff;

    SolverAztec _solver;
    //SolverPETSC _solver;

    Real _time;
};


template<typename Mesh> void NavierStokesSolverIP<Mesh>::
eval(Vector& fx0,Vector& gx0,Vector x0, int status) {

    iterate(0.0);
    for (UInt i=0; i < nDimensions*_dim_u ; ++i)
        fx0[i] = _u[i];

}

//
//                                         IMPLEMENTATION
//
template<typename Mesh> NavierStokesSolverIP<Mesh>::
NavierStokesSolverIP(const GetPot& data_file, const RefFE& refFE,
                     const QuadRule& Qr, const QuadRule& bdQr,
                     BC_Handler& BCh):
    NavierStokesHandler<Mesh>(data_file,refFE,refFE,Qr,bdQr,Qr,bdQr,BCh),
    _pattM_u_block(_dof_u),
    _pattM_u(_pattM_u_block,"diag"),
    _pattC(_dof_u,_mesh,nDimensions+1),
    _M_u(_pattM_u),
    _CStokes(_pattC),
    _C(_pattC),
    _elmatC(_fe_u.nbNode,nDimensions,nDimensions),
    _elmatM_u(_fe_u.nbNode,nDimensions,nDimensions),
    _elmatD(_fe_u.nbNode,nDimensions+1,nDimensions),
    _elmatDtr(_fe_u.nbNode,nDimensions,nDimensions+1),
    _elmatP(_fe_u.nbNode,nDimensions+1,nDimensions+1),
    _elvec(_fe_u.nbNode,nDimensions),
    _f_u(_dim_u),
    _b((nDimensions+1)*_dim_u),
    _x((nDimensions+1)*_dim_u),
    _diff((nDimensions+1)*_dim_u) 
{


    _solver.setOptionsFromGetPot(data_file, "fluid/aztec");
    //_solver.setOptionsFromGetPot(data_file, "fluid/petsc");

    std::cout << std::endl;
    std::cout << "O-  Pressure unknowns: " << _dim_p     << std::endl;
    std::cout << "O-  Velocity unknowns: " << _dim_u     << std::endl
              << std::endl;
    std::cout << "O-  Computing mass and Stokes matrices... " << std::flush;

    Chrono chrono;
    chrono.start();

    // Matrices initialization
    _M_u.zeros();
    _CStokes.zeros();

    // Number of velocity components
    UInt nc_u=_u.nbcomp();

    //inverse of dt:
    // Real dti=1./_dt;




    Real gamma;

    // Elementary computation and matrix assembling
    // Loop on elements
    for(UInt i = 1; i <= _mesh.numVolumes(); i++){

        _fe_u.updateFirstDeriv(_mesh.volumeList(i));

        _elmatC.zero();
        _elmatM_u.zero();
        _elmatD.zero();
        _elmatDtr.zero();

        // stiffness strain
        stiff_strain(2.0*_mu,_elmatC,_fe_u);
        //stiff_div(0.5*_fe_u.diameter(),_elmatC,_fe_u);

        // mass
        //mass(_rho*dti,_elmatM_u,_fe_u,0,0,nDimensions);
        // _elmatC.mat() += _elmatM_u.mat();


        for(UInt ic=0;ic<nc_u;ic++){
            for(UInt jc=0;jc<nc_u;jc++) {
                // stiffness
                assemb_mat(_CStokes,_elmatC,_fe_u,_dof_u,ic,jc);
            }
            // mass
            //assemb_mat(_M_u,_elmatM_u,_fe_u,_dof_u,ic,ic);

            // div
            grad(ic,1.0,_elmatDtr,_fe_u,_fe_u,ic,nDimensions);
            div( ic,-1.0,_elmatD  ,_fe_u,_fe_u,nc_u,ic);

            // assembling
            assemb_mat(_CStokes,_elmatDtr,_fe_u,_dof_u,ic,nc_u);
            assemb_mat(_CStokes,_elmatD,_fe_u,_dof_u,nc_u,ic);

        }
    }


    std::cout << std::endl;

    UInt iElAd1, iElAd2;
    CurrentFE fe1(_refFE_u,getGeoMap(_mesh),_Qr_u);
    CurrentFE fe2(_refFE_u,getGeoMap(_mesh),_Qr_u);

    _elmatP.zero();

    // Elementary computation and matrix assembling
    // Loop on interior faces
    for (UInt i=_mesh.numBFaces()+1; i<=_mesh.numFaces(); ++i){


        // Updating face staff
        _feBd_u.updateMeas( _mesh.face(i) );

        gamma = _feBd_u.measure()/32.0;   // P1
        //  gamma = _feBd_u.measure()/128.0; // P2
        // gamma = _feBd_u.measure()*sqrt(_feBd_u.measure())/8.0; // P1 p non smooth

        iElAd1 = _mesh.face(i).ad_first();
        iElAd2 = _mesh.face(i).ad_second();

        fe1.updateFirstDeriv(_mesh.volumeList( iElAd1 ) );
        fe2.updateFirstDeriv(_mesh.volumeList( iElAd2 ) );

        ipstab_grad(gamma,_elmatP, fe1, fe1, _feBd_u,nDimensions,nDimensions);
        assemb_mat(_CStokes,_elmatP,fe1,_dof_u,nDimensions,nDimensions);

        ipstab_grad(gamma,_elmatP, fe2, fe2, _feBd_u,nDimensions,nDimensions);
        assemb_mat(_CStokes,_elmatP,fe2,_dof_u,nDimensions,nDimensions);

        ipstab_grad(-gamma,_elmatP, fe1, fe2, _feBd_u,nDimensions,nDimensions);
        assemb_mat(_CStokes,_elmatP,fe1,fe2,_dof_u,nDimensions,nDimensions);

        ipstab_grad(-gamma,_elmatP, fe2, fe1, _feBd_u,nDimensions,nDimensions);
        assemb_mat(_CStokes,_elmatP,fe2,fe1,_dof_u,nDimensions,nDimensions);

    }

    _x = 0.0;

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
    // _CStokes.spy("CS.m");


}

template<typename Mesh>
void NavierStokesSolverIP<Mesh>::
timeAdvance(const Function source, const Real& time) {


    _time = time;

    std::cout << std::endl;
    std::cout << "O== Now we are at time "<< _time << " s." << std::endl;

    // Number of velocity components
    UInt nc_u=_u.nbcomp();

    std::cout << "  o-  Updating mass term on right hand side... "
              << std::flush;

    Chrono chrono;
    chrono.start();

    // Right hand side for the velocity at time
    _f_u=0.0;

    // loop on volumes: assembling source term
    for(UInt i=1; i<=_mesh.numVolumes(); ++i){
        _elvec.zero();
        _fe_u.updateJacQuadPt(_mesh.volumeList(i));

        for (UInt ic=0; ic<nc_u; ++ic) {
            // compute local vector
            compute_vec(source,_elvec,_fe_u,_time,ic);
            
            // assemble local vector into global one
            assemb_vec(_f_u,_elvec,_fe_u,_dof_u,ic);
        }
    }

    //  For the moment steady: without mass term
    //
    //  _f_u += _M_u * _u;

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
}




template<typename Mesh>
void NavierStokesSolverIP<Mesh>::iterate(const Real& time) {

    Chrono  chrono;


    _time = time;
    chrono.start();
    _C=_CStokes;

    chrono.stop();

    std::cout << "  o-  Stokes matrix was copied in " << chrono.diff() << "s."
              << std::endl;

    std::cout << "  o-  Updating convective term... " << std::flush;

    // Number of velocity components
    UInt nc_u=_u.nbcomp();


    chrono.start();

    // loop on volumes
    for(UInt i=1; i<=_mesh.numVolumes(); ++i){

        _fe_u.updateFirstDeriv(_mesh.volumeList(i)); // as updateFirstDer

        _elmatC.zero();

        UInt eleID = _fe_u.currentId();
        // Non linear term, Semi-implicit approach
        // ULoc contains the velocity values in the nodes
        for (UInt k=0 ; k<(UInt)_fe_u.nbNode ; k++){
            UInt  iloc = _fe_u.patternFirst(k);
            for (UInt ic=0; ic<nc_u; ++ic){
                UInt ig=_dof_u.localToGlobal(eleID,iloc+1)-1+ic*_dim_u;
                _elvec.vec()[iloc+ic*_fe_u.nbNode] = _rho * _u(ig);
            }
        }

        // Stabilising term: div u^n u v
        //    mass_divw(0.5*_rho,_elvec,_elmatC,_fe_u,0,0,nc_u);

        // loop on components
        for (UInt ic=0; ic<nc_u; ++ic){
            // compute local convective term and assembling
            grad(0,_elvec,_elmatC,_fe_u,_fe_u,ic,ic);
            grad(1,_elvec,_elmatC,_fe_u,_fe_u,ic,ic);
            grad(2,_elvec,_elmatC,_fe_u,_fe_u,ic,ic);
            assemb_mat(_C,_elmatC,_fe_u,_dof_u,ic,ic);
        }
    }

    UInt iElAd1, iElAd2, ig, iFaEl,iloc;
    CurrentFE fe1(_refFE_u,getGeoMap(_mesh),_Qr_u);
    CurrentFE fe2(_refFE_u,getGeoMap(_mesh),_Qr_u);
    Real gamma_b,gamma_u ,bmax,bmax_u;

    ElemVec beta(_feBd_u.nbNode,nDimensions); // local trace of the velocity

    typedef ID (*FTOP)(ID const _localFace, ID const _point);

    FTOP ftop=0;

    switch( _fe_u.nbNode ) {
        case 4:
            ftop = LinearTetra::fToP;
            break;
        case 10:
            ftop = QuadraticTetra::fToP;
            break;
        case 8:
            ftop = LinearHexa::fToP;
            break;
        case 20:
            ftop = QuadraticHexa::fToP;
            break;
        default:
            ERROR_MSG("This refFE is not allowed with IP stabilisation");
            break;
    }

    _elmatC.zero();

    // Elementary computation and matrix assembling
    // Loop on interior faces
    for (UInt i=_mesh.numBFaces()+1; i<=_mesh.numFaces(); ++i){


        // Updating face staff
        _feBd_u.updateMeas( _mesh.face(i) );


        iElAd1 = _mesh.face(i).ad_first();
        iElAd2 = _mesh.face(i).ad_second();

        fe1.updateFirstDeriv(_mesh.volumeList( iElAd1 ));
        fe2.updateFirstDeriv(_mesh.volumeList( iElAd2 ));
        
        // local id of the face in its adjacent element
        iFaEl = _mesh.face(i).pos_first();

        for(int in=0; in < _feBd_u.nbNode; ++in) {
            iloc = ftop(iFaEl,in+1);
            for(int ic=0; ic < fe1.nbCoor; ++ic) {
                ig=_dof_u.localToGlobal(iElAd1,iloc+1)-1+ic*_dim_u;
                beta.vec()[ic*_feBd_u.nbNode + in]= _u(ig);
            }
        }

        bmax = fabs(beta.vec()[0]);
        for (int l=1; l < int(fe1.nbCoor*_feBd_u.nbNode); ++l) {
            if ( bmax < fabs(beta.vec()[l]) )
                bmax = fabs(beta.vec()[l]);
        }

        bmax_u = bmax;
        if ( bmax_u < _feBd_u.measure() )
            bmax_u = _feBd_u.measure();

        gamma_b = 0.125*_feBd_u.measure()/bmax_u;
        gamma_u = 0.125*_feBd_u.measure()*bmax;
        //gamma_u = 0.125*sqrt(_feBd_u.measure())*bmax;


        _elmatC.zero();
        //ipstab_grad(gamma_u,_elmatC, fe1, fe1, _feBd_u,0,0,nDimensions);
        ipstab_bgrad(gamma_b,_elmatC, fe1, fe1,beta,_feBd_u,0,0,nDimensions);
        ipstab_div(gamma_u,_elmatC, fe1, fe1, _feBd_u);
        for (UInt ic=0; ic<nDimensions; ++ic)
            for (UInt jc=0; jc<nDimensions; ++jc)
                assemb_mat(_C,_elmatC,fe1,_dof_u,ic,jc);

        _elmatC.zero();
        //ipstab_grad(gamma_u,_elmatC, fe2, fe2, _feBd_u,0,0,nDimensions);
        ipstab_bgrad(gamma_b,_elmatC, fe2, fe2,beta,_feBd_u,0,0,nDimensions);
        ipstab_div(gamma_u,_elmatC, fe2, fe2, _feBd_u);
        for (UInt ic=0; ic<nDimensions; ++ic)
            for (UInt jc=0; jc<nDimensions; ++jc)
                assemb_mat(_C,_elmatC,fe2,_dof_u,ic,jc);

        _elmatC.zero();
        //ipstab_grad(-gamma_u,_elmatC, fe1, fe2, _feBd_u,0,0,nDimensions);
        ipstab_bgrad(-gamma_b,_elmatC,fe1,fe2,beta,_feBd_u,0,0,nDimensions);
        ipstab_div(  -gamma_u,_elmatC,fe1,fe2,_feBd_u);
        for (UInt ic=0; ic<nDimensions; ++ic)
            for (UInt jc=0; jc<nDimensions; ++jc)
                assemb_mat(_C,_elmatC,fe1,fe2,_dof_u,ic,jc);

        _elmatC.zero();
        //ipstab_grad(-gamma_u,_elmatC, fe2, fe1, _feBd_u,0,0,nDimensions);
        ipstab_bgrad(-gamma_b,_elmatC, fe2, fe1,beta,_feBd_u,0,0,nDimensions);
        ipstab_div(-gamma_u,_elmatC, fe2, fe1, _feBd_u);
        for (UInt ic=0; ic<nDimensions; ++ic)
            for (UInt jc=0; jc<nDimensions; ++jc)
                assemb_mat(_C,_elmatC,fe2,fe1,_dof_u,ic,jc);

    }


    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;



    // for BC treatment (done at each time-step)
    std::cout << "  o-  Applying boundary conditions... " << std::flush;
    chrono.start();



    _b = 0.0;
    for (UInt i=0; i<nDimensions*_dim_u; ++i) {
        _b[i]=_f_u[i];
    }


    // BC manage for the velocity
    if ( !_BCh_u.bdUpdateDone() )
        _BCh_u.bdUpdate(_mesh, _feBd_u, _dof_u);
    bc_manage(_C, _b, _mesh, _dof_u, _BCh_u, _feBd_u, 1.0, _time);


    //if (_BCh_u.fullEssential())
    _C.diagonalize(nDimensions*_dim_u, 1.0, _b, pexact(_mesh.point(1).x(),
                                             _mesh.point(1).y() ,
                                             _mesh.point(1).z()) );
    //  _C.diagonalize_row(nDimensions*_dim_u, 1.0);
    //_b[nDimensions*_dim_u]= pexact(_mesh.point(1).x(),
    //                     _mesh.point(1).y(),
    //                     _mesh.point(1).z());

    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    _solver.setMatrix(_C);



    // ---------------
    // C * V = F
    // ---------------

    for (UInt i=0; i<nDimensions*_dim_u; ++i)
        _x[i]=_u[i];

    _diff = _x;

    std::cout << "  o-  Solving system...  " << std::flush;
    chrono.start();
    _solver.solve(_x, _b);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    for (UInt i=0; i<nDimensions*_dim_u; ++i) {
        _u[i]=_x[i];
        //   std::cout << i+1 << ": " << _x[i] << std::endl;
    }

    for (UInt i=0; i<_dim_u; ++i) {
        _p[i]=_x[i+nDimensions*_dim_u];
        // std::cout << i+1 << ": " << _x[i+nDimensions*_dim_u] << std::endl;
    }

    Real norm_p=0.0;
    Real norm_u=0.0;

    for(UInt i = 1; i <= _mesh.numVolumes(); i++){

        _fe_u.updateFirstDeriv(_mesh.volumeList(i));

        norm_p += elem_L2_diff_2(_p,pexact,_fe_u,_dof_u);
        norm_u += elem_L2_diff_2(_u,uexact,_fe_u,_dof_u,0.0,int(nc_u));

    }
    std::cout << std::endl;
    std::cout << " - L2 pressure error = " << sqrt(norm_p) << std::endl;
    std::cout << " - L2 velocity error = " << sqrt(norm_u) << std::endl;
    std::cout << std::endl;

}
}
#endif
