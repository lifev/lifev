/* parabolicSolver.hpp

ParabolicSolver class
This class solves parabolic equations, that
are as Quarteroni says in his book, equations
that really I've a bit generalized to this form
nu(t,x,y,z,u)\frac{\partial u}{\partial t}-\nabla\cdot(mu(t,x,y,z,u)\nabla u)+
sigma(t,x,y,z,t)u=f(t,x,y,z,u)
with boundary conditions depending on u too, see testsuite/lifesolver

This works by iteration, so for big problems it can take a while.
The convergence criteria can be another one if you want, I've chosed
a simple |u^{last}-u^{previous}|/|u^last| < epsilon
where you enter epsilon in the call to automaticSolver.
*/


#include <life/lifesolver/timeSolver.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefilters/readMesh3D.hpp>    //if someone moves saveNetgenSolution from here, plz update
#include <stdio.h>    //for NULL



namespace LifeV
{


/* solves common equations of the form:
   \frac{\partial u}{\partial t} + Lu=f
   where
   L=\mu\Delta+b\dot \Nabla+\sigma
   (Nabla not still implemented with variable coeff(x,y,z,t,u))
   with boundary conditions depending on solution,
   like in boundary convection in heat transfer problems
*/
template <typename MESHTYPE>
class ParabolicSolver
{
public:
    typedef Real (*function_type)(Real t,Real x,Real y,Real z,Real u);
    typedef ScalUnknown<Vector>    vector_type;
    typedef MSRMatr<Real> matrix_type;
    typedef MESHTYPE mesh_type;

    explicit ParabolicSolver(function_type _nu,
                             function_type _mu, function_type _sigma,
                             const mesh_type& _mesh,CurrentFE& _fe,
                             CurrentBdFE& _feBd,
                             const Dof& _dof,
                             const BCHandler& _BCh,
                             function_type _m_fct,
                             vector_type _U0,
                             Real _timeStep=1.0,UInt _nSteps=4);
    explicit ParabolicSolver(function_type _nu,
                             function_type _mu,
                             const mesh_type& _mesh,CurrentFE& _fe,
                             CurrentBdFE& _feBd,
                             const Dof& _dof,
                             const BCHandler& _BCh,
                             function_type _m_fct,
                             vector_type _U0,
                             Real _timeStep=1.0,UInt _nSteps=4);
    ParabolicSolver();
    ~ParabolicSolver();

    void assembleM();
    void assembleAt();
    void assembleBt();
    void assemble();
    void tuneSolver();
    vector_type solve();
    UInt timeInc();
    UInt timeDec();
    void reset();
    UInt getStep();

    void automaticSolver(Real epsilon=0.1,UInt maxiter=10,UInt steps=0);
  
    vector_type getU(UInt _nStep);
    void setUin(const vector_type& _Uin);
    void saveSolution(std::string dir="sol");
private:
    void init();


private:
    UInt _m_dim;
    function_type _m_nu,_m_mu,_m_sigma,_m_fct;
    const mesh_type& _m_mesh;
    /* warning: this objects haven't a good copy constructor,
       so I keep references to them, but please don't destroy 
       them from the outside while still used by ParabolicSolver */
    CurrentFE& _m_fe;
    CurrentBdFE& _m_feBd;
    const BCHandler& _m_BCh;
    const Dof& _m_dof;

    Real _m_time0;
    Real _m_timeStep;
    Real _m_timeCur;
    UInt _m_nSteps;
    UInt _m_nStepCur;
    vector_type    _m_Uin;
    TimeSolver    _m_solver;

    MSRPatt _m_msrPatt;
    matrix_type _m_Mt;
    matrix_type _m_Mt_1;
    matrix_type _m_At_1;
    matrix_type _m_At;
    vector_type _m_bt_1;
    vector_type _m_bt;
public:
    std::vector<vector_type> _m_Ut;    //solution

};

template <typename MESHTYPE>
ParabolicSolver<MESHTYPE>::ParabolicSolver(){
    std::cout<<"use another constructor"<<std::endl;
    ABORT();
}
template <typename MESHTYPE>
ParabolicSolver<MESHTYPE>::~ParabolicSolver(){}

template <typename MESHTYPE>
ParabolicSolver<MESHTYPE>::ParabolicSolver(function_type _nu,
                                           function_type _mu,
                                           function_type _sigma,
                                           const mesh_type& _mesh,CurrentFE& _fe,
                                           CurrentBdFE& _feBd,
                                           const Dof& _dof,
                                           const BCHandler& _BCh,
                                           function_type _fct,
                                           vector_type _U0,
                                           Real _timeStep,UInt _nSteps):
    _m_dim(_dof.numTotalDof()),
    _m_nu(_nu),
    _m_mu(_mu),_m_sigma(_sigma),
    _m_fct(_fct),
    _m_mesh(_mesh),
    _m_fe(_fe),
    _m_feBd(_feBd),
    _m_BCh(_BCh),
    _m_dof(_dof),
    _m_time0(0),
    _m_timeStep(_timeStep),
    _m_timeCur(0),
    _m_nSteps(_nSteps),
    _m_nStepCur(0),
    _m_Uin(_U0),
    _m_solver(_m_dim),
    _m_Mt(),_m_Mt_1(),_m_At_1(),_m_At(),
    _m_bt_1(_m_dim),_m_bt(_m_dim),
    _m_Ut(_nSteps,_U0) //ZeroVector(_m_dim)),
{
    init();
}
template <typename MESHTYPE>
ParabolicSolver<MESHTYPE>::ParabolicSolver(function_type _nu,
                                           function_type _mu, 
                                           const mesh_type& _mesh,CurrentFE& _fe,
                                           CurrentBdFE& _feBd,
                                           const Dof& _dof,
                                           const BCHandler& _BCh,
                                           function_type _fct,
                                           vector_type _U0,
                                           Real _timeStep,UInt _nSteps):
    _m_dim(_dof.numTotalDof()),
    _m_nu(_nu),
    _m_mu(_mu),_m_sigma(NULL),
    _m_fct(_fct),
    _m_mesh(_mesh),
    _m_fe(_fe),
    _m_feBd(_feBd),
    _m_BCh(_BCh),
    _m_dof(_dof),
    _m_time0(0),
    _m_timeStep(_timeStep),
    _m_timeCur(0),
    _m_nSteps(_nSteps),
    _m_nStepCur(0),
    _m_Uin(_U0),
    _m_solver(_m_dim),
    _m_Mt(),_m_Mt_1(),_m_At_1(),_m_At(),
    _m_bt_1(_m_dim),_m_bt(_m_dim),
    _m_Ut(_nSteps,_U0)    //ZeroVector(_m_dim)),
{
    init();
}


template <typename MESHTYPE>
void
ParabolicSolver<MESHTYPE>::init()
{
    _m_dim=_m_dof.numTotalDof();
    _m_Ut[0]=_m_Uin;
    for(UInt i=1;i<_m_nSteps;i++){
        _m_Ut[i]=ZeroVector(_m_dim);    //do you like more initial value as default?
    }
    _m_bt_1=ZeroVector(_m_dim);
    _m_bt=ZeroVector(_m_dim);
    _m_msrPatt=MSRPatt(_m_dof);
    _m_Mt=matrix_type(_m_msrPatt);
    _m_Mt_1=matrix_type(_m_msrPatt);
    _m_At_1=matrix_type(_m_msrPatt);
    _m_At=matrix_type(_m_msrPatt);
}
template <typename MESHTYPE>
UInt
ParabolicSolver<MESHTYPE>::timeInc()
{
    _m_At_1=_m_At;
    _m_bt_1=_m_bt;
    _m_Mt_1=_m_Mt;
    _m_timeCur+=_m_timeStep;
    return ++_m_nStepCur;
}
/* you cannot solve after this call */
template <typename MESHTYPE>
UInt
ParabolicSolver<MESHTYPE>::timeDec()
{
    _m_At=_m_At_1;
    _m_bt=_m_bt_1;
    _m_Mt=_m_Mt_1;
    _m_timeCur-=_m_timeStep;
    return --_m_nStepCur;
}

template <typename MESHTYPE>
void
ParabolicSolver<MESHTYPE>::reset()
{
    _m_timeCur=_m_time0;
    _m_nStepCur=0;
}
template <typename MESHTYPE>
UInt
ParabolicSolver<MESHTYPE>::getStep()
{
    return _m_nStepCur;
}
/* you must timeInc after automaticSolver if you want to solve again */
template <typename MESHTYPE>
void
ParabolicSolver<MESHTYPE>::automaticSolver(Real epsilon,UInt maxiter,UInt steps)
{
    Real norm_u2,norm_diff;
    vector_type u1(_m_dim),u2(_m_dim);
    if(steps==0)steps=_m_nSteps;
    UInt iter,laststep=Min(_m_nSteps,_m_nStepCur+steps);
    while(_m_nStepCur<laststep)
        {
            iter=0;
            do
                {
                    u1=solve();
                    u2=solve();
                    norm_u2=sqrt(dot(u2,u2));
                    norm_diff=sqrt(dot((u1-u2),(u1-u2)));
                    std::cout<<"automaticSolver: iter "<<_m_nStepCur<<"->"<<iter<<" norm_u2="
                             <<norm_u2<<" norm_diff="<<norm_diff<<std::endl; 
                } while(norm_diff>epsilon*norm_u2 && ++iter<maxiter);
            if(iter==maxiter)std::cout<<"ParabolicSolver: warning "<<iter<<" iterations"
                                 " not enought"<<std::endl;
            timeInc();
        }
    timeDec();
}

template <typename MESHTYPE>
void
ParabolicSolver<MESHTYPE>::assembleM()
{
    ElemMat elmat(_m_fe.nbNode,1,1);
    for(UInt i = 1; i<=_m_mesh.numVolumes(); i++){
        _m_fe.updateFirstDerivQuadPt(_m_mesh.volumeList(i));
        elmat.zero();
        mass(_m_nu,elmat,_m_fe,_m_dof,_m_Uin,_m_timeCur);
        //mass(1.,elmat,_m_fe);
        assemb_mat(_m_Mt,elmat,_m_fe,_m_dof);
    }
    bcManageMtimeUDep( _m_Mt, _m_dof, _m_BCh, 1.0);
}
template <typename MESHTYPE>
void
ParabolicSolver<MESHTYPE>::assembleAt()
{
    ElemMat elmat(_m_fe.nbNode,1,1);
    if(_m_mu!=NULL){
        for(UInt i = 1; i<=_m_mesh.numVolumes(); i++){
            _m_fe.updateFirstDerivQuadPt(_m_mesh.volumeList(i));
            elmat.zero();
            stiff(_m_mu,elmat,_m_fe,_m_dof,_m_Uin,_m_timeCur);
            assemb_mat(_m_At,elmat,_m_fe,_m_dof);
        }
    }
    if(_m_sigma!=NULL){
        for(UInt i = 1; i<=_m_mesh.numVolumes(); i++){
            _m_fe.updateFirstDerivQuadPt(_m_mesh.volumeList(i));
            elmat.zero();
            mass(_m_sigma,elmat,_m_fe,_m_dof,_m_Uin,_m_timeCur);
            assemb_mat(_m_At,elmat,_m_fe,_m_dof);
        }
    }//todo: add grad here and implement it too as UDep func
}
template <typename MESHTYPE>
void
ParabolicSolver<MESHTYPE>::assembleBt()
{
    ElemVec elvec(_m_fe.nbNode,1);
    for(UInt i = 1; i<=_m_mesh.numVolumes(); i++){
        _m_fe.updateFirstDerivQuadPt(_m_mesh.volumeList(i));
        elvec.zero();
        source(_m_fct,elvec,_m_fe,_m_dof,_m_Uin,_m_timeCur);
        assemb_vec(_m_bt,elvec,_m_fe,_m_dof);
    }
}


template <typename MESHTYPE>
void
ParabolicSolver<MESHTYPE>::assemble()
{
    if(_m_nStepCur==0){    //for time0 there is no iteration as U0 is known
        assembleM();
        assembleAt();
        assembleBt();
        bcManage(_m_mu,_m_At,_m_bt,_m_mesh,_m_dof,_m_BCh,_m_feBd,1.0,_m_timeCur,_m_Uin);
        timeInc();
    }
    assembleM();    
    assembleAt();
    assembleBt();
    bcManage(_m_mu,_m_At,_m_bt,_m_mesh,_m_dof,_m_BCh,_m_feBd,1.0,_m_timeCur,_m_Uin);
}
template <typename MESHTYPE>
void
ParabolicSolver<MESHTYPE>::tuneSolver()
{
    _m_solver.setM_k(_m_Mt);
    _m_solver.setM_k_1(_m_Mt_1);
    _m_solver.setA_k(_m_At);
    _m_solver.setA_k_1(_m_At_1);
    _m_solver.setb_k(_m_bt);
    _m_solver.setb_k_1(_m_bt_1);
    _m_solver.settheta(0.5);
    _m_solver.setu_k_1(_m_Uin);
    _m_solver.setdt(_m_timeStep);
}
template <typename MESHTYPE>
typename ParabolicSolver<MESHTYPE>::vector_type
ParabolicSolver<MESHTYPE>::solve()
{
    assemble();
    tuneSolver();
    return _m_Ut[_m_nStepCur]=_m_Uin=_m_solver.solve();
}

template <typename MESHTYPE>
typename ParabolicSolver<MESHTYPE>::vector_type 
ParabolicSolver<MESHTYPE>::getU(UInt _nStep)
{
    if(_nStep>_m_nStepCur)
        ERROR_MSG("ParabolicSolver::getU: out of band index\n");
    return _m_Ut[_nStep];
}

template <typename MESHTYPE>
void
ParabolicSolver<MESHTYPE>::setUin(const vector_type& _Uin)
{
    _m_Uin=_Uin;
}


/* I can think only of netgen format as netgen is
   the only free program that I know
*/
template <typename MESHTYPE>
void
ParabolicSolver<MESHTYPE>::saveSolution(std::string dir)
{
    UInt i;
    for(i=0;i<=_m_nStepCur;i++)
        {
            saveNetgenSolution(dir+"/u"+i+".sol",_m_Ut[i],std::string("u")+i);
        }
}

}
