// Main programm for coupled calculation of mass transport in two domains
// (the arterial lumen and the arterial wall)

// author:M. Prosi                                august/04
#include <life/lifecore/life.hpp>
#include <life/lifesolver/convDiffReactSolverPC.hpp>
#include <life/lifefem/dofInterface3Dto3D.hpp>
#include <life/lifecore/chrono.hpp>
#include "ud_functions.hpp"
#include <life/lifecore/GetPot.hpp>
#include <life/lifefilters/ensight7Writer.hpp>

int main(int argc, char** argv)
{
 using namespace LifeV;
 using namespace std;
 // ********** Reading from data files *****************************************

  GetPot command_line(argc,argv);
  const char* data_file_name1 = command_line.follow("data-lumen", 2, "-f","--file");
  GetPot data_file1(data_file_name1);
  const char* data_file_name2 = command_line.follow("data-wall", 2, "-f","--file");
  GetPot data_file2(data_file_name2);

// ***** Parameters describing the transport properties of the wall layers *****

  Real u_filt = 1.78E-6;
  Real U_0 = 15.7;
  Real L_0 = 0.67;
  Real Klag = 0.1164;
  Real epsilon = 0.15;

  Real P = 2.0e-8;
  Real kappa = 0.6176;
  Real s = 2.134e-3;

// ********** Boundary conditions definitions for mass transport **************
  BCFunctionBase c_inflow(c1);
  BCFunctionMixte c_wall(alpha,beta);   // Permeability boundary condition
  BCHandler BCh_cl(3);

  BCFunctionBase c_adv(cadv);
  BCHandler BCh_cw(2);

// ****** Concentration class lumen: cdrlumen **********************************
  ConvDiffReactSolverPC< RegionMesh3D<LinearHexa> >
     cdrlumen(data_file1, feHexaQ1, quadRuleHexa8pt, quadRuleQuad4pt, BCh_cl);
  cdrlumen.showMe();

// ******** Concentration class wall: cdrwall **********************************
  ConvDiffReactSolverPC< RegionMesh3D<LinearHexa> >
     cdrwall(data_file2, feHexaQ1, quadRuleHexa8pt, quadRuleQuad4pt, BCh_cw);
  cdrwall.showMe();

  UInt dim_cl = cdrlumen.mesh().numVertices();
  UInt dim_cw = cdrwall.mesh().numVertices();
  Real tol = 1.0e-6;

  Vector cw_interface(dim_cw);
  Vector cl_interface(dim_cl);

  Vector cw_mixte(dim_cw);
  Vector cl_mixte(dim_cl);

  Vector Pw_var(dim_cw);
  Vector Pl_var(dim_cl);

// Produce a "hole" in the permeability of the testgeometry (10 times higher permeability)

  for(UInt ii=0; ii < dim_cw; ii++){
     Pw_var(ii)=P;}
  for(UInt ii=3; ii < 4; ii++){
     for(UInt jj=749; jj < 768; jj++){
    Pw_var(jj+ii*777)=10*P;}}

  for(UInt ii=0; ii < dim_cl; ii++){
     Pl_var(ii)=P;}
  for(UInt ii=5; ii < 6; ii++){
     Pl_var(44+ii*992)=10*P;
     for(UInt jj=717; jj < 725; jj++){
    Pl_var(jj+ii*992)=10*P;}
     Pl_var(45+ii*992)=10*P;
     for(UInt jj=806; jj < 814; jj++){
    Pl_var(jj+ii*992)=10*P;}
     Pl_var(46+ii*992)=10*P;}

//========================================================================================
//  DATA INTERFACING BETWEEN BOTH SOLVERS
//========================================================================================
//
// Passing data from the lumen to the wall
  boost::shared_ptr<DofInterface3Dto3D> dofLumentoWall( new DofInterface3Dto3D(feHexaQ1,
                                                                               cdrwall.cDof(),
                                                                               feHexaQ1,
                                                                               cdrlumen.cDof()) );
  dofLumentoWall->update(cdrwall.mesh(), 3, cdrlumen.mesh(), 6, tol);
  BCVectorInterface cl_coupling(cl_interface, dim_cl, dofLumentoWall );
//  cl_coupling.setMixteCoef((-1.0/epsilon)*(u_filt*(((s*kappa)/2)-Klag)-P));
  for(UInt ii=0; ii < dim_cl; ii++){
     cl_mixte(ii)=(-1.0/epsilon)*(u_filt*(((s*kappa)/2)-Klag)-Pl_var(ii));}
  cl_coupling.setMixteVec( cl_mixte );

// Passing data from the wall to the lumen
  boost::shared_ptr<DofInterface3Dto3D> dofWalltoLumen( new DofInterface3Dto3D(feHexaQ1,
                                                                               cdrlumen.cDof(),
                                                                               feHexaQ1,
                                                                               cdrwall.cDof()) );
  dofWalltoLumen->update(cdrlumen.mesh(), 6, cdrwall.mesh(), 3, tol);
  BCVectorInterface cw_coupling(cw_interface, dim_cw, dofWalltoLumen );
//  cw_coupling.setMixteCoef(P-u_filt*(1.0-((s*kappa)/2)));
  for(UInt ii=0; ii < dim_cw; ii++){
     cw_mixte(ii)=Pw_var(ii)-u_filt*(1.0-((s*kappa)/2));}
  cw_coupling.setMixteVec( cw_mixte );

  BCh_cl.addBC("coupling part",   6, Mixte, Scalar, cw_coupling);// Interface to wall part
  BCh_cl.addBC("C-Endothel",   2, Mixte, Scalar, c_wall);        // Permeability boundary condition
  BCh_cl.addBC("C-InFlow", 1, Essential,  Scalar, c_inflow);     // Concentration profile

  BCh_cw.addBC("C-Endothel",   3, Mixte, Scalar, cl_coupling);   // Interface to lumen
  BCh_cw.addBC("C-Adventitia", 1, Essential,  Scalar, c_adv);    // Concentration at the adventitia

// *********** Initialization **************************************************
// this are properties of the fluid classes, they are not included here (analytical
// solution of the flow) -> hardcoded timestep, endtime and initialisation time
  Real dt = 1.0;
  Real T  = 1.0;
  Real startT = 0.0;

  if(startT > 0.0){
     ostringstream indexin;
     std::string cinname;
     indexin << (startT*100);
     cinname = "concentration-lumen.res"+indexin.str();
     cdrlumen.initialize(cinname);
     cinname = "concentration-wall.res"+indexin.str();
     cdrwall.initialize(cinname);
  }
  else{
     std::cout << "initialize concentrations with cl_0 and cw_0" << std::endl;
     cdrlumen.initialize(cl_0,0.0,dt);
     cdrwall.initialize(cw_0,0.0,dt);
  }

// ********* Analytical solution of the velocity in the lumen ******************
// straight tube, Reynoldsnumber in the lumen Re=300
  for(UInt i=0; i < dim_cl; i++){
     Real xp = cdrlumen.mesh().point(i+1).x();
     Real yp = cdrlumen.mesh().point(i+1).y();
     Real zp = cdrlumen.mesh().point(i+1).z();
     cdrlumen.u_c()(i)=(2.0/L_0)*u_filt*(2.0-(4.0/(L_0*L_0))*(xp*xp+yp*yp))*xp;
     cdrlumen.u_c()(i+cdrlumen.mesh().numVertices())=(2.0/L_0)*u_filt*(2.0-(4.0/(L_0*L_0))*(xp*xp+yp*yp))*yp;
     cdrlumen.u_c()(i+2*cdrlumen.mesh().numVertices())=2.0*U_0*(1.0-(4.0/(L_0*L_0))*(xp*xp+yp*yp))*(1.0-(4.0*u_filt/(U_0*L_0))*zp);}
// ********* Analytical solution of the velocity in the wall *******************
  for(UInt i=0; i < dim_cw; i++){
     Real xp = cdrwall.mesh().point(i+1).x();
     Real yp = cdrwall.mesh().point(i+1).y();
     cdrwall.u_c()(i)=(Klag/epsilon)*u_filt*(L_0/2.0)*xp/(xp*xp+yp*yp);
     cdrwall.u_c()(i+cdrwall.mesh().numVertices())=(Klag/epsilon)*u_filt*(L_0/2.0)*yp/(xp*xp+yp*yp);
     cdrwall.u_c()(i+2*cdrwall.mesh().numVertices())=0.0;}

  Vector cw_old(dim_cw);
  Vector cl_old(dim_cl);
  Vector cw_delta(dim_cw);
  Vector cl_delta(dim_cl);

  cw_old=cdrwall.c();
  cl_old=cdrlumen.c();

// ************ Temporal loop **************************************************
  for (Real time=startT+dt ; time <= T; time+=dt) {

     cdrlumen.timeAdvance(fc,time);
     cdrwall.timeAdvance(fc,time);

// ****** Coupling of concentration solvers between the wall and the lumen *****
// ****** inner iterations (till now fix number of iterations)

     for(Int iter=0; iter < 4; iter ++){

    std::cout << iter << std::endl;

    for(UInt ii=0; ii < dim_cw; ii++){
       cw_interface(ii)=cdrwall.c()(ii)*(1/epsilon)*(Pw_var(ii)-u_filt*((s*kappa)/2));}

    cdrlumen.iterate(time);

    for(UInt ii=0; ii < dim_cl; ii++){
       cl_interface(ii)=cdrlumen.c()(ii)*(Pl_var(ii)+u_filt*((s*kappa)/2));}

    cdrwall.iterate(time);

    cw_delta=cdrwall.c()-cw_old;
    cl_delta=cdrlumen.c()-cl_old;

    std::cout << "difference in the wall: " << norm_inf(cw_delta) << std::endl;
    std::cout << "difference in the lumen: " << norm_inf(cl_delta) << std::endl;

    cw_old=cdrwall.c();
    cl_old=cdrlumen.c();
     }


    cout <<  sqrt(cdrlumen.mesh().point(4964).x()*cdrlumen.mesh().point(4964).x()+
                     cdrlumen.mesh().point(4964).y()*cdrlumen.mesh().point(4964).y())
            << ", " << cdrlumen.c()(4964) << endl;
    cout <<  sqrt(cdrlumen.mesh().point(4971).x()*cdrlumen.mesh().point(4971).x()+
                     cdrlumen.mesh().point(4971).y()*cdrlumen.mesh().point(4971).y())
            << ", " << cdrlumen.c()(4971) << endl;
    cout <<  sqrt(cdrlumen.mesh().point(4978).x()*cdrlumen.mesh().point(4978).x()+
                     cdrlumen.mesh().point(4978).y()*cdrlumen.mesh().point(4978).y())
            << ", " << cdrlumen.c()(4978) << endl;
    cout <<  sqrt(cdrlumen.mesh().point(4985).x()*cdrlumen.mesh().point(4985).x()+
                     cdrlumen.mesh().point(4985).y()*cdrlumen.mesh().point(4985).y())
            << ", " << cdrlumen.c()(4985) << endl;
    cout <<  sqrt(cdrlumen.mesh().point(4992).x()*cdrlumen.mesh().point(4992).x()+
                     cdrlumen.mesh().point(4992).y()*cdrlumen.mesh().point(4992).y())
            << ", " << cdrlumen.c()(4992) << endl;
    for(int i=0; i<4;i++){
    cout <<  sqrt(cdrlumen.mesh().point(5234+i*9).x()*cdrlumen.mesh().point(5234+i*9).x()+
                     cdrlumen.mesh().point(5234+i*9).y()*cdrlumen.mesh().point(5234+i*9).y())
     << ", " << cdrlumen.c()(5234+i*9) << endl;}
    cout <<  sqrt(cdrlumen.mesh().point(4999).x()*cdrlumen.mesh().point(4999).x()+
                     cdrlumen.mesh().point(4999).y()*cdrlumen.mesh().point(4999).y())
            << ", " << cdrlumen.c()(4999) << endl;
    for(int i=0; i<9;i++){
    cout <<  sqrt(cdrlumen.mesh().point(5597+i*9).x()*cdrlumen.mesh().point(5597+i*9).x()+
                     cdrlumen.mesh().point(5597+i*9).y()*cdrlumen.mesh().point(5597+i*9).y())
     << ", " << cdrlumen.c()(5597+i*9) << endl;}
    cout <<  sqrt(cdrlumen.mesh().point(5006).x()*cdrlumen.mesh().point(5006).x()+
                     cdrlumen.mesh().point(5006).y()*cdrlumen.mesh().point(5006).y())
            << ", " << cdrlumen.c()(5006) << endl;


    for(int i=0; i<21;i++){
       cout <<  sqrt(cdrwall.mesh().point(19+777*3+37*i).x()*cdrwall.mesh().point(19+777*3+37*i).x()+
                     cdrwall.mesh().point(19+777*3+37*i).y()*cdrwall.mesh().point(19+777*3+37*i).y())
            << ", " << cdrwall.c()(19+777*3+37*i) << endl;}



    outensight7Mesh3D(cdrwall.mesh(), cdrwall.u_c(), cdrwall.c(),time);

// ************* saving result on file *****************************************
    ostringstream indexout;
    indexout << (time*100);
    std::string voutname;
    voutname = "concentration-lumen.res"+indexout.str();
    std::fstream Resfile_cl(voutname.c_str(),ios::out | ios::binary);
    Resfile_cl.write((char*)&cdrlumen.c()(1),cdrlumen.c().size()*sizeof(double));
    Resfile_cl.close();
    voutname = "concentration-wall.res"+indexout.str();
    std::fstream Resfile_cw(voutname.c_str(),ios::out | ios::binary);
    Resfile_cw.write((char*)&cdrwall.c()(1),cdrwall.c().size()*sizeof(double));
    Resfile_cw.close();

// ************* creating Ensight output file **********************************
// TODO: creatig outputfiles for two domains
  }

  return 0;
}
