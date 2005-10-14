namespace LifeV
{

  Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    return 0.;
  }
  Real u0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    return 0.;
  }
  
  Real p(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    return 0;
  }
  
  Real fZero(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    return 0.0;
  }
   
  Real g(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    if ( t < 0.005 )
      return 1.3320e4;
    return 0.0;
  } 

  Real g_anev(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    if ( t < 0.005 )
      return 1.0e4;
    return 0.0;
  } 

  Real d0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {  
    return 0.;
  }
  Real w0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {  
    return 0.;
  }



  Real bdDisp(const Real& t, const Real& x, const Real& y, const Real& z,
              const ID& i) {
    
    Real R=sqrt(x*x+y*y);
    Real omega = acos(-1.0)*400;
    
    switch(i){
    case 1:
      return sin(t*omega)*0.02*x/R;
      break;
    case 2:
      return sin(t*omega)*0.02*y/R;
    case 3:
      return 0.0;
      break;
      
    }
    return 0.0;
  }


  
  Real qND(const Real& t)
  {
    if (t<1.0/7.0) {
      return
	10.9506162830641998*t-193.580197870147003*t*t*t;
    }
    else if ((t>=1.0/7.0)&&(t<2.0/7.0)) {
      return
	1.12874750944692-.901232566128460033*t-82.9629419443489980*(t-1./7.)*(t-1./7.)+281.900989350737973*(t-1./7.)*(t-1./7.)*(t-1./7.);
    }
    else if ((t>=2.0/7.0)&&(t<3.0/7.0)) {
      return
	2.09876743387152-7.34568601855032987*t+37.8517677773960983*(t-2./7.)*(t-2./7.)-21.3007595328068007*(t-2./7.)*(t-2./7.)*(t-2./7.);
    }
    else if ((t>=3.0/7.0)&&(t<4.0/7.0)) {
      return
	-1.26684713156991+2.16497664032978986*t+28.7228708347645991*(t-3./7.)*(t-3./7.)-139.416951219511986*(t-3./7.)*(t-3./7.)*(t-3./7.);
    }
    else if ((t>=4.0/7.0)&&(t<5.0/7.0)) {
      return
	-.899016832703540+1.83577945723119007*t-31.0272511164547993*(t-4./7.)*(t-4./7.)+94.3095644108550033*(t-4./7.)*(t-4./7.)*(t-4./7.);
    }
    else if ((t>=5.0/7.0)&&(t<6.0/7.0)) {
     return
       .950496049467536-1.25509446925455004*t+9.39113363105458988*(t-5./7.)*(t-5./7.)-22.7603064239091992*(t-5./7.)*(t-5./7.)*(t-5./7.);
    }
    else if ((t>=6.0/7.0)&&(t<1.)) {
      return
     -0.296557883888689e-1+0.345984197870137026e-1*t-.363283407763643984*(t-6./7.)*(t-6./7.)+.847661284781833002*(t-6./7.)*(t-6./7.)*(t-6./7.);
    }
  }


Real g_caro(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  Real T=0.84;                  //Period (s)
  Real pi=acos(-1.0);
  Real t_loc=(t+0.33)/T-int((t+0.33)/T);

  return (qND(t_loc)+0.4)/1.4*1.3e4;
}




  Real u_in(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    Real R=0.85; 		//Tube radius (cm)
    Real T=0.84;			//Period (s)
    Real pi=acos(-1.0);
    Real t_loc=t/T-int(t/T);
    
    switch(i) {
    case 1:
    case 2:
      return 0.0;
      break;
    case 3:  	// profile in tanh
    //			std::cout << t_loc << "  " <<qND(t_loc)  << std::endl;
      return 
	126.18*qND(t_loc)/pi/(R*R)*20/18*(1-(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R);
    break;
    }
  return 0.0;
  
  }




}
