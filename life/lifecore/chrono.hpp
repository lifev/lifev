#ifndef _CHRONO_H_INCLUDE
#define _CHRONO_H_INCLUDE
#include <time.h>
class Chrono{
  clock_t _t1,_t2,_dt;
public:
  Chrono()
  {
	_dt=0;
  }
  void start(){ _t1=clock();};
  void stop(){ _t2=clock(); _dt += _t2 - _t1;};
  double diff(){ return (1.*( _t2 - _t1))/CLOCKS_PER_SEC;};
  double diff_cumul(){ return (1.*_dt/CLOCKS_PER_SEC);};
};
#endif



// $Id: chrono.hpp,v 1.1 2004-02-08 09:09:24 prudhomm Exp $
