#include <climits>
#include "markers_base.hpp"

//  ***********************************************************************************************************
//                                           IMPLEMENTATION
//  ***********************************************************************************************************
//MarkerTraits_Base
const MarkerTraits_Base::EntityFlag MarkerTraits_Base::NULLFLAG=LONG_MIN;

MarkerTraits_Base::EntityFlag MarkerTraits_Base::strongerFlag(EntityFlag const & a,EntityFlag const & b)
{
  if (a == NULLFLAG| b==NULLFLAG) return NULLFLAG;
  return a > b ? a : b ;
}

MarkerTraits_Base::EntityFlag MarkerTraits_Base::weakerFlag(EntityFlag const & a,EntityFlag const & b)
{
  if (a == NULLFLAG| b==NULLFLAG) return NULLFLAG;
  return a < b ? a : b ;
}

//$Id: markers_base.cpp,v 1.1 2004-02-08 09:09:22 prudhomm Exp $
//$Log: markers_base.cpp,v $
//Revision 1.1  2004-02-08 09:09:22  prudhomm
//finally added the new life libraries layout
//
//life/lifecore core library
//life/lifemesh mesh library
//life/lifefem fem library
//
//more to come in the future
//
//Revision 1.1  2002/11/15 13:03:32  forma
//Modifications related to the new marker classes and the new mesh
//checking routines.
//
