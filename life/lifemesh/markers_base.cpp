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

