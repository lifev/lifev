#ifndef HH_MARKERS_HH_
#define HH_MARKERS_HH_
#include "markers_base.hpp"  
/*! \file markers.h
  \brief A Simple implementation of Markers
  \author Luca Formaggia
  \Version $Revision: 1.2 $
  
  This is the simplest implementation of the markers, which just adopts the
  basis marker classes defined in marker_base.h.

  Specialised markers can be implemented using this "template" as
  reference.
*/

//! The simples MArkerCommon: uses all defaults
typedef MarkerCommon<MarkerTraits_Base> DefMarkerCommon;

//!Expose EntityFlag
typedef MarkerTraits_Base::EntityFlag EntityFlag;

//!Expose NULLFLAG
static const EntityFlag NULLFLAG=MarkerTraits_Base::NULLFLAG;
#endif

