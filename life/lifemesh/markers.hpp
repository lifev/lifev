/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003 LifeV Team
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.
  
  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef HH_MARKERS_HH_
#define HH_MARKERS_HH_
#include "markers_base.hpp"  
/*! \file markers.h
  \brief A Simple implementation of Markers
  \author Luca Formaggia
  \Version $Revision: 1.3 $
  
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

