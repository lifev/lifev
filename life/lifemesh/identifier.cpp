/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
  
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
/*!
  \file identifier.cc
  \brief Implementations for identifier.h
  \version 1.0
  \author M.A. Fernandez
  \date 07/2002

*/
#include "identifier.hpp"

Identifier_Natural::Identifier_Natural(const ID& id, const SimpleVect<ID>& bdltg):Identifier_Base(id) {
  _bdltg.reserve(bdltg.size());
  _bdltg.insert(_bdltg.end(),bdltg.begin(),bdltg.end());
}

  
Identifier_Natural::Identifier_Natural(const ID& id):Identifier_Base(id) {
}
