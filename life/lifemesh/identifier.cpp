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
