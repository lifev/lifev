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
#include "vecUnknown.hpp"

//! the case of VectorBlock type
PhysVectUnknown<VectorBlock>::
PhysVectUnknown(UInt const Ndof):_vec(Ndof,nDimensions),_size(Ndof){}

  //! the case of VectorBlock type
GenericVecHdl<VectorBlock>::
GenericVecHdl(const ScalUnknown<VectorBlock> &RhScU1,
	      const ScalUnknown<VectorBlock> &RhScU2):
  _vec(RhScU1.size()+RhScU2.size()),_size(RhScU1.size()+RhScU2.size()),
  _nbcomp(RhScU1.nbcomp()+RhScU2.nbcomp())
{
  for (UInt i=0; i<RhScU1.size(); ++i)
    _vec.numBlock(i)= RhScU1.vec().numBlock(i);
  for (UInt i=0; i<RhScU2.size(); ++i)
    _vec.numBlock(i+RhScU1.size())= RhScU2.vec().numBlock(i);
}

//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::
GenericVecHdl(const PhysVectUnknown<VectorBlock> &RhPhVU,
	      const ScalUnknown<VectorBlock> &RhScU):
  _vec(RhPhVU.size()+RhScU.size()),_size(RhPhVU.size()+RhScU.size()),
  _nbcomp(RhPhVU.nbcomp()+RhScU.nbcomp())
{
  for (UInt i=0; i<RhPhVU.size(); ++i)
    _vec.numBlock(i)=RhPhVU.vec().numBlock(i);
  for (UInt i=0; i<RhScU.size(); ++i)
    _vec.numBlock(RhPhVU.size()+i)=RhScU.vec().numBlock(i);
}

//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::
GenericVecHdl(const ScalUnknown<VectorBlock> &RhScU,
	      const PhysVectUnknown<VectorBlock> &RhPhVU):
  _vec(RhPhVU.size()+RhScU.size()),_size(RhPhVU.size()+RhScU.size()),
  _nbcomp(RhPhVU.nbcomp()+RhScU.nbcomp())
{
  for (UInt i=0; i<RhScU.size(); ++i)
    _vec.numBlock(i)=RhScU.vec().numBlock(i);
  for (UInt i=0; i<RhPhVU.size(); ++i)
    _vec.numBlock(RhScU.size()+i)=RhPhVU.vec().numBlock(i);
}

//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::
GenericVecHdl(const PhysVectUnknown<VectorBlock> &RhPhVU,
	      const GenericVecHdl<VectorBlock> &RhGenVec):
  _vec(RhPhVU.size()+RhGenVec.size()),
  _size(RhPhVU.size()+RhGenVec.size()),
  _nbcomp(RhPhVU.nbcomp()+RhGenVec.nbcomp())
{
  for (UInt i=0;i<RhPhVU.size();++i)
    _vec.numBlock(i)=RhPhVU.vec().numBlock(i);
  for (UInt i=0;i<RhGenVec.size();++i)
    _vec.numBlock(RhPhVU.size()+i)=RhGenVec.vec().numBlock(i);
}
//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::
GenericVecHdl(const GenericVecHdl<VectorBlock> &RhGenVec,
	      const PhysVectUnknown<VectorBlock> &RhPhVU):
  _vec(RhPhVU.size()+RhGenVec.size()),
  _size(RhPhVU.size()+RhGenVec.size()),
  _nbcomp(RhPhVU.nbcomp()+RhGenVec.nbcomp())
{
  for (UInt i=0;i<RhGenVec.size();++i)
    _vec.numBlock(i)=RhGenVec.vec().numBlock(i);
  for (UInt i=0;i<RhPhVU.size();++i)
    _vec.numBlock(i+RhGenVec.size())=RhPhVU.vec().numBlock(i);
}
