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
#ifndef ELEMVEC_H
#define ELEMVEC_H

#include <vector>

#include "lifeV.hpp"
#include <iomanip>
#include "tab.hpp"

namespace LifeV
{
class ElemVec
    :
    public Tab1d
{
    int _nBlockRow; // number of block rows
    std::vector<int> _nRow; // _nRow[i]=nb of rows in the i-th block row
    std::vector<int> _firstRow;//_firstRow[i]=index of first row of i-th block row
public:
    typedef Tab1d super;

    ElemVec(int nNode1,int nbr1);
    ElemVec(int nNode1,int nbr1,
	    int nNode2,int nbr2);
    ElemVec(int nNode1,int nbr1,
	    int nNode2,int nbr2,
	    int nNode3,int nbr3);

    ElemVec& operator=( super const& __v )
	{
	    if ( this == &__v )
		return *this;
	    super& __super = (super&)*this;
	    __super = super(*this);
	    return *this;
	}

    Tab1d& vec(){return *this;};
    const Tab1d& vec() const {return *this;};
    int nBlockRow()const{return _nBlockRow;}
    Tab1dView block(int i)
	{
	    return (*this)(SubArray(_nRow[i],_firstRow[i]));
	}
    //  inline void zero(){_vec=Tab1d(_nBlockRow*_vec.N(),0.0);};
    inline void zero()
	{
	    super& __super = (super&)*this;
	    __super = super(this->N(),0.0);
	}
    void showMe(std::ostream& c=std::cout);
};
}


#endif

