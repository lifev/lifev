#include "dof.hpp"

 //! Constructor
Dof::Dof(const LocalDofPattern& _fe, UInt off):fe(_fe),_offset(off),_totalDof(0),
					       _nEl(0),nlv(0),nle(0),nlf(0),_ltg()
{
  for (UInt i=0; i<5; ++i)_ncount[i]=0;
}

//! Copy constructor
Dof::Dof(const Dof & dof2 ):fe(dof2.fe),_offset(dof2._offset),
			    _totalDof(dof2._totalDof),_nEl(dof2._nEl),
			    nlv(dof2.nlv), nle(dof2.nle), nlf(dof2.nlf),
			    _ltg(dof2._ltg)
{
  if (&dof2 == this) return;
  
  //fe=dof2.fe;
  //  _offset=dof2._offset;
  //  _ltg=dof2._ltg;
  //  _totalDof=dof2._totalDof;
    for (UInt i=0; i<5; ++i)_ncount[i]=dof2._ncount[i];
}

//! Ouput
void Dof::showMe(ostream  & out, bool verbose) const{
  out<< " Degree of Freedom (Dof) Object"<<endl;
  out<< " Total Dof Stored             "<<_totalDof<<endl;
  out<< " With offset (min. Dof Id) =  "<<_offset<<endl;
  out<< " Dof's on Vertices  from "<< _ncount[0]<<" , to:"<<_ncount[1]-1<<endl;
  out<< " Dof's on Edges     from "<< _ncount[1]<<" , to:"<<_ncount[2]-1<<endl;
  out<< " Dof's on Faces     from "<< _ncount[2]<<" , to:"<<_ncount[3]-1<<endl;
  out<< " Dof's on Volumes   from "<< _ncount[3]<<" , to:"<<_ncount[4]-1<<endl;
  if (verbose){
    out<<"************************************************************"<<endl;
    out<<"           Local to Global DOF table"<<endl;
    out<<"************************************************************"<<endl;
    out<<"Element Id   Loc. N.   Global N.  #  Element Id  Loc. N. Global N. "<<endl;

      
    for (UInt i=0; i<_nEl;++i){
      for (UInt j=0; j<numLocalDof();++j)
	{
	  out.width(10);
	  out<<i+1;
	  out.width(10);
	  out<<j+1;
	  out.width(10);
	  out<<localToGlobal(i+1,j+1);
	  out<<" # ";
	  if(j % 2!=0)out<<endl;
	}
      
    }
    
  }
  
}
