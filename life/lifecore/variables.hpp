#ifndef _VARIABLES_H
#define _VARIABLES_H
#include<vector>
#include<string>
#include<typeinfo>
#include "tab.hpp"


template<typename UnkType,int size,int nblock>
class Unknown {
 
 public:
  /*  Unknown(string name):_name(name),_size(size),_nb(nblock) 
   {if (size==1) {unk=0}
    else (unk(size,nblock));
    return;};

  Unknown():_name(" "),_size(size) 
   {if (size==1) {unk=0;}
    else unk(size,nblock);
    return;};
  */

  int size() {return _size;}; //physical dimension of the unknown
  int ndof() {return _ndof;}; //degrees of freedom for the unknown
  int nb() {return _nb;};
  string name() {return _name;};

  vector<UnkType> unk; 
  void set_vect(int ndof){unk(ndof);};
  

  //bc_selection & c.;


 private:
 string _name;
 int _size;
 int _nb;
 int _ndof;
}



typedef Unknown<double,1,1> ScaUnk; // Generic Scalar unknown

typedef Unknown<ElemVec,nDimension,1> PhyUnk; // Physical Vector unknown (e.g. velocity in NS problem)

#endif

// Porzione di codice:
/*
ScaUnk p("pressure");
PhyUnk u("velocity");

*/
