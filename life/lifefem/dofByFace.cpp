#include "dofByFace.hpp"

namespace LifeV
{

void DofByFace::showMe(std::ostream& out, bool verbose) const{
  out << "--------------------------------------------------------------------------------" << std::endl;
  out << " Degree of freedom by face object " << std::endl;
  out << "--------------------------------------------------------------------------------" << std::endl;

  out << " Offset (min Dof Id) = " << _offset << std::endl;
  out << " Number of internal faces = " << numIFaces() << std::endl;
  out << " Number of local dof per face = " << _numLocalDofByFace << std::endl;

  if(verbose){
    out << "*********************************************************************************" << std::endl;
    out << " Local-to-global DOF table (DOF grouped by internal face)" << std::endl;
    out << "*********************************************************************************" << std::endl;
    out << "=================================================================================" << std::endl;
    out << "Face ID     Local DOF   Global DOF  " << std::endl;
    out << "=================================================================================" << std::endl;
    
    for(UInt i = 0; i < _numIFaces; ++i){
      for(UInt j = 0; j < _numLocalDofByFace; ++j){
	out.width(12);
	out << i + 1;
	out.width(12);
	out << j + 1;
	out.width(12);
	out << localToGlobal(i+1, j+1);
	out << " # ";
	if(j % 2 != 0) out << std::endl;
      } // for j
    } //for i
  } // if verbose
}
}
