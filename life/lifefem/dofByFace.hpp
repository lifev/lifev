#ifndef _DOFBYFACE_HH
#define _DOFBYFACE_HH

#include "life.hpp"
#include "SimpleVect.hpp"
#include "localDofPattern.hpp"
#include "dof.hpp"
#include <algorithm>

/*!
  \file DofByFace.h
  \brief local-to-global table with DOF grouped by (internal) face
*/
namespace LifeV
{
/*!
  \class DofByFace
  \brief The class for local-to-global table with DOF grouped by (internal) face
  \author D. A. Di Pietro
  \date 01/2004
*/

class DofByFace{
 public:
  typedef SimpleArray<ID> Container;

  const LocalDofPattern& fe;
  /*!
    This constructors simply inspects the second part of mesh.faceList
    (from mesh.numBFaces() to mesh.numFaces(), where internal faces are
    stored) and copies the DOF of both adjacent and opposite element of
    the j-th face consecutively in the j-th column of the local-to-global
    table with DOF grouped by (internal) face.
  */

  DofByFace(const LocalDofPattern & _fe, UInt offSet = 1):
    fe(_fe),
    _offset(offSet)
    {
    }

  template<typename Mesh, typename DOF>
    void update(Mesh& mesh, DOF& dof);

  /*!
    Returns the global numbering of a DOF, given an internal face and the
    local numbering
    \param faceId the internal face ID
    \param localDOF the local DOF numbering (starting from 1)
    \return The global numbering of the DOF
  */
  inline ID localToGlobal(const ID faceId, const ID localDOF) const {
    return _ltg(localDOF, faceId);
  }

  inline UInt numLocalDofByFace() const {
    return _numLocalDofByFace;
  }

  inline UInt numTotalDof() const {
    return _numTotalDof;
  }

  //! Number of internal faces in the mesh
  UInt numIFaces() const {return _numIFaces;};

  //! Output
  void showMe(std::ostream& out = std::cout, bool verbose = false) const;

 private:
  UInt _offset;
  UInt _numTotalDof;
  UInt _numLocalDofByFace;
  UInt _numIFaces;
  Container _ltg;
};

//******************************************************************************
// Implementations
//******************************************************************************

template<typename Mesh, typename DOF>
void DofByFace::update(Mesh& mesh, DOF& dof){

  _numTotalDof = dof.numTotalDof();
  _numLocalDofByFace = 2 * dof.numLocalDof();

  UInt numBFaces = mesh.numBFaces();
  UInt numFaces = mesh.numFaces();

  _numIFaces = numFaces - numBFaces;

  UInt numLocalDof = dof.numLocalDof();

  ID iAd, iOp;

  _ltg.reshape(_numLocalDofByFace, _numIFaces);

  if(numFaces == numBFaces){
    std::cout << "Internal faces not stored" << std::endl;
    return;
  }else{
    for(UInt j = numBFaces + 1; j <= numFaces; j++){
      iAd = (ID)mesh.faceList(j).ad_first(); //! ID of the adjacent element
      iOp = (ID)mesh.faceList(j).ad_second(); //! ID of the opposite element

      for(UInt i = 1; i <= numLocalDof; i++){
	_ltg(i, j - numBFaces) = dof.localToGlobal(iAd, i);
	_ltg(i + numLocalDof, j - numBFaces) = dof.localToGlobal(iOp, i);

      }
    }
  }
}
}
#endif
