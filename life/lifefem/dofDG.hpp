/*!
  \file dof.h
  \brief Degrees of freedom, the class that provides the localtoglobal table
  \version 1.0
  \author D. A. Di Pietro
  \date 10/2004
*/

#ifndef _DOFDG_HH
#define _DOFDG_HH

#include "life.hpp"
#include "SimpleVect.hpp"
#include "localDofPattern.hpp"
#include <algorithm>

namespace LifeV
{
/*! Local-to-global table

 This class provides the localtoglobal table that relates the local DOF of
 a finite element to its global numbering. It needs a LocalDofPattern in
 order to obtain all the necessary information about the local pattern. In
 fact it stores a copy of it so to make the local pattern available, if
 needed.

 It is useless until is has not been set up on a specific RegionMesh. This is accomplished either by
 passing the mesh to the constructor, or calling the method Dof::update().

 \note The methods bulds the table for ALL degrees of freedom, i.e. it does not handle any essential
 boundary condition.
*/
class DofDG {
public:
  //! Type for the localToGlobal table.
  typedef SimpleArray<UInt> Container;

  //! The pattern of the local degrees of freedom.
  /*! It is exposed so that is is possible to interrogate it directly.*/
  const LocalDofPattern& fe; // be careful : use fe.nbLocalDof (fe.nbDof does not exist !)

  /*! The minimal constructor
    \param _fe is the LocalDofPattern on which the ref FE is built
    \param Offset: the smallest Dof numbering. It might be used if we want the
    degrees of freedom numbering start from a specific value.
  */
  DofDG(const LocalDofPattern& _fe, UInt offSet=1);

  DofDG(const DofDG & dof2);

  //! Constructor accepting a mesh as parameter
  /*!
    \param mesh a RegionMesh3D
    \param _fe is the LocalDofPattern on which the ref FE is built
    \param Offset: the smalest Dof numbering. It might be used if we want the
    degrees of freedom numbering start from a specific value.
  */
  template<typename Mesh>
    DofDG(Mesh& mesh, const LocalDofPattern& _fe, UInt offSet=1);

  //! Build the localToGlobal table
  /*!
    \param mesh A RegionMesh3D
    Updates the LocaltoGlobal array
  */
  template <typename Mesh>
    void update(Mesh &);

  //! The total number of Dof
  inline UInt numTotalDof() const {return _totalDof;}

  //! The number of local Dof (nodes) in the finite element
  inline UInt numLocalDof() const {return fe.nbLocalDof;}

  //! Return the specified entries of the localToGlobal table
  /*!
    Returns the global numbering of a DOF, given an element and the local numbering
    \param ELId the element ID
    \param localNode the local DOF numbering (starting from 1)
    \return The numbering of the DOF
  */
  inline ID localToGlobal(const ID ElId, const ID localNode) const {
    return _ltg(localNode,ElId);
  }

  //! Number of elements in mesh
  UInt numElements() const {return _nEl;}

  //! Number of local vertices (in a elment)
  UInt numLocalVertices() const {return nlv;}

  //! Number of local edges (in a elment)
  UInt numLocalEdges() const {return nle;}

  //! Number of local faces (in a elment)
  UInt numLocalFaces() const {return nlf;}

  //! Ouput
  void showMe(std::ostream  & out=std::cout, bool verbose=false) const;

private:
  UInt _offset;
  UInt _totalDof;
  UInt _nEl;
  UInt nlv;
  UInt nle;
  UInt nlf;
  Container _ltg;
  UInt _ncount[5];
};



/********************************************************************
                      IMPLEMENTATIONS
********************************************************************/

//! Constructor that builds the localToglobal table
template <typename Mesh>
DofDG::DofDG(Mesh& mesh, const LocalDofPattern& _fe, UInt off):DofDG(_fe,off)
{update(mesh);}


//! Build the localToGlobal table
template<typename Mesh>
void DofDG::update(Mesh& M){

  typedef  typename Mesh::VolumeShape GeoShape;

  // Some useful local variables, to save some typing
  UInt nldpe=fe.nbDofPerEdge;
  UInt nldpv=fe.nbDofPerVertex;
  UInt nldpf=fe.nbDofPerFace;
  UInt nldpV=fe.nbDofPerVolume;

  nlv = GeoShape::numVertices;
  nle = GeoShape::numEdges;
  nlf = GeoShape::numFaces;

  _nEl = M.numElements();

  UInt nV = M.numVolumes();

  UInt i,l,ie;

  UInt nldof = nldpV + nldpe * nle + nldpv * nlv + nldpf * nlf;

  ASSERT_PRE( nldof == UInt(fe.nbLocalDof), "Something wrong in FE specification") ;

  _totalDof = _nEl * (nldpV + nle * nldpe + nlv * nldpv + nlf * nldpf);
  _ltg.reshape(nldof,nV);

  // Make sure the mesh has everything needed
  bool update_edges(nldpe !=0 && ! M.hasLocalEdges());
  bool update_faces(nldpf !=0 && ! M.hasLocalFaces());

  if (update_edges) M.updateElementEdges();
  if (update_faces) M.updateElementFaces();
  //  ASSERT_PRE( !(nldpe !=0 && M.hasLocalEdges()) , "Element edges stuff have not been updated") ;
  //  ASSERT_PRE( !(nldpf !=0 && M.hasLocalFaces()) , "Element faces stuff have not been updated") ;
  //ASSERT_PRE( (nldpe == 0 || M.hasLocalEdges()) , "Element edges stuff have not been updated") ;
  //ASSERT_PRE( (nldpf == 0 || M.hasLocalFaces()) , "Element faces stuff have not been updated") ;


  unsigned int gcount(_offset);
  unsigned int lcount = 0;
  unsigned int lc;

  // Vertex Based Dof
  _ncount[0]=gcount;

  if (nldpv >0 )
    for (ie=1; ie<= _nEl; ++ie){
      lc=0;
      for (i=1; i<=nlv; ++i)
	for (l=0; l<nldpv; ++l)
	  // Discontinuous Elements require that the degrees of freedom of geometrically corresponding
	  // nodes belonging to adjacent elements be decoupled
	  _ltg(++lc,ie) = gcount + (ie - 1) * nlv * nldpv + (i - 1) * nldpv + l;
    }

  // Edge Based Dof

  gcount += _nEl * nldpv * nlv;
  lcount += nldpv * nlv;
  _ncount[1]=gcount;

  if (nldpe >0 )
    for (ie=1; ie<= _nEl; ++ie){
      lc=lcount;
      for (i=1; i<=nle; ++i)
	for (l=0; l<nldpe; ++l)
	  _ltg(++lc,ie) = gcount + (ie - 1) * nle * nldpe + (i - 1) * nldpe + l;
    }

  // Face  Based Dof

  gcount += _nEl * nldpe* nle;
  lcount += nldpe * nle;

  _ncount[2]=gcount;
  if (nldpf >0 )
    for (ie=1; ie<= _nEl; ++ie){
      lc=lcount;
      for (i=1; i<=nlf; ++i)
	for (l=0; l<nldpf; ++l)
	  _ltg(++lc,ie) = gcount + (ie - 1) * nlf * nldpf + (i - 1) * nldpf + l;
    }

  // Volume  Based Dof
  gcount += _nEl * nldpf * nlf;
  lcount += nldpf * nlf;

  _ncount[3]=gcount;
  if (nldpV >0 )
    for (ie=1; ie<= _nEl; ++ie) {
      lc=lcount;
      for (l=0; l<nldpV; ++l)
	_ltg(++lc,ie)=gcount+(ie-1)*nldpV +l;
    }
  gcount += _nEl * nldpV;

  _ncount[4]=gcount;
  ASSERT_POS(gcount -_offset == _totalDof , "Something wrong in Dof Setup " << " gcount = " << gcount << " offset = " <<_offset <<" totalDof = " <<_totalDof) ;
  if (update_edges) M.cleanElementEdges();
  if (update_faces) M.cleanElementFaces();
}
}
#endif
