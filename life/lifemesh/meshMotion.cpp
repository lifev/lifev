#include "meshMotion.hpp"

// This method gives a reference to the computed harmonic extension.
Vector& HarmonicExtension
::getDisplacement() { 
  return _disp.vec();
}

const Dof& HarmonicExtension::
dofMesh() const {
  return _dof_mesh;
}
