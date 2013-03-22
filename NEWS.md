LifeV 3.6.1:
============

Added Multiscale module. Cf documentation.

Increase the way the Darcy solvers handle the xml parameter lists.

Fixed bug in FSI:

* Removing spikes in the harmonic extension (within the implicit fsi_monolithic).
* Fixing a problem related to the assembly of the shape derivative
* Adding some assert to avoid wrong sets of parameters
* It now possible to choose P2 finite element for the harmonic extension (the mesh remains a P1 mesh ⇒ only the P1 nodes are moved! In this case a warning menton this fact).
* Removing a few warnings.



LifeV 3.4.4:
============
Bug fixes in ExporterVTK:

* the pvdu and vtu files names should be stored in the pvd (and pvdu)
  file(s) without the postDir. Otherwise paraview will look for them
  in postDir/postDir (path is relative to the leading file for paraview)
* the update of the CurrentFE (during the points map computation) should
  be done with the newer method, otherwise we don't store correctly the
  z-coordinate of 2D entities which do not lie on the xy plane (and
  similarly for 1D entities which do not lie on the x axis).

LifeV 3.4.3:
============

Removed old files from autotools.

LifeV 3.4.2:
============

Update of the Lifev-tutorial (FSI part unchanged)


LifeV 3.4.1:
============

* Rewriting the content of the LinearSolver class to be cleaner
* Bug fix (in some circumstances there were seg. faults)
* Adding a new "preconditioner" which is basically using a LinearSolver object
  as a preconditioner (to perform subiterations)
* Introduction of the concept of matrices and vectors structure… you can now
  store the structure as an object. This is useful to work with
  MatrixEpetraStructured.
* New functions in MatrixEpetraStructuredUtility.hpp to be able to exploit block
  structure using EpetraMatrix (typically if you do not have a
  MatrixEpetraStructured).
* New wrapper for operator (ConfinedOperator). You can apply an operator only on
   a chosen section of a vector.


LifeV 3.4.0:
============

Implementation of a physical solver for structural mechanics problems using an
abstract class to handle the constitutive laws that are currently available.
Moreover, these laws have been inserted in the FSI framework (partitioned and
monolithic) together with the TimeAdvance classes to have high order in time
discretization in the case of FSI problems.


LifeV 3.2.0:
============
* New modular structure and configuration using CMake similar to Trilinos

* Distribution through GitHub

New features are not described here for brevity.

LifeV  2.2.0:
=============

* This is a major release, the last one in the 2.x.x series

* For the next release, LifeV will be modularized, allower finer grained control over what parts of the library to install

* This release allows configuration either with the configure script, either with CMake. The CMake based configuration (undergoing a refactoring) will be the prefered configuration method and starting with the next release, the autoconf based configuration will be deprecated

* Check out the new features in LifeV: the monolithic FSI solver, the one dimensional model and solver, an experimental zero dimensional solver and others


LifeV  2.0.0:
=============

* Major release, the source code has been refactored in order to follow the naming and coding guidelines. Please consult the updated doxygen documentation for the new class and method names.

* This version of LifeV includes zero based numbering for all entities (user defined functions, boundary conditions, mesh entities etc.)

The following is some advice for porting applications to the new numbering scheme:

    In most cases, only a few things need to be changed when moving an application from one-based to zero-based numbering:
    1. user defined functions: now the component index (ic) starts from zero.
    2. adding boundary conditions: if you use "Component" mode, make sure that the component vector starts from zero.
    3. Loops, in particular the loop over mesh elements (volumes) start from zero.
    4. localToGlobal(...) function must be replaced with localToGlobalMap(..), which has zero-based entities both in input and output.

    If, after these changes, there are +1 or -1 here and there in the code, it is likely (but there are many exceptions) that there is something wrong.



LifeV  1.2.3:
=============

heart modeling has some errors that are fixed on the head.
It is now removed from 1.2.3 and will be in release 1.4.X

LifeV  1.2.2:
=============

This is the first release candidate of the parallel version of LifeV.

Recent new features includes:
* better doxygen documentation (although a lot of work is still needed in this direction)
* gcc 4.4.2 compatibility
* fixed distribution files


LifeV  1.1.1:
=============

New Features
------------

Distribution 1.2.x will the first public distribution
of the lifev-parallel.
The serial version is continuing on a separate cvs module, lifev-serial.


Obsolete (Branch OBSOLETE)
--------------------------
This are no longer in the head, nor distributed.
Developers can access to these files by cvs co -r OBSOLETE filename .
If the files are not obsolete anymore, they can be merged to the head.

 life/lifealg/PressureMatrixSolver.hpp  \
 life/lifefilters/medit.hpp  \
 life/lifesolver/darcyHandler.hpp  \
 life/lifesolver/darcySolverBase.hpp  \
 life/lifesolver/darcySolver.{hpp,cpp}  \
 life/lifesolver/dataDarcy.hpp  \
 life/lifefem/assembDG.hpp  \
 life/lifefem/refEleDG.{hpp,cpp}  \
 life/lifefem/refFEDG.{hpp,cpp}  \
 life/lifefem/dofDG.{hpp,cpp}  \
 life/lifefem/geoMapDG.{hpp,cpp}  \
 life/lifefem/refHdivFE.{hpp,cpp}  \
 life/lifefem/refHybridFE.{hpp,cpp}  \
 life/lifealg/cblas.hpp  \
 life/lifealg/clapack.hpp  \
 life/lifefem/currentBdDG.{hpp,cpp}  \
 life/lifefem/currentFEDG.{hpp,cpp}  \
 life/lifefem/currentHdivFE.{hpp,cpp}  \
 life/lifefem/currentIFDG.{hpp,cpp} \
 life/lifefem/defQuadRuleHybridFE.cpp


Removed (tagged REMOVED)
------------------------
These files are removed. The last tag that has them is LIFEV_1_1_0.

 life/lifearray/pattern.hpp  \
 life/lifearray/pattern.cpp  \
 life/lifealg/algebraic_facto.hpp  \
 life/lifealg/dataAztec.hpp  \
 life/lifealg/givens.hpp  \
 life/lifealg/gmres.hpp  \
 life/lifealg/inexactLU.hpp  \
 life/lifealg/iteration.hpp  \
 life/lifealg/preconditioner.hpp  \
 life/lifealg/saddlePointCG.hpp  \
 life/lifealg/SolverAztec.hpp  \
 life/lifealg/SolverBase.hpp  \
 life/lifealg/SolverPETSC.{hpp,cpp}  \
 life/lifealg/SolverUMFPACK.{hpp,cpp}  \
 life/lifealg/triDiagCholesky.hpp  \
 life/lifealg/triDiagLU.hpp  \
 life/lifefilters/ensight7Writer.hpp  \
 life/lifefilters/gmsh_wrtrs.hpp  \
 life/lifefilters/gmv_wrtrs.hpp  \
 life/lifefilters/openDX_wrtrs.{hpp,cpp}  \
 life/lifefilters/vtk_wrtrs.{hpp,cpp}  \
 life/lifearray/blockMatrix.hpp  \
 life/lifearray/boostmatrix.hpp  \
 life/lifearray/connectivity.hpp  \
 life/lifearray/CSRMatrix.hpp  \
 life/lifearray/MixedMatrix.hpp  \
 life/lifearray/MSRMatrix.hpp  \
 life/lifearray/sparseArray.{hpp,cpp}  \
 life/lifearray/VBRMatrix.hpp  \
 life/lifefem/assembGeneric.hpp  \
 life/lifefem/dataTransient.{hpp,cpp}  \
 life/lifefem/dofByFace.{hpp,cpp}  \
 life/lifefem/elemOper_ext.hpp  \
 life/lifefem/interpolate.hpp  \
 life/lifefem/Operator01.hpp  \
 life/lifefem/v2elemVec.{hpp,cpp}  \
 life/lifefem/values.hpp  \
 life/lifemesh/basicOneDMesh.{hpp,cpp}  \
 life/lifemesh/subDomainNeighbors.{hpp,cpp}  \
 life/lifesolver/AFSolvers.hpp  \
 life/lifesolver/convDiffReactHandler.hpp  \
 life/lifesolver/convDiffReactSolverPC.hpp  \
 life/lifesolver/dataConvDiffReact.hpp  \
 life/lifesolver/NavierStokesAleHandler.hpp  \
 life/lifesolver/NavierStokesAleSolverCT.hpp  \
 life/lifesolver/NavierStokesAleSolver.hpp  \
 life/lifesolver/NavierStokesAleSolverPC.hpp  \
 life/lifesolver/NavierStokesWithFlux.hpp  \
 life/lifesolver/nsip.hpp  \
 life/lifesolver/simplifiedStructure.hpp  \
 life/lifesolver/parabolicSolver.hpp  \
 life/lifesolver/timeSolver.{hpp,cpp}


Bug Fixes
---------
